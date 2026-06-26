import datetime as dt
import logging
import re
from pathlib import Path

import numpy as np

from antigen.io import get_fits_file_time, get_fits_files_in_path
from antigen import config

logger = logging.getLogger('antigen.datasets')


def _match_user_label(user_label: str, any_contains=None, must_contain=None, must_not_contain=None) -> bool:
    """Return True if user_label satisfies provided substring constraints.

    - any_contains: pass if empty or any term appears (case-insensitive)
    - must_contain: pass if empty or all terms appear (case-insensitive)
    - must_not_contain: fail if any term appears (case-insensitive)
    """
    ul = (user_label or '').lower()
    any_contains = [s.lower() for s in (any_contains or []) if s]
    must_contain = [s.lower() for s in (must_contain or []) if s]
    must_not_contain = [s.lower() for s in (must_not_contain or []) if s]

    if any_contains:
        if not any(s in ul for s in any_contains):
            return False
    if must_contain:
        if not all(s in ul for s in must_contain):
            return False
    if must_not_contain:
        if any(s in ul for s in must_not_contain):
            return False
    return True


def select_files_by_calibration_config(unit_records, selection_cfg):
    """Select calibration files using config-driven rules.

    Args:
        unit_records (list[dict]): Records as returned by parse_fits_file_tree.
        selection_cfg (dict): calibration_selection section from YAML.

    Returns:
        (cal_dict, sci_files): where cal_dict has keys bias/flat/arc (if present), and sci_files is list of filenames
                               inferred as non-calibration frames.
    """
    if not selection_cfg or not isinstance(selection_cfg, dict):
        return None, None

    # Normalize keys and pre-extract lists
    cal_keys = [k for k in selection_cfg.keys() if k in ('bias', 'flat', 'arc', 'dark', 'twilight', 'twilight_flat')]

    # Initialize calibration dict including optional 'dark'
    cal_dict = {k: [] for k in ('bias', 'flat', 'arc', 'dark')}
    selected_set = set()

    for cal_type in cal_keys:
        rules = selection_cfg.get(cal_type, {}) or {}
        ft_any = [s.lower() for s in rules.get('frametype_any', [])]
        any_contains = rules.get('userlabel_any_contains', []) or rules.get('userlabel_must_contain', []) or []
        must_contain = rules.get('userlabel_must_contain_all', [])  # optional, if provided
        must_not = rules.get('userlabel_must_not_contain', [])

        # Build list for this type
        matched = []
        for rec in unit_records:
            ftype = (rec.get('frame_type') or '').lower()
            ulabel = rec.get('user_label') or ''
            if ft_any and ftype not in ft_any:
                continue
            if not _match_user_label(ulabel, any_contains=any_contains, must_contain=must_contain, must_not_contain=must_not):
                continue
            matched.append(rec['filename'])

        # Map cal_type aliases to manifest keys
        key_map = {
            'bias': 'bias',
            'flat': 'flat',
            'twilight': 'flat',
            'twilight_flat': 'flat',
            'arc': 'arc',
            'dark': 'dark',  # include optional darks in manifest
        }
        manifest_key = key_map.get(cal_type)
        if manifest_key:
            cal_dict[manifest_key].extend(matched)
        selected_set.update(matched)

        # Optional warnings
        if rules.get('warn_on_empty_selection') and len(matched) == 0:
            logger.warning(f"Config selection for '{cal_type}' returned 0 files.")

    # Science files = those not selected as calibrations
    filenames = [rec['filename'] for rec in unit_records]
    sci_files = [f for f in filenames if f not in selected_set]

    return cal_dict, sci_files


def parse_reduce_file_name(fits_filename):
    """
    Parse a FITS file name into components.

    Example file names:
        reduction_HD128_dither_2_20240606T051235_gcms_VP1B_multi.fits
        reduction_HD141_20240606T051235_virus2_D3G_multi.fits
    Args:
        fits_filename (pathlike): full-path filename of reduced fits file (pathlib.Path object)
    Returns:
        info (dict): keys = ['target', 'dither_num', 'obs_date', 'instrument', 'element', file_name']
    """

    pattern = (
        r"reduction_(?P<target>[A-Za-z0-9+_\-]+?)"
        r"(?:_dither_(?P<dither>\d+))?"
        r"_(?P<date>\d{8}T\d{6})"
        r"_(?P<instrument>[a-zA-Z0-9]+)"
        r"_(?P<element>[a-zA-Z0-9]+)_multi"
    )

    path = Path(fits_filename).expanduser()
    file_name_stem = path.stem

    match = re.match(pattern, file_name_stem)
    if not match:
        return None

    info = match.groupdict()
    d = info.pop("dither", None)
    info["dither_number"] = int(d) if d else 1
    info["file_name"] = Path(fits_filename).name
    return info

def parse_fits_file_name(fits_filename, expected_prefix_parts=8, expected_extension='.fits'):
    """
    Purpose: Parse FITS filenames written by VIRUS2 exposure code
             Expects a suffix of '.fits' extension,
             Expects exactly 8 words in the stem,
             Expects and supports either dunder _ or dot . delimiters between stem words

    Note: Example expected filename pattern is
          ROOT_PATH/VIRUS2/20250618/0000001/D3G/VIRUS2_20250618_0000005_test_D3G_exp01_20250619T003023.0_test.fits
          <instrument>_<obsdate>_<obsid>_<frametype>_<specid>_exp<exposureindex>_<utctime>_<userlabel>.fits

    Args:
        fits_filename (str, pathlike): full-path filename of FITS file containing VIRUS2 obs data, string or pathlib.Path object
        expected_prefix_parts (int): Number of parts or words expected to be parsed from the filename stem/base, after stripping the ".fits" extension
        expected_extension (str): e.g. default = ".fits"

    Returns:
        filename_metadata (dict): keys = ['filename', 'instrument',
            'obs_date', 'obs_id', 'frame_type', spec_id', exp_index', 'utc_str', 'user_label']
    """
    path = Path(fits_filename).expanduser()

    # Validate extension
    if path.suffix.lower() != expected_extension:
        raise ValueError(f"ERROR: Expected extension {expected_extension}, but got {path.suffix}, for fits_filename={fits_filename}")

    # Strip off the file name extension
    file_name_stem = path.stem

    # Split the file name stem into at most 8 parts, allows ONLY dunder delimiters,
    # but with the maxsplit, this prevents splitting the last word which is a user-word that can be literally anything
    # including the delimiters WITHIN the word.
    number_of_delimiters_to_split = expected_prefix_parts - 1
    filename_words = re.split('_', file_name_stem, maxsplit=number_of_delimiters_to_split)

    if len(filename_words) != expected_prefix_parts:
        logger.warning(f'WARNING: Cannot parse filename, returning None; '
                       f' Expected pattern of 8 words delimited by underscores. '
                       f' file_name_stem={file_name_stem}')
        filename_metadata = None
    else:
        filename_metadata = dict()
        filename_metadata['filename']   = fits_filename
        filename_metadata['instrument'] = filename_words[0]
        filename_metadata['obs_date']   = filename_words[1]
        filename_metadata['obs_id']     = filename_words[2]
        filename_metadata['frame_type'] = filename_words[3]
        filename_metadata['spec_id']    = filename_words[4]
        filename_metadata['exp_index']  = filename_words[5]
        filename_metadata['utc_str']    = filename_words[6]
        filename_metadata['user_label'] = filename_words[7]

        utc_str = filename_metadata['utc_str']
        try:
            obs_datetime = dt.datetime.strptime(utc_str, "%Y%m%dT%H%M%S.%f")
            filename_metadata['utc_str_date'] = obs_datetime.strftime("%Y-%m-%d")
            filename_metadata['utc_str_time'] = obs_datetime.strftime("%H:%M:%S")
        except ValueError:
            print(f'WARNING: could not parse UTC datetime string: {utc_str}')
            filename_metadata['utc_str_date'] = None
            filename_metadata['utc_str_time'] = None

    return filename_metadata


def get_fits_filenames(root_path, instrument='VIRUS2', date=None, verbose=False):
    """
    Purpose: Search file-tree below given root_path, find all FITS filenames that match VIRUS2 name pattern
    Note: Expected filename pattern is
        ROOT_PATH/VIRUS2/20250618/0000001/D3G/VIRUS2_20250618_0000005_test_D3G_exp01_20250619T003023.0_test.fits
        VIRUS2_<obsdate>_<obsid>_<frametype>_<specid>_exp<exposureindex>_<utctime>_<userlabel>.fits

    Args:
        root_path (Path): top-level of file-tree containing VIRUS2 exposure FITS files
        instrument (str, optional): e.g. VIRUS2, defaults to VIRUS2
        date (str, optional): date string of format 'YYYYMMDD', defaults to None
        verbose (bool, optional): if True, print more info to console, defaults to False

    Returns:
        file_names (list(str)): List of file name strings
    """

    logger.info(f'Searching for FITS files under in root_path={root_path} for date={date} ...')

    # Construct the path components
    virus_root_path = Path(root_path) / instrument
    date_dir = virus_root_path / (date if date else '*')

    # find all files matching this file-tree glob pattern
    fits_filenames = sorted(date_dir.glob('*/*/*.fits'))

    # Converting list of Paths to list of strings
    fits_filenames = [str(filename) for filename in fits_filenames]
    if verbose:
        num_files = len(fits_filenames)
        logger.info(f'VERBOSE: Found {num_files} files under.')
        if num_files < 1:
            raise FileNotFoundError(
                f'ERROR: found no files matching pattern: {date_dir}/<obsid>/<specid>/*.fits. Exiting...')
    return fits_filenames


def parse_fits_file_tree(root_path, instrument='VIRUS2', date=None, verbose=False):
    """
    Purpose: Parse all FITS filenames found in file-tree below given root_path
    Note: Expected filename pattern is
        ROOT_PATH/VIRUS2/20250618/0000001/D3G/VIRUS2_20250618_0000005_test_D3G_exp01_20250619T003023.0_test.fits
        VIRUS2_<obsdate>_<obsid>_<frametype>_<specid>_exp<exposureindex>_<utctime>_<userlabel>.fits

    Args:
        root_path (Path): top-level of file-tree containing VIRUS2 exposure FITS files
        instrument (str, optional): e.g. VIRUS2, defaults to VIRUS2
        date (str, optional): date string of format 'YYYYMMDD', defaults to None
        verbose (bool, optional): if True, print more info to console, defaults to False

    Returns:
        metadata_records (list(dict)): List of dictionaries with keys =  ['filename', 'instrument',
                 'obs_date', 'obs_id', 'frame_type', spec_id', exp_index', 'utc_str', 'user_label']
    """
    fits_filenames = get_fits_filenames(root_path, instrument, date, verbose)

    metadata_records = list()
    for filename in fits_filenames:
        filename_metadata = parse_fits_file_name(filename)
        if filename_metadata:
            metadata_records.append(filename_metadata)

    return metadata_records


def get_elements_within_time_radius(element_list, time_list, time_center, time_radius):
    """
    Purpose: Given a list of elements and a list of MJD times for those elements,
    find the elements that fall within a time_window centered on a time_center

    Args:
        element_list (list): list of elements corresponding to the times in time_list
        time_list (list(numeric)): list of MJD times corresponding to the times for the elements in element_list
        time_center (numeric): MJD time
        time_radius (numeric): a delta which will be used to search for elements that have a time inside [time_center - dt, time_center + dt]

    Returns:
        elements_inside: elements within the time_radius of the time_center
    """
    elements_inside = [
        element for element, time in zip(element_list, time_list)
        if abs(time - time_center) <= time_radius
    ]
    return elements_inside


def get_matching_filenames(file_name_list, type_list, match_keywords):
    """
    Purpose: Finds filenames that match a list of keywords by checking if any of the keywords
    appear in the associated types.

    Args:
        file_name_list (list(str)): List of filenames to check.
        type_list (list(str)): List of types or descriptions corresponding to the filenames.
        match_keywords (list(str)): List of keywords to search for within the types/descriptions.

    Returns:
        matched_filenames (list(str)): list of filenames from file_name_list that match any of the match_keywords.
    """
    matched_filenames = []
    for file_name, type_name in zip(file_name_list, type_list):
        for word in match_keywords:
            if word.lower() in str(type_name).lower():
                matched_filenames.append(file_name)
    return matched_filenames


def check_file_count(label, filenames, minimum_count, root_data_path, unit):
    """
    Check whether a list of files meets the minimum required count, and log a warning if not.

    Args:
        label (str): Descriptive label for the file type (e.g., 'bias_label=BIAS').
        filenames (list): List of filenames to check.
        minimum_count (int): Minimum number of files required.
        root_data_path (Path or str): Path where the files were searched.
        unit (str): Identifier for the current instrument or observation unit.

    Returns:
        bool: True if the number of files is below the required minimum (i.e., check failed), False otherwise.
    """
    num_files = len(filenames)
    if num_files < minimum_count:
        logger.warning(
            f"Searched {root_data_path}, unit={unit}, {label} "
            f"found {num_files}, needed >= {minimum_count}"
        )
        return True
    return False


def build_calibration_groups(calibrations, time_radius):
    """
    Groups calibration frames (bias, arc, flat) that are within `time_radius` MJD of each other.

    Args:
        calibrations (list of dict): Each dict has 'name', 'mjd', and 'type'
        time_radius (float): time window (in MJD days) for grouping

    Returns:
        list of dict: calibration groups with lists of biases, arcs, and flats
    """
    remaining = calibrations.copy()
    groups = []

    while remaining:
        # Start with the first remaining time
        first = remaining[0]
        first_time = first['mjd']
        grouped = get_elements_within_time_radius(remaining, [cal['mjd'] for cal in remaining],
                                                  first_time, time_radius)

        # Iterate with the new average time (once should be enough, but you could iterate until no change occurs)
        new_average_time = np.mean([cal['mjd'] for cal in grouped])
        grouped = get_elements_within_time_radius(remaining, [cal['mjd'] for cal in remaining],
                                                  new_average_time, time_radius)

        final_average_time = np.mean([cal['mjd'] for cal in grouped])
        # Organize by type (including optional 'dark')
        group = {'bias': [], 'arc': [], 'flat': [], 'dark': []}
        for cal in grouped:
            typ = cal['type']
            if typ not in group:
                # In case of unexpected type, initialize dynamically
                group[typ] = []
            group[typ].append(cal['name'])

        # Keep the final average time of the cal group
        group['mjd'] = final_average_time

        # Remove used calibrations
        used_names = {cal['name'] for cal in grouped}
        remaining = [cal for cal in remaining if cal['name'] not in used_names]
        groups.append(group)

    return groups


def assign_science_to_groups(science_frames, cal_groups):
    """
    Assign science frames to the closest matching calibration group.

    Args:
        science_frames (list of dict): Each has 'name' and 'mjd'
        cal_groups (list of dict): Each group has 'bias', 'arc', 'flat' lists and 'mjd' of the central time

    Returns:
        cal_groups (list of dict): Each dict contains a calibration group and a list of matched science frames.
    """
    for group in cal_groups:
        group['sci'] = []

    for sci in science_frames:
        sci_time = sci['mjd']

        # Find group center times
        group_times = [group['mjd'] for group in cal_groups]

        closest_idx = np.argmin(np.abs(np.array(group_times) - sci_time))
        cal_groups[closest_idx]['sci'].append(sci['name'])

    return cal_groups

def get_calibration_and_science_files(unit_records, obs_name,
                                      bias_label, flat_label, twilight_flat_label,
                                      arc_label, dark_label, selection_cfg=None):
    """
    Categorize files into calibration and science types.

    Args:
        unit_records (list): Records for a unit.
        obs_name (str): Science target name.
        bias_label, flat_label, twilight_flat_label, arc_label, dark_label (str): Frame type labels.

    Returns:
        cal_files_dict (dict): calibration files dictionary
        sci_files (list): science filenames
    """
    # If a selection config is provided, use it
    if selection_cfg:
        cal_by_cfg, sci_by_cfg = select_files_by_calibration_config(unit_records, selection_cfg)
        if cal_by_cfg is not None:
            logger.info('Using calibration_selection rules from config to choose files.')
            return cal_by_cfg, sci_by_cfg

    # Fallback: legacy behavior using label keywords
    filenames = [rec['filename'] for rec in unit_records]
    types = [rec['frame_type'] for rec in unit_records]

    bias = get_matching_filenames(filenames, types, [bias_label])
    twi = get_matching_filenames(filenames, types, [twilight_flat_label])
    dome = get_matching_filenames(filenames, types, [flat_label])
    arc = get_matching_filenames(filenames, types, [arc_label])
    dark = get_matching_filenames(filenames, types, [dark_label])

    flats = twi if twi else dome

    if obs_name:
        sci_files = get_matching_filenames(filenames, types, [obs_name])
    else:
        all_cal = set(bias + flats + arc + dark)
        sci_files = [filename for filename in filenames if filename not in all_cal]

    cal_files_dict = {'bias': bias, 'flat': flats, 'arc': arc, 'dark': dark}
    return cal_files_dict, sci_files


def validate_calibration_counts(cal_files, unit, root_data_path, minimum=1):
    """
    Ensure each calibration type has at least the minimum required files.

    Args:
        cal_files (dict): Dictionary of calibration file lists.
        unit (str): Unit name.
        root_data_path (Path): Root directory.
        minimum (int): Minimum required files per type.

    Returns:
        bool: True if all checks pass.
    """
    fail = False
    for key, files in cal_files.items():
        logger.info(f'Checking {key} files for unit={unit} in {root_data_path}, found {len(files)} files')
        # Dark frames are optional; do not fail if none are present
        if key == 'dark':
            continue
        if len(files) < minimum:
            logger.warning(f'Searched {root_data_path}, unit={unit}, found {len(files)} < {minimum} {key} files')
            fail = True
    return not fail


def prepare_group_input_dicts(cal_files, sci_filenames):
    """
    Build calibration and science input lists for grouping.

    Args:
        cal_files (dict): Calibration file dict.
        sci_filenames (list): Science file names.

    Returns:
        cal_dict_list (list of dicts): Calibration list of dictionaries with type, filename, and mjd.
        sci_dict_list (list of dicts): Science file list of dictionaries.
    """
    cal_dict_list = []
    for typ in ['bias', 'flat', 'arc', 'dark']:
        if typ in cal_files:
            times = [get_fits_file_time(f) for f in cal_files.get(typ, [])]
            cal_dict_list += [{'type': typ, 'mjd': t, 'name': f} for t, f in zip(times, cal_files.get(typ, []))]

    sci_times = [get_fits_file_time(f) for f in sci_filenames]
    sci_dict_list = [{'type': 'sci', 'mjd': t, 'name': f} for t, f in zip(sci_times, sci_filenames)]

    return cal_dict_list, sci_dict_list


def build_dataset_records(cal_groups, instrument, unit, obs_date, obs_name):
    """
    Build manifest-style dataset records.

    Args:
        cal_groups (list): Groups of calibration + science files.
        instrument (str): Instrument name.
        unit (str): Unit identifier.
        obs_date (str): Observation date.
        obs_name (str): Science target name.

    Returns:
        list of dict: Dataset records.
    """
    records = []
    now_string = dt.datetime.now().strftime('%Y%m%d_%H%M%S')

    for group in cal_groups:
        record = {
            'reduction_name': f'antigen_manifest_{now_string}',
            'unit_date': 'unknown',
            'unit_instrument': instrument,
            'unit_id': unit,
            'obs_date': obs_date,
            'obs_name': obs_name,
            'in_folder': './',
            'observation_files': group['sci'],
            'calibration_files': {
                'bias': group['bias'],
                'flat': group['flat'],
                'arc': group['arc'],
                'dark': group.get('dark', [])
            }
        }
        records.append(record)

    return records


def build_dataset_from_reduced_files(file_directory):
    """
    Build a dataset record grouped by target
    Args:
        file_directory: Directory path of the reduced fits files.

    Returns:
        list of dict: Dataset records.
    """
    file_list = get_fits_files_in_path(file_directory)  # List of path objects
    now_string = dt.datetime.now().strftime('%Y%m%d_%H%M%S')
    reduction_name = f"antigen_advanced_manifest_{now_string}"

    records_by_target = {}
    for file_name in file_list:
        info = parse_reduce_file_name(file_name)

        if info is None:
            continue

        target = info["target"]
        if target not in records_by_target:
            records_by_target[target] = {
                "reduction_name": reduction_name,
                "target": target,
                "unit_instrument": info["instrument"],
                "unit_id": info["element"],
                "in_folder": file_directory,
                "ndithers": 0,  # will update later
                "reduced_files": [],
                "dither_number": [],
                # Will set 'dither_file' after ndithers is known
            }

        record = records_by_target[target]
        record["reduced_files"].append(info["file_name"])
        record["dither_number"].append(info["dither_number"])
        # Ensure ndithers reflects unique dithers
        record["ndithers"] = len(set(record["dither_number"]))

    # After collecting all, compute dither_file paths per target using same logic as fiber.load_fiber_positions
    base_path = config.get_base_config_path()
    for target, rec in records_by_target.items():
        instr = str(rec.get("unit_instrument", "")).lower()
        ndith = int(rec.get("ndithers", 1) or 1)
        # Build default path; if ndithers==1 we still record the 1pt file if present
        try:
            dither_path = (Path(base_path) / instr / f"{instr}_dither_{ndith}pt.lis").expanduser().resolve()
            rec["dither_file"] = str(dither_path)
        except Exception:
            rec["dither_file"] = None

    logger.info(f"Built dataset with {len(records_by_target)} targets and {len(file_list)} files.")
    return list(records_by_target.values())


def find_datasets(in_folder, obs_date, obs_name, time_radius,
                  bias_label, arc_label, dark_label, flat_label, twilight_flat_label,
                  instrument='VIRUS2'):
    """
    Search FITS file tree and generate groupings of calibration and science files by FITS header times.

    Args:
        in_folder (str): Root path to search.
        obs_date (str): Observation date in YYYYMMDD.
        obs_name (str): Object name to search for in science frames.
        time_radius (float): MJD time radius for grouping science with calibration.
        bias_label (str): Label for bias frames.
        arc_label (str): Label for arc frames.
        dark_label (str): Label for dark frames.
        flat_label (str): Label for dome flat frames.
        twilight_flat_label (str): Label for twilight flat frames.
        instrument (str): Instrument name (default 'VIRUS2').

    Returns:
        dataset_records (list of dict): Dataset records for manifest file
    """
    root_data_path = Path(in_folder).expanduser().resolve()
    if not root_data_path.is_dir():
        raise NotADirectoryError(f'Input path does not exist: {root_data_path}')

    metadata_records = parse_fits_file_tree(root_data_path, instrument=instrument, date=obs_date, verbose=True)
    unique_units_found = list(set([record['spec_id'] for record in metadata_records]))
    dataset_records = []

    for unit in unique_units_found:
        unit_records = [record for record in metadata_records if record['spec_id'] == unit]

        # Try to load calibration selection from config for this instrument/element
        selection_cfg = None
        try:
            cfg_yaml = config.read_config_yaml(instrument.lower(), unit.upper(), validate=False)
            selection_cfg = cfg_yaml.get('calibration_selection') if isinstance(cfg_yaml, dict) else None
            if selection_cfg:
                logger.info(f"Loaded calibration_selection from config for {instrument}/{unit}.")
        except Exception as e:
            logger.debug(f"Could not load calibration_selection for {instrument}/{unit}: {e}")

        cal_files, sci_filenames = get_calibration_and_science_files(
            unit_records, obs_name,
            bias_label, flat_label, twilight_flat_label, arc_label, dark_label,
            selection_cfg=selection_cfg)

        if not sci_filenames:
            logger.warning(f'unit={unit}, no matching science files')

        if not validate_calibration_counts(cal_files, unit, root_data_path):
            logger.warning(f'unit={unit}, calibration files insufficient. Skipping.')

        cal_dict_list, sci_dict_list = prepare_group_input_dicts(cal_files, sci_filenames)

        cal_groups = build_calibration_groups(cal_dict_list, time_radius)
        cal_groups = assign_science_to_groups(sci_dict_list, cal_groups)

        dataset_records += build_dataset_records(cal_groups, instrument, unit, obs_date, obs_name)
    return dataset_records