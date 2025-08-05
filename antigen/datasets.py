import datetime as dt
import logging
import re
from pathlib import Path

import numpy as np

from antigen import io
from antigen.io import get_fits_file_time

logger = logging.getLogger('antigen.datasets')


def parse_fits_file_name(fits_filename, expected_prefix_parts=8, expected_extension='.fits'):
    """
    Purpose: Parse FITS filenames written by VIRUS2 exposure code
             Expects a suffix of '.fits' extension,
             Expects exactly 8 words in the stem,
             Expects and supports either dunder _ or dot . delimiters between stem words

    Note: Example expected filename pattern is
          ROOT_PATH/VIRUS2/20250618/0000001/D3G/VIRUS2_20250618_0000005_test_D3G_exp01_20250619T003023.0_test.fits
          VIRUS2_<obsdate>_<obsid>_<frametype>_<specid>_exp<exposureindex>_<utctime>_<userlabel>.fits

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
        file_names (list(dict)): List of files
    """

    logger.info(f'Searching for FITS files under in root_path={root_path} for date={date} ...')

    # Construct the path components
    virus_root_path = Path(root_path) / instrument
    date_dir = virus_root_path / (date if date else '*')

    # find all files matching this file-tree glob pattern
    fits_filenames = sorted(date_dir.glob('*/*/*.fits'))

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


def parse_file_dir_obs_id(file_name):
    """
    Purpose: Extract the obs ID string e.g. 0000001 from expected FITS obs filename pattern like the following:
        ROOT_PATH/VIRUS2/20250618/0000001/D3G/VIRUS2_20250618_0000005_test_D3G_exp01_20250619T003023.0_test.fits

    Args:
        file_name (str): VIRUS2 obs FITS filename

    Returns:
        dir_obs_id (int): e.g. 1 for directory name string '0000001'
    """
    dir_name = Path(file_name).parent.parent.name
    try:
        dir_obs_id = int(dir_name)
    except ValueError:
        dir_obs_id = None
    return dir_obs_id


def get_fits_file_obs_ids(fits_file_names):
    """
    Purpose: Extract the obs ID string e.g. 0000001 from ALL fits_file_names passed in, assumed to be found in a VIRUS2 obs FITS file-tree
    Note: Expected filename pattern is for each filename in the input list:
        ROOT_PATH/VIRUS2/20250618/0000001/D3G/VIRUS2_20250618_0000005_test_D3G_exp01_20250619T003023.0_test.fits
        VIRUS2_<obsdate>_<obsid>_<frametype>_<specid>_exp<exposureindex>_<utctime>_<userlabel>.fits

    Args:
        fits_file_names (list(str)): list of VIRUS2 obs FITS filenames

    Returns:
        obs_id_list (list(str)): list of string IDs for directory name string '0000001'

    Note: replaces lines like [int(os.path.basename(os.path.dirname(os.path.dirname(fn)))) for fn in bias_filenames]
    """
    obs_id_list = list()
    for filename in fits_file_names:
        filename_metadata = parse_fits_file_name(filename)
        if filename_metadata:
            obs_id_list.append(filename_metadata['obs_id'])
    return obs_id_list


def get_file_block_break_indices(fits_file_names):
    """
    Purpose: for a given list of VIRUS2 FITS files, which have integer obs-id in their full path names,
    find all the obs-ids, and find the location of the last ID before a "jump" or "break" larger than 1, e.g.
        In [1]: obs_id_list = [1,2,3,10]
        In [2]: np.where(np.diff(obs_id_list) > 1)[0]
        Out[3]: array([2])

    Args:
        fits_file_names (list(str)): list of FITS filepaths, assumed to contain a substring like "/00001/"

    Returns:
        obs_id_break_inds (np.ndarray): a numpy array of the indices of the right-edge/last-element in each contiguous block for filenames
            e.g. if obs_ids for the files are [1,2,3,11,12,13,21,22,23] then obs_id_break_inds = [2,5]
    """
    obs_id_list = get_fits_file_obs_ids(fits_file_names)
    obs_id_ints = [int(thing) for thing in obs_id_list]
    if len(obs_id_ints) > 1:
        breaks = np.where(np.diff(obs_id_ints) > 1)
        if len(breaks) > 0:
            obs_id_break_inds = breaks[0]
        else:
            obs_id_break_inds = None
    elif len(obs_id_ints) == 1:
        obs_id_break_inds = None
    else:
        raise ValueError(f'Failed to compute gaps in obs_id_list={obs_id_list}')
    return obs_id_break_inds


def get_element_with_closest_time(element_list, time_list, target_time):
    """
    Purpose: Given a list of elements and a list of MJD times for those elements,
    find the element that has a time closest to the target_time

    Args:
        element_list (list): list of elements corresponding to the times in time_list
        time_list (list(numeric)): list of MJD times corresponding to the times for the elements in element_list
        target_time (numeric): MJD time

    Returns:
        closest_element: element for closest_time from time_list compared to target_time
    """
    # Find index of time closest to target_time
    index = np.argmin(np.abs(np.array(time_list) - target_time))

    # Use the index directly in the original file list
    closest_element = element_list[index]
    return closest_element


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
    time_distances = np.abs(np.array(time_list) - time_center)
    mask_inside = (time_distances <= time_radius)
    elements_inside = np.array(element_list)[mask_inside].tolist()
    return elements_inside


def get_file_breakind_times(filenames, breakind):
    """
    Creates a list of master calibration images and corresponding times
    by splitting the input list of filenames at given indices.

    Args:
        filenames (list(list(str))): List of FITS file paths containing calibration data.
        breakind (list(int)): List of indices to split the filenames into different chunks.

    Returns
        chunked_file_names (list(str)): list of files, for each chunk.
        chunked_file_times (list(float)): mean observation time (MJD float) for each file chunk.
    """

    # Define break points for splitting the filenames into chunks
    breakind1 = np.hstack([0, breakind])  # Start indices for chunks
    breakind2 = np.hstack([breakind, len(filenames)+1])  # End indices for chunks

    chunked_file_names = []
    chunked_file_times = []

    # Iterate over the file-list chunks defined by breakind1 and breakind2
    for bk1, bk2 in zip(breakind1, breakind2):
        # Collect and preprocess frames within the current chunk
        chunk_files = [f for cnt, f in enumerate(filenames)
                       if ((cnt > bk1) * (cnt < bk2))]  # Only include files in the current file-list-chunk

        # Extract observation times (MJD) for frames in the current chunk
        chunk_times = [io.get_fits_file_time(filename) for filename in chunk_files]

        # Append the median frame and the mean time for the current chunk
        chunked_file_names.append(chunk_files)
        chunked_file_times.append(np.mean(chunk_times))

    return chunked_file_names, chunked_file_times


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


def find_datasets(in_folder, obs_date, obs_name, reduce_all, time_radius,
                  bias_label, arc_label, dark_label, flat_label, twilight_flat_label,
                  instrument='VIRUS2'):
    """
    Purpose: Search FITS file tree and generate groupings of calibration and science files by FITS header times,
    organized into dataset records to then later help build a dataset file manifest

    Args:
        in_folder (str): Root path where reduction input file tree is located
        obs_date (str): Observation calendar date string formatted as YYYYMMDD
        obs_name (str): Observation object/target name, e.g. from FITS header card
        reduce_all (bool): Reduce all files found under infolder file tree
        time_radius (float) will group calibration files with times that fall within this distance from a given science file time
        bias_label (str): string label from FITS file header for bias frames
        arc_label (str): string label from FITS file header for arc frames
        dark_label (str): string label from FITS file header for dark frames
        flat_label (str): string label from FITS file header for flat frames
        twilight_flat_label (str): string label from FITS file header for twilight frames
        instrument (str, optional): instrument to use, e.g. VIRUS2, defaults to VIRUS2

    Returns:
        dataset_records (list(dict)): list of dicts/records, where each dict has the following keys:
            record['obs_time']: Obs time as MJD float
            record[observation_files']: filename of obs file to be reduced
            record['calibration_files']['bias']: list of bias FITS file names to use for calibration
            record['calibration_files']['flat']: list of flat FITS file names to use for calibration
            record['calibration_files']['arc']: list of arc FITS file names to use for calibration
    """

    root_data_path = Path(in_folder).expanduser().resolve()
    if not root_data_path.is_dir():
        raise NotADirectoryError(f'ERROR: user-specified input path does not exist: {root_data_path}')

    # time filtering: Get FITS files in file-tree that match the obs_date
    metadata_records = parse_fits_file_tree(root_data_path, instrument=instrument, date=obs_date, verbose=True)
    units_found = [record['spec_id'] for record in metadata_records]
    unique_units_found = list(set(units_found))

    dataset_records = []
    for unit in unique_units_found:
        logger.info(f'Search for subset of file records with matching unit={unit} ...')
        unit_records = [record for record in metadata_records if record['spec_id'] == unit]

        # Group subsets of filenames into lists based on the frame types
        unit_frame_types = [record['frame_type'] for record in unit_records]
        unit_filenames = [record['filename'] for record in unit_records]

        bias_filenames    = get_matching_filenames(unit_filenames, unit_frame_types, [bias_label])
        twiflt_filenames  = get_matching_filenames(unit_filenames, unit_frame_types, [twilight_flat_label])
        domeflt_filenames = get_matching_filenames(unit_filenames, unit_frame_types, [flat_label])
        arc_filenames     = get_matching_filenames(unit_filenames, unit_frame_types, [arc_label])
        dark_filenames    = get_matching_filenames(unit_filenames, unit_frame_types, [dark_label])

        if obs_name is not None:
            sci_filenames = get_matching_filenames(unit_filenames, unit_frame_types, [obs_name])
        else:
            sci_filenames = []

        if len(sci_filenames) == 0 and not reduce_all:
            logger.warning(f'unit={unit}, found ZERO matching files for obs_name={obs_name}. Continuing to next unit...')
            continue
        else:
            logger.info(f'unit={unit}, Found len(sci_filenames)={len(sci_filenames)} matching files for obs_name={obs_name}')

        if reduce_all:
            calibration_files = set(bias_filenames) | set(twiflt_filenames) | set(domeflt_filenames) | set(arc_filenames) | set(dark_filenames)
            sci_filenames = [name for name in unit_filenames if name not in calibration_files]
            logger.info(f'reduce_all==True, expanding dataset to reduce all non_calibration_files')

        if len(twiflt_filenames) > 0:
            flt_filenames = twiflt_filenames
        else:
            flt_filenames = domeflt_filenames

        # =============================================================================
        # Validate number of files found before attempting to use diffs on file obs_ids
        # Use exceptions to exit process if needed frame-types file counts were not found
        # =============================================================================
        fail_bias = False
        fail_flat = False
        fail_arc = False
        minimum_file_count_for_break = 1

        bias_minimum_count = minimum_file_count_for_break
        num_bias_files = len(bias_filenames)
        if num_bias_files < bias_minimum_count:
            fail_bias = True
            logger.warning(f'Searched {root_data_path}, unit={unit}, found bias_label={bias_label}, '
                           f'found {num_bias_files}, needed >= {bias_minimum_count}')

        flat_minimum_count = minimum_file_count_for_break
        num_flt_files = len(flt_filenames)
        if num_flt_files < flat_minimum_count:
            fail_flat = True
            logger.warning(f'Searched {root_data_path}, unit={unit}, flat_label={flat_label}, '
                           f'found {num_flt_files}, needed >= {flat_minimum_count}')

        arc_minimum_count = minimum_file_count_for_break
        num_arc_files = len(arc_filenames)
        if num_arc_files < arc_minimum_count:
            fail_arc = False
            logger.warning(f'Searched {root_data_path}, unit={unit}, arc_label={arc_label}, '
                           f'found {num_arc_files}, needed >= {arc_minimum_count}')

        if fail_bias or fail_flat or fail_arc:
            logger.warning(f'Did not find enough calibration files to process unit={unit}. Continuing to next unit...')
            continue
        # =============================================================================
        # Use the filename obs_id numbers for grouping/splitting contiguous blocks of observations files
        # =============================================================================
        bias_times = [get_fits_file_time(name) for name in bias_filenames]
        flt_times  = [get_fits_file_time(name) for name in flt_filenames]
        arc_times  = [get_fits_file_time(name) for name in arc_filenames]

        # Generate manifest file for each science file
        logger.info(f'Found {len(sci_filenames)} science files. '
                    f'Attempting to find calibrations files withing time_radius of each ...')
        for sci_file in sci_filenames:
            fail_bias = False
            fail_flt = False
            fail_arc = False
            # For each science file, find the calibration files within a time radius
            obs_time_mjd = get_fits_file_time(sci_file)
            time_center = obs_time_mjd
            bias_files = get_elements_within_time_radius(bias_filenames, bias_times, time_center, time_radius)
            flt_files = get_elements_within_time_radius(flt_filenames, flt_times, time_center, time_radius)
            arc_files = get_elements_within_time_radius(arc_filenames, arc_times, time_center, time_radius)

            if len(bias_files) == 0:
                fail_bias = True
            if len(flt_files) == 0:
                fail_flt = True
            if len(arc_files) == 0:
                fail_arc = True

            if fail_bias or fail_arc or fail_flt:
                logger.warning(f'FAIL: found ZERO calibration files for sci_file={sci_file}: '
                               f'time_center={time_center}, time_radius={time_radius}')
                continue

            num_cals = len(bias_files) + len(flt_files) + len(arc_files)
            logger.info(f'PASS: found {num_cals} calibration files for sci_file={sci_file}: '
                        f'time_center={time_center}, time_radius={time_radius}')

            record = dict()
            now_string = dt.datetime.now().strftime('%Y%m%d_%H%M%S')
            record['reduction_name'] = f'antigen_manifest_{now_string}'
            record['unit_date'] = 'unknown'
            record['unit_instrument'] = 'VIRUS2'
            record['unit_id'] = unit
            record['obs_date'] = obs_date
            record['obs_name'] = obs_name
            record['in_folder'] = './'
            record['observation_files'] = sci_filenames
            record['calibration_files'] = dict()
            record['calibration_files']['bias'] = bias_filenames
            record['calibration_files']['flat'] = flt_filenames
            record['calibration_files']['arc'] = arc_filenames

            dataset_records.append(record)

    return dataset_records
