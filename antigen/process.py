import datetime
import os
import warnings

from antigen import utils
from antigen.datasets import parse_fits_file_tree, get_matching_filenames
from antigen.datasets import get_elements_within_time_radius
from antigen.io import get_fits_file_time
from antigen.reduce_virus2 import process_unit


def build_manifest_records(infolder, obs_date, obs_name, reduce_all, time_radius,
                           bias_label, arc_label, dark_label, flat_label, twilight_flat_label):
    """
    Purpose: Search FITS file tree and generate a list of manifest records,
    one for each "obs id break block" in the files

    Args:
        infolder (str): Root path where reduction input file tree is located
        obs_date (str): Observation calendar date string formatted as YYYYMMDD
        obs_name (str): Observation object/target name, e.g. from FITS header card
        reduce_all (bool): Reduce all files found under infolder file tree
        time_radius (float) will group calibration files with times that fall within this distance from a given science file time
        bias_label (str): string label from FITS file header for bias frames
        arc_label (str): string label from FITS file header for arc frames
        dark_label (str): string label from FITS file header for dark frames
        flat_label (str): string label from FITS file header for flat frames
        twilight_flat_label (str): string label from FITS file header for twilight frames

    Returns:
        manifest_records (list(dict)): list of dicts/records, where each dict has the following keys:
            record['obs_time']: Obs time as MJD float
            record[observation_files']: filename of obs file to be reduced
            record['calibration_files']['bias']: list of bias FITS file names to use for calibration
            record['calibration_files']['flat']: list of flat FITS file names to use for calibration
            record['calibration_files']['arc']: list of arc FITS file names to use for calibration
    """

    ROOT_DATA_PATH = os.path.abspath(infolder)
    if not os.path.isdir(ROOT_DATA_PATH):
        raise NotADirectoryError(f'ERROR: user-specified input path does not exist: {ROOT_DATA_PATH}')

    # time filtering: Get FITS files in file-tree that match the obs_date
    metadata_records = parse_fits_file_tree(ROOT_DATA_PATH, date=obs_date, verbose=True)
    units_found = [record['spec_id'] for record in metadata_records]
    unique_units_found = list(set(units_found))

    manifest_records = []
    for unit in unique_units_found:
        # find subset of file records with matching spec_id, e.g. 'D3G' and get their filenames and frame types
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
            warnings.warn(f'unit={unit}, found ZERO matching files for obs_name={obs_name}. Continuing to next unit...')
            continue

        if reduce_all:
            calibration_files = set(bias_filenames) | set(twiflt_filenames) | set(domeflt_filenames) | set(arc_filenames) | set(dark_filenames)
            sci_filenames = [name for name in unit_filenames if name not in calibration_files]

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
            warnings.warn(f'WARNING: unit={unit}, Searched {ROOT_DATA_PATH}, found BIAS label = {bias_label}, '
                          f'found {num_bias_files}, needed >= {bias_minimum_count}')

        flat_minimum_count = minimum_file_count_for_break
        num_flt_files = len(flt_filenames)
        if num_flt_files < flat_minimum_count:
            fail_flat = True
            warnings.warn(f'WARNING: unit={unit}, Searched {ROOT_DATA_PATH}, FLAT label = {flat_label}, '
                          f'found {num_flt_files}, needed >= {flat_minimum_count}')

        arc_minimum_count = minimum_file_count_for_break
        num_arc_files = len(arc_filenames)
        if num_arc_files < arc_minimum_count:
            fail_arc = False
            warnings.warn(f'WARNING: unit={unit}, Searched {ROOT_DATA_PATH}, ARC label = {arc_label}, '
                          f'found {num_arc_files}, needed >= {arc_minimum_count}')

        if fail_bias or fail_flat or fail_arc:
            warnings.warn(f'WARNING: did not find enough calibration files to process unit={unit}. Continuing to next unit...')
            continue
        # =============================================================================
        # Use the filename obs_id numbers for grouping/splitting contiguous blocks of observations files
        # =============================================================================
        bias_times = [get_fits_file_time(name) for name in bias_filenames]
        flt_times  = [get_fits_file_time(name) for name in flt_filenames]
        arc_times  = [get_fits_file_time(name) for name in arc_filenames]

        # Generate manifest file for each science file
        print(f'Found {len(sci_filenames)} science files. '
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
                print(f'FAIL: found ZERO calibration files for sci_file={sci_file}: '
                      f'time_center={time_center}, time_radius={time_radius}')
                continue

            num_cals = len(bias_files) + len(flt_files) + len(arc_files)
            print(f'PASS: found {num_cals} calibration files for sci_file={sci_file}: '
                  f'time_center={time_center}, time_radius={time_radius}')

            record = dict()
            now_string = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
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
            manifest_records.append(record)

    return manifest_records


def process(infolder, outfolder, obs_date, obs_name, reduce_all, time_radius,
            bias_label, arc_label, dark_label, flat_label, twilight_flat_label):
    """
    Purpose: data reduction pipeline to process VIRUS2 observation files

    Args:
        infolder (str): Root path where reduction input file tree is located
        outfolder (str): Path where reduction output files will be written
        obs_date (str): Observation calendar date string formatted as YYYYMMDD
        obs_name (str): Observation object/target name, e.g. from FITS header card
        reduce_all (bool): Reduce all files found under infolder file tree
        bias_label (str): string label from FITS file header for bias frames
        arc_label (str): string label from FITS file header for arc frames
        dark_label (str): string label from FITS file header for dark frames
        flat_label (str): string label from FITS file header for flat frames
        twilight_flat_label (str): string label from FITS file header for twilight frames

    Returns:
        None
    """
    log = utils.setup_logging('virus2_reductions')

    manifest_records = build_manifest_records(infolder, obs_date, obs_name, reduce_all, time_radius,
                                              bias_label, arc_label, dark_label, flat_label, twilight_flat_label)

    os.makedirs(outfolder, exist_ok=True)

    for record in manifest_records:
        unit_id = record['unit_id']
        log.info(f'Processing reduction for unit = {unit_id}:')
        try:
            output_fits_filename = process_unit(record, outfolder)
            log.error(f'Processing reduction for unit = {unit_id}: PASS: wrote reduction to FITS file {output_fits_filename}')
        except Exception as error:
            log.error(f'Processing reduction for unit = {unit_id}: FAILED: {error}')

    return None
