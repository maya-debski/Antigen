# -*- coding: utf-8 -*-

import datetime as dt
import glob
import os


def parse_fits_file_name(fits_filename):
    """
    Purpose: Parse FITS filenames written by VIRUS2 exposure code
    Note: Expected filename pattern is
        ROOT_PATH/VIRUS2/20250618/0000001/D3G/VIRUS2_20250618_0000005_test_D3G_exp01_20250619T003023.0_test.fits
        VIRUS2_<obsdate>_<obsid>_<frametype>_<specid>_exp<exposureindex>_<utctime>_<userlabel>.fits

    Args:
        fits_filename (str): full-path filename of FITS file containing VIRUS2 obs data

    Returns:
        filename_metadata (dict): keys = ['filename', 'instrument',
            'obs_date', 'obs_id', 'frame_type', spec_id', exp_index', 'utc_str', 'user_label']
    """
    filename_words = fits_filename.split('_')
    if len(filename_words) != 8:
        print(f'WARNING: Cannot parse filename, returning None; '
              f' Expected pattern of 8 words delimited by underscores. '
              f' filename={fits_filename}')
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


def parse_fits_file_tree(root_path, date=None, verbose=False):
    """
    Purpose: Parse all FITS filenames found in file-tree below given root_path
    Note: Expected filename pattern is
        ROOT_PATH/VIRUS2/20250618/0000001/D3G/VIRUS2_20250618_0000005_test_D3G_exp01_20250619T003023.0_test.fits
        VIRUS2_<obsdate>_<obsid>_<frametype>_<specid>_exp<exposureindex>_<utctime>_<userlabel>.fits

    Args:
        root_path (str): top-level of file-tree containing VIRUS2 exposure FITS files
        date (str): date string of format 'YYYYMMDD'
        verbose (bool): if True, print more info to console

    Returns:
        metadata_records (list(dict)): List of dictionaries with keys =  ['filename', 'instrument',
                 'obs_date', 'obs_id', 'frame_type', spec_id', exp_index', 'utc_str', 'user_label']
    """

    if verbose:
        print(f'VERBOSE: Searching for FITS files under in root_path={root_path} for date={date} ...')

    if date:
        fits_filenames = sorted(glob.glob(os.path.join(root_path, 'VIRUS2', date, '*', '*', '*.fits')))
    else:
        fits_filenames = sorted(glob.glob(os.path.join(root_path, 'VIRUS2', '*', '*', '*', '*.fits')))

    if verbose:
        num_files = len(fits_filenames)
        print(f'VERBOSE: Found {num_files} files under.')
        if num_files < 1:
            raise FileNotFoundError(f'ERROR: found no files to process. Exiting...')

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
    parent_dir = os.path.dirname(os.path.dirname(file_name))  # Go up two levels
    dir_name = os.path.basename(parent_dir)  # Get the directory name
    try:
        dir_obs_id = int(dir_name)
    except ValueError:
        dir_obs_id = None  # or raise, or log error
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
    """
    obs_id_list = list()
    for filename in fits_file_names:
        filename_metadata = parse_fits_file_name(filename)
        if filename_metadata:
            obs_id_list.append(filename_metadata['obs_id'])
    return obs_id_list
