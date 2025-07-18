import datetime as dt
import glob
import os
from pathlib import Path
import re

import numpy as np
from astropy.io import fits
from astropy.time import Time

from antigen import config


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
        fits_filename (str): full-path filename of FITS file containing VIRUS2 obs data
        expected_prefix_parts (int): Number of parts or words expected to be parsed from the filename stem/base, after stipping the ".fits" extension
        expected_extension (str): e.g. default = ".fits"

    Returns:
        filename_metadata (dict): keys = ['filename', 'instrument',
            'obs_date', 'obs_id', 'frame_type', spec_id', exp_index', 'utc_str', 'user_label']
    """
    path = Path(fits_filename)

    # Validate extension
    if path.suffix.lower() != expected_extension:
        raise ValueError(f"ERROR: Expected extension {expected_extension}, but got {path.suffix}, for fits_filename={fits_filename}")

    # Strip off the file name extension
    file_name_stem = path.stem

    # Split the file name stem into at most 8 parts, allows ONLY dunder delimiters,
    # but with the maxsplit, this prevents splitting the last word which is a user-word that can be literally anything
    # including the delimiters WITHIN the word.
    filename_words = re.split('_', file_name_stem, maxsplit=expected_prefix_parts)

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


def get_fits_filenames(root_path, date=None, verbose=False):
    """
    Purpose: Search file-tree below given root_path, find all FITS filenames that match VIRUS2 name pattern
    Note: Expected filename pattern is
        ROOT_PATH/VIRUS2/20250618/0000001/D3G/VIRUS2_20250618_0000005_test_D3G_exp01_20250619T003023.0_test.fits
        VIRUS2_<obsdate>_<obsid>_<frametype>_<specid>_exp<exposureindex>_<utctime>_<userlabel>.fits

    Args:
        root_path (str): top-level of file-tree containing VIRUS2 exposure FITS files
        date (str): date string of format 'YYYYMMDD'
        verbose (bool): if True, print more info to console

    Returns:
        file_names (list(dict)): List of files
    """

    if verbose:
        print(f'VERBOSE: Searching for FITS files under in root_path={root_path} for date={date} ...')

    # Construct the path components
    virus_root_path = os.path.join(root_path, 'VIRUS2')
    if date:
        date_dir = os.path.join(virus_root_path, date)
    else:
        date_dir = os.path.join(virus_root_path, '*')

    # Add wildcard subdirectories
    obs_id_dir = os.path.join(date_dir, '*')
    spec_id_dir = os.path.join(obs_id_dir, '*')
    fits_glob_pattern = os.path.join(spec_id_dir, '*.fits')

    # find all files matching this file-tree glob pattern
    fits_filenames = sorted(glob.glob(fits_glob_pattern))

    if verbose:
        num_files = len(fits_filenames)
        print(f'VERBOSE: Found {num_files} files under.')
        if num_files < 1:
            raise FileNotFoundError(f'ERROR: found no files matching fits_glob_pattern={fits_glob_pattern}. Exiting...')

    return fits_filenames


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
    fits_filenames = get_fits_filenames(root_path, date, verbose)

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
    if len(obs_id_list) > 1:
        obs_id_break_inds = np.where(np.diff(obs_id_list) > 1)[0]
    elif len(obs_id_list) == 1:
        obs_id_break_inds = obs_id_list[-1]
    else:
        raise ValueError(f'Failed to compute gaps in obs_id_list={obs_id_list}')
    return obs_id_break_inds


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


def write_fits(skysubrect_adv, skysubrect, specrect, errorrect, header, channel, outfolder):
    """
    Purpose: Writes the sky-subtracted, rectified spectra and error data to a FITS file,
    preserving the header information and adding necessary meta-information.

    Args:
        skysubrect_adv (np.ndarray): 2D numpy array, The advanced sky-subtracted spectrum.
        skysubrect (np.ndarray): 2D numpy array, The basic sky-subtracted spectrum.
        specrect (np.ndarray): 2D numpy array, The rectified spectrum.
        errorrect (np.ndarray): 2D numpy array, The error associated with the rectified spectrum.
        header (fits.Header): The header information to be preserved in the output FITS file.
        channel (str): name of frequency channel, options = ['g', 'b', 'r', 'd']
        outfolder (str): existing path of directory to write the FITS file to

    Returns:
        None
    """
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder, exist_ok=True)

    hdulist = []  # List to store HDU objects for the FITS file

    # Loop through the data arrays and create HDUs for each
    for image, ftp in zip([skysubrect_adv, skysubrect, specrect, errorrect],
                          [fits.PrimaryHDU, fits.ImageHDU, fits.ImageHDU, fits.ImageHDU]):

        # Create an HDU object from each image, setting it to 'float32' type
        hdu = ftp(np.array(image, dtype='float32'))

        # Remove any conflicting CD matrix elements first
        for key in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CDELT1', 'CDELT2']:
            if key in hdu.header:
                del hdu.header[key]

        # Define your wavelength solution
        _, CHANNEL_DEF_WAVE = config.get_channel_config_virus2()
        def_wave = CHANNEL_DEF_WAVE[channel]
        wavelength_step = def_wave[1] - def_wave[0]  # Compute wavelength step

        # Set WCS parameters correctly
        hdu.header['CRVAL1'] = def_wave[0]  # First wavelength (Angstroms)
        hdu.header['CRPIX1'] = 1  # First pixel corresponds to first wavelength
        hdu.header['CD1_1'] = wavelength_step  # Set CD1_1 to match wavelength step
        hdu.header['CTYPE1'] = 'WAVE'  # Spectral axis label

        # Set fiber axis metadata
        hdu.header['CRVAL2'] = 1  # Reference value for fiber index
        hdu.header['CRPIX2'] = 1  # First pixel for fiber axis
        hdu.header['CD2_2'] = 1  # Step of 1 fiber per index
        hdu.header['CTYPE2'] = 'FIBER'  # Labeling fiber axis

        # Copy relevant keys from the input header, avoiding duplicates
        for key in header.keys():
            if key in hdu.header:
                continue
            if ('CCDSEC' in key) or ('DATASEC' in key):  # Exclude CCDSEC and DATASEC keys
                continue
            if ('BSCALE' in key) or ('BZERO' in key):  # Exclude BSCALE and BZERO keys
                continue
            try:
                hdu.header[key] = header[key]  # Copy header data to the new HDU
            except:
                continue

        # Format the output filename using observation date and object name
        t = Time(header['DATE-OBS'] + 'T' + header['UT'])
        objname = '_'.join(header['OBJECT'].split())
        iname = '_'.join([objname, t.strftime('%Y%m%dT%H%M%S'), 'multi'])  # Generate filename

        # Append the HDU to the list
        hdulist.append(hdu)

    # TODO: iname, accidental scope spill from for loop to function scope; replace with intentional assignment
    # Write the HDU list to the output file, overwriting if necessary
    fits.HDUList(hdulist).writeto(os.path.join(outfolder, iname + '.fits'), overwrite=True)

    return None
