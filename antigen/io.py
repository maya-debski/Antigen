import datetime as dt
import glob
import os
from pathlib import Path
import psutil
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
        print(f'WARNING: Cannot parse filename, returning None; '
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


def get_fits_file_time(fits_file_name):
    """
    Purpose: Read FITS header, get header cards 'DATE-OBS' and 'UT', construct a time stamp, convert it to MJD

    Args:
        fits_file_name (str): file name for FITS file
    Returns:
        time_mjd (float): MJD time

    """
    with fits.open(fits_file_name) as hdul:
        header = hdul[0].header
    timestamp = header['DATE-OBS'] + 'T' + header['UT']
    obs_time = Time(timestamp, format='isot', scale='utc')
    mjd = obs_time.mjd
    return mjd


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
        chunk_times = [Time(fits.open(filename)[0].header['DATE-OBS']).mjd for filename in chunk_files]

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

    # Format the output filename using observation date and object name
    obj_time_string = Time(header['DATE-OBS'] + 'T' + header['UT']).strftime('%Y%m%dT%H%M%S')
    header_card_object = header['OBJECT']
    if len(header_card_object.strip()) > 0:
        obj_name_string = '_'.join(header['OBJECT'].split())
    else:
        obj_name_string = 'ObjectedCardEmpty'
    image_name_stem = f'reduction_{obj_name_string}_{obj_time_string}_multi'

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

        # Append the HDU to the list
        hdulist.append(hdu)

    # Write the HDU list to the output file, overwriting if necessary
    output_filename = os.path.abspath(os.path.join(outfolder, image_name_stem + '.fits'))
    fits.HDUList(hdulist).writeto(output_filename, overwrite=True)

    return output_filename


def load_fits(fits_filename):
    """
    Purpose: Open the FITS file and extract header cards needed to construct observation MJD time
    Reads only the zeroth element returned by astropy.io.fits.open()

    Args:
        fits_filename (str): full-path file name for FITS file to be read
    Returns:
        obs_data (fits.Header): astropy FITS Header object
        obs_header (fits.Data): astropy FITS Data object
    """
    with fits.open(fits_filename) as fob:
        obs_data = fob[0].data
        obs_header = fob[0].header
    return obs_data, obs_header


def get_fits_header_mjd(fits_header):
    """
    Purpose: Open the FITS file and extract header cards needed to construct observation MJD time

    Args:
        fits_header (astropy.io.fits.Header): FITS header read from FITS obs file, expected cards 'DATE-OBS' and 'UT'
    Returns:
        obs_mjd (float): Modified Julian Date (JD - 2400000.5)
    """
    datetime_string = fits_header['DATE-OBS'] + 'T' + fits_header['UT']
    astro_time = Time(datetime_string)
    obs_mjd = astro_time.mjd
    return obs_mjd


def read_fits(file_name, read_data=False, use_memmap=False):
    """
	Purpose: Load image data and header from FITS file.
	Note: Assumes only a single HDU (Header/Data Unit) in the FITS file, e.g. len(hdu_list) ==1

	Args:
        file_name (str): FITS file name to be read
        read_data (bool): do NOT load array data unless read_data==True
        use_memmap (bool): `memmap` for astropy.io.fits.open() read of FITS file, for large files
	Returns:
        data_array (list(np.ndarray), None): numpy array, or None if read_data==False
        fits_header (astropy.io.fits.Header): dict-like object with header cards, use convert_fits_header_to_dict() if needed
	"""
    if not os.path.isfile(file_name):
        raise FileNotFoundError(file_name)

    if not use_memmap:
        # check it anyway!
        use_memmap = is_fits_memory_map_needed(file_name)
        if use_memmap:
            print(f'WARNING: file too large compared to system memory, read will be SLOW: forcing use of memmap ...')

    if not read_data:
        print(f'WARNING: only reading FITS header, since input read_data={read_data}')

    with fits.open(file_name, memmap=use_memmap) as hdul:
        fits_header = hdul[0].header
        if not read_data:
            data_array = None
        else:
            data_array = hdul[0].data

        # Manually apply scaling if BZERO, BSCALE are present
        if read_data and use_memmap:
            bzero = fits_header.get('BZERO', 0)
            bscale = fits_header.get('BSCALE', 1)
            if bscale != 1 or bzero != 0:
                data_array = (data_array * bscale) + bzero

    return data_array, fits_header


def is_fits_memory_map_needed(file_name=None, file_size_bytes=None):
    """
    Purpose: Checks available system memory and suggests whether to
             use memory-mapping (memmap) based on file size.
             Either file_name or file_size must be provided.

    Args:
        file_name (str): The name of the file to check (optional).
        file_size_bytes (int): The size of the file in bytes (optional).
    Returns:
        memory_mapping_needed (bool): True if memmap is recommended, False if not.
    Raises:
        ValueError if neither file_name nor file_size is provided.
    """
    if file_name is None and file_size_bytes is None:
        raise ValueError("Either file_name or file_size_bytes must be provided.")

    if file_size_bytes is None and file_name is not None:
        file_size_bytes = os.path.getsize(file_name)

    # Get available memory (in bytes)
    available_memory_bytes = psutil.virtual_memory().available
    available_memory_MB = available_memory_bytes / (1024 ** 2)
    file_size_MB = file_size_bytes / (1024 ** 2)
    print(f"Available memory: {available_memory_MB:.2f} MB")
    print(f"File size: {file_size_MB:.2f} MB")

    # Suggest using memmap if file size exceeds 50% of available memory
    if file_size_MB > (0.5 * available_memory_MB):
        print("Large file detected, consider using FITS memmap.")
        memory_mapping_needed = True
    else:
        print("Sufficient memory available, FITS memmap not needed.")
        memory_mapping_needed = False
    return memory_mapping_needed


def convert_fits_header_to_dict(fits_header, keep_card_comments=False):
    """
    Purpose: convert FITS Header object to dict

    Args:
        fits_header (astropy.io.fits.Header): FITS header object, collection of Card objects
        keep_card_comments (bool): If True, preserves values as nested dicts required to recreate FITS header object.
                                   If False, just a flat dict, ideal for print/write to JSON.
    Returns:
        fits_header_dict (dict): See comments about keep_card_comments input kwarg.
    """
    fits_header_dict = {}
    for card in fits_header.cards:
        keyword = card.keyword
        value = card.value
        comment = card.comment
        # append string values for "keys" that can appear more than once (e.g. COMMENT) in a FITS header
        if keyword in ['COMMENT', 'HISTORY'] and keyword in fits_header_dict.keys():
            value = fits_header_dict[keyword] + value
        if keep_card_comments:
            fits_header_dict[keyword] = {'value': value, 'comment': comment}
        else:
            fits_header_dict[keyword] = value

    return fits_header_dict


def write_fits_header_txt(fits_header, file_name):
    """
    Purpose: write astropy.io.fits Header object to a plain text file

    Args:
        fits_header (astropy.io.fits.Header): FITS Header object
        file_name (str): name of TXT file to write PLAIN TXT-converted FITS header
    Returns:
        fits_header_str (str): plain-text version of FITS header written to file_name
    """
    if os.path.isfile(file_name):
        raise FileExistsError()
    with open(file_name, 'w') as fob:
        fits_header_str = fits_header.tostring(sep='\n', endcard=False)
        fob.write(fits_header_str)
    return fits_header_str
