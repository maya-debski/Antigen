import os
import psutil

import numpy as np
from astropy.io import fits
from astropy.time import Time

from antigen import config


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


def get_fits_file_time(fits_file_name, instrument='VIRUS2'):
    """
    Purpose: Read FITS header, get header cards 'DATE-OBS' and 'UT', construct a time stamp, convert it to MJD

    Note: VIRUS2 uses both 'DATE-OBS' card to contain YYYYmmdd and 'UT' card to contain HH:MM:SS
    Note: GCMS and VIRUS-W use 'DATE-OBS' to contain YYYYmmddTHH:MM:SS

    Args:
        fits_file_name (str): file name for FITS file
        instrument (str): choices = ('VIRUS2', 'GCMS')
    Returns:
        obs_time_mjd (float): Modified Julian Date, float decimal days
    """
    VALID_INSTRUMENTS = ('VIRUS2', 'VIRUSW', 'GCMS')
    if instrument not in VALID_INSTRUMENTS:
        raise ValueError(f'ERROR: input instrument={instrument} not in VALID_INSTRUMENTS={VALID_INSTRUMENTS}')

    with fits.open(fits_file_name) as hdul:
        header = hdul[0].header

    if instrument == 'VIRUS2':
        time_stamp = header['DATE-OBS'] + 'T' + header['UT']
        obs_time = Time(time_stamp, format='isot', scale='utc')
        obs_time_mjd = obs_time.mjd
    elif instrument in ('VIRUSW', 'GCMS'):
        time_stamp = header['DATE-OBS']
        obs_time = Time(time_stamp, format='isot', scale='utc')
        obs_time_mjd = obs_time.mjd
    else:
        obs_time_mjd = None
    return obs_time_mjd
