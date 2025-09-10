import logging
import os
from pathlib import Path
import psutil

import numpy as np
from astropy.io import fits
from astropy.table import Table
from astropy.time import Time
from getCalspec.getCalspec import Calspec, is_calspec
from getCalspec.rebuild import rebuild_tables, rebuild_cache

from antigen import config

logger = logging.getLogger('antigen.io')


def get_fits_extension_map():
    """
    Return a consistent mapping of extension names for reduction FITS files.

    Returns:
        dict: Mapping of human-readable keys to FITS EXTNAME values.
    """
    return {
        "skysubrect_adv": "SKYSUB_PCA",
        "skysubrect": "SKYSUB",
        "specrect": "NOSKYSUB",
        "errorrect": "ERROR",
    }

def write_fits(skysubrect_adv, skysubrect, specrect, errorrect, header, config_dict, outfolder):
    """
    Purpose: Writes the sky-subtracted, rectified spectra and error data to a FITS file,
    preserving the header information and adding necessary meta-information.

    Args:
        skysubrect_adv (np.ndarray): 2D numpy array, The advanced sky-subtracted spectrum.
        skysubrect (np.ndarray): 2D numpy array, The basic sky-subtracted spectrum.
        specrect (np.ndarray): 2D numpy array, The rectified spectrum.
        errorrect (np.ndarray): 2D numpy array, The error associated with the rectified spectrum.
        header (fits.Header): The header information to be preserved in the output FITS file.
        config_dict (dict): Dictionary containing the configuration information.
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
    instrument = config_dict['instrument']
    instrument_element = config_dict['instrument_element']
    image_name_stem = f'reduction_{obj_name_string}_{obj_time_string}_{instrument}_{instrument_element}_multi'

    extmap = get_fits_extension_map()

    # Loop through the data arrays and create HDUs for each
    for image, ftp, extname in zip([skysubrect_adv, skysubrect, specrect, errorrect],
                                   [fits.PrimaryHDU, fits.ImageHDU, fits.ImageHDU, fits.ImageHDU],
                                   extmap.values()):

        # Create an HDU object from each image, setting it to 'float32' type
        hdu = ftp(np.array(image, dtype='float32'))
        hdu.header["EXTNAME"] = extname

        # Remove any conflicting CD matrix elements first
        for key in ['CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'CDELT1', 'CDELT2']:
            if key in hdu.header:
                del hdu.header[key]

        # Define your wavelength solution
        def_wave = np.linspace(config_dict['start_wavelength'], config_dict['end_wavelength'],
                               config_dict['detector_dimensions']['X'])
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

def get_airmass(header):
    """Extract airmass value from FITS header.

    Args:
        header (dict): FITS header dictionary containing 'AIRMASS' key

    Returns:
        float: Airmass value from the header

    Raises:
        KeyError: If AIRMASS key is not present in header
        ValueError: If AIRMASS value cannot be converted to float
    """
    return float(header['AIRMASS'])

def load_reduced_data(base_folder, filenames, extname='SKYSUB_PCA'):
    """
    Load fiber bundle data and error arrays from target reduced products
    made by Antigen at an earlier step

    Args:
        base_folder (str): base folder path
        filenames (list): list of file paths to reduced frames for target
        extname (str, optional): Extension name to use. Defaults to 'SKYSUB_PCA'.

    Returns:
        data (array): sky subtracted frames
        error (array): error frames
        header (fits header): header of first file.
    """
    data_list, err_list = [], []
    for filename in filenames:
        with fits.open(Path(base_folder) / filename, memmap=False) as hdul:
            try:
                data_list.append(hdul[extname].data)
                err_list.append(hdul["ERROR"].data)
            except KeyError as e:
                raise KeyError(f"Missing extension in {filename}: {e}")

    data = np.vstack(data_list)
    error = np.vstack(err_list)
    header = fits.getheader(Path(base_folder) /  filenames[0], 0)

    return data, error, header


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

def validate_fits_header(header, required_cards=None):
    """
    Check that required FITS header keywords exist.

    Args:
        header (fits.Header): FITS header to check.
        required_cards (list, optional): Keywords to validate. Default: ['OBJECT', 'DATE-OBS', 'UT'].

    Returns:
        is_valid (bool): True if all required keywords are present in header.
        missing (list): list of keywords not present in header.
    """
    if required_cards is None:
        required_cards = ['OBJECT', 'DATE-OBS', 'UT']

    missing = [key for key in required_cards if key not in header]
    is_valid = len(missing) == 0
    return is_valid, missing


def load_fits_header(fits_filename, validate=True, strict=False):
    """
    Load FITS header and optionally validate required keywords.

    Args:
        fits_filename (str): Full path to FITS file.
        validate (bool): Whether to validate required keywords.
        strict (bool): If True, raise error if validation fails.

    Returns:
        header (fits.Header): FITS Header object
        header_is_valid (bool or None)): True if header is valid, False otherwise. None if no validation.
    """
    if not os.path.isfile(fits_filename):
        raise FileNotFoundError(f"File not found: {fits_filename}")

    try:
        with fits.open(fits_filename) as fob:
            if len(fob) == 0:  # Empty FITS file
                raise OSError(f"Empty FITS file: {fits_filename}")
            header = fob[0].header
    except OSError as e:
        if strict:
            raise OSError(f"Failed to read FITS file '{fits_filename}': {e}")
        else:
            logger.warning(f"Failed to read FITS file '{fits_filename}': {e}")
            return None, None

    header_is_valid = None
    if validate:
        header_is_valid, missing = validate_fits_header(header)
        if not header_is_valid:
            msg = f"{fits_filename} missing required cards: {', '.join(missing)}"
            if strict:
                raise ValueError(msg)
            else:
                print(msg)

    return header, header_is_valid


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


def load_calspec_spectrum(name, spec_type="stis", date="latest", check_cache=False):
    """Load a CALSPEC spectrum as an Astropy Table.

    This function verifies that the requested source is part of the CALSPEC
    library, optionally updates the local CALSPEC cache, and loads the
    requested spectrum into an Astropy Table.

    Args:
        name (str): Name of the CALSPEC source (case-insensitive).
        spec_type (str, optional): Type of spectrum to load (e.g., "stis").
            Defaults to "stis".
        date (str, optional): Date of the version to load. Defaults to "latest".
        check_cache (bool, optional): If True, update the CALSPEC cache and tables
            before loading. Defaults to False.

    Returns:
        astropy.table.Table: Spectrum table containing wavelength and flux columns.

    Raises:
        ValueError: If the requested source is not in the CALSPEC library.
        RuntimeError: If the spectrum cannot be retrieved.
    """
    # Optional cache update
    if check_cache:
        logger.info("Updating CALSPEC cache and tables...")
        rebuild_tables()
        rebuild_cache()

    # Verify the name
    if not is_calspec(name):
        raise ValueError(f"Source '{name}' not found in CALSPEC library.")

    # Get the spectrum
    try:
        calspec_obj = Calspec(name)
        spectrum = calspec_obj.get_spectrum_table()
    except Exception as e:
        raise RuntimeError(f"Failed to retrieve CALSPEC spectrum for {name}") from e

    # Handle both upper and lower case column names
    wave_col = next((col for col in ['wavelength', 'WAVELENGTH']
                     if col in spectrum.colnames), None)
    flux_col = next((col for col in ['flux', 'FLUX']
                     if col in spectrum.colnames), None)

    if wave_col is None or flux_col is None:
        raise ValueError("Cannot find wavelength/flux columns in calibration spectrum table")

    # Get the column values, handling potential dictionary case
    wave_data = spectrum[wave_col].value
    flux_data = spectrum[flux_col].value
    return wave_data, flux_data


def read_extinction_table(file_path=None):
    """Read the McDonald Observatory extinction curve.

    This function reads the extinction curve from a fixed-width, two-column ASCII file.
    The file contains wavelength and extinction (magnitudes per airmass).

    Args:
        file_path (str or Path, optional): Path to the extinction data file.
            Defaults to ``antigen/config_files/mcdonald_extinction.dat``.

    Returns:
        astropy.table.Table: A table with two columns:
            - wavelength (float): Wavelength in Angstroms.
            - mags_per_airmass (float): Extinction in magnitudes per airmass.

    Raises:
        FileNotFoundError: If the extinction file does not exist.
        OSError: If the file cannot be parsed into two numeric columns.
    """
    base_path = config.get_base_config_path()
    default_path = base_path /  "extinction" / "mcdonald_extinction.dat"

    file_path = Path(file_path) if file_path else default_path

    if not file_path.exists():
        raise FileNotFoundError(f"Extinction file not found: {file_path}")

    try:
        table = Table.read(
            file_path,
            format="ascii.fixed_width_two_line",
            names=("wavelength", "mags_per_airmass")
        )
    except Exception as e:
        raise OSError(f"Could not parse extinction file: {e}")

    return table


def get_fits_files_in_path(input_path, pattern='*.fits'):
    """
    Find and return all FITS files in the given directory.

    Args:
        input_path (str): Path to the directory to search.
        pattern (str, optional): Glob pattern to match files. Defaults to '*.fits'.

    Returns:
        list (Path): A sorted list of Path objects matching the pattern.
    """
    return sorted(Path(input_path).glob(pattern))


def write_cube(filename, cube, wavelength, header, x_grid, y_grid, pixel_size, overwrite=True):
    """Write a spectral cube to a FITS file with relevant metadata.

    The function saves a datacube (lambda, y, x) to a FITS file, carrying
    over telescope/instrument metadata from the input header and adding
    WCS-like keywords for the cube dimensions.

    Args:
        filename (str): Path to the output FITS file.
        cube (ndarray): 3D datacube with shape (N_lambda, Ny, Nx).
        wavelength (ndarray): 1D Array of wavelength values corresponding to the first axis of the cube.
        header (fits.Header): Header from a fiber frame containing telescope/instrument information to be propagated.
        x_grid (ndarray): 2D array of x positions.
        y_grid (ndarray): 2D array of y positions.
        pixel_size (float): Size of pixel in arcsec.
        overwrite (bool, optional): If True, overwrite an existing file. Default is True.

    Returns:
        None

    Notes:
        - The output FITS file has the datacube stored in the primary HDU.
        - The header is updated to include:
            * `CRPIX1/2/3` : reference pixels
            * `CRVAL1/2/3` : reference values (x, y, wavelength)
            * `CDELT1/2/3` : increments (pixel size, wavelength step)
            * `CTYPE1/2/3` : axis type labels
        - `CDELT3` is computed from the median step in the wavelength array.

    """
    # Copy header to avoid modifying input
    hdr = header.copy()

    # Get cube shape
    n_lambda, ny, nx = cube.shape

    # Spectral step
    if len(wavelength) > 1:
        delta_lambda = np.median(np.diff(wavelength))
    else:
        delta_lambda = 1.0  # fallback if only one wavelength

    # Update WCS-like metadata
    hdr["NAXIS"]  = 3
    hdr["NAXIS1"] = nx
    hdr["NAXIS2"] = ny
    hdr["NAXIS3"] = n_lambda

    hdr["CTYPE1"] = "X"
    hdr["CTYPE2"] = "Y"
    hdr["CTYPE3"] = "WAVE"

    hdr["CDELT1"] = pixel_size
    hdr["CDELT2"] = pixel_size
    hdr["CD1_1"] = pixel_size
    hdr["CD2_2"] = pixel_size
    hdr["CDELT3"] = delta_lambda

    hdr["CRPIX1"] = 1
    hdr["CRPIX2"] = 1
    hdr["CRPIX3"] = 1

    hdr["CRVAL1"] = x_grid[0][0]
    hdr["CRVAL2"] = y_grid[0][0]
    hdr["CRVAL3"] = wavelength[0]

    # Save FITS
    hdu = fits.PrimaryHDU(data=cube.astype("float32"), header=hdr)
    hdu.writeto(filename, overwrite=overwrite)
