import logging

from astropy.stats import biweight_location as biweight
from astropy.table import Table
import numpy as np
from scipy.interpolate import interp1d

from antigen import config

logger = logging.getLogger('antigen.fiber')

def get_fiber_to_fiber(spectrum, n_chunks=100):
    """
    Computes the fiber-to-fiber correction by normalizing each fiber's spectrum
    to the average spectrum across all fibers, then smooths the correction factor
    using interpolation.

    Parameters:
        spectrum (2D array): Array of spectra from multiple fibers (fibers x wavelength).
        n_chunks (int, optional): Number of chunks to split the wavelength range into
                                  for smoothing. Default is 100.

    Returns:
        initial_ftf (2D array): Initial fiber-to-fiber correction factors.
        ftf (2D array): Smoothed fiber-to-fiber correction factors.
    """
    # Compute the average spectrum across all fibers using a robust biweight statistic
    average = biweight(spectrum, axis=0, ignore_nan=True)

    # Calculate the initial fiber-to-fiber correction by dividing each fiber by the average spectrum
    initial_ftf = spectrum / average[np.newaxis, :]

    # Create a wavelength grid and divide it into chunks for smoothing
    columns = np.arange(spectrum.shape[1])
    chunked_columns = np.array([np.mean(chunk) for chunk in np.array_split(columns, n_chunks)])

    # Initialize the smoothed correction array
    ftf = spectrum * 0.

    # Loop through each fiber to compute the smoothed correction factor
    for i in np.arange(len(spectrum)):
        # Compute the biweight statistic for each chunk of the initial correction factor
        chunked_ftf = np.array([biweight(chunk, ignore_nan=True) for chunk in np.array_split(initial_ftf[i], n_chunks)])

        # Select valid (finite) values for interpolation
        sel = np.isfinite(chunked_ftf)
        if sel.sum() == 0.:
            continue
        # Interpolate the correction factor using quadratic interpolation
        interp_func = interp1d(chunked_columns[sel], chunked_ftf[sel], kind='quadratic', bounds_error=False,
                     fill_value='extrapolate')

        # Apply the interpolation to the full wavelength range
        ftf[i] = interp_func(columns)

    # Return both the initial and smoothed fiber-to-fiber correction factors
    return initial_ftf, ftf


def load_fiber_positions(instrument, ndithers, config_dict):
    """
    Load fiber base positions and apply dither offsets.

    Args:
        instrument (str): Instrument name.
        ndithers (int): Number of dithers.
        config_dict (dict): Config dict containing 'ifu_x' and 'ifu_y'.

    Returns:
        fiber_x (array): X-positions of fibers including dither pattern
        fiber_y (array): Y-positions of fibers including dither pattern
    """
    # --- Input validation ---
    if "ifu_x" not in config_dict or "ifu_y" not in config_dict:
        raise ValueError("config_dict must contain 'ifu_x' and 'ifu_y' arrays.")

    fiber_x_base, fiber_y_base = config_dict["ifu_x"], config_dict["ifu_y"]

    base_path = config.get_base_config_path()

    dither_file = base_path / instrument / f"{instrument}_dither_{ndithers}pt.lis"

    # --- Handle single-dither (no offset) case ---
    if ndithers == 1:
        logger.info(f"No dither pattern needed for {instrument}, ndithers=1")
        return fiber_x_base.copy(), fiber_y_base.copy()

    # --- Load dither offsets ---
    if not dither_file.exists():
        raise FileNotFoundError(f"Dither file not found: {dither_file}")

    try:
        dither_table = Table.read(dither_file, format="ascii")
        dither_offsets = np.array([dither_table[col].data for col in dither_table.colnames]).T
    except Exception as e:
        raise RuntimeError(f"Failed to read dither file {dither_file}: {e}")

    # Expect 2 columns: dx, dy
    if dither_offsets.shape[1] != 2:
        raise ValueError(f"Dither file {dither_file} must have 2 columns (dx, dy).")

    if dither_offsets.shape[0] != ndithers:
        logger.warning(
            f"Dither file has {dither_offsets.shape[0]} dithers but ndithers={ndithers}. "
            "Proceeding with file contents."
        )

    # --- Apply dithers ---
    fiber_x = np.hstack([fiber_x_base - dx for dx, dy in dither_offsets])
    fiber_y = np.hstack([fiber_y_base - dy for dx, dy in dither_offsets])

    logger.info(f"Loaded fiber positions for {instrument} with {len(dither_offsets)} dithers.")

    return fiber_x, fiber_y
