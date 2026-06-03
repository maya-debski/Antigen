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
    # Guard against zeros/NaNs in the average to avoid divide-by-zero/invalid warnings
    den = average[np.newaxis, :]
    initial_ftf = np.full_like(spectrum, np.nan)
    valid = np.isfinite(den) & (den != 0)
    np.divide(spectrum, den, out=initial_ftf, where=valid)

    # Create a wavelength grid and divide it into chunks for smoothing
    columns = np.arange(spectrum.shape[1])
    chunked_columns = np.array([np.mean(chunk) for chunk in np.array_split(columns, n_chunks)])

    # Initialize the smoothed correction array
    ftf = np.zeros_like(spectrum)

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


def get_fiber_bounds(fiber_x, fiber_y):
    """
    Computes robust bounds [x_min, x_max, y_min, y_max] for fiber positions.

    - Ignores non-finite entries (NaN/Inf).
    - If all entries are non-finite, returns a small default box [-5, 5, -5, 5].
    - Ensures non-degenerate bounds by expanding zero-width ranges slightly.

    Args:
        fiber_x (ndarray): Fiber x coordinates.
        fiber_y (ndarray): Fiber y coordinates.

    Returns:
        list[float]: [x_min, x_max, y_min, y_max]
    """
    fx = np.asarray(fiber_x, dtype=float).ravel()
    fy = np.asarray(fiber_y, dtype=float).ravel()
    m = np.isfinite(fx) & np.isfinite(fy)
    if not np.any(m):
        logger.warning("All fiber positions are non-finite; using default bounds [-5,5,-5,5].")
        return [-5.0, 5.0, -5.0, 5.0]

    fx = fx[m]
    fy = fy[m]
    try:
        x_min = float(np.nanmin(fx))
        x_max = float(np.nanmax(fx))
        y_min = float(np.nanmin(fy))
        y_max = float(np.nanmax(fy))
    except Exception:
        logger.warning("Failed to compute finite bounds; using default [-5,5,-5,5].")
        return [-5.0, 5.0, -5.0, 5.0]

    # Ensure non-degenerate ranges
    if not np.isfinite(x_min) or not np.isfinite(x_max) or x_min == x_max:
        mu = np.nanmean(fx) if np.isfinite(np.nanmean(fx)) else 0.0
        x_min, x_max = mu - 5.0, mu + 5.0
        logger.warning("Degenerate/invalid X bounds; expanded to [%.2f, %.2f] around mean.", x_min, x_max)
    if not np.isfinite(y_min) or not np.isfinite(y_max) or y_min == y_max:
        mu = np.nanmean(fy) if np.isfinite(np.nanmean(fy)) else 0.0
        y_min, y_max = mu - 5.0, mu + 5.0
        logger.warning("Degenerate/invalid Y bounds; expanded to [%.2f, %.2f] around mean.", y_min, y_max)

    return [x_min, x_max, y_min, y_max]

def load_fiber_positions(instrument, ndithers, dither_numbers, fiber_x_base, fiber_y_base):
    """
    Load the fiber positions for a given instrument and apply the appropriate dither offsets.

    This function retrieves the base fiber coordinates from the instrument configuration
    (typically the focal-plane positions of each fiber) and then applies telescope dither
    offsets to produce the full set of fiber positions for the requested dither pattern.
    The dither offsets are read from an instrument-specific dither file, which encodes
    the step sizes (dx, dy) for each dither position. If no dithering is used
    (ndithers == 1), the base positions are returned unchanged.

    This is typically used at the start of a reduction to construct the fiber footprint
    on the sky (or detector) for all exposures in a dataset.

    Args:
        instrument (str): Instrument name (used to locate the correct dither file).
        ndithers (int): Number of dithers in the observing sequence.
        dither_numbers (list of int): Indices of the dithers to apply (1-based).
        fiber_x_base (ndarray): X positions of fibers in arcseconds
        fiber_y_base (ndarray): Y positions of fibers in arcseconds

    Returns:
        fiber_x (ndarray): X positions of fibers with dither offsets applied.
        fiber_y (ndarray): Y positions of fibers with dither offsets applied.
    """
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
    # Convert dither count to 1...N and not 0...N-1.
    if np.any(dither_numbers == 0):
        dither_numbers += 1

    dithers_observed = [dither_offsets[int(dither_number - 1)] for dither_number in dither_numbers]
    fiber_x = np.hstack([fiber_x_base - dx for dx, dy in dithers_observed])
    fiber_y = np.hstack([fiber_y_base - dy for dx, dy in dithers_observed])

    logger.info(f"Loaded positions for {instrument} with {len(dithers_observed)} dithers and {dither_file} pattern.")

    return fiber_x, fiber_y
