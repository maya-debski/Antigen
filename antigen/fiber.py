import numpy as np
from astropy.stats import biweight_location as biweight
import logging
from scipy.interpolate import interp1d

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
