import numpy as np
from photutils.detection import DAOStarFinder
from photutils.background import MMMBackground
from astropy.table import Table

def detect_sources(original_image, fwhm=2.0, threshold_sigma=5.0, brightest_only=False, nbright=1):
    """
    Detect potential sources in an image using DAOStarFinder.

    The background is calculated using a mode estimator of the form (3 * median) - (2 * mean)

    Args:
        original_image (ndarray): The synthetic image reconstructed from fiber spectra.
        fwhm (float): Approximate FWHM of the PSF in pixels.
        threshold_sigma (float): Detection threshold in units of the background RMS.
        brightest_only (bool): If True, return only the brightest source(s).
        nbright (int): Number of brightest sources to return if brightest_only=True.

    Returns:
        sources (Table): Astropy Table with detected sources (id, xcentroid, ycentroid, flux, etc.).
    """
    # Correct NaNs to zero
    image = original_image.copy()
    image[np.isnan(image)] = 0.0

    # Estimate background and noise
    bkg_estimator = MMMBackground()
    bkg_value = bkg_estimator(image)
    data_sub = image - bkg_value

    # Use robust std estimate for threshold
    sigma = np.std(data_sub)
    threshold = threshold_sigma * sigma

    # Run DAOStarFinder
    daofind = DAOStarFinder(fwhm=fwhm, threshold=threshold)
    sources = daofind(data_sub)

    if sources is None:
        return Table()  # return empty if nothing detected

    # Optionally select only the brightest
    if brightest_only:
        sources.sort('flux')
        sources = sources[::-1]  # brightest first
        sources = sources[:nbright]

    return sources
