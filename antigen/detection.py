import logging

import numpy as np
from photutils.detection import DAOStarFinder
from photutils.background import MMMBackground
from astropy.table import Table
from astropy.stats import mad_std

from antigen import fiber
from antigen import cube

logger = logging.getLogger('antigen.detection')

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
    sigma = mad_std(data_sub)
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

def detect_brightest_source(fiber_x, fiber_y, reduced_spectra, fiber_area):
    """
    Detect the brightest source in a collapsed fiber image.

    This function collapses the reduced spectra across wavelength to create a
    synthetic "white-light" image of the field, projects the fiber fluxes into
    image space, and runs source detection to locate the brightest object.
    It returns both the object catalog and the coordinates of the detected
    source in image units.

    Args:
        fiber_x (ndarray): X positions of fibers (1D array of length Nfibers).
        fiber_y (ndarray): Y positions of fibers (1D array of length Nfibers).
        reduced_spectra (ndarray): Reduced spectra with shape (Nfibers, Nlambda).
        fiber_area (float): Effective area of each fiber (used in flux projection).

    Returns:
        sources (Table): Source catalog from detection routine.
        x_coord (float): X coordinate of the brightest source centroid.
        y_coord (float): Y coordinate of the brightest source centroid.
        X (ndarray): Grid of X coordinates corresponding to the detection image.
        Y (ndarray): Grid of Y coordinates corresponding to the detection image.

    Raises:
        RuntimeError: If no sources are detected.
    """
    # Collapse flux over wavelength to make detection image
    collapsed_fiber_flux = np.nanmedian(reduced_spectra, axis=1)
    bounds = fiber.get_fiber_bounds(fiber_x, fiber_y)

    detection_image, X, Y = cube.fibers_to_image(
        fiber_x, fiber_y, collapsed_fiber_flux, fiber_area, bounds=bounds, method="gdw", k=5, sigma=2.0
    )

    sources = detect_sources(detection_image, brightest_only=True)
    if len(sources) == 0:
        raise RuntimeError("No sources detected in the collapsed fiber image.")

    j, i = (int(sources['xcentroid']), int(sources['ycentroid']))
    x_coord = X[i, j]
    y_coord = Y[i, j]

    logger.info("Detected brightest source near %.1f, %.1f", x_coord, y_coord)

    return sources, x_coord, y_coord, X, Y