import logging

import numpy as np
from astropy.stats import biweight_location as biweight
from astropy.modeling.functional_models import Moffat2D
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import least_squares

logger = logging.getLogger('antigen.psf')

def build_psf_interpolator(r, seeing, scale=0.2, alpha=3.5, fiber_radius=2.1):
    """
    Create an interpolator for Moffat PSF flux fractions within circular fiber apertures.

    Args:
        r (ndarray): 1D array of fiber radii (arcsec).
        seeing (ndarray): 1D array of FWHM values (arcsec).
        scale (float, optional): Pixel grid scale in arcsec. Default is 0.05.
        alpha (float, optional): Moffat shape parameter. Default is 3.5.
        fiber_radius (float, optional): Radius of the IFU fibers (arcsecs)

    Returns:
        LinearNDInterpolator: Interpolator to compute enclosed PSF flux.
    """
    boxsize = np.max(r) * 2.0
    x = np.arange(-boxsize / 2., boxsize / 2. + scale, scale)
    y = np.arange(-boxsize / 2., boxsize / 2. + scale, scale)
    xgrid, ygrid = np.meshgrid(x, y)

    V = np.zeros((len(seeing), len(r)))

    for j, fwhm in enumerate(seeing):
        moffat = Moffat2D()
        moffat.alpha.value = alpha
        moffat.gamma.value = 0.5 * fwhm / np.sqrt(2 ** (1. / alpha) - 1.)
        Z = moffat(xgrid, ygrid)
        Z /= Z.sum()

        for i, ri in enumerate(r):
            d = np.sqrt((xgrid - ri) ** 2 + (ygrid - 0.) ** 2)
            sel = d <= fiber_radius
            adj = np.pi * fiber_radius ** 2 / (sel.sum() * scale ** 2)
            V[j, i] = np.sum(Z[sel]) * adj

    R, S = np.meshgrid(r, seeing)
    return LinearNDInterpolator(list(zip(R.ravel(), S.ravel())), V.ravel())


def model_psf(params, fiber_x, fiber_y, interp):
    """
    Compute PSF model flux in each fiber.

    Args:
        params (ndarray): PSF parameters [x0, y0, FWHM, total_flux].
        fiber_x (ndarray): Fiber X coordinates.
        fiber_y (ndarray): Fiber Y coordinates.
        interp (LinearNDInterpolator): PSF flux interpolator.

    Returns:
        array: Model flux per fiber.
    """
    x0, y0, fwhm, flux_total = params
    r = np.sqrt((fiber_x - x0) ** 2 + (fiber_y - y0) ** 2)
    fractions = interp(r, fwhm)
    return flux_total * fractions


def _psf_residuals(params, fiber_x, fiber_y, fiber_flux, fiber_error, interp):
    """
    Calculate residuals between observed and modeled fiber flux.

    Args:
        params (ndarray): PSF parameters.
        fiber_x (ndarray): Fiber X coordinates.
        fiber_y (ndarray): Fiber Y coordinates.
        fiber_flux (ndarray): Observed fiber fluxes.
        fiber_error (ndarray): Fiber flux uncertainties.
        interp (LinearNDInterpolator): PSF flux interpolator.

    Returns:
        residuals (array): chi residuals
    """
    model = model_psf(params, fiber_x, fiber_y, interp)
    residuals = (fiber_flux - model) / fiber_error
    return residuals


def fit_psf(data, error, fiber_x, fiber_y, interp, initial_x, initial_y,
            extraction_radius=20., Nchunks=20):
    """
    Fit the PSF to fiber data in chunks.

    Args:
        data (ndarray): Fiber flux data.
        error (ndarray): Fiber error data.
        fiber_x (ndarray): Fiber X coordinates.
        fiber_y (ndarray): Fiber Y coordinates.
        interp (LinearNDInterpolator): PSF flux interpolator.
        initial_x (float): initial guess for star x coordinates.
        initial_y (float): initial guess for star y coordinates.
        extraction_radius (float, optional): Extraction radius (arcsec).
        Nchunks (int, optional): Number of chunks. Default is 20.

    Returns:
        source_x (ndarray): location in x of source as a function of wave
        source_y (ndarray): location in y of source as a function of wave
        source_fwhm (ndarray): PSF fwhm of source as a function of wave
    """
    datachunks = np.array_split(data, Nchunks, axis=1)
    errorchunks = np.array_split(error, Nchunks, axis=1)

    source_x = np.zeros(Nchunks)
    source_y = np.zeros_like(source_x)
    source_fwhm = np.zeros_like(source_x)
    indices = np.arange(data.shape[1])
    little_inds = [np.mean(xi) for xi in np.array_split(indices, Nchunks)]

    for idx, (dchunk, echunk) in enumerate(zip(datachunks, errorchunks)):
        fiber_flux = biweight(dchunk, ignore_nan=True, axis=1)
        fiber_error = np.sqrt(np.nansum(echunk ** 2, axis=1)) / np.sqrt(dchunk.shape[1])

        # Use brightest fiber as initial guess
        # Add an option for a guess??
        fiber_dist = np.sqrt((fiber_x - initial_x) ** 2 + (fiber_y - initial_y) ** 2)
        sel = np.isfinite(fiber_flux) & (fiber_dist < extraction_radius)

        # fiber-weighted initial x,y positions
        x_init = np.nansum(fiber_x[sel] * fiber_flux[sel]) / np.nansum(fiber_flux[sel])
        y_init = np.nansum(fiber_y[sel] * fiber_flux[sel]) / np.nansum(fiber_flux[sel])

        # Initial PSF params in the frame of x_init, y_init
        initial = np.array([0.0, 0.0, 1.5, np.nansum(fiber_flux[sel])])

        # Fit the PSF restricted to the fibers within extraction radius / 2.
        result = least_squares(
            _psf_residuals, x0=initial,
            args=(fiber_x[sel] - initial_x, fiber_y[sel] - initial_y,
                  fiber_flux[sel], fiber_error[sel], interp)
        )

        params = result.x
        params[0] += initial_x
        params[1] += initial_y
        source_x[idx] = params[0]
        source_y[idx] = params[1]
        source_fwhm[idx] = params[2]
    source_x = np.polyval(np.polyfit(little_inds, source_x, 2), indices)
    source_y = np.polyval(np.polyfit(little_inds, source_y, 2), indices)
    source_fwhm = np.polyval(np.polyfit(little_inds, source_fwhm, 2), indices)

    return source_x, source_y, source_fwhm


def build_psf_weights(source_x, source_y, source_fwhm, fiber_x, fiber_y,
                      interp):
    '''
    Build the PSF weights (covering fractions) for each fiber as function
    of wavelength

    Args:
        source_x (ndarray): X centroid positions as a function of wavelength
        source_y (ndarray): Y centroid positions as a function of wavelength
        source_fwhm (ndarray): FWHM value as a function of wavelength
        fiber_x (ndarray): Fiber X coordinates.
        fiber_y (ndarray): Fiber Y coordinates.

    Returns:
        weights (ndarray): PSF weights for each fiber as a function of wavelength

    '''
    r = np.sqrt((fiber_x[:, np.newaxis] - source_x[np.newaxis, :]) ** 2 +
                (fiber_y[:, np.newaxis] - source_y[np.newaxis, :]) ** 2)
    fractions = interp(r, source_fwhm)
    weights = fractions / np.nansum(fractions, axis=0)[np.newaxis, :]
    weights[np.isnan(weights)] = 0.0
    return weights