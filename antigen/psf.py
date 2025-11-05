import logging

import numpy as np
from astropy.stats import biweight_location as biweight
from astropy.modeling.functional_models import Moffat2D
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import least_squares

from antigen import dar
from antigen.plot import plot_psf_fit_diagnostics, plot_dar_models_comparison

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


def fit_psf(data, error, fiber_x, fiber_y, interp, initial_x, initial_y, wavelength,
            extraction_radius=20., Nchunks=20, save_diagnostics=True, output_dir="."):
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
        wavelength (ndarray): wavelength array.
        extraction_radius (float, optional): Extraction radius (arcsec).
        Nchunks (int, optional): Number of chunks. Default is 20.
        save_diagnostics (bool, optional): Save diagnostic plots. Default is True.
        output_dir (str, optional): Directory to save diagnostics. Default is ".".

    Returns:
        source_x (ndarray): location in x of source as a function of wave
        source_y (ndarray): location in y of source as a function of wave
        source_fwhm (ndarray): PSF fwhm of source as a function of wave
        poly_x (ndarray): polynomial coefficients for x position
        poly_y (ndarray): polynomial coefficients for y position
        poly_fwhm (ndarray): polynomial coefficients for fwhm
        fit_diagnostics (dict): Dictionary containing fit data for diagnostics
    """
    import matplotlib.pyplot as plt
    from pathlib import Path

    datachunks = np.array_split(data, Nchunks, axis=1)
    errorchunks = np.array_split(error, Nchunks, axis=1)

    source_x = np.zeros(Nchunks)
    source_y = np.zeros_like(source_x)
    source_fwhm = np.zeros_like(source_x)
    indices = np.arange(data.shape[1])
    little_inds = [int(np.mean(xi)) for xi in np.array_split(indices, Nchunks)]
    little_waves = wavelength[little_inds]

    # Store all fit data for diagnostics
    fit_diagnostics = {
        'chunk_data': [],
        'wavelengths': little_waves
    }

    for idx, (dchunk, echunk) in enumerate(zip(datachunks, errorchunks)):
        logger.debug(f"Processing PSF fit chunk {idx+1}/{Nchunks}")
        
        fiber_flux = biweight(dchunk, ignore_nan=True, axis=1)
        fiber_error = np.sqrt(np.nansum(echunk ** 2, axis=1)) / np.sqrt(dchunk.shape[1])

        # Use brightest fiber as initial guess
        fiber_dist = np.sqrt((fiber_x - initial_x) ** 2 + (fiber_y - initial_y) ** 2)
        sel = np.isfinite(fiber_flux) & (fiber_dist < extraction_radius)
        
        logger.debug(f"Chunk {idx}: Selected {np.sum(sel)} fibers within extraction radius")
        
        if np.sum(sel) < 4:  # Need minimum fibers for fit
            logger.warning(f"Chunk {idx}: Insufficient valid fibers ({np.sum(sel)}) for PSF fit, skipping")
            source_x[idx] = initial_x if idx == 0 else source_x[idx-1] 
            source_y[idx] = initial_y if idx == 0 else source_y[idx-1]
            source_fwhm[idx] = 1.5 if idx == 0 else source_fwhm[idx-1]
            continue

        # fiber-weighted initial x,y positions
        weights = fiber_flux[sel]
        if np.sum(weights) == 0:
            logger.warning(f"Chunk {idx}: Zero total flux, using initial positions")
            x_init = initial_x
            y_init = initial_y
        else:
            x_init = np.nansum(fiber_x[sel] * weights) / np.nansum(weights)
            y_init = np.nansum(fiber_y[sel] * weights) / np.nansum(weights)

        # Initial PSF params in the frame of x_init, y_init
        total_flux = np.nansum(fiber_flux[sel])
        initial = np.array([0.0, 0.0, 1.5, total_flux])
        
        logger.debug(f"Chunk {idx}: Initial params: x_off={initial[0]:.3f}, y_off={initial[1]:.3f}, "
                    f"fwhm={initial[2]:.3f}, flux={initial[3]:.1f}")
        
        # Validate initial parameters
        if not np.all(np.isfinite(initial)):
            logger.error(f"Chunk {idx}: Non-finite initial parameters: {initial}")
            raise ValueError(f"PSF fit chunk {idx} has non-finite initial parameters")
        
        try:
            # Fit the PSF restricted to the fibers within extraction radius
            result = least_squares(
                _psf_residuals, x0=initial,
                args=(fiber_x[sel] - initial_x, fiber_y[sel] - initial_y,
                      fiber_flux[sel], fiber_error[sel], interp)
            )
            
            if not result.success:
                logger.warning(f"Chunk {idx}: PSF fit did not converge: {result.message}")
            
            params = result.x
            
            # Validate fitted parameters
            if not np.all(np.isfinite(params)):
                logger.error(f"Chunk {idx}: Non-finite fitted parameters: {params}")
                raise ValueError(f"PSF fit chunk {idx} produced non-finite parameters")
            
            source_x[idx] = params[0] + initial_x
            source_y[idx] = params[1] + initial_y  
            source_fwhm[idx] = params[2]
            total_flux = params[3]
            
            # Calculate model values at fiber positions for diagnostics
            model_flux = model_psf(params, fiber_x[sel] - initial_x, fiber_y[sel] - initial_y, interp)
            
            # Store diagnostic data for this chunk
            chunk_diagnostics = {
                'chunk_idx': idx,
                'wavelength': little_waves[idx],
                'x_rel': fiber_x[sel] - initial_x,  # Relative to initial guess
                'y_rel': fiber_y[sel] - initial_y,
                'x_abs': fiber_x[sel],  # Absolute positions
                'y_abs': fiber_y[sel],
                'fiber_flux': fiber_flux[sel] / total_flux,
                'fiber_error': fiber_error[sel] / total_flux,
                'model_flux': model_flux / total_flux,
                'fitted_x': source_x[idx],
                'fitted_y': source_y[idx],
                'fitted_fwhm': source_fwhm[idx],
                'r': np.sqrt((fiber_x[sel] - source_x[idx])**2 + (fiber_y[sel] - source_y[idx])**2),
                'flux_model_ratio': fiber_flux[sel] / model_flux
            }
            fit_diagnostics['chunk_data'].append(chunk_diagnostics)
            
            logger.debug(f"Chunk {idx}: Fitted params: x={source_x[idx]:.3f}, y={source_y[idx]:.3f}, "
                        f"fwhm={source_fwhm[idx]:.3f}")
            
        except Exception as e:
            logger.error(f"Chunk {idx}: PSF fit failed with error: {str(e)}")
            logger.error(f"Debug info - sel_count: {np.sum(sel)}, x_init: {x_init:.3f}, y_init: {y_init:.3f}")
            logger.error(f"fiber_flux[sel] range: [{np.nanmin(fiber_flux[sel]):.3f}, {np.nanmax(fiber_flux[sel]):.3f}]")
            logger.error(f"fiber_error[sel] range: [{np.nanmin(fiber_error[sel]):.3f}, {np.nanmax(fiber_error[sel]):.3f}]")
            raise
    
    # Fit polynomials to the results
    poly_x = np.polyfit(little_waves, source_x[int(Nchunks/2)] - source_x, 2)
    poly_y = np.polyfit(little_waves, source_y[int(Nchunks/2)] - source_y, 2)
    poly_fwhm = np.polyfit(little_waves, source_fwhm, 2)

    # Create diagnostic plots
    if save_diagnostics and fit_diagnostics['chunk_data']:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
        # Save diagnostics under a subfolder within the provided output_dir
        plot_psf_fit_diagnostics(fit_diagnostics, output_dir=output_path / "psf_diagnostics")


    return source_x, source_y, source_fwhm, poly_x, poly_y, poly_fwhm


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


def fit_psf_and_build_dar_model(reduced_spectra, reduced_error, fiber_x, fiber_y,
                                psf_interp, x_coord, y_coord, def_wave,
                                extraction_radius=2.0, nchunks=20, output_dir=".", header=None):
    """
    Fit the PSF to the reduced spectra and construct a DAR model.

    Args:
        reduced_spectra (ndarray): Flux array of shape (Nfibers, Nlambda).
        reduced_error (ndarray): Error array of shape (Nfibers, Nlambda).
        fiber_x (ndarray): Fiber X positions.
        fiber_y (ndarray): Fiber Y positions.
        psf_interp (callable): Interpolator for PSF profile.
        x_coord (float): X coordinate of detected source centroid.
        y_coord (float): Y coordinate of detected source centroid.
        def_wave (ndarray): Wavelength grid.
        extraction_radius (float, optional): Extraction radius in arcsec. Defaults to 2.0.
        nchunks (int, optional): Number of wavelength chunks for fitting. Defaults to 20.
        output_dir (str, optional): Directory to save diagnostics. Defaults to ".".
    Returns:
        tuple:
            dar_model (DARModel): Differential atmospheric refraction model.
            measured_fwhm (ndarray): FWHM of PSF per chunk.
    """
    params = fit_psf(
        reduced_spectra, reduced_error, fiber_x, fiber_y, psf_interp,
        x_coord, y_coord, def_wave,
        extraction_radius=extraction_radius, Nchunks=nchunks,
        save_diagnostics=True, output_dir=output_dir
    )
    _, _, measured_fwhm, coeff_x, coeff_y, _ = params

    # Build DARModel. If polynomial coeffs are available from PSF fit, use them
    if coeff_x is not None and coeff_y is not None:
        dar_model = dar.DARModel(wave=def_wave, header=header)
        # If header is provided, also build a header-based DAR model and compare
        if header is not None:
            try:
                # Evaluate both models across the wavelength grid
                dx_psf_instr = np.polyval(coeff_x, def_wave)
                dy_psf_instr = np.polyval(coeff_y, def_wave)
                # Convert PSF-fitted instrument-axis shifts to East/North for comparison
                # Convention: +X (instrument) = West => East = -X; +Y (instrument) = North
                dx_psf = -dx_psf_instr  # East
                dy_psf = dy_psf_instr   # North
                dx_hdr = np.polyval(dar_model.coeffs_x, def_wave)
                dy_hdr = np.polyval(dar_model.coeffs_y, def_wave)
                # Plot comparison (align both models to a common reference wavelength for fair comparison)
                plot_dar_models_comparison(
                    def_wave, dx_psf, dy_psf, dx_hdr, dy_hdr,
                    output_dir=output_dir, align_to_ref=True, ref_wave=5000.0
                )
                logger.info("DAR model comparison plot saved under dar_diagnostics (aligned to common ref wavelength).")
            except Exception as e:
                logger.warning(f"Failed to build/compare header-based DAR model: {e}")
    else:
        dar_model = dar.DARModel(wave=def_wave, header=header)
    return dar_model, measured_fwhm