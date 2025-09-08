import logging
import warnings

import numpy as np
from scipy.signal import savgol_filter

from antigen import detection
from antigen import extinction
from antigen import io
from antigen import psf
from antigen import spectra

warnings.filterwarnings("ignore")
logger = logging.getLogger('antigen.calibrate')


def measure_response(obs_wave, obs_flux, std_wave, std_flux, window=251):
    """Compute the instrument response function from a standard star.

    Response is defined as:

        R(λ) = F_std(λ) / F_obs(λ)

    Optionally, a simple moving-average smoothing can be applied.

    Args:
        obs_wave (np.ndarray): Wavelength grid of the observed (extinction-corrected) spectrum.
        obs_flux (np.ndarray): Observed (corrected) flux values.
        std_wave (np.ndarray): Wavelength grid of the reference standard star spectrum.
        std_flux (np.ndarray): Reference flux of the standard star
        window (int, optional): Window size for moving-average smoothing (default 251).

    Returns:
        response (np.ndarray): Instrument response values.
    """
    std_on_obs = np.interp(obs_wave, std_wave, std_flux)
    with np.errstate(divide="ignore", invalid="ignore"):
        response = np.where(obs_flux > 0, std_on_obs / obs_flux, np.nan)
    finite_values = np.isfinite(response)
    response = np.interp(obs_wave, obs_wave[finite_values], response[finite_values],
                        left=response[finite_values][0],  right=response[finite_values][-1])
    response = savgol_filter(response, window, 3)
    return response


def load_calibration_data(standard_name):
    """
    Load standard star spectrum and extinction curve for calibration.

    Args:
        standard_name (str): Name of standard star.

    Returns:
        tuple:
            calspectrum_table (Table): CALSPEC standard star spectrum.
            extinction_table (Table): Atmospheric extinction table.
    """
    calspectrum_table = io.load_calspec_spectrum(standard_name)
    extinction_table = io.read_extinction_table()
    return calspectrum_table, extinction_table


def extract_optimal_spectrum(reduced_spectra, reduced_error, dar_model,
                             sources, X, Y, measured_fwhm, fiber_x, fiber_y,
                             psf_interp, def_wave):
    """
    Extract an optimal 1D spectrum using PSF-weighted fiber fluxes.

    Args:
        reduced_spectra (np.ndarray): Flux array of shape (Nfibers, Nlambda).
        reduced_error (np.ndarray): Error array of shape (Nfibers, Nlambda).
        dar_model (DARModel): DAR model for source position correction.
        sources (Table): Detected source catalog.
        X (np.ndarray): Grid of X coordinates.
        Y (np.ndarray): Grid of Y coordinates.
        measured_fwhm (np.ndarray): FWHM array from PSF fit.
        fiber_x (np.ndarray): Fiber X positions.
        fiber_y (np.ndarray): Fiber Y positions.
        psf_interp (callable): PSF interpolator.
        def_wave (np.ndarray): Wavelength grid.

    Returns:
        tuple:
            spectrum (np.ndarray): Extracted flux spectrum.
            spectrum_error (np.ndarray): Extracted error spectrum.
    """
    source_x, source_y = dar_model(
        def_wave,
        sources["xcentroid"] + X[0, 0],
        sources["ycentroid"] + Y[0, 0]
    )
    source_fwhm = np.median(measured_fwhm)

    weights = psf.build_psf_weights(source_x, source_y, source_fwhm,
                                    fiber_x, fiber_y, psf_interp)
    spectrum, spectrum_error = spectra.get_optimal_spectrum(
        reduced_spectra, reduced_error, weights
    )

    # Mask bad edges
    spectrum[:4] = np.nan
    spectrum_error[:4] = np.nan

    return spectrum, spectrum_error


def apply_extinction(def_wave, spectrum, spectrum_error, extinction_table, airmass):
    """
    Apply atmospheric extinction correction to extracted spectrum.

    Args:
        def_wave (np.ndarray): Wavelength grid.
        spectrum (np.ndarray): Flux spectrum.
        spectrum_error (np.ndarray): Error spectrum.
        extinction_table (Table): Table with 'wavelength' and 'mags_per_airmass'.
        airmass (float): Observed airmass.

    Returns:
        tuple:
            spectrum (np.ndarray): Extinction-corrected spectrum.
            spectrum_error (np.ndarray): Extinction-corrected errors.
    """
    return extinction.apply_extinction_correction(
        def_wave, spectrum, spectrum_error,
        extinction_table["wavelength"],
        extinction_table["mags_per_airmass"],
        airmass
    )


def build_psf_and_dar(fiber_x, fiber_y, def_wave, reduced_spectra, reduced_error,
                      extraction_radius, fiber_radius, psf_seeing_grid=None):
    """Create the PSF interpolator, detect a bright source, and fit the DAR model.

    This function:
      1) Builds a radially symmetric PSF interpolator on a radius grid sized
         relative to the extraction radius.
      2) Detects the brightest source in a collapsed frame for PSF anchoring.
      3) Fits the PSF as a function of wavelength to derive a DAR model and
         measures the seeing (FWHM) as a function of wavelength.

    Args:
        fiber_x (np.ndarray): Fiber x-positions after dithering adjustments.
        fiber_y (np.ndarray): Fiber y-positions after dithering adjustments.
        def_wave (np.ndarray): Rectified wavelength grid.
        reduced_spectra (np.ndarray): Stacked reduced spectra.
        reduced_error (np.ndarray): Stacked error estimates.
        extraction_radius (float): Extraction radius taken from args.
        fiber_radius (float): Fiber radius determined from config_dict['fiber_radius']
        psf_seeing_grid (np.ndarray, optional): 1D array of seeing values (arcsec)
            to tabulate the PSF interpolator. If not provided, uses
            ``np.linspace(0.5, 5.0, 45)``.

    Returns:
        dar_model (callable): Function mapping (wavelength, fiber_x, fiber_y)
          to refraction-corrected positions; suitable for downstream use.
        measured_fwhm (np.ndarray): Estimated FWHM (arcsec) vs wavelength.

    Raises:
        RuntimeError: If source detection or PSF/DAR fitting fails.
    """
    fiber_area = np.pi * (fiber_radius ** 2)

    # PSF interpolator grid
    if psf_seeing_grid is None:
        psf_seeing_grid = np.linspace(0.5, 5.0, 45)
    interpolation_radius = 1.5 * extraction_radius
    r = np.linspace(0.0, interpolation_radius, 101)

    # 1) PSF interpolator
    psf_interp = psf.build_psf_interpolator(
        r,
        seeing=psf_seeing_grid,
        fiber_radius=fiber_radius
    )

    # 2) Brightest source for PSF anchoring
    sources, x_coord, y_coord, X, Y = detection.detect_brightest_source(
        fiber_x, fiber_y, reduced_spectra, fiber_area
    )

    # 3) Fit PSF across wavelength to build DAR and seeing model
    dar_model, measured_fwhm = psf.fit_psf_and_build_dar_model(
        reduced_spectra, reduced_error, fiber_x, fiber_y,
        psf_interp, x_coord, y_coord, def_wave,
        extraction_radius=extraction_radius, nchunks=20
    )

    return {
        "dar_model": dar_model,
        "measured_fwhm": measured_fwhm,
        "psf_interp": psf_interp,
        "sources": sources,
        "X": X,
        "Y": Y,
    }

