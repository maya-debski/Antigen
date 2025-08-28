import logging
import warnings

import numpy as np

from antigen import config
from antigen import detection
from antigen import extinction
from antigen import fiber
from antigen import io
from antigen import psf
from antigen import spectra
from antigen import utils
from antigen import wavelength


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
    from scipy.signal import savgol_filter
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


def measure_and_apply_response(def_wave, spectrum, spectrum_error,
                               cal_spectrum_table, output_folder, window=251):
    """
    Measure the instrument response and apply it to the extracted spectrum.

    Args:
        def_wave (np.ndarray): Wavelength grid.
        spectrum (np.ndarray): Flux spectrum.
        spectrum_error (np.ndarray): Error spectrum.
        cal_spectrum_table (Table): CALSPEC standard star spectrum.
        output_folder (str): Output folder for plots.
        window (int, optional): Smoothing window for response. Defaults to 251.

    Returns:
        tuple:
            spectrum (np.ndarray): Response-corrected spectrum.
            spectrum_error (np.ndarray): Response-corrected errors.
            response (np.ndarray): Derived instrument response function.
    """
    response = measure_response(
        def_wave, spectrum,
        cal_spectrum_table["wavelength"], cal_spectrum_table["flux"],
        window=window
    )

    spectrum *= response
    spectrum_error *= response
    from antigen import plot  # import here, not at module top
    plot.plot_spectrum_with_standard(
        spectrum, spectrum_error, def_wave,
        cal_spectrum_table["wavelength"], cal_spectrum_table["flux"],
        1. / response, outfolder=output_folder
    )

    return spectrum, spectrum_error, response

def prepare_response_inputs(dataset_manifest, args):
    """Prepare all inputs needed for building the instrument response.

    This function:
      1) Builds the instrument configuration dictionary.
      2) Validates the dataset manifest, CLI args, and configuration.
      3) Loads fiber positions (accounting for dithers).
      4) Retrieves the rectified wavelength grid.
      5) Loads reduced spectra, errors, and header.
      6) Computes convenience values used downstream.

    Args:
        dataset_manifest (dict): Dataset manifest returned from
            ``build_dataset_from_reduced_files()``. Must include:
            - 'unit_instrument' (str)
            - 'unit_id' (str)
            - 'ndithers' (int)
            - 'dither_number' (int)
            - 'in_folder' (str)
            - 'reduced_files' (list[str])
        args (argparse.Namespace): Command line arguments. Must include:
            - extraction_radius (float)

    Returns:
        dict: A dictionary with prepared inputs:
            - 'config_dict' (dict): Instrument configuration.
            - 'fiber_x' (np.ndarray): Fiber x-positions after dithering adjustments.
            - 'fiber_y' (np.ndarray): Fiber y-positions after dithering adjustments.
            - 'def_wave' (np.ndarray): Rectified wavelength grid.
            - 'reduced_spectra' (np.ndarray): Stacked reduced spectra (n_fiber × n_wave).
            - 'reduced_error' (np.ndarray): Stacked error estimates (n_fiber × n_wave).
            - 'header' (fits.Header): Representative FITS header from the reduced data.
            - 'fiber_area' (float): Geometric fiber area in square arcseconds.
            - 'extraction_radius' (float): Extraction radius taken from args.

    Raises:
        ValueError: If required inputs are missing or invalid.
    """
    # 1) Build config for this instrument/element
    config_dict = config.build_config_for_element(
        dataset_manifest['unit_instrument'].lower(),
        dataset_manifest['unit_id'].upper()
    )

    # 2) Validate everything early
    utils.validate_inputs(dataset_manifest, args, config_dict)

    # 3) Fiber positions (adjusted for dithers)
    fiber_x, fiber_y = fiber.load_fiber_positions(
        dataset_manifest['unit_instrument'],
        dataset_manifest['ndithers'],
        dataset_manifest['dither_number'],
        config_dict
    )

    # 4) Wavelength grid
    def_wave = wavelength.get_rectified_wavelength(config_dict)

    # 5) Reduced data
    reduced_spectra, reduced_error, header = io.load_reduced_data(
        dataset_manifest['in_folder'],
        dataset_manifest['reduced_files']
    )

    # 6) Convenience scalars
    fiber_area = np.pi * (config_dict['fiber_radius'] ** 2)
    extraction_radius = args.extraction_radius

    return {
        'config_dict': config_dict,
        'fiber_x': fiber_x,
        'fiber_y': fiber_y,
        'def_wave': def_wave,
        'reduced_spectra': reduced_spectra,
        'reduced_error': reduced_error,
        'header': header,
        'fiber_area': fiber_area,
        'extraction_radius': extraction_radius,
    }


def build_psf_and_dar(prep, psf_seeing_grid=None):
    """Create the PSF interpolator, detect a bright source, and fit the DAR model.

    This function:
      1) Builds a radially symmetric PSF interpolator on a radius grid sized
         relative to the extraction radius.
      2) Detects the brightest source in a collapsed frame for PSF anchoring.
      3) Fits the PSF as a function of wavelength to derive a DAR model and
         measures the seeing (FWHM) as a function of wavelength.

    Args:
        prep (dict): Prepared inputs returned by :func:`prepare_response_inputs`.
            Must contain:
            - 'config_dict', 'fiber_x', 'fiber_y', 'def_wave',
              'reduced_spectra', 'reduced_error', 'fiber_area',
              'extraction_radius'
        psf_seeing_grid (np.ndarray, optional): 1D array of seeing values (arcsec)
            to tabulate the PSF interpolator. If not provided, uses
            ``np.linspace(0.5, 5.0, 45)``.

    Returns:
        tuple:
            - dar_model (callable): Function mapping (wavelength, fiber_x, fiber_y)
              to refraction-corrected positions; suitable for downstream use.
            - measured_fwhm (np.ndarray): Estimated FWHM (arcsec) vs wavelength.

    Raises:
        RuntimeError: If source detection or PSF/DAR fitting fails.
    """
    config_dict = prep['config_dict']
    fiber_x = prep['fiber_x']
    fiber_y = prep['fiber_y']
    def_wave = prep['def_wave']
    reduced_spectra = prep['reduced_spectra']
    reduced_error = prep['reduced_error']
    fiber_area = prep['fiber_area']
    extraction_radius = prep['extraction_radius']

    # PSF interpolator grid
    if psf_seeing_grid is None:
        psf_seeing_grid = np.linspace(0.5, 5.0, 45)
    interpolation_radius = 1.5 * extraction_radius
    r = np.linspace(0.0, interpolation_radius, 101)

    # 1) PSF interpolator
    psf_interp = psf.build_psf_interpolator(
        r,
        seeing=psf_seeing_grid,
        fiber_radius=config_dict['fiber_radius']
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

    return dar_model, measured_fwhm, psf_interp, sources, X, Y



def build_response(dataset_manifest, args):
    """Build the instrument response function.

    Args:
        dataset_manifest (dict): dataset manifest dictionary returned from build_dataset_from_reduced_files()
        output_path (str): output file path to which this method will write a reduced FITS file
        args (argparse.Namespace): command line arguments
    """

    # Loading the data inputs
    prep = prepare_response_inputs(dataset_manifest, args)

    # Construct the PSF and DAR model
    dar_model, measured_fwhm, psf_interp, sources, X, Y = build_psf_and_dar(prep)

    # Load calibration data
    cal_spectrum_table, extinction_table = load_calibration_data(args.standard_name)
    fiber_x = prep['fiber_x']
    fiber_y = prep['fiber_y']
    def_wave = prep['def_wave']
    reduced_spectra = prep['reduced_spectra']
    reduced_error = prep['reduced_error']
    header = prep['header']

    # Extract spectrum
    spectrum, spectrum_error = extract_optimal_spectrum(
        reduced_spectra, reduced_error, dar_model,
        sources, X, Y, measured_fwhm,
        fiber_x, fiber_y, psf_interp, def_wave
    )

    # Apply extinction correction
    spectrum, spectrum_error = apply_extinction(
        def_wave, spectrum, spectrum_error, extinction_table, float(header["AIRMASS"])
    )

    # Measure and apply response
    spectrum, spectrum_error, response = measure_and_apply_response(
        def_wave, spectrum, spectrum_error,
        cal_spectrum_table, args.output_folder, window=251
    )

