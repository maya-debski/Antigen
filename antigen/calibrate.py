import logging
import warnings

import numpy as np
from numpy import ndarray, dtype
from scipy.signal import savgol_filter

from antigen import config
from antigen import cube
from antigen import dar
from antigen import detection
from antigen import extinction
from antigen import fiber
from antigen import io
from antigen import plot
from antigen import psf
from antigen import spectra
from antigen import wavelength


warnings.filterwarnings("ignore")
logger = logging.getLogger('antigen.calibrate')


def measure_response(obs_wave, obs_flux, std_wave, std_flux, window=51):
    """Compute the instrument response function from a standard star.

    Response is defined as:

        R(λ) = F_std(λ) / F_obs(λ)

    Optionally, a simple moving-average smoothing can be applied.

    Args:
        obs_wave (ndarray): Wavelength grid of the observed (extinction-corrected) spectrum.
        obs_flux (ndarray): Observed (corrected) flux values.
        std_wave (ndarray): Wavelength grid of the reference standard star spectrum.
        std_flux (ndarray): Reference flux of the standard star
        window (int, optional): Window size for moving-average smoothing (default 11).

    Returns:
        response (ndarray): Instrument response values.
    """
    std_on_obs = np.interp(obs_wave, std_wave, std_flux)
    with np.errstate(divide="ignore", invalid="ignore"):
        response = np.where(obs_flux > 0, std_on_obs / obs_flux, np.nan)
    finite_values = np.isfinite(response)
    response = np.interp(obs_wave, obs_wave[finite_values], response[finite_values],
                        left=response[finite_values][0],  right=response[finite_values][-1])
    response = savgol_filter(response, window, 3)
    return response

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
        fiber_x, fiber_y, collapsed_fiber_flux, fiber_area, bounds=bounds
    )

    sources = detection.detect_sources(detection_image, brightest_only=True)
    if len(sources) == 0:
        raise RuntimeError("No sources detected in the collapsed fiber image.")

    j, i = (int(sources['xcentroid']), int(sources['ycentroid']))
    x_coord = X[i, j]
    y_coord = Y[i, j]

    logger.info("Detected brightest source near %.1f, %.1f", x_coord, y_coord)

    return sources, x_coord, y_coord, X, Y

def fit_psf_and_build_dar_model(reduced_spectra, reduced_error, fiber_x, fiber_y,
                                psf_interp, x_coord, y_coord, def_wave,
                                extraction_radius=2.0, nchunks=20):
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

    Returns:
        tuple:
            dar_model (DARModel): Differential atmospheric refraction model.
            measured_fwhm (ndarray): FWHM of PSF per chunk.
    """
    params = psf.fit_psf(
        reduced_spectra, reduced_error, fiber_x, fiber_y, psf_interp,
        x_coord, y_coord, def_wave,
        extraction_radius=extraction_radius, Nchunks=nchunks
    )
    _, _, measured_fwhm, coeff_x, coeff_y, _ = params
    dar_model = dar.DARModel(coeff_x, coeff_y)
    return dar_model, measured_fwhm


def build_and_write_cube(def_wave, reduced_spectra, fiber_x, fiber_y, fiber_area,
                         header, dar_model, pixel_size, output_file="test.fits"):
    """
    Build a data cube from fiber spectra and write it to disk.

    Args:
        def_wave (ndarray): Wavelength grid.
        reduced_spectra (ndarray): Flux array of shape (Nfibers, Nlambda).
        fiber_x (ndarray): Fiber X positions.
        fiber_y (ndarray): Fiber Y positions.
        fiber_area (float): Effective area of fibers.
        header (fits.Header): Header metadata for FITS file.
        dar_model (DARModel): DAR model for correcting positions.
        pixel_size (float): Cube pixel size in arcsec.
        output_file (str, optional): Path for output FITS file. Defaults to "test.fits".

    Returns:
        tuple:
            data_cube (ndarray): Constructed 3D data cube.
            x_grid (ndarray): Cube X coordinate grid.
            y_grid (ndarray): Cube Y coordinate grid.
    """
    logger.info("Building Cube")
    data_cube, x_grid, y_grid = cube.make_cube(
        def_wave, reduced_spectra, fiber_x, fiber_y, fiber_area,
        pixel_size=pixel_size, method="gdw", k=15,
        dar_model=dar_model, sigma=2.0
    )
    io.write_cube(output_file, data_cube, def_wave, header,
                  x_grid, y_grid, pixel_size, overwrite=True)
    return data_cube, x_grid, y_grid


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
        reduced_spectra (ndarray): Flux array of shape (Nfibers, Nlambda).
        reduced_error (ndarray): Error array of shape (Nfibers, Nlambda).
        dar_model (DARModel): DAR model for source position correction.
        sources (Table): Detected source catalog.
        X (ndarray): Grid of X coordinates.
        Y (ndarray): Grid of Y coordinates.
        measured_fwhm (ndarray): FWHM array from PSF fit.
        fiber_x (ndarray): Fiber X positions.
        fiber_y (ndarray): Fiber Y positions.
        psf_interp (callable): PSF interpolator.
        def_wave (ndarray): Wavelength grid.

    Returns:
        tuple:
            spectrum (ndarray): Extracted flux spectrum.
            spectrum_error (ndarray): Extracted error spectrum.
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
        def_wave (ndarray): Wavelength grid.
        spectrum (ndarray): Flux spectrum.
        spectrum_error (ndarray): Error spectrum.
        extinction_table (Table): Table with 'wavelength' and 'mags_per_airmass'.
        airmass (float): Observed airmass.

    Returns:
        tuple:
            spectrum (ndarray): Extinction-corrected spectrum.
            spectrum_error (ndarray): Extinction-corrected errors.
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
        def_wave (ndarray): Wavelength grid.
        spectrum (ndarray): Flux spectrum.
        spectrum_error (ndarray): Error spectrum.
        cal_spectrum_table (Table): CALSPEC standard star spectrum.
        output_folder (str): Output folder for plots.
        window (int, optional): Smoothing window for response. Defaults to 251.

    Returns:
        tuple:
            spectrum (ndarray): Response-corrected spectrum.
            spectrum_error (ndarray): Response-corrected errors.
            response (ndarray): Derived instrument response function.
    """
    response = measure_response(
        def_wave, spectrum,
        cal_spectrum_table["wavelength"], cal_spectrum_table["flux"],
        window=window
    )

    spectrum *= response
    spectrum_error *= response

    plot.plot_spectrum_with_standard(
        spectrum, spectrum_error, def_wave,
        cal_spectrum_table["wavelength"], cal_spectrum_table["flux"],
        1. / response, outfolder=output_folder
    )

    return spectrum, spectrum_error, response



def build_response(dataset_manifest, args):
    """Build the instrument response function.

    Args:
        dataset_manifest (dict): dataset manifest dictionary returned from build_dataset_from_reduced_files()
        output_path (str): output file path to which this method will write a reduced FITS file
        args (argparse.Namespace): command line arguments
    """
    # Loading the config dictionary
    config_dict = config.build_config_for_element(dataset_manifest['unit_instrument'].lower(),
                                                  dataset_manifest['unit_id'].upper())

    # Adjust the fiber positions for the dithering of the observation
    fiber_x, fiber_y = fiber.load_fiber_positions(dataset_manifest['unit_instrument'],
                                                  dataset_manifest['ndithers'],
                                                  dataset_manifest['dither_number'],
                                                  config_dict)

    # Grab the rectified wavelength from parameters in config_dict
    def_wave = wavelength.get_rectified_wavelength(config_dict)

    # Load the reduced frames and vertically stack them together
    reduced_spectra, reduced_error, header = io.load_reduced_data(dataset_manifest['in_folder'],
                                                                  dataset_manifest['reduced_files'])

    # A pre-requisite for PSF extraction is the area of a fiber
    fiber_area = np.pi * config_dict['fiber_radius'] **2

    # Interpolation radius should be larger than extraction radius or an error will occur
    extraction_radius = args.extraction_radius
    interpolation_radius = extraction_radius * 1.5

    # Construct the PSF interpolator
    r = np.linspace(0, interpolation_radius, 101)
    psf_interp = psf.build_psf_interpolator(r, seeing=np.linspace(0.5, 5.0, 45),
                                            fiber_radius=config_dict['fiber_radius'])

    # Detect the brightest source in a collapsed frame
    sources, x_coord, y_coord, X, Y = detect_brightest_source(fiber_x, fiber_y, reduced_spectra, fiber_area)

    # Fit PSF and build DAR model
    dar_model, measured_fwhm = fit_psf_and_build_dar_model(
        reduced_spectra, reduced_error, fiber_x, fiber_y,
        psf_interp, x_coord, y_coord, def_wave,
        extraction_radius=extraction_radius, nchunks=20
    )

    # Build and save cube
    data_cube, x_grid, y_grid = build_and_write_cube(
        def_wave, reduced_spectra, fiber_x, fiber_y, fiber_area,
        header, dar_model, args.pixel_size, output_file="test.fits"
    )

    # Load calibration data
    cal_spectrum_table, extinction_table = load_calibration_data(args.standard_name)

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

