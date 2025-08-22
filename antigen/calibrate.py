from typing import Any
import logging

import numpy as np
from numpy import ndarray, dtype
from scipy.signal import savgol_filter

from antigen import config
from antigen import cube
from antigen import fiber
from antigen import io
from antigen import wavelength

logger = logging.getLogger('antigen.calibrate')

def apply_extinction_correction(wave, flux, flux_err, ext_wave, A_lambda_airmass, airmass):
    """Apply extinction correction to an observed spectrum.

    The correction uses the relation:

        F_true = F_obs * 10^(0.4 * A_lambda * airmass)

    where A_lambda is the extinction (in magnitudes) per unit airmass.

    Args:
        wave (ndarray): 1-D wavelength array of the observed spectrum
        flux (ndarray): 1-D or 2-D observed flux values (erg/s/Å/cm^2).
        flux_err (ndarray or None): 1-D or 2-D Observed flux uncertainties.
        ext_wave (ndarray): 1-D Wavelength grid of the extinction table.
        A_lambda_airmass (ndarray): 1-D Extinction in magnitudes per  airmass at each `ext_wave`.
        airmass (float): Airmass at which to apply the extinction correction.

    Returns:
        flux_corr (ndarray): Extinction-corrected flux on the observed wavelength grid.
        flux_err_corr (ndarray): Extinction-corrected flux  uncertainty on the same grid.
    """
    A_interp = np.interp(wave, ext_wave, A_lambda_airmass, left=A_lambda_airmass[0],  right=A_lambda_airmass[-1])
    factor = 10.0 ** (0.4 * A_interp * airmass)
    # Ensure factor can broadcast against flux (1D or 2D)
    if flux.ndim == 2:
        factor = factor[np.newaxis, :]  # shape (1, Nwave)
    flux_corr = flux * factor
    flux_err_corr = flux_err * factor
    return flux_corr, flux_err_corr


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


def build_response(dataset_manifest, output_path, standard_name, pixel_size):
    """Build the instrument response function.

    Args:
        dataset_manifest (dict): dataset manifest dictionary returned from build_dataset_from_reduced_files()
        output_path (str): output file path to which this method will write a reduced FITS file
        pixel_size (float): pixel size in arcsec.
    """
    config_dict = config.build_config_for_element(dataset_manifest['unit_instrument'].lower(),
                                                  dataset_manifest['unit_id'].upper())
    fiber_x, fiber_y = fiber.load_fiber_positions(dataset_manifest['unit_instrument'],
                                                  dataset_manifest['ndithers'],
                                                  dataset_manifest['dither_number'],
                                                  config_dict)

    def_wave = wavelength.get_rectified_wavelength(config_dict)

    reduced_spectra, reduced_error, header = io.load_reduced_data(dataset_manifest['in_folder'],
                                                                  dataset_manifest['reduced_files'])

    fiber_area = np.pi * config_dict['fiber_radius'] **2
    logger.info('Building Cube')
    data_cube, x_grid, y_grid = cube.make_cube(def_wave, reduced_spectra, fiber_x, fiber_y, fiber_area,
                                               pixel_size=pixel_size, method='gdw', k=7)
    io.write_cube('test.fits', data_cube, def_wave, header, x_grid, y_grid, pixel_size,
                  overwrite=True)
    calspectrum_table = io.load_calspec_spectrum(standard_name)
    # 2) load_data of standard star: Open reduced fiber spectra from a directory.
    # 3) build_psf_interpolator    : Build a Moffat PSF interpolator over radius.
    # 4) fit_psf                   : Fit source X/Y and FWHM vs wavelength.
    # 5) build_psf_weights         : Turn PSF model into fiber weights.
    # 6) get_spectrum              : Extract an optimal 1D spectrum.
    # 7) correction extinction     : Correct with an extinction curve + airmass.
    # 8) grab standard star spectrum: Load reference standard star SED.
    # 9) measure_response          : Compute the instrument response function.
