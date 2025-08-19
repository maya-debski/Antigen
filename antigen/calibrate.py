from typing import Any

import numpy as np
from numpy import ndarray, dtype
from scipy.signal import savgol_filter

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
