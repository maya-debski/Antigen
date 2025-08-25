import logging
from scipy.interpolate import interp1d

logger = logging.getLogger('antigen.extinction')

def compute_extinction_factor(obs_wave, ext_wave, A_lambda_airmass, airmass):
    """Compute extinction scaling factor for an observed spectrum.

    Args:
        obs_wave (ndarray): 1-D array of observed wavelengths.
        ext_wave (ndarray): 1-D wavelength grid of the extinction table.
        A_lambda_airmass (ndarray): 1-D extinction in magnitudes per airmass
            at each `ext_wave`.
        airmass (float): Airmass at which to apply the extinction correction.

    Returns:
        ndarray: Extinction correction factor on the observed wavelength grid.
    """
    interp = interp1d(ext_wave, A_lambda_airmass, bounds_error=False, fill_value="extrapolate")
    A_lambda = interp(obs_wave)
    factor = 10 ** (0.4 * A_lambda * airmass)
    return factor


def apply_extinction_correction(obs_wave, flux, flux_err, ext_wave, A_lambda_airmass, airmass):
    """Apply extinction correction to an observed spectrum.

    Args:
        obs_wave (ndarray): 1-D wavelength array of the observed spectrum.
        flux (ndarray): 1-D or 2-D observed flux values.
        flux_err (ndarray or None): 1-D or 2-D observed flux uncertainties.
        ext_wave (ndarray): 1-D wavelength grid of the extinction table.
        A_lambda_airmass (ndarray): 1-D extinction in magnitudes per airmass
            at each `ext_wave`.
        airmass (float): Airmass at which to apply the extinction correction.

    Returns:
         flux_corr (ndarray): Extinction-corrected flux values.
         flux_err_corr (ndarray): Extinction-corrected flux uncertainties.
    """
    factor = compute_extinction_factor(obs_wave, ext_wave, A_lambda_airmass, airmass)
    flux_corr = flux * factor
    flux_err_corr = flux_err * factor

    logger.info("Applied extinction correction at airmass %.2f", airmass)
    return flux_corr, flux_err_corr
