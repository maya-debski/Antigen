import numpy as np
from astropy.stats import mad_std
from astropy.convolution import Gaussian1DKernel, convolve
import logging

from antigen.sky import identify_sky_pixels
from antigen.trace import _get_peak_positions

logger = logging.getLogger('antigen.wavelength')

def get_wavelength(spectra, trace_positions, valid_fibers, arc_pixel_guesses, arc_wavelengths,
                   use_kernel=True, peak_threshold=None, reference_fiber_index=130):
    """
    Computes the wavelength solution for each fiber by identifying arc lamp emission line positions
    in each fiber's spectrum, fits a smooth curve across the fiber direction, and derives a polynomial
    mapping from pixel position to physical wavelength.

    Args:
        spectra (ndarray): 2D array of flux values; each row is a fiber spectrum.
        trace_positions (ndarray): 2D array with trace column positions for each fiber.
        valid_fibers (ndarray): Boolean array indicating which fibers have valid data.
        arc_pixel_guesses (ndarray): Initial pixel location guesses for each arc line.
        arc_wavelengths (ndarray): Known wavelengths of arc lines.
        use_kernel (bool): Whether to apply kernel smoothing when detecting peaks. Default is True.
        peak_threshold (float, optional): Minimum peak height to be considered a valid arc line. Default is None.
        reference_fiber_index (int): Index of the fiber used as the starting point. Default is 130.

    Returns:
        wavelength_solution (ndarray): Wavelength array for each fiber.
        residuals (ndarray): Residuals from trace fitting per line.
        fitted_trace_positions (ndarray): Smoothed trace-space positions of arc lines.
        arc_pixel_locations (ndarray): Raw arc line pixel positions per fiber.
    """
    logger.info("Starting wavelength solution fit for %d out of %d fibers",
                valid_fibers.sum(), spectra.shape[0])

    dispersion_axis = np.arange(trace_positions.shape[1])

    arc_pixel_locations = _compute_arc_positions(spectra, valid_fibers, reference_fiber_index, arc_pixel_guesses,
                                                 arc_wavelengths, use_kernel, peak_threshold)

    fitted_arc_line_positions, residuals = _fit_arc_line_positions(arc_pixel_locations, trace_positions, valid_fibers)

    wavelength_solution = _fit_wavelength_polynomials(fitted_arc_line_positions, arc_pixel_locations, arc_wavelengths,
                                                      valid_fibers, dispersion_axis)

    num_good_lines = np.sum(residuals < 0.2)
    logger.info("Wavelength calibration complete using %d out of %d arc lines",
                num_good_lines, len(residuals))

    return wavelength_solution, residuals, fitted_arc_line_positions, arc_pixel_locations


def _compute_arc_positions(spectra, valid_fibers, reference_index,
                            arc_pixel_guesses, arc_wavelengths, use_kernel, peak_threshold=None):
    """
    Detects arc line pixel positions for each fiber, starting from a reference and spreading out.

    Args:
        spectra (ndarray): 2D array of flux values; each row is a fiber spectrum.
        valid_fibers (ndarray): Boolean array indicating which fibers have valid data.
        reference_index (int): Index of the fiber used as the starting point. Default is 130.
        arc_pixel_guesses (ndarray): Initial pixel location guesses for each arc line.
        arc_wavelengths (ndarray): Known wavelengths of arc lines.
        use_kernel (bool): Whether to apply kernel smoothing when detecting peaks. Default is True.
        peak_threshold (float, optional): Minimum peak height to be considered a valid arc line. Default is None.

    Returns:
        arc_pixels (ndarray): Pixel locations of arc lines for each fiber.
    """
    num_fibers = spectra.shape[0]
    num_lines = len(arc_wavelengths)
    arc_pixels = np.zeros((num_fibers, num_lines))

    logger.debug("Starting arc line detection from reference fiber %d", reference_index)
    ref_flux = spectra[reference_index]
    mask, continuum = identify_sky_pixels(ref_flux, per=5)
    ref_flux_corrected = ref_flux - continuum
    if peak_threshold is None:
        noise = mad_std(ref_flux_corrected[~mask], ignore_nan=True)
        peak_threshold = 10. * noise
    arc_pixels[reference_index] = _get_arclines_fiber(ref_flux_corrected, arc_pixel_guesses,
                                                     limit=peak_threshold, use_kernel=use_kernel)

    for direction in [+1, -1]:
        for step in range(1, num_fibers):
            i = reference_index + direction * step
            if not (0 <= i <num_fibers):
                continue
            if not valid_fibers[i]:
                continue
            _, cont = identify_sky_pixels(spectra[i], per=5)
            flux_corrected = spectra[i] - cont
            last_fit = arc_pixels[i - direction]
            fitted_result = _get_arclines_fiber(flux_corrected, last_fit,
                                           limit=peak_threshold, use_kernel=use_kernel)
            arc_pixels[i] = fitted_result

    return arc_pixels


def _fit_arc_line_positions(arc_pixels, trace_positions, valid_fibers):
    """
    Fits 4th-order polynomials across fibers to model trace-space location of each arc line.

    Returns:
        fitted_arc_positions (ndarray): Fitted pixel locations per fiber and arc line.
        residuals (ndarray): MAD standard deviation of fit errors for each arc line.
    """
    num_fibers, num_lines = arc_pixels.shape
    dispersion_axis = np.arange(trace_positions.shape[1])
    fitted_arc_positions = np.zeros_like(arc_pixels)
    residuals = np.zeros(num_lines)

    logger.debug("Fitting arc lines across trace positions")

    for line_idx in range(num_lines):
        interpolated_x = np.zeros(num_fibers)

        # Fill bad fiber positions using nearest good fiber
        bad_indices = np.where(~valid_fibers)[0]
        good_indices = np.where(valid_fibers)[0]
        for b in bad_indices:
            arc_pixels[b, line_idx] = arc_pixels[
                good_indices[np.argmin(np.abs(b - good_indices))], line_idx
            ]

        # Interpolate each fiber’s trace-space position
        for fiber_idx in range(num_fibers):
            interpolated_x[fiber_idx] = np.interp(
                arc_pixels[fiber_idx, line_idx], dispersion_axis, trace_positions[fiber_idx]
            )

        valid = (arc_pixels[:, line_idx] > 0) & np.isfinite(interpolated_x)
        if np.sum(valid) < 5:
            logger.warning("Too few valid points to fit arc line %d", line_idx)
            continue

        poly_coeff = np.polyfit(interpolated_x[valid], arc_pixels[valid, line_idx], 4)
        fitted_arc_positions[:, line_idx] = np.polyval(poly_coeff, interpolated_x)
        residuals[line_idx] = mad_std(fitted_arc_positions[:, line_idx] - arc_pixels[:, line_idx],
                                      ignore_nan=True)

    return fitted_arc_positions, residuals


def _fit_wavelength_polynomials(fitted_trace_positions, arc_pixels, arc_wavelengths,
                                valid_fibers, dispersion_axis):
    """
    Computes final wavelength solution by fitting 5th-order polynomials per fiber.

    Returns:
        wavelength_solution (ndarray): 2D wavelength array per fiber.
    """
    num_fibers = fitted_trace_positions.shape[0]
    residuals = mad_std(fitted_trace_positions - arc_pixels, axis=0, ignore_nan=True)
    good_arc_lines = residuals < 0.2

    if np.sum(good_arc_lines) < 3:
        logger.warning("Fewer than 3 good arc lines detected — wavelength fit may be unstable.")

    wavelength_solution = np.zeros((num_fibers, len(dispersion_axis)))

    for fiber_idx in range(num_fibers):
        if not valid_fibers[fiber_idx]:
            continue
        coeffs = np.polyfit(fitted_trace_positions[fiber_idx][good_arc_lines],
                            arc_wavelengths[good_arc_lines], 5)
        wavelength_solution[fiber_idx] = np.polyval(coeffs, dispersion_axis)

    return wavelength_solution


def _get_arclines_fiber(spectrum, init_loc=None, limit=100, use_kernel=True):
    """
    Identifies arc line positions in a given spectrum by detecting peaks.

    Args:
        spectrum (ndarray): 1D array representing the spectrum of a fiber.
        init_loc (ndarray, optional): Initial guess locations for arc lines. Default is None.
        limit (float): Minimum peak value to consider a valid arc line. Default is 1000.
        use_kernel (bool): Whether to apply a box kernel convolution to smooth the spectrum. Default is True.

    Returns:
        ndarray: Array of arc line positions (pixel indices) in the spectrum.
    """

    # Apply box kernel convolution to smooth the spectrum if use_kernel is True
    if use_kernel:
        kernel = Gaussian1DKernel(1.0)
        convolved_spectrum = convolve(spectrum, kernel)
    else:
        convolved_spectrum = spectrum.copy()

    # Identify peaks in the spectrum by finding zero-crossings in the first derivative
    diff_array = convolved_spectrum[1:] - convolved_spectrum[:-1]
    loc = np.where((diff_array[:-1] > 0) & (diff_array[1:] < 0))[0]

    # Filter peaks based on the limit threshold
    peaks = convolved_spectrum[loc + 1]
    loc = loc[peaks > limit] + 1

    # Refine peak positions using quadratic interpolation
    loc = _get_peak_positions(convolved_spectrum, loc)

    # Match refined peak positions with initial guess locations, if provided
    if init_loc is not None:
        final_loc = []
        for i in init_loc:
            final_loc.append(loc[np.argmin(np.abs(np.array(loc) - i))])
        loc = final_loc

    return loc