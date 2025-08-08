import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.stats import sigma_clip, mad_std, biweight_location as biweight
import logging
from scipy.interpolate import interp1d
from scipy.ndimage import percentile_filter

logger = logging.getLogger('antigen.fiber')

def identify_sky_pixels(sky, per=50, size=50):
    """
    Identifies sky pixels by applying a percentile filter and sigma-clipping.

    Parameters:
        sky (array-like): Input sky intensity values.
        per (int, optional): Percentile value for the filter. Default is 50 (median).
        size (int, optional): Size of the filter window. Default is 50.

    Returns:
        tuple: A boolean mask array indicating sky pixels and the filtered continuum array.
    """
    # Apply a percentile filter to smooth the sky data and estimate the continuum
    cont = percentile_filter(sky, per, size=size)

    try:
        # Apply sigma-clipping to identify outliers (sky pixels)
        # Use MAD-based standard deviation for robust statistics
        mask = sigma_clip(sky - cont, masked=True, maxiters=None,
                          stdfunc=mad_std, sigma=5)
    except:
        # Fallback for older versions of sigma_clip where maxiters was iters
        mask = sigma_clip(sky - cont, iters=None, stdfunc=mad_std, sigma=5)

    # Return the mask (True for sky pixels) and the filtered continuum
    return mask.mask, cont


def rectify(science_spectra, error_spectra, wave_all, def_wave):
    """
    Rectifies scientific and error spectra by interpolating them onto a common wavelength grid.

    Parameters:
        science_spectra (2D array): Array of scientific spectra to be rectified.
        error_spectra (2D array): Corresponding error spectra for each scientific spectrum.
        wave_all (2D array): Wavelength grids corresponding to each input spectrum.
        def_wave (1D array): Target wavelength grid for interpolation.

    Returns:
        tuple:
            - science_rect (2D array): Rectified scientific spectra on the target wavelength grid.
            - error_rect (2D array): Rectified error spectra on the target wavelength grid.
    """
    # Initialize arrays to store rectified scientific spectra and errors
    science_rect = np.zeros((science_spectra.shape[0], len(def_wave)))
    error_rect = np.zeros((science_spectra.shape[0], len(def_wave)))

    # Loop through each spectrum to interpolate onto the target wavelength grid
    for i in np.arange(science_spectra.shape[0]):
        # Compute wavelength bin sizes for flux normalization
        diff_wave = np.diff(wave_all[i])
        diff_wave = np.hstack([diff_wave[0], diff_wave])  # Ensure length matches the wavelength array

        # Interpolate the scientific spectrum, normalizing by wavelength bin size
        science_rect[i] = np.interp(def_wave, wave_all[i], science_spectra[i] / diff_wave,
                               left=np.nan, right=np.nan)

        # Interpolate the error spectrum, normalizing by wavelength bin size
        error_rect[i] = np.interp(def_wave, wave_all[i], error_spectra[i] / diff_wave,
                                 left=np.nan, right=np.nan)

    return science_rect, error_rect


def get_fiber_to_fiber(spectrum, n_chunks=100):
    """
    Computes the fiber-to-fiber correction by normalizing each fiber's spectrum
    to the average spectrum across all fibers, then smooths the correction factor
    using interpolation.

    Parameters:
        spectrum (2D array): Array of spectra from multiple fibers (fibers x wavelength).
        n_chunks (int, optional): Number of chunks to split the wavelength range into
                                  for smoothing. Default is 100.

    Returns:
        initial_ftf (2D array): Initial fiber-to-fiber correction factors.
        ftf (2D array): Smoothed fiber-to-fiber correction factors.
    """
    # Compute the average spectrum across all fibers using a robust biweight statistic
    average = biweight(spectrum, axis=0, ignore_nan=True)

    # Calculate the initial fiber-to-fiber correction by dividing each fiber by the average spectrum
    initial_ftf = spectrum / average[np.newaxis, :]

    # Create a wavelength grid and divide it into chunks for smoothing
    columns = np.arange(spectrum.shape[1])
    chunked_columns = np.array([np.mean(chunk) for chunk in np.array_split(columns, n_chunks)])

    # Initialize the smoothed correction array
    ftf = spectrum * 0.

    # Loop through each fiber to compute the smoothed correction factor
    for i in np.arange(len(spectrum)):
        # Compute the biweight statistic for each chunk of the initial correction factor
        chunked_ftf = np.array([biweight(chunk, ignore_nan=True) for chunk in np.array_split(initial_ftf[i], n_chunks)])

        # Select valid (finite) values for interpolation
        sel = np.isfinite(chunked_ftf)
        if sel.sum() == 0.:
            continue
        # Interpolate the correction factor using quadratic interpolation
        interp_func = interp1d(chunked_columns[sel], chunked_ftf[sel], kind='quadratic', bounds_error=False,
                     fill_value='extrapolate')

        # Apply the interpolation to the full wavelength range
        ftf[i] = interp_func(columns)

    # Return both the initial and smoothed fiber-to-fiber correction factors
    return initial_ftf, ftf


def get_wavelength(spectra, trace_positions, valid_fibers, arc_pixel_guesses, arc_wavelengths,
                   use_kernel=True, peak_threshold=100, reference_fiber_index=130):
    """
    Computes the wavelength solution for each fiber in a spectrograph based on trace and spectral data.

    Args:
        spectra (ndarray): 2D array of flux values; each row is a fiber spectrum.
        trace_positions (ndarray): 2D array with trace column positions for each fiber.
        valid_fibers (ndarray): Boolean array indicating which fibers have valid data.
        arc_pixel_guesses (ndarray): Initial pixel location guesses for each arc line.
        arc_wavelengths (ndarray): Known wavelengths of arc lines.
        use_kernel (bool): Whether to apply kernel smoothing when detecting peaks. Default is True.
        peak_threshold (float): Minimum peak height to be considered a valid arc line. Default is 100.
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
                            arc_pixel_guesses, arc_wavelengths, use_kernel, peak_threshold):
    """
    Detects arc line pixel positions for each fiber, starting from a reference and spreading out.

    Args:
        spectra (ndarray): 2D array of flux values; each row is a fiber spectrum.
        valid_fibers (ndarray): Boolean array indicating which fibers have valid data.
        reference_index (int): Index of the fiber used as the starting point. Default is 130.
        arc_pixel_guesses (ndarray): Initial pixel location guesses for each arc line.
        arc_wavelengths (ndarray): Known wavelengths of arc lines.
        use_kernel (bool): Whether to apply kernel smoothing when detecting peaks. Default is True.
        peak_threshold (float): Minimum peak height to be considered a valid arc line. Default is 100.

    Returns:
        arc_pixels (ndarray): Pixel locations of arc lines for each fiber.
    """
    num_fibers = spectra.shape[0]
    num_lines = len(arc_wavelengths)
    arc_pixels = np.zeros((num_fibers, num_lines))

    logger.debug("Starting arc line detection from reference fiber %d", reference_index)
    ref_flux = spectra[reference_index]
    _, continuum = identify_sky_pixels(ref_flux, per=5)
    ref_flux_corrected = ref_flux - continuum
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
        B = Gaussian1DKernel(1.0)
        y1 = convolve(spectrum, B)
    else:
        y1 = spectrum.copy()

    # Identify peaks in the spectrum by finding zero-crossings in the first derivative
    diff_array = y1[1:] - y1[:-1]
    loc = np.where((diff_array[:-1] > 0) & (diff_array[1:] < 0))[0]

    # Filter peaks based on the limit threshold
    peaks = y1[loc + 1]
    loc = loc[peaks > limit] + 1

    # Refine peak positions using quadratic interpolation
    loc = _get_peak_positions(y1, loc)

    # Match refined peak positions with initial guess locations, if provided
    if init_loc is not None:
        final_loc = []
        for i in init_loc:
            final_loc.append(loc[np.argmin(np.abs(np.array(loc) - i))])
        loc = final_loc

    return loc

def get_spectra(array_flt, array_trace, npix=5):
    """
    Extract spectra by dividing the flat field and averaging the central pixels.

    Parameters
    ----------
    array_flt : 2D numpy array
        Twilight image.
    array_trace : 2D numpy array
        Trace for each fiber.
    npix : int, optional
        Number of pixels to extract around the trace center. Default is 5.

    Returns
    -------
    spec : 2D numpy array
        Extracted and rectified spectrum for each fiber.
    """

    # Initialize the output spectrum array
    spec = np.zeros((array_trace.shape[0], array_trace.shape[1]))

    # Get the number of rows in the flat field image
    N = array_flt.shape[0]

    # Create an array of x-axis pixel indices
    x = np.arange(array_flt.shape[1])

    # Calculate the lower and upper bounds for pixel extraction
    LB = int((npix + 1) / 2)  # Lower bound
    HB = -LB + npix + 1       # Upper bound

    # Iterate through each fiber
    for fiber in np.arange(array_trace.shape[0]):

        # Skip fibers with trace positions too close to the image edges
        if np.round(array_trace[fiber]).min() < LB:
            continue
        if np.round(array_trace[fiber]).max() >= (N - LB):
            continue

        # Convert trace positions to integer indices
        indv = np.round(array_trace[fiber]).astype(int)

        # Iterate through pixels around the trace center
        for j in np.arange(-LB, HB):

            # Calculate weight for the lower boundary pixel
            if j == -LB:
                w = indv + j + 1 - (array_trace[fiber] - npix / 2.)

            # Calculate weight for the upper boundary pixel
            elif j == HB - 1:
                w = (npix / 2. + array_trace[fiber]) - (indv + j)

            # Assign weight 1 for central pixels
            else:
                w = 1.

            # Add the weighted pixel values to the spectrum
            spec[fiber] += array_flt[indv + j, x] * w

    # Normalize the spectrum by the number of extracted pixels
    return spec


def get_spectra_error(array_flt, array_trace, npix=5):
    """
    Extract spectra by dividing the flat field and averaging the central
    two pixels

    Parameters
    ----------
    array_flt : 2d numpy array
        Twilight image
    array_trace : 2d numpy array
        Trace for each fiber
    npix : int, optional
        Number of pixels for averaging (default is 5)

    Returns
    -------
    twi_spectrum : 2d numpy array
        Rectified twilight spectrum for each fiber
    """

    # Initialize spectrum array to store extracted spectra
    spec = np.zeros((array_trace.shape[0], array_trace.shape[1]))

    # Get number of rows in the flat field image
    N = array_flt.shape[0]

    # Create an array of x-coordinates for the flat field image
    x = np.arange(array_flt.shape[1])

    # Calculate bounds for pixel averaging
    LB = int((npix + 1) / 2)
    HB = -LB + npix + 1

    # Iterate over each fiber to extract its spectrum
    for fiber in np.arange(array_trace.shape[0]):
        # Skip fibers with traces too close to image edges
        if np.round(array_trace[fiber]).min() < LB:
            continue
        if np.round(array_trace[fiber]).max() >= (N - LB):
            continue

        # Convert trace positions to integer indices
        indv = np.round(array_trace[fiber]).astype(int)

        # Loop over neighboring pixels for averaging
        for j in np.arange(-LB, HB):
            if j == -LB:
                # Calculate weight for lower boundary pixels
                w = indv + j + 1 - (array_trace[fiber] - npix / 2.)
            elif j == HB - 1:
                # Calculate weight for upper boundary pixels
                w = (npix / 2. + array_trace[fiber]) - (indv + j)
            else:
                # Set weight to 1 for central pixels
                w = 1.

            # Accumulate weighted sum of squared values from the flat field
            spec[fiber] += array_flt[indv + j, x] ** 2 * w

    # Return the root mean square error normalized by npix
    return np.sqrt(spec)


def get_spectra_chi2(array_flt, array_sci, array_err, array_trace, npix=5):
    """
    Extract spectra by dividing the flat field and averaging the central
    two pixels

    Parameters
    ----------
    array_flt : 2d numpy array
        Twilight image
    array_sci : 2d numpy array
        Science image
    array_err : 2d numpy array
        Error estimate for each pixel
    array_trace : 2d numpy array
        Trace for each fiber
    npix : int, optional
        Number of pixels for averaging (default is 5)

    Returns
    -------
    spec : 2d numpy array
        Chi-squared spectra for each fiber
    """

    # Initialize spectrum array to hold chi-squared values
    spec = np.zeros((array_trace.shape[0], array_trace.shape[1]))

    # Get the number of rows in the flat field image
    N = array_flt.shape[0]

    # Create an array of x-coordinates for the images
    x = np.arange(array_flt.shape[1])

    # Calculate bounds for pixel averaging
    LB = int((npix + 1) / 2)
    HB = -LB + npix + 1

    # Iterate over each fiber to extract its chi-squared spectrum
    for fiber in np.arange(array_trace.shape[0]):
        # Initialize a chi-squared array with shape (npix+1, 3, len(x))
        chi2 = np.zeros((npix + 1, 3, len(x)))

        # Skip fibers with traces too close to the image edges
        if np.round(array_trace[fiber]).min() < LB:
            continue
        if np.round(array_trace[fiber]).max() >= (N - LB):
            continue

        # Convert trace positions to integer indices
        indv = np.round(array_trace[fiber]).astype(int)

        # Loop over neighboring pixels for averaging
        for j in np.arange(-LB, HB):
            # Calculate weights for boundary pixels
            if j == -LB:
                w = indv + j + 1 - (array_trace[fiber] - npix / 2.)
            elif j == HB - 1:
                w = (npix / 2. + array_trace[fiber]) - (indv + j)
            else:
                # Use a weight of 1 for central pixels
                w = 1.

            # Apply weights to science, flat field, and error images
            chi2[j + LB, 0] = array_sci[indv + j, x] * w
            chi2[j + LB, 1] = array_flt[indv + j, x] * w
            chi2[j + LB, 2] = array_err[indv + j, x] * w

        # Compute the normalization factor for the flux
        norm = chi2[:, 0].sum(axis=0) / chi2[:, 1].sum(axis=0)

        # Calculate the chi-squared numerator: (data - model)^2
        num = (chi2[:, 0] - chi2[:, 1] * norm[np.newaxis, :]) ** 2

        # Calculate the denominator: (error + regularization term)^2
        denom = (chi2[:, 2] + 0.01 * chi2[:, 0].sum(axis=0)[np.newaxis, :]) ** 2

        # Compute the chi-squared value for each fiber
        spec[fiber] = 1. / (1. + npix) * np.sum(num / denom, axis=0)

    # Return the final chi-squared spectrum array
    return spec


def _get_peak_positions(column_profile, peak_x_positions):
    """
    Refine the vertical (y) positions of fiber peaks using a 3-point quadratic interpolation.

    Args:
        column_profile (np.ndarray): 1D array representing a collapsed spatial profile (e.g., mean over chunk).
        peak_x_positions (np.ndarray): Approximate x-positions of detected peaks in the profile.

    Returns:
        refined_peaks (np.ndarray): Subpixel-accurate y-positions of fiber peaks.
    """
    y_range = np.arange(column_profile.shape[0])
    offsets = np.array([
        peak_x_positions - 1,
        peak_x_positions,
        peak_x_positions + 1
    ], dtype=int)

    # Remove any positions that fall out of bounds
    valid = (offsets >= 0) & (offsets < column_profile.shape[0])
    valid = valid.all(axis=0)
    offsets = offsets[:, valid]

    refined_peaks = y_range[offsets[1]] - (
        (column_profile[offsets[2]] - column_profile[offsets[0]]) /
        (2. * (column_profile[offsets[2]] - 2. * column_profile[offsets[1]] + column_profile[offsets[0]]))
    )
    return refined_peaks


def _estimate_offset(nominal, detected, max_sep=3.0, offset_range=(-10, 10)):
    """
    Estimate a global vertical offset between nominal and detected peak positions.

    This function tests a range of offsets and selects the one that
    maximizes the cross-correlation between the first derivatives of the nominal and detected positions

    Args:
        nominal (np.ndarray): 1D array of expected fiber y-positions.
        detected (np.ndarray): 1D array of detected peak y-positions.
        max_sep (float): Maximum allowed distance (in pixels) for a peak to be considered a match.
        offset_range (tuple): Tuple of (min_offset, max_offset) to search (inclusive of min, exclusive of max).

    Returns:
        best_offset (float): The offset (in pixels) that gives the best alignment between nominal and detected positions.
    """
    nominal = np.asarray(nominal)
    detected = np.asarray(detected)

    # Interpolate spacing (first differences) onto a common grid
    nominal_diff = np.diff(nominal)
    detected_diff = np.diff(detected)

    # Use y-positions of centers between fibers as the x-axis for diffs
    nominal_x = 0.5 * (nominal[1:] + nominal[:-1])
    detected_x = 0.5 * (detected[1:] + detected[:-1])

    # Define a common y-grid for interpolation
    x_min = int(min(nominal_x[0], detected_x[0])) - 10
    x_max = int(max(nominal_x[-1], detected_x[-1])) + 10
    x_grid = np.arange(x_min, x_max + 1)

    # Interpolate the diffs onto the common grid
    nominal_interp = np.interp(x_grid, nominal_x, nominal_diff, left=0, right=0)
    detected_interp = np.interp(x_grid, detected_x, detected_diff, left=0, right=0)

    # Cross-correlate over range of integer offsets
    scores = []
    test_range = np.linspace(offset_range[0], offset_range[1], 100)
    for offset in test_range:
        shifted_detected = np.interp(x_grid + offset, x_grid, detected_interp)
        score = np.sum(nominal_interp * shifted_detected)
        scores.append(score)

    best_offset = test_range[np.argmax(scores)]
    return best_offset


def _match_peaks_to_nominal(nominal, detected, offset, max_sep=3.0):
    """
    Match detected peaks to nominal positions using a global offset and nearest-neighbor search.

    Each nominal position is shifted by the given offset, and then the nearest
    unmatched detected peak within max_sep is assigned to that nominal fiber.

    Args:
        nominal (np.ndarray): 1D array of expected fiber y-positions.
        detected (np.ndarray): 1D array of detected peak y-positions.
        offset (float): Global vertical offset to apply to nominal positions before matching.
        max_sep (float): Maximum allowed distance (in pixels) to consider a peak a valid match.

    Returns:
        matches (np.ndarray): Array of matched peak positions (same shape as nominal), with np.nan for unmatched fibers.
    """
    corrected_nominal = nominal + offset
    matches = np.full_like(nominal, fill_value=np.nan, dtype=np.float64)
    used = set()

    for i, nom_y in enumerate(corrected_nominal):
        diffs = np.abs(detected - nom_y)
        candidates = [j for j in np.argsort(diffs) if diffs[j] <= max_sep and j not in used]
        if candidates:
            j = candidates[0]
            matches[i] = detected[j]
            used.add(j)
    return matches

def _evaluate_trace_chunk(flat_column, x_center, chunk_index, nominal_positions, good_fiber_indices, exclude_fiber_mask,
                          total_fibers, max_sep=4.0):
    """
    Process a single image chunk: detect peaks, estimate offset, match to nominal,
    and extrapolate excluded fiber positions.

    Args:
        flat_column (np.ndarray): 1D median profile of the current chunk.
        x_center (float): X position (center) of this chunk.
        chunk_index (int): Index of the chunk.
        nominal_positions (np.ndarray): Expected positions for each fiber.
        good_fiber_indices (np.ndarray): Indices of fibers not excluded.
        exclude_fiber_mask (np.ndarray): Boolean mask for excluded fibers.
        total_fibers (int): Total number of fibers.
        max_sep (float): Maximum matching distance.

    Returns:
        fiber_trace: 1D array of matched/refined fiber positions (length = total_fibers).
    """
    logger.debug(f"Processing chunk {chunk_index+1} at x={x_center:.1f}")

    diff = np.diff(flat_column)
    peak_indices = np.where((diff[:-1] > 0) & (diff[1:] < 0))[0]
    peak_indices = peak_indices[peak_indices > 2]

    if len(peak_indices) == 0:
        logger.warning(f"No initial peaks found at chunk {x_center:.1f}")
        return np.zeros(total_fibers)

    peak_values = flat_column[peak_indices + 1]
    strong_peaks = peak_indices[peak_values > 0.3 * np.median(peak_values)] + 1

    if len(strong_peaks) == 0:
        logger.warning(f"No strong peaks found at chunk {x_center:.1f}")
        return np.zeros(total_fibers)

    refined_peaks = _get_peak_positions(flat_column, strong_peaks)
    nominal_good = nominal_positions[good_fiber_indices]

    offset = _estimate_offset(nominal_good, refined_peaks, max_sep=max_sep)
    matched_peaks = _match_peaks_to_nominal(nominal_good, refined_peaks, offset, max_sep=max_sep)

    fiber_trace = np.zeros(total_fibers)
    fiber_trace[good_fiber_indices] = np.where(np.isnan(matched_peaks), 0.0, matched_peaks)

    for excluded_idx in np.where(exclude_fiber_mask)[0]:
        nearest_good = np.argmin(np.abs(excluded_idx - good_fiber_indices))
        ref_idx = good_fiber_indices[nearest_good]
        if fiber_trace[ref_idx] != 0.0:
            fiber_trace[excluded_idx] = (
                fiber_trace[ref_idx] +
                nominal_positions[excluded_idx] -
                nominal_positions[ref_idx]
            )

    return fiber_trace


def get_trace(flat_image, nominal_trace_positions, exclude_fiber_mask, num_chunks=80):
    """
    Extract smoothed fiber traces from a twilight flat field image.

    Args:
        flat_image (np.ndarray): 2D array representing the flat field image.
        nominal_trace_positions (np.ndarray): 1D array of expected y-positions for each fiber.
        exclude_fiber_mask (np.ndarray): 1D boolean array marking fibers to exclude from tracing.
        num_chunks (int): Number of chunks to bin image in dispersion direction

    Returns:
        trace_array (np.ndarray): 2D array (fibers × x-pixels) of smoothed fiber traces.
        good_fiber_mask (np.ndarray): 1D boolean array indicating usable fibers.
        raw_trace_matrix (np.ndarray): 2D array (fibers × chunks) of raw trace positions.
        x_chunk_centers (np.ndarray): 1D array of x-axis chunk centers used for trace extraction.
    """
    total_fibers = len(nominal_trace_positions)
    good_fiber_indices = np.where(~exclude_fiber_mask)[0]

    logger.info("Beginning fiber trace extraction")
    image = flat_image

    chunk_indices = np.array_split(np.arange(image.shape[1]), num_chunks)
    x_chunk_centers = np.array([np.mean(chunk) for chunk in chunk_indices])
    chunked_image = np.array_split(image, num_chunks, axis=1)

    # Discard edges to avoid poor flat data
    chunked_image = chunked_image[1:-1]
    x_chunk_centers = x_chunk_centers[1:-1]

    flattened_chunks = [np.median(chunk, axis=1) for chunk in chunked_image]
    raw_trace_matrix = np.zeros((total_fibers, len(flattened_chunks)))

    center_idx = len(flattened_chunks) // 2
    center_column = flattened_chunks[center_idx]
    x_center = x_chunk_centers[center_idx]

    # Use nominal positions as-is
    center_trace = _evaluate_trace_chunk(center_column, x_center, center_idx, nominal_trace_positions,
                                         good_fiber_indices,  exclude_fiber_mask, total_fibers)
    raw_trace_matrix[:, center_idx] = center_trace
    initial_offset = np.nanmedian(center_trace - nominal_trace_positions)
    logger.info(f"Best fiber position offset for central chunk: {initial_offset:.2f}")

    # Estimate global offset once using center

    for direction in [+1, -1]:
        for step in range(1, len(flattened_chunks)):
            k = center_idx + direction * step
            if not (0 <= k < len(flattened_chunks)):
                continue

            flat_column = flattened_chunks[k]
            x_center = x_chunk_centers[k]

            nominal_estimate = raw_trace_matrix[:, k - direction]

            trace = _evaluate_trace_chunk(flat_column, x_center, k, nominal_estimate, good_fiber_indices,
                                          exclude_fiber_mask, total_fibers)
            raw_trace_matrix[:, k] = trace

    logger.info("Smoothing fiber traces with 7th-order polynomial")
    trace_array = np.zeros((total_fibers, image.shape[1]))
    x_full = np.arange(image.shape[1])

    for fiber_index in range(total_fibers):
        valid = raw_trace_matrix[fiber_index] != 0.0
        if np.any(valid):
            poly_coeffs = np.polyfit(x_chunk_centers[valid], raw_trace_matrix[fiber_index, valid], deg=7)
            trace_array[fiber_index] = np.polyval(poly_coeffs, x_full)
        else:
            logger.debug(f"Fiber {fiber_index} has no valid trace values.")
    good_fiber_mask = ~exclude_fiber_mask
    logger.info("Trace extraction complete")

    from astropy.io import fits
    fits.PrimaryHDU(raw_trace_matrix).writeto('test.fits', overwrite=True)
    return trace_array, good_fiber_mask, raw_trace_matrix, x_chunk_centers
