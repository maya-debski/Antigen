import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.io import fits
from astropy.stats import sigma_clip, mad_std, biweight_location as biweight
import logging
from scipy.interpolate import interp1d
from scipy.ndimage import percentile_filter

from antigen import config

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

    arc_pixel_locations = _compute_arc_positions(
        spectra, trace_positions, valid_fibers, reference_fiber_index,
        arc_pixel_guesses, arc_wavelengths, use_kernel, peak_threshold
    )

    fitted_arc_line_positions, residuals = _fit_arc_line_positions(
        arc_pixel_locations, trace_positions, valid_fibers
    )

    wavelength_solution = _fit_wavelength_polynomials(
        fitted_arc_line_positions, arc_pixel_locations, arc_wavelengths,
        valid_fibers, dispersion_axis
    )

    num_good_lines = np.sum(residuals < 0.2)
    logger.info("Wavelength calibration complete using %d out of %d arc lines",
                num_good_lines, len(residuals))

    return wavelength_solution, residuals, fitted_arc_line_positions, arc_pixel_locations


def _compute_arc_positions(spectra, trace_positions, valid_fibers, reference_index,
                            arc_pixel_guesses, arc_wavelengths, use_kernel, peak_threshold):
    """
    Detects arc line pixel positions for each fiber, starting from a reference and spreading out.

    Args:
        spectra (ndarray): 2D array of flux values; each row is a fiber spectrum.
        trace_positions (ndarray): 2D array with trace column positions for each fiber.
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

    # Backward pass
    last_fit = arc_pixel_guesses.copy()
    for i in range(reference_index - 1, -1, -1):
        if not valid_fibers[i]:
            continue
        _, cont = identify_sky_pixels(spectra[i], per=5)
        flux_corrected = spectra[i] - cont
        last_fit = _get_arclines_fiber(flux_corrected, last_fit,
                                      limit=peak_threshold, use_kernel=use_kernel)
        arc_pixels[i] = last_fit
        logger.debug("Backward fiber %d fitted", i)

    # Forward pass
    last_fit = arc_pixel_guesses.copy()
    for i in range(reference_index + 1, num_fibers):
        if not valid_fibers[i]:
            continue
        _, cont = identify_sky_pixels(spectra[i])
        flux_corrected = spectra[i] - cont
        last_fit = _get_arclines_fiber(flux_corrected, last_fit,
                                      limit=peak_threshold, use_kernel=use_kernel)
        arc_pixels[i] = last_fit
        logger.debug("Forward fiber %d fitted", i)

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
    peaks = y1[loc]

    # Helper function to refine peak positions using quadratic interpolation
    def get_trace_chunk(flat, XN):
        YM = np.arange(flat.shape[0])
        inds = np.zeros((3, len(XN)))
        inds[0] = XN - 1
        inds[1] = XN
        inds[2] = XN + 1
        inds = inds.astype(int)

        # Quadratic interpolation to refine peak positions
        Trace = (YM[inds[1]] - (flat[inds[2]] - flat[inds[0]]) /
                 (2. * (flat[inds[2]] - 2. * flat[inds[1]] + flat[inds[0]])))
        return Trace

    # Refine peak positions using quadratic interpolation
    loc = get_trace_chunk(y1, loc)

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


def get_trace(twilight, trace_rows, exclude_fibers):
    """
    Extract fiber traces from a twilight flat field image.

    Args:
        twilight (np.ndarray): 2D array representing the twilight flat field image used to determine fiber locations.
        trace_rows (np.ndarray): 1D array of nominal y-positions for each fiber's trace in pixels.
        exclude_fibers (np.ndarray):  array float specifying fibers to exclude.

    Returns:
        trace (np.ndarray): 2D array of the calculated trace positions for each fiber across the image.
        good_fiber_mask (np.ndarray): 1D boolean array indicating which fibers are valid (not excluded).
        Trace (np.ndarray): Intermediate 2D array of trace positions before polynomial smoothing.
        xchunks (np.ndarray): 1D array of x-axis chunk centers used during trace extraction.
    """

    total_fibers = len(trace_rows)

    exclude_mask = np.array(exclude_fibers, dtype=bool)

    good = np.where(~exclude_mask)[0]
    N1 = len(good)

    def get_trace_chunk(flat, XN):
        YM = np.arange(flat.shape[0])
        inds = np.zeros((3, len(XN)))
        inds[0] = XN - 1.
        inds[1] = XN + 0.
        inds[2] = XN + 1.
        inds = np.array(inds, dtype=int)
        Trace = (YM[inds[1]] - (flat[inds[2]] - flat[inds[0]]) /
                 (2. * (flat[inds[2]] - 2. * flat[inds[1]] + flat[inds[0]])))
        return Trace

    image = twilight
    N = 80
    xchunks = np.array([np.mean(x)
                        for x in np.array_split(np.arange(image.shape[1]), N)])
    chunks = np.array_split(image, N, axis=1)

    # Remove first and last chunks (border cleanup)
    chunks = chunks[1:-1]
    xchunks = xchunks[1:-1]
    flats = [np.mean(chunk, axis=1) for chunk in chunks]

    Trace = np.zeros((total_fibers, len(chunks)))
    k = 0
    P = []

    for flat, x in zip(flats, xchunks):
        diff_array = flat[1:] - flat[:-1]
        loc = np.where((diff_array[:-1] > 0.) & (diff_array[1:] < 0.))[0]
        loc = loc[loc > 2]

        peaks = flat[loc + 1]
        loc = loc[peaks > 0.3 * np.median(peaks)] + 1
        P.append(loc)

        trace = get_trace_chunk(flat, loc)
        T = np.zeros(total_fibers)

        if len(trace) > N1:
            trace = trace[-N1:]

        if len(trace) == N1:
            T[good] = trace
            for missing in np.where(exclude_mask)[0]:
                gind = np.argmin(np.abs(missing - good))
                T[missing] = T[good[gind]] + trace_rows[missing] - trace_rows[good[gind]]

        if len(trace) == total_fibers:
            T = trace

        Trace[:, k] = T
        k += 1

    x = np.arange(twilight.shape[1])
    trace = np.zeros((Trace.shape[0], twilight.shape[1]))
    for i in range(Trace.shape[0]):
        sel = Trace[i, :] != 0.
        if np.any(sel):
            trace[i] = np.polyval(np.polyfit(xchunks[sel], Trace[i, sel], 7), x)

    good_fiber_mask = ~exclude_mask
    return trace, good_fiber_mask, Trace, xchunks
