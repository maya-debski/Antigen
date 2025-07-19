import numpy as np
from astropy.convolution import Gaussian1DKernel, convolve
from astropy.io import fits
from astropy.stats import sigma_clip, mad_std, biweight_location as biweight
from scipy.interpolate import interp1d
from scipy.ndimage import percentile_filter

from antigen import config


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


def rectify(scispectra, errspectra, wave_all, def_wave):
    """
    Rectifies scientific and error spectra by interpolating them onto a common wavelength grid.

    Parameters:
        scispectra (2D array): Array of scientific spectra to be rectified.
        errspectra (2D array): Corresponding error spectra for each scientific spectrum.
        wave_all (2D array): Wavelength grids corresponding to each input spectrum.
        def_wave (1D array): Target wavelength grid for interpolation.

    Returns:
        tuple:
            - scirect (2D array): Rectified scientific spectra on the target wavelength grid.
            - errorrect (2D array): Rectified error spectra on the target wavelength grid.
    """
    # Initialize arrays to store rectified scientific spectra and errors
    scirect = np.zeros((scispectra.shape[0], len(def_wave)))
    errorrect = np.zeros((scispectra.shape[0], len(def_wave)))

    # Placeholder for wavelength grid indices (used in optional refinement)
    indices1 = np.ones(scirect.shape, dtype=int)

    # Loop through each spectrum to interpolate onto the target wavelength grid
    for i in np.arange(scispectra.shape[0]):
        # Compute wavelength bin sizes for flux normalization
        dw = np.diff(wave_all[i])
        dw = np.hstack([dw[0], dw])  # Ensure length matches the wavelength array

        # Interpolate the scientific spectrum, normalizing by wavelength bin size
        scirect[i] = np.interp(def_wave, wave_all[i], scispectra[i] / dw,
                               left=np.nan, right=np.nan)

        # Interpolate the error spectrum, normalizing by wavelength bin size
        errorrect[i] = np.interp(def_wave, wave_all[i], errspectra[i] / dw,
                                 left=np.nan, right=np.nan)

        # Store indices for possible further refinement (optional block below)
        indices1[i] = np.searchsorted(wave_all[i], def_wave) + i * 1032

    # Optional: Weighted error interpolation for more accurate error propagation
    # Uncomment the following block if needed for error refinement:
    #
    # x_var = (def_wave[np.newaxis, :] * np.ones((scirect.shape[0], 1))).ravel()
    # x_fix = wave_all.ravel()
    # indices1 = indices1.ravel()
    # indices2 = indices1 - 1
    # indices2[indices2 < 0] = 0
    # indices1[indices1 >= len(x_fix)] = len(x_fix) - 1
    #
    # distances1 = np.abs(x_fix[indices1] - x_var)
    # distances2 = np.abs(x_fix[indices2] - x_var)
    # total_distance = distances1 + distances2
    # weight1 = distances1 / total_distance
    # weight2 = distances2 / total_distance
    # errorrect = (weight2**2 * errspectra.ravel()[indices1]**2 +
    #              weight1**2 * errspectra.ravel()[indices2]**2)
    # errorrect = np.sqrt(errorrect)
    # errorrect = np.reshape(errorrect, scirect.shape)
    # errorrect[np.isnan(scirect)] = np.nan

    # Return the rectified scientific and error spectra
    return scirect, errorrect


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
        tuple:
            - initial_ftf (2D array): Initial fiber-to-fiber correction factors.
            - ftf (2D array): Smoothed fiber-to-fiber correction factors.
    """
    # Compute the average spectrum across all fibers using a robust biweight statistic
    average = biweight(spectrum, axis=0, ignore_nan=True)

    # Calculate the initial fiber-to-fiber correction by dividing each fiber by the average spectrum
    initial_ftf = spectrum / average[np.newaxis, :]

    # Create a wavelength grid and divide it into chunks for smoothing
    X = np.arange(spectrum.shape[1])
    x = np.array([np.mean(chunk) for chunk in np.array_split(X, n_chunks)])

    # Initialize the smoothed correction array
    ftf = spectrum * 0.

    # Loop through each fiber to compute the smoothed correction factor
    for i in np.arange(len(spectrum)):
        # Compute the biweight statistic for each chunk of the initial correction factor
        y = np.array([biweight(chunk, ignore_nan=True) for chunk in np.array_split(initial_ftf[i], n_chunks)])

        # Select valid (finite) values for interpolation
        sel = np.isfinite(y)
        if sel.sum() == 0.:
            continue
        # Interpolate the correction factor using quadratic interpolation
        I = interp1d(x[sel], y[sel], kind='quadratic', bounds_error=False, fill_value='extrapolate')

        # Apply the interpolation to the full wavelength range
        ftf[i] = I(X)

    # Return both the initial and smoothed fiber-to-fiber correction factors
    return initial_ftf, ftf


def get_wavelength(spectrum, trace, good, xref, lines, use_kernel=True, limit=100, fiberref=config.VIRUS2_FIBER_REF):
    """
    Computes the wavelength solution for each fiber in a spectrograph based on trace and spectral data.

    Args:
        spectrum (ndarray): 2D array of spectra, each row corresponding to a fiber.
        trace (ndarray): 2D array with trace positions for each fiber.
        good (ndarray): Boolean array indicating which fibers have valid data.
        use_kernel (bool): Whether to apply kernel smoothing when identifying arc lines. Default is True.
        limit (float): Limit on how far to search for matching arc lines. Default is 100.
        fiberref (int): default = config.VIRUS2_FIBER_REF=130

    Returns:
        tuple: (wavelength, res, X, W)
            - wavelength (ndarray): Wavelength solution for each fiber.
            - res (ndarray): Residuals from the biweight mean calculation.
            - X (ndarray): Adjusted positions in trace space for arc lines.
            - W (ndarray): Arc line positions for each fiber.
    """

    # Initialize wavelength array and starting fiber position
    init_fiber = fiberref
    wavelength = np.zeros_like(spectrum)
    loc = xref.copy()

    # W will store arc line positions for each fiber
    W = np.zeros((trace.shape[0], len(lines)))
    mask, cont = identify_sky_pixels(spectrum[init_fiber], per=5)  # Identify sky lines
    y = spectrum[init_fiber] - cont  # Subtract continuum
    W[init_fiber] = get_arclines_fiber(y, loc, limit=limit, use_kernel=use_kernel)

    # Process fibers before the reference fiber in reverse order
    for i in np.arange(init_fiber)[::-1]:
        mask, cont = identify_sky_pixels(spectrum[i], per=5)  # Identify sky lines
        y = spectrum[i] - cont  # Subtract continuum
        if good[i]:  # Only process if the fiber is marked as good
            loc = get_arclines_fiber(y, loc, limit=limit,
                                     use_kernel=use_kernel)
            W[i] = loc

    # Reset location and process fibers after the reference fiber
    loc = xref.copy()
    for i in np.arange(init_fiber + 1, spectrum.shape[0]):
        mask, cont = identify_sky_pixels(spectrum[i])
        y = spectrum[i] - cont
        if good[i]:
            loc = get_arclines_fiber(y, loc, limit=limit,
                                     use_kernel=use_kernel)
            W[i] = loc

    # Initialize X (adjusted trace positions) and residuals array
    X = np.zeros_like(W)
    xall = np.arange(trace.shape[1])
    res = np.zeros(W.shape[1])

    # Interpolate missing values and fit polynomial to each arc line
    for i in range(W.shape[1]):
        x = np.zeros(W.shape[0])
        bad = np.where(~good)[0]
        gind = np.where(good)[0]

        # Fill in missing arc line positions for bad fibers
        for b in bad:
            W[b, i] = W[gind[np.argmin(np.abs(b - gind))], i]

        # Interpolate positions in trace space
        for j in range(W.shape[0]):
            x[j] = np.interp(W[j, i], xall, trace[j])

        # Fit a 4th-order polynomial to arc line positions
        sel = (W[:, i] > 0) * np.isfinite(x)

        X[:, i] = np.polyval(np.polyfit(x[sel], W[sel, i], 4), x)

        # Compute residuals using biweight mean
        res[i] = mad_std(X[:, i] - W[:, i], ignore_nan=True)
    fits.HDUList([fits.PrimaryHDU(X), fits.ImageHDU(lines)]).writeto('test.fits', overwrite=True)
    goodlines = res < 0.2
    # Compute final wavelength solution for each fiber
    for j in range(W.shape[0]):
        if good[j]:
            wavelength[j] = np.polyval(np.polyfit(X[j][goodlines],
                                                  lines[goodlines], 5), xall)

    return wavelength, res, X, W


def get_arclines_fiber(spectrum, init_loc=None, limit=100, use_kernel=True):
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
    return spec / npix


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
    return np.sqrt(spec) / npix


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


def get_trace(twilight, ref):
    """
    Extract fiber traces from a twilight flat field image.

    Parameters
    ----------
    twilight : 2d numpy array
        Twilight flat field image used to determine fiber locations.

    Returns
    -------
    trace : 2d numpy array
        The calculated trace positions for each fiber across the image.
    good : 1d numpy array (boolean)
        Boolean mask indicating which fibers are valid (non-missing).
    """

    # Determine the number of valid (good) fibers
    N1 = np.isfinite(ref['px']).sum()
    good = np.where(np.isfinite(ref['px']))[0]  # Indices of good fibers

    # Helper function to calculate trace positions for a chunk of the image
    def get_trace_chunk(flat, XN):
        # YM represents the y-axis pixel coordinates
        YM = np.arange(flat.shape[0])

        # Create a 3-row array for XN-1, XN, and XN+1 indices
        inds = np.zeros((3, len(XN)))
        inds[0] = XN - 1.
        inds[1] = XN + 0.
        inds[2] = XN + 1.
        inds = np.array(inds, dtype=int)

        # Calculate the trace using a second-order derivative method
        Trace = (YM[inds[1]] - (flat[inds[2]] - flat[inds[0]]) /
                 (2. * (flat[inds[2]] - 2. * flat[inds[1]] + flat[inds[0]])))
        return Trace

    # Assign the input image to a variable
    image = twilight

    # Determine the number of chunks based on whether the image is binned
    N = 80

    # Split the x-axis into chunks and calculate the mean x-position for each chunk
    xchunks = np.array([np.mean(x)
                        for x in np.array_split(np.arange(image.shape[1]), N)])

    # Split the image into vertical chunks
    chunks = np.array_split(image, N, axis=1)

    # Calculate the mean flat field for each chunk
    flats = [np.mean(chunk, axis=1) for chunk in chunks]

    # Initialize an array to hold the trace positions for each fiber
    Trace = np.zeros((len(ref), len(chunks)))

    # Initialize a counter and a list to store peak positions
    k = 0
    P = []

    # Iterate over each chunk to calculate the fiber traces
    for flat, x in zip(flats, xchunks):
        # Calculate the difference between adjacent pixels
        diff_array = flat[1:] - flat[:-1]

        # Identify peaks by finding zero-crossings in the difference array
        loc = np.where((diff_array[:-1] > 0.) & (diff_array[1:] < 0.))[0]
        loc = loc[loc > 2]  # Ignore peaks near the image edges

        # Filter out weak peaks
        peaks = flat[loc + 1]
        loc = loc[peaks > 0.3 * np.median(peaks)] + 1

        # Store the detected peak positions
        P.append(loc)

        # Get the trace positions for the detected peaks
        trace = get_trace_chunk(flat, loc)

        # Initialize an array to hold the trace for this chunk
        T = np.zeros((len(ref)))

        # If the number of detected peaks exceeds the number of good fibers, trim the excess
        if len(trace) > N1:
            trace = trace[-N1:]

        # If the number of detected peaks matches the number of good fibers
        if len(trace) == N1:
            T[good] = trace
            # Interpolate missing fibers based on nearby good fibers
            for missing in np.where(np.isnan(ref['px']))[0]:
                gind = np.argmin(np.abs(missing - good))
                T[missing] = T[good[gind]] + ref['px'][missing] - ref['px'][good[gind]]

        # If the number of detected peaks matches the total number of fibers
        if len(trace) == len(ref):
            T = trace

        # Store the calculated trace for this chunk
        Trace[:, k] = T
        k += 1

    # Fit a 7th-order polynomial to smooth the traces across the x-axis
    x = np.arange(twilight.shape[1])
    trace = np.zeros((Trace.shape[0], twilight.shape[1]))
    for i in np.arange(Trace.shape[0]):
        sel = Trace[i, :] != 0.
        trace[i] = np.polyval(np.polyfit(xchunks[sel], Trace[i, sel], 7), x)

    # Return the final trace array and the good fiber mask
    return trace, np.isfinite(ref['px']), Trace, xchunks
