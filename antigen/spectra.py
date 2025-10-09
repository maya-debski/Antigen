import numpy as np
import logging
from scipy.sparse import csr_matrix


logger = logging.getLogger('antigen.spectra')

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

    # Calculate old wavelength widths using gradient for each fiber
    old_widths = np.gradient(wave_all, axis=1)

    # Convert e- per pixel to e-/A before resampling
    scispectra /= old_widths
    errspectra /= old_widths

    # Calculate new wavelength widths from def_wave
    new_widths = np.gradient(def_wave)

    # Convert error to variance for proper propagation
    variance = errspectra**2
    
    # Use resample_flux for flux-conserving resampling on full 2D arrays
    # For 2D input, we need to pass the wavelength centers for each fiber
    # Since wave_all is 2D and def_wave is 1D, we need to handle this properly
    
    # Initialize output arrays
    scirect = np.zeros((scispectra.shape[0], len(def_wave)))
    errorrect = np.zeros((scispectra.shape[0], len(def_wave)))
    
    # Process each fiber (resample_flux expects consistent old_centers for 2D)
    # We'll use the approach where each fiber has its own wavelength solution
    for i in range(scispectra.shape[0]):
        flux_out, var_out = resample_flux(
            flux=scispectra[i],
            old_centers=wave_all[i],
            old_widths=old_widths[i],
            new_centers=def_wave,
            new_widths=new_widths,
            variance=variance[i],
            return_density=True
        )
        
        scirect[i] = flux_out
        errorrect[i] = np.sqrt(var_out)

    return scirect, errorrect

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

    # Return the root mean square error normalized
    return np.sqrt(spec)

def extract_fiber_spectra(array_flt, array_trace, array_error=None, npix=5):
    """Extract per-fiber spectra and corresponding uncertainties.

    This function wraps two lower-level routines, ``get_spectra`` and
    ``get_spectra_error``. If ``array_error`` is not provided, a synthetic
    per-pixel uncertainty image is constructed assuming Poisson statistics,
    i.e., :math:`\\sigma = \\sqrt{\\max(\\mathrm{counts}, 0)}`.

    Args:
        array_flt (np.ndarray): 2D image of counts to extract from. Assumed to be
            in units where Poisson variance equals the value (e.g., electrons).
        array_trace (np.ndarray): Trace description consumed by ``get_spectra`` and
            ``get_spectra_error`` (e.g., fiber center positions vs. pixel row/col).
        array_error (np.ndarray | None, optional): 2D per-pixel uncertainty image.
            If ``None``, synthetic Poisson errors are computed from ``array_flt``.
            Defaults to ``None``.
        npix (int, optional): Half-width or aperture parameter passed through to the
            underlying extraction functions. Defaults to ``5``.

    Returns:
        spectra: 2D array with shape (nfiber, nwave) (or as defined by
          ``get_spectra``) containing extracted fluxes.
        spectra_err: 2D array of the same shape as ``spectra`` containing
          propagated uncertainties.

    Notes:
        - Negative values in ``array_flt`` are clipped to zero when computing
          synthetic Poisson errors to avoid invalid square roots.
        - This function assumes that any gain conversion has already been applied
          to ``array_flt`` if needed.

    See Also:
        get_spectra, get_spectra_error
    """
    # 1) Flux extraction
    spectra = get_spectra(array_flt, array_trace, npix=npix)

    # 2) Error image (use provided, or synthesize from Poisson statistics)
    if array_error is None:
        # Clip negatives to 0 to avoid NaNs from sqrt on background-subtracted images
        poisson_var = np.clip(np.asarray(array_flt, dtype=np.float32), a_min=0.0, a_max=None)
        error_image = np.sqrt(poisson_var)
    else:
        error_image = np.asarray(array_error, dtype=np.float32)

    # 3) Error extraction, mirroring the flux path
    spectra_err = get_spectra_error(error_image, array_trace, npix=npix)

    return spectra, spectra_err


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


def get_optimal_spectrum(data, error, weights):
    """
    Extract a weighted 1D spectrum from a fiber data

    Args:
        data (ndarray): 2D array of fiber spectra within aperture
        error (ndarray): 2D array of fiber spectral errors
        weights (ndarray): 2D array of fiber weights

    Returns:
        spectrum (ndarray): Weighted 1D spectrum.
        spectrum_error (ndarray): Weighted 1D error spectrum.
    """

    var = error ** 2
    spectrum = (np.nansum(data / var * weights, axis=0) /
                np.nansum(weights ** 2 / var, axis=0))
    spectrum_error = np.sqrt(np.nansum(weights, axis=0) /
                             np.nansum(weights ** 2 / var, axis=0))

    return spectrum, spectrum_error


def convert_spectral_units(spectrum, error, wave, exposure_time, area=5.81e4):
    """
    Incoming units are e- / A

    Outgoing units are ergs / S / A/ cm^2

    Args:
        spectrum (ndarray): 2D array of fiber spectra (units: e- / A)
        error (ndarray): 2D array of fiber spectral errors (units: e- / A)
        wave (ndarray): 1D array of fiber wavelengths
        exposure_time (ndarray): 1D array of fiber exposure times
        area (float): Area of telescope mirror in cm^2

    Returns:
        spectrum (array): 2D array of fiber in (units: ergs / S / A/ cm^2)
    """
    flux_factor = (6.626e-27 * 2.99792e18 / wave)
    total_factor = flux_factor / area / exposure_time
    converted_flux = spectrum * total_factor
    converted_error = error * total_factor
    return converted_flux, converted_error

def correct_fiber_to_fiber(rectified_spectra, rectified_error, rectified_ftf):
    """
    Corrects fiber-to-fiber variations in spectroscopic data.

    Args:
        rectified_spectra: The array representing the rectified spectra before fiber-to-fiber
            correction is applied.
        rectified_error: The array representing the errors associated with the rectified
            spectra before fiber-to-fiber correction.
        rectified_ftf: The array representing the fiber-to-fiber response function values
            that will be used for the correction.

    Returns:
        corrected_spectra (np.ndarray): The array representing the corrected spectra after
        corrected_error (np.ndarray): The array representing the errors associated with the corrected
    """
    corrected_spectra = rectified_spectra / rectified_ftf
    corrected_error = rectified_error / rectified_ftf
    return corrected_spectra, corrected_error


def edges_from_centers(centers, widths=None):
    """Compute bin edges from pixel centers and optional per-pixel widths.

    If widths is None, edges are midpoints between adjacent centers, with
    the first/last edges extrapolated.

    Args:
        centers (array-like): Monotonic array of wavelength centers, shape (N,).
        widths (array-like | None): Optional bin widths (same length as centers).

    Returns:
        np.ndarray: Bin edges, shape (N+1,). Strictly increasing.
    """
    c = np.asarray(centers, dtype=float)
    if widths is not None:
        w = np.asarray(widths, dtype=float)
        if w.shape != c.shape:
            raise ValueError("widths must match centers shape.")
        edges = np.empty(c.size + 1, dtype=float)
        edges[:-1] = c - 0.5 * w
        edges[1:]  = c + 0.5 * w
        # Guard against numerical non-monotonicity
        if not np.all(edges[1:] > edges[:-1]):
            raise ValueError("Computed edges are not strictly increasing.")
        return edges

    # Widths unknown: use midpoints
    mid = 0.5 * (c[1:] + c[:-1])
    left = c[0] - (mid[0] - c[0])
    right = c[-1] + (c[-1] - mid[-1])
    edges = np.concatenate([[left], mid, [right]])
    if not np.all(edges[1:] > edges[:-1]):
        raise ValueError("Centers must be strictly monotonic to infer edges.")
    return edges


def build_overlap_matrix(old_edges, new_edges, return_density=True):
    """Build a sparse overlap matrix that conserves flux between binnings.

    Each output bin receives the integral of input flux *over the overlap width*.
    If `return_density=True`, the integral is divided by the new bin width so
    the output is flux density again.

    Args:
        old_edges (array-like): Input bin edges, shape (N_old+1,).
        new_edges (array-like): Output bin edges, shape (N_new+1,), typically uniform.
        return_density (bool): If True, divide by new widths to return density.
                               If False, return integrated flux per output bin.

    Returns:
        csr_matrix: W with shape (N_new, N_old).
                    Output_flux = W @ input_flux   (for 1D arrays)
                    For 2D arrays F[n_spec, N_old]: F_out = F @ W.T
    """
    old_edges = np.asarray(old_edges, dtype=float)
    new_edges = np.asarray(new_edges, dtype=float)

    if not (np.all(old_edges[1:] > old_edges[:-1]) and np.all(new_edges[1:] > new_edges[:-1])):
        raise ValueError("Edges must be strictly increasing.")

    N_old = old_edges.size - 1
    N_new = new_edges.size - 1
    old_left, old_right = old_edges[:-1], old_edges[1:]
    new_left, new_right = new_edges[:-1], new_edges[1:]
    new_widths = new_right - new_left

    # For each new bin, find candidate old bins that can overlap
    # Use searchsorted to bracket quickly, then compute exact overlaps.
    # We accumulate COO pieces then build CSR.
    data, rows, cols = [], [], []

    # Precompute for speed
    old_min, old_max = old_left[0], old_right[-1]

    # Map from new bin index to old index range that might overlap
    # Start with where new_left/new_right would fall in old_edges
    start_idx = np.clip(np.searchsorted(old_right, new_left, side="right") - 1, 0, N_old-1)
    end_idx   = np.clip(np.searchsorted(old_left, new_right, side="left"),        0, N_old)

    for i in range(N_new):
        if new_right[i] <= old_min or new_left[i] >= old_max:
            continue  # completely outside range

        j0 = start_idx[i]
        j1 = end_idx[i]  # this is one past the last possible overlapping bin
        if j1 <= j0:
            # Edge case; still try computing overlap just with j0
            j_candidates = [j0]
        else:
            j_candidates = range(j0, j1)

        for j in j_candidates:
            # Overlap length between [new_left[i], new_right[i]] and [old_left[j], old_right[j]]
            overlap = max(0.0, min(new_right[i], old_right[j]) - max(new_left[i], old_left[j]))
            if overlap <= 0.0:
                continue

            w = overlap
            if return_density:
                # Convert integrated overlap to density by 1 / new bin width
                w = w / new_widths[i]

            data.append(w)
            rows.append(i)
            cols.append(j)

    W = csr_matrix((data, (rows, cols)), shape=(N_new, N_old))
    return W


def resample_flux(flux, old_centers=None, old_widths=None, old_edges=None,
    new_centers=None, new_widths=None, new_edges=None,
    variance=None, return_density=True):
    """Flux-conserving resample from irregular to uniform (or any) grid.

    You may pass centers+widths or explicit edges for old/new grids.
    If passing centers without widths, edges are inferred by midpoints.

    Args:
        flux (ndarray): Input flux. Shape (N_old,)
                        Interpreted as flux density if `return_density=True`,
                        otherwise as integrated flux per bin.
        old_centers (ndarray|None): Input bin centers, shape (N_old,).
        old_widths (ndarray|None): Input bin widths, shape (N_old,).
        old_edges (ndarray|None): Input edges, shape (N_old+1,). Overrides centers/widths if provided.
        new_centers (ndarray|None): Target bin centers, shape (N_new,).
        new_widths (ndarray|None): Target bin widths, shape (N_new,).
        new_edges (ndarray|None): Target edges, shape (N_new+1,). Overrides centers/widths if provided.
        variance (ndarray|None): Variance array matching `flux` shape.
        return_density (bool): If True, return flux density. If False, return integrated flux.

    Returns:
        (flux_out, var_out, W):
            flux_out (ndarray): Resampled flux, shape (N_new,) or (N_spec, N_new).
            var_out (ndarray|None): Resampled variance (same shape) if `variance` provided, else None.
            W (csr_matrix): Overlap matrix with shape (N_new, N_old).
    """
    # Build old/new edges
    if old_edges is None:
        if old_centers is None:
            raise ValueError("Provide old_edges or old_centers (with optional old_widths).")
        old_edges = edges_from_centers(old_centers, old_widths)
    if new_edges is None:
        if new_centers is None:
            raise ValueError("Provide new_edges or new_centers (with optional new_widths).")
        new_edges = edges_from_centers(new_centers, new_widths)

    # Build flux-conserving overlap matrix
    W = build_overlap_matrix(old_edges, new_edges, return_density=return_density)

    f = np.asarray(flux, dtype=float)
    if f.ndim == 1:
        f_out = W @ f
        v_out = None
        if variance is not None:
            var = np.asarray(variance, dtype=float)
            if var.shape != f.shape:
                raise ValueError("variance must match flux shape.")
            # Variance propagation: Var(Ax) = A (Var) A^T for diagonal Var -> sum(w^2 * var)
            # With sparse W, compute (W.power(2) @ var)
            v_out = (W.power(2) @ var)
        return f_out, v_out
    else:
        raise ValueError("flux must be 1D.")