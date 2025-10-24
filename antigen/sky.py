import numpy as np
from astropy.stats import sigma_clip, mad_std, biweight_location as biweight
from scipy.interpolate import interp1d
from scipy.ndimage import percentile_filter
from sklearn.decomposition import PCA

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

    mask = sigma_clip(sky - cont, masked=True, maxiters=None, stdfunc=mad_std, sigma=5)

    # Return the mask (True for sky pixels) and the filtered continuum
    return mask.mask, cont

def get_skymask(sky, per=50, size=50, niter=3):
    """
    Iteratively identifies and masks sky pixels in an input array by applying
    a percentile filter and sigma-clipping.

    Parameters:
        sky (array-like): Input sky intensity values.
        per (int, optional): Percentile value for the filter. Default is 50 (median).
        size (int, optional): Size of the filter window. Default is 50.
        niter (int, optional): Number of iterations for refining the sky mask. Default is 3.

    Returns:
        tuple: A boolean mask array indicating sky pixels and the final filtered continuum array.
    """
    # Keep a copy of the original sky data for final comparison
    sky_orig = sky * 1.0

    # Iteratively refine the sky mask
    for i in np.arange(niter):
        # Apply a percentile filter to estimate the continuum
        cont = percentile_filter(sky, per, size=size)

        try:
            # Apply sigma-clipping to identify sky pixels using robust statistics
            mask = sigma_clip(sky - cont, masked=True, maxiters=None,
                              stdfunc=mad_std, sigma=5, sigma_lower=500)
        except:
            # Fallback for older versions of sigma_clip
            mask = sigma_clip(sky - cont, iters=None, stdfunc=mad_std,
                              sigma=5, sigma_lower=500)

        # Update the sky values for masked pixels using the continuum
        sky[mask.mask] = cont[mask.mask]

    # Perform a final sigma-clipping pass using the original sky data
    try:
        mask = sigma_clip(sky_orig - cont, masked=True, maxiters=None,
                          stdfunc=mad_std, sigma=5, sigma_lower=500)
    except:
        mask = sigma_clip(sky_orig - cont, iters=None, stdfunc=mad_std,
                          sigma=5, sigma_lower=500)

    # Return the final mask and continuum array
    return mask.mask, cont


def subtract_sky(spectra, good):
    """
    Subtracts sky contributions from input spectra using biweight calculations.

    This function processes the input spectra to calculate and subtract the sky
    contribution, ensuring accurate data for further analysis. It identifies sky
    pixels, computes the biweighted sky spectrum, and subtracts the sky spectrum
    from the original input spectra.

    Args:
        spectra: ndarray
            A 2D array where rows represent fibers and columns represent wavelength
            bins.
        good: ndarray
            A 1D boolean array indicating whether a fiber is good (True) or not
            (False).

    Returns:
        sky_subtracted_frame: The sky-subtracted spectra.
        init_sky: The biweighted sky spectrum for each fiber.
    """

    # Get the number of fibers and number of wavelength bins
    nfibs, N = spectra.shape

    # Define range for biweight calculation (middle third of the data)
    n1 = int(1. / 3. * N)
    n2 = int(2. / 3. * N)

    # Calculate the biweight of spectra over the middle third of each fiber's data
    biweighted_spectrum = biweight(spectra[:, n1:n2], axis=1, ignore_nan=True)

    # Identify sky pixels based on the biweighted data and apply a mask
    mask, cont = identify_sky_pixels(biweighted_spectrum[good], size=15)

    # Create a mask for fibers that are not good and are sky fibers
    m1 = ~good
    m1[good] = mask
    skyfibers = ~m1

    # Compute the biweighted sky spectrum based on sky fibers
    init_sky = biweight(spectra[skyfibers], axis=0, ignore_nan=True)

    # Subtract the sky spectrum from the original spectra

    sky_subtracted_frame = spectra - init_sky[np.newaxis, :]
    return sky_subtracted_frame, init_sky


def get_pca_sky_residuals(data, ncomponents=config.DEFAULT_PCA_COMPONENTS_SKY):
    """
    Perform PCA on the input data to extract the principal components and
    reconstruct the data using a specified number of components.

    Parameters
    ----------
    data : 2d numpy array
        Input data where rows represent samples (e.g., spectra) and columns represent features (e.g., wavelengths).
    ncomponents : int, optional
        The number of principal components to retain for reconstruction. Default is 5.

    Returns
    -------
    pca : PCA object
        The fitted PCA model.
    A : 2d numpy array
        The reconstructed data using the first `ncomponents` principal components.
    """
    # TODO: fix magic number default for components, use CONFIG params?

    # Initialize PCA with the specified number of components
    pca = PCA(n_components=ncomponents)

    # Fit the PCA model and transform the data into the principal components space
    H = pca.fit_transform(data)

    # Reconstruct the data using the principal components
    A = np.dot(H, pca.components_)

    # Return the PCA model and the reconstructed data
    return pca, A


def get_residual_map(data, pca):
    """
    Compute the residual map by subtracting a PCA-based model from the input data.

    Parameters
    ----------
    data : 2d numpy array
        Input data where rows represent samples (e.g., spectra) and columns represent features.
    pca : PCA object
        Fitted PCA model used for reconstruction of the data.

    Returns
    -------
    res : 2d numpy array
        The residual map, which is the difference between the input data and the PCA model.
    """

    # Initialize the residual map with zeros without propagating NaNs
    # Avoid data * 0. which can emit RuntimeWarning when data contains NaNs
    res = np.zeros_like(data)

    # Loop over each column (feature) in the input data
    for i in np.arange(data.shape[1]):

        # Compute the absolute deviation from the median for each column
        A = np.abs(data[:, i] - np.nanmedian(data[:, i]))

        # Identify the "good" data points (those within 3 times the median deviation)
        good = A < (3. * np.nanmedian(A))

        # Select finite values and good points for the model fitting
        sel = np.isfinite(data[:, i]) * good

        # Compute the PCA coefficients for the selected data points
        coeff = np.dot(data[sel, i], pca.components_.T[sel])

        # Reconstruct the model based on the computed coefficients
        model = np.dot(coeff, pca.components_)

        # Store the residual (difference) for this feature in the residual map
        res[:, i] = model

    # Return the computed residual map
    return res


def get_arc_pca(arcskysub, good, mask, components=config.DEFAULT_PCA_COMPONENTS_ARC):
    """
    Perform PCA on the arc-sky-subtracted data with preprocessing to remove
    bad data points and mask the non-relevant pixels.

    Parameters
    ----------
    arcskysub : 2d numpy array
        The arc-sky-subtracted data to be processed.
    good : 1d numpy array of bools
        Boolean mask indicating which fibers are good (non-bad).
    mask : 1d numpy array of bools
        Mask indicating the relevant pixels (e.g., fiber locations).
    components : int, optional
        Number of PCA components to retain. Default is 15.

    Returns
    -------
    pca : PCA object
        The fitted PCA model.
    """
    # TODO: fix magic number default for components, use CONFIG params?

    # Initialize X as the arc-sky-subtracted data
    X = arcskysub

    # Set values outside the mask to 0 and apply the "good" fiber mask
    X[:, ~mask] = 0.
    X[~good] = 0.

    # Transpose X to have samples in rows and features in columns
    X = X.swapaxes(0, 1)

    # Calculate the mean and standard deviation along the rows (fiber axis)
    M = np.nanmean(X, axis=1)
    Var = np.nanstd(X, axis=1)

    # Prevent division by zero by setting variance of zero to 1
    Var[Var == 0.] = 1.

    # Normalize the data by subtracting the mean and dividing by the variance
    X = (X - M[:, np.newaxis]) / Var[:, np.newaxis]

    # Replace NaN values in the normalized data with 0
    X[np.isnan(X)] = 0.

    # Perform PCA on the normalized data
    pca, A = get_pca_sky_residuals(X, ncomponents=components)

    # Return the fitted PCA model
    return pca


def get_continuum(skysub, masksky, nbins=config.DEFAULT_SKY_CONTINUUM_BINS):
    """
    Compute the continuum by interpolating the biweighted median of
    sky-subtracted data, with masked values excluded.

    Parameters
    ----------
    skysub : 2d numpy array
        Sky-subtracted data where each row represents a spectrum.
    masksky : 1d numpy array of bools
        Mask indicating which sky pixels should be excluded in the continuum calculation.
    nbins : int, optional
        Number of bins to divide the spectrum into for the biweight calculation. Default is 50.

    Returns
    -------
    bigcont : 2d numpy array
        The continuum for each spectrum in the sky-subtracted data.
    """
    # TODO: fix magic number default for nbins, use CONFIG params?

    # Initialize the output array for the continuum with zeros without propagating NaNs
    # Avoid skysub * 0. which can emit RuntimeWarning when skysub contains NaNs
    bigcont = np.zeros_like(skysub)

    # Loop over each row (spectrum) in the sky-subtracted data
    for j in np.arange(skysub.shape[0]):
        # Copy the current spectrum and mask the sky pixels
        y = skysub[j] * 1.
        y[masksky] = np.nan

        # Divide the spectrum into bins and calculate the mean of each bin
        x = np.array([np.mean(chunk) for chunk in np.array_split(np.arange(len(y)), nbins)])

        # Calculate the biweighted median for each bin
        # Prefilter NaNs/Infs to avoid RuntimeWarnings inside astropy stats
        z = np.array([
            (biweight(chunk[np.isfinite(chunk)], ignore_nan=True)
             if np.isfinite(chunk).sum() > 0 else np.nan)
            for chunk in np.array_split(y, nbins)
        ])

        # Select bins with finite values for interpolation
        sel = np.isfinite(z)

        # If there are enough valid bins, perform quadratic interpolation
        if sel.sum() > 5:
            I = interp1d(x[sel], z[sel], kind='quadratic', bounds_error=False,
                         fill_value='extrapolate')

            # Store the interpolated continuum for the current spectrum
            bigcont[j] = I(np.arange(len(y)))

    # Return the computed continuum for all spectra
    return bigcont

def advanced_sky_subtraction(spectra, good_fiber_mask, pca=None, 
                            skymask_size=25, apply_pca_residuals=True,
                            continuum_bins=50):
    """
    Perform advanced sky subtraction with optional PCA residual correction.
    
    This function provides a comprehensive sky subtraction pipeline that can
    optionally apply PCA-based residual correction for improved sky removal.
    
    Args:
        spectra (numpy.ndarray): Input spectra (fibers x wavelength)
        good_fiber_mask (numpy.ndarray): Boolean mask of good (non-sky) fibers
        pca (sklearn.decomposition.PCA, optional): Pre-fitted PCA model for
            residual correction. If None, no PCA correction is applied.
        skymask_size (int, optional): Size for sky mask generation. Defaults to 25.
        apply_pca_residuals (bool, optional): Whether to apply PCA residual
            correction. Defaults to True.
        continuum_bins (int, optional): Number of bins for continuum estimation.
            Defaults to 50.
            
    Returns:
        numpy.ndarray or tuple: Sky-subtracted spectra. If return_intermediates
            is True, returns tuple of (sky_subtracted, residuals, continuum).
            
    Example:
        >>> # Basic sky subtraction
        >>> sky_sub = sky.advanced_sky_subtraction(spectra, good_fibers)
        >>> 
        >>> # Advanced with PCA correction
        >>> sky_sub = sky.advanced_sky_subtraction(
        ...     spectra, good_fibers, pca=pca_model, skymask=mask
        ... )
        >>> 
        >>> # With diagnostics
        >>> sky_sub, residuals, cont = sky.advanced_sky_subtraction(
        ...     spectra, good_fibers, pca=pca_model, skymask=mask,
        ...     return_intermediates=True
        ... )
    """
    # Perform basic sky subtraction
    sky_subtracted_basic, biweighted_spectrum = subtract_sky(spectra, good_fiber_mask)
    
    residuals = None
    continuum = None
    sky_subtracted_advanced = sky_subtracted_basic
    if apply_pca_residuals and pca is not None:
        # Expand sky mask for better continuum estimation
        skymask, continuum = get_skymask(biweighted_spectrum, size=skymask_size)
        expanded_skymask = skymask.copy()
        expanded_skymask[1:] += skymask[:-1]  # Add previous pixel
        expanded_skymask[:-1] += skymask[1:]  # Add next pixel
        
        # Estimate continuum
        continuum = get_continuum(sky_subtracted_basic, expanded_skymask,
                                nbins=continuum_bins)
        
        # Calculate residuals after continuum subtraction
        residual_input = sky_subtracted_basic - continuum
        residuals = get_residual_map(residual_input, pca)
        
        # Mask residuals where sky mask is invalid
        residuals[:, ~expanded_skymask] = 0.0
        
        # Apply residual correction
        sky_subtracted_advanced = sky_subtracted_basic- residuals
    
    return sky_subtracted_advanced, sky_subtracted_basic, residuals, continuum


def generate_pca_from_arc_spectra(arc_spectra, good_fiber_mask, ftf_correction, 
                                  skymask_size=25, pca_components=15):
    """
    Generate PCA model from rectified arc lamp spectra for sky subtraction.
    
    This function processes arc spectra through the same sky processing pipeline
    used for science data to create a PCA model that can be applied for 
    residual sky subtraction in science observations.
    
    Args:
        arc_spectra (numpy.ndarray): Rectified arc lamp spectra (fibers x wavelength)
        good_fiber_mask (numpy.ndarray): Boolean mask of good fibers
        ftf_correction (numpy.ndarray): Fiber-to-fiber correction array
        skymask_size (int, optional): Size for sky mask generation. Defaults to 25.
        pca_components (int, optional): Number of PCA components. Defaults to 15.
            
    Returns:
        tuple: A tuple containing:
            - pca (sklearn.decomposition.PCA): Fitted PCA model
            - skymask (numpy.ndarray): Sky mask used for processing
            - continuum (numpy.ndarray): Continuum model
            - biweighted_spectrum (numpy.ndarray): Biweighted spectrum
    """
    # Apply fiber-to-fiber correction to arc spectra
    arc_corrected = arc_spectra / ftf_correction

    # Generate biweighted spectrum for sky mask creation
    biweighted_spectrum = biweight(arc_corrected, axis=0, ignore_nan=True)

    # Generate sky mask and continuum
    skymask, continuum = get_skymask(biweighted_spectrum, size=skymask_size)

    # Perform sky subtraction on arc spectra  
    arc_sky_subtracted, arc_biweight_spectrum = subtract_sky(arc_corrected, good_fiber_mask)
        
    pca = get_arc_pca(arc_sky_subtracted, good_fiber_mask, skymask, 
                      components=pca_components)
    
    return pca, skymask, continuum, biweighted_spectrum