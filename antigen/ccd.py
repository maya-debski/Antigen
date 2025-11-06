import logging

import numpy as np
from astropy.stats import biweight_location as biweight
from astropy.time import Time
from astropy.modeling import models, fitting
from astropy.stats import biweight_location, sigma_clip, mad_std

from antigen.io import load_fits
from antigen.plot import plot_scattered_light_fit

logger = logging.getLogger('antigen.ccd')


def prep_image(image, overscan_length, flip_x=False, flip_y=False, rotate=False, add_rows=0):
    """
    Measure the bias in the overscan region, trim the overscan and subtract the bias
    Perform any flips that are necessary

    Args:
        image (np.ndarray): FITS 2D image data array
        overscan_length (int): Length of overscan region in pixels
        flip_x (bool, optional): Flag to flip image in x direction. Defaults to False
        flip_y (bool, optional): Flag to flip image in y direction. Defaults to False
        rotate (bool, optional): Flag to rotate image 90 degrees. Defaults to False
        add_rows (int, optional): Number of rows to add as padding. Defaults to 0

    Returns:
        image (np.ndarray): Oriented fits 2D image data array, corrected for what amplifier it comes from.
    """
    if overscan_length == 0:
        # virus-w specific bias subtraction and image preparation
        # Split into two regions and  bias from each
        left_region = image[:, :1024] - biweight(image[-43:, :1024])
        right_region = image[:, 1124:] - biweight(image[-43:, 1124:])

        # Combine the two regions
        image = np.hstack([left_region, right_region])
    else:
        bias_value = biweight(image[:, -(overscan_length-2):])
        image = image[:, :-overscan_length] - bias_value
    if rotate:
        image = np.rot90(image)
    if flip_x:
        image = np.flip(image, axis=1)
    if flip_y:
        image = np.flip(image, axis=0)

    if add_rows:
        new_image = np.zeros((len(image) + add_rows, image.shape[1]))
        new_image[add_rows:, :] = image
        image = new_image

    if overscan_length == 0:
        image = image[:, 50:] # Trim the first 50 columns

    return image


def base_reduction(data, master_bias, overscan_length, gain, read_noise, flip_x, flip_y, rotate, add_rows,
                   header):
    """
    Performs basic reduction on raw astronomical image data. This includes preprocessing,
    bias subtraction, gain correction, and error estimation. The function takes into account
    overscan regions, flipping and rotation of the image, as well as the addition of padding rows.

    Args:
        data: The raw 2D image array that needs processing.
        master_bias: The 2D master bias array to subtract from the raw image.
        overscan_length: The length (in pixels) of the overscan region to crop before processing.
        gain: The gain factor to convert image pixel values from counts to electrons.
        read_noise: The read noise level (in electrons) for the image sensor.
        flip_x: A boolean to specify if the image needs flipping along the x-axis.
        flip_y: A boolean to specify if the image needs flipping along the y-axis.
        rotate: The integer number of 90-degree counter-clockwise rotations to apply.
        add_rows: The number of rows to add as padding to the processed image.
        header: The metadata associated with the image, usually in the form of a dictionary.

    Returns:
        image: The processed 2D image array with bias subtracted and gain applied.
        error_estimate: The estimated error for each pixel in the image.
        exposure_time: The exposure time of the image in seconds.

    Raises:
        ValueError: If any input dimensions mismatch or if invalid alignment parameters
        are provided.
        TypeError: If inputs are not of the expected data types.
    """
    # TODO: update docstring to match input args, function signature

    # Preprocess the raw image (e.g., background subtraction, padding)
    image = prep_image(data, overscan_length, flip_x, flip_y, rotate, add_rows)

    # Subtract the master bias from the image
    image[:] -= master_bias

    # Apply gain correction to convert counts to electrons
    image[:] *= gain

    # Calculate the error estimate (read noise + photon noise)
    error_estimate = np.sqrt(read_noise**2 + np.where(image > 0., image, 0.))

    # Return the reduced image and the error estimate

    try:
        exposure_time = header['EXPTIME']
    except KeyError:
        exposure_time = 0.0
    return image, error_estimate, exposure_time


def make_master_cal(filenames, frame_type, overscan_length, flip_x, flip_y, rotate, add_rows):
    """
    Load a list of calibration FITS files, apply basic CCD preprocessing, and create a master calibration frame.

    Args:
        filenames (list(str)): list of filenames
        frame_type (str): type of frame (e.g., arc, bias, flat)
        overscan_length (int): Length of overscan region in pixels
        flip_x (bool, optional): Flag to flip image in x direction. Defaults to False
        flip_y (bool, optional): Flag to flip image in y direction. Defaults to False
        rotate (bool, optional): Flag to rotate image 90 degrees. Defaults to False
        add_rows (int, optional): Number of rows to add as padding. Defaults to 0

    Returns:
        master_cal (np.ndarray): median stacked calibration frame
        master_cal_time (float): average MJD of the frames that were stacked
    """
    # Extract from the files, re-oriented by prep_image()
    # Extract observation times (MJD) for frames in the current chunk
    frames, times = ([], [])
    for filename in filenames:
        frame, header = load_fits(filename)
        prepped_frame = prep_image(frame, overscan_length, flip_x, flip_y, rotate, add_rows)
        frames.append(prepped_frame)
        mjd = Time(header['DATE-OBS']).mjd
        times.append(mjd)

    # Compute median frame and the mean time for the current chunk
    master_cal      = np.nanmedian(frames, axis=0)  # maybe biweight() as an alternate method
    master_cal_time = np.mean(times)
    return master_cal, master_cal_time, f"master_{frame_type}.fits"


def make_mask_for_trace(image, trace, fiber_profile_mask_size=7):
    """
    Creates a boolean mask to exclude regions near the trace from further processing.

    Args:
        image (np.ndarray): 2D numpy array image for which the mask is created.
        trace (np.ndarray): 2D numpy array describing the trace position as a function of fiber
        fiber_profile_mask_size (int, optional): Vertical size (in pixels) of the mask region
            above and below the trace. Default is 11.

    Returns:
        mask: 2D binary mask with the same shape as the input image. The region
        around the trace is set to 1 (or `wave` values, if provided), indicating
        where the mask covers; the rest is set to 0.
    """
    mask = np.zeros(image.shape, dtype=bool)
    columns = np.arange(image.shape[1])
    lower_bound = -int(fiber_profile_mask_size / 2)
    for fiber_trace in trace:
        for column in columns:
            if np.isnan(fiber_trace[column]):
                continue
            bottom = int(fiber_trace[column]) + lower_bound
            top = bottom + fiber_profile_mask_size
            mask[bottom:top, column] = True
    return mask


def get_scattered_light(image, trace, fiber_profile_mask_size=11, poly_order=2,
                       group_size=40, sigma_clip_threshold=3.0, outfolder=None):
    """
    Estimates scattered light background using grouped pixel averaging and 2D polynomial fitting.
    
    This function masks fiber trace regions, groups neighboring unmasked pixels into blocks,
    computes robust statistics for each group, then fits an Astropy Polynomial2D model
    to the resulting sparse but high-quality background measurements.
    
    Args:
        image (np.ndarray): 2D numpy image array.
        trace (np.ndarray): 2D numpy array describing trace position as function of fiber
        fiber_profile_mask_size (int, optional): Vertical size (in pixels) of mask region
            above and below traces. Default is 11.
        poly_order (int, optional): Order of 2D polynomial fit. Default is 3.
        group_size (int, optional): Size of pixel groups (group_size x group_size blocks).
            Default is 40.
        sigma_clip_threshold (float, optional): Sigma threshold for outlier rejection within
            each group. Default is 3.0.
        outfolder (str, optional): Output folder for diagnostic plots. If None, no plots saved.

    Returns:
        background_subtracted (np.ndarray): Image with scattered light subtracted
        scattered_light (np.ndarray): 2D scattered light model
            
    Raises:
        ValueError: If polynomial order is too high for available groups
        ImportError: If astropy modeling is not available
        
    Example:
        >>> # Standard usage with 4x4 pixel groups and 3rd order polynomial
        >>> bg_sub, bg_model = get_scattered_light(image, trace)
        >>> 
        >>> # Higher resolution with 2x2 groups and 5th order polynomial  
        >>> bg_sub, bg_model = get_scattered_light(image, trace, group_size=2, poly_order=5)
    """

    
    logger.info(f"Estimating scattered light using grouped pixel averaging "
               f"(group_size={group_size}x{group_size}, poly_order={poly_order})")
    
    # Create mask for fiber traces
    mask = make_mask_for_trace(image, trace, fiber_profile_mask_size=fiber_profile_mask_size)
    
    # Get coordinates of all pixels
    rows, cols = np.mgrid[0:image.shape[0], 0:image.shape[1]]
    
    # Find groups of neighboring unmasked pixels
    logger.info("Grouping neighboring unmasked pixels...")
    group_x, group_y, group_values = _group_unmasked_pixels(image, mask, rows, cols, group_size, sigma_clip_threshold)
    
    n_groups = len(group_x)
    logger.info(f"Created {n_groups:,} pixel groups from {np.sum(~mask):,} unmasked pixels")
    
    if n_groups == 0:
        raise ValueError("No background pixel groups available - fiber_profile_mask_size may be too large")
    
    # Check if we have enough groups for the polynomial order
    n_coeffs = (poly_order + 1) * (poly_order + 2) // 2  # Number of coefficients in 2D polynomial
    if n_groups < n_coeffs:
        raise ValueError(
            f"Insufficient pixel groups for 2D polynomial fitting. "
            f"Need ≥{n_coeffs} groups for poly_order={poly_order}, "
            f"but only have {n_groups} groups. "
            f"Reduce poly_order to ≤{_max_poly_order_for_groups(n_groups)} "
            f"or decrease group_size/fiber_profile_mask_size."
        )
    
    # Normalize coordinates to [0,1] for numerical stability
    x_norm = group_x / (image.shape[1] - 1)
    y_norm = group_y / (image.shape[0] - 1)
    
    # Create and fit Astropy Polynomial2D model
    logger.info(f"Fitting Polynomial2D model (degree={poly_order}) to {n_groups} groups")
    
    # Initialize the 2D polynomial model
    poly_model = models.Polynomial2D(degree=poly_order)
    
    # Choose fitter (LevMarLSQFitter is robust and handles outliers well)
    fitter = fitting.LevMarLSQFitter()
    
    try:
        # Fit the model
        fitted_model = fitter(poly_model, x_norm, y_norm, group_values)
        
        # Check fit quality
        residuals = group_values - fitted_model(x_norm, y_norm)
        rms_residual = mad_std(residuals)
        
        logger.info(f"Polynomial fit complete: RMS residual = {rms_residual:.2f}")
        
    except Exception as e:
        raise ValueError(f"Polynomial2D fitting failed: {str(e)}. "
                        "Try reducing poly_order or increasing group_size.")

    # Generate diagnostic plot if output folder is provided
    if outfolder is not None:
        logger.info("Creating scattered light fit diagnostic plot")
        try:
            plot_scattered_light_fit(x_norm, y_norm, group_values, residuals,
                                   poly_order, rms_residual, n_groups, outfolder)
            logger.info("Diagnostic plot saved to scattered_light_fit_diagnostics.png")
        except Exception as e:
            logger.warning(f"Failed to create diagnostic plot: {str(e)}")

    # Generate scattered light model for entire image
    logger.info("Generating scattered light model for full image")
    
    # Create normalized coordinate grids for full image
    cols_full, rows_full = np.mgrid[0:image.shape[1], 0:image.shape[0]]
    x_full_norm = cols_full / (image.shape[1] - 1)
    y_full_norm = rows_full / (image.shape[0] - 1)
    
    # Evaluate model over entire image (note: Polynomial2D expects (x,y) order)
    scattered_light = fitted_model(x_full_norm, y_full_norm).T
    
    # Subtract scattered light
    background_subtracted = image - scattered_light
    
    logger.info("Scattered light background estimation complete")
    
    return background_subtracted, scattered_light


def _group_unmasked_pixels(image, mask, rows, cols, group_size, sigma_clip_threshold):
    """
    Groups neighboring unmasked pixels and computes robust statistics for each group.
    
    Args:
        image (np.ndarray): Input image
        mask (np.ndarray): Boolean mask (True = masked pixels)
        rows, cols (np.ndarray): Row and column coordinate arrays
        group_size (int): Size of pixel groups (group_size x group_size)
        sigma_clip_threshold (float): Sigma threshold for clipping within groups
        
    Returns:
        tuple: (group_x, group_y, group_values) arrays
    """
    group_x = []
    group_y = []
    group_values = []
    
    # Iterate over image in group_size x group_size blocks
    for row_start in range(0, image.shape[0], group_size):
        for col_start in range(0, image.shape[1], group_size):
            # Define group boundaries
            row_end = min(row_start + group_size, image.shape[0])
            col_end = min(col_start + group_size, image.shape[1])
            
            # Extract group data
            group_mask = mask[row_start:row_end, col_start:col_end]
            group_image = image[row_start:row_end, col_start:col_end]
            group_rows = rows[row_start:row_end, col_start:col_end]
            group_cols = cols[row_start:row_end, col_start:col_end]
            
            # Only process groups with at least some unmasked pixels
            unmasked_in_group = ~group_mask
            if not np.any(unmasked_in_group):
                continue
                
            # Extract unmasked values from this group
            unmasked_values = group_image[unmasked_in_group]
            unmasked_rows = group_rows[unmasked_in_group]
            unmasked_cols = group_cols[unmasked_in_group]
            
            # Skip groups with too few pixels
            if len(unmasked_values) < 3:
                continue
            
            # Apply sigma clipping to remove outliers within the group
            if sigma_clip_threshold > 0:
                clipped_values = sigma_clip(unmasked_values, sigma=sigma_clip_threshold, masked=True)
                if np.ma.count(clipped_values) < 2:  # Need at least 2 values
                    continue
                final_values = clipped_values[~clipped_values.mask]
                # Get corresponding coordinates for non-clipped values
                valid_indices = ~clipped_values.mask
                final_rows = unmasked_rows[valid_indices]
                final_cols = unmasked_cols[valid_indices]
            else:
                final_values = unmasked_values
                final_rows = unmasked_rows
                final_cols = unmasked_cols
            
            # Compute robust statistic for this group
            if len(final_values) >= 3:
                group_value = biweight_location(final_values)
            else:
                continue
            
            # Use mean position of contributing pixels as group center
            group_center_row = np.mean(final_rows)
            group_center_col = np.mean(final_cols)
            
            group_x.append(group_center_col)
            group_y.append(group_center_row)
            group_values.append(group_value)
    
    return np.array(group_x), np.array(group_y), np.array(group_values)


def _max_poly_order_for_groups(n_groups):
    """
    Calculate maximum polynomial order possible for given number of groups.
    
    For a 2D polynomial of degree n, the number of coefficients is (n+1)(n+2)/2.
    This function finds the largest n such that (n+1)(n+2)/2 <= n_groups.
    """
    for order in range(1, 20):  # Reasonable upper limit
        n_coeffs = (order + 1) * (order + 2) // 2
        if n_coeffs > n_groups:
            return order - 1
    return 19  # Fallback maximum


def detect_cosmic_rays(image, error=None, gain=1.0, readnoise=6.5, sigclip=4.5, 
                      sigfrac=0.3, objlim=5.0, satlevel=65535.0, niter=4, 
                      sepmed=True, cleantype='meanmask', fsmode='median'):
    """
    Detect and flag cosmic rays in CCD images using the LA Cosmic algorithm.
    
    This function uses the astroscrappy implementation of the LA Cosmic algorithm
    to identify cosmic ray hits in astronomical images. The algorithm works by
    detecting sharp, positive deviations from the local background that are
    characteristic of cosmic ray events.
    
    Args:
        image (numpy.ndarray): Input image data array. Should include the sky 
            background or a mean background level.
        error (numpy.ndarray, optional): Error/uncertainty array corresponding
            to the input image. If None, error will be estimated from image
            statistics. Defaults to None.
        gain (float, optional): Detector gain in e-/ADU. Used to properly 
            estimate photon noise. Defaults to 1.0.
        readnoise (float, optional): Detector read noise in electrons. Used
            for error estimation. Defaults to 6.5.
        sigclip (float, optional): Laplacian-to-noise limit for cosmic ray
            detection. Lower values detect fainter cosmic rays but may flag
            more false positives. Defaults to 4.5.
        sigfrac (float, optional): Fractional detection limit for neighboring
            pixels. Helps distinguish cosmic rays from compact astronomical
            sources. Defaults to 0.3.
        objlim (float, optional): Minimum contrast between LA Cosmic kernel
            and fine structure noise. Higher values preserve more real
            astronomical structure. Defaults to 5.0.
        satlevel (float, optional): Saturation level of the detector in ADU.
            Pixels above this level are flagged as saturated. Defaults to 65535.0.
        niter (int, optional): Number of iterations of the LA Cosmic algorithm.
            More iterations can detect additional cosmic rays but increase
            computation time. Defaults to 4.
        sepmed (bool, optional): Use separable median filter instead of full
            2D median for improved speed. Defaults to True.
        cleantype (str, optional): Type of cleaning to apply. Options are
            'median', 'medmask', 'meanmask', or 'idw'. Defaults to 'meanmask'.
        fsmode (str, optional): Method for removing large-scale structure.
            Options are 'median' or 'convolve'. Defaults to 'median'.
            
    Returns:
        tuple: A tuple containing:
            - cleaned_image (numpy.ndarray): Input image with cosmic rays
              replaced by interpolated values.
            - cosmic_ray_mask (numpy.ndarray): Boolean mask where True
              indicates pixels flagged as cosmic rays.
              
    Raises:
        ImportError: If astroscrappy package is not installed.
        ValueError: If input image is not a valid numpy array.
        
    Example:
        >>> import numpy as np
        >>> from antigen import ccd
        >>> # Load your image data
        >>> image_data = np.random.normal(1000, 50, (1024, 1024))  # Example image
        >>> # Detect cosmic rays
        >>> clean_image, cr_mask = ccd.detect_cosmic_rays(
        ...     image_data, gain=2.3, readnoise=4.5, sigclip=5.0
        ... )
        >>> # Report number of cosmic rays found
        >>> num_cosmic_rays = np.sum(cr_mask)
        >>> print(f"Detected {num_cosmic_rays} cosmic ray pixels")

    References:
        van Dokkum, P. G. 2001, PASP, 113, 1420 [[1]](https://astroscrappy.readthedocs.io/en/latest/api/astroscrappy.detect_cosmics.html)
        (Original LA Cosmic algorithm paper)
    """
    try:
        from astroscrappy import detect_cosmics
    except ImportError:
        raise ImportError(
            "astroscrappy package is required for cosmic ray detection. "
            "Install it with: conda install astroscrappy or pip install astroscrappy."
        )


    # Convert image to float32 and ensure C-contiguity
    image_float32 = np.ascontiguousarray(image, dtype=np.float32)

    # Convert error array to float32 and ensure C-contiguity if provided
    error_float32 = None
    if error is not None:
        if not isinstance(error, np.ndarray):
            raise ValueError("Error array must be a numpy array")
        error_float32 = np.ascontiguousarray(error, dtype=np.float32)


    # Call detect_cosmics with specified parameters
    cosmic_ray_mask, cleaned_image = detect_cosmics(
        image_float32,
        inmask=None,
        invar=error_float32**2 if error_float32 is not None else None,
        gain=gain,
        readnoise=readnoise,
        sigclip=sigclip,
        sigfrac=sigfrac,
        objlim=objlim,
        satlevel=satlevel,
        niter=niter,
        sepmed=sepmed,
        cleantype=cleantype,
        fsmode=fsmode,
        verbose=False
    )
    
    return cleaned_image, cosmic_ray_mask