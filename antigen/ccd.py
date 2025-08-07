import numpy as np
from astropy.stats import biweight_location as biweight
from astropy.time import Time

from scipy.interpolate import LSQUnivariateSpline

from antigen.io import load_fits


def prep_image(image, config_dict):
    """
    Measure the bias in the overscan region, trim the overscan and subtract the bias
    Perform any flips that are necessary

    Args:
        image (np.ndarray): FITS 2D image data array
        config_dict (dict): Dictionary of configuration parameters

    Returns:
        image (np.ndarray): Oriented fits 2D image data array, corrected for what amplifier it comes from.
    """
    # TODO: update docstring to match input args, function signature
    # TODO: the amplifier flip and over-scan should be read from a detector CONFIG
    # Should be read by config rather than hardcoded
    overscan_length = config_dict['overscan_length']
    flip_x = config_dict['flip_x']
    flip_y = config_dict['flip_y']
    rotate = config_dict['rotate']
    add_rows = config_dict['add_rows']

    bias_value = biweight(image[:, -(overscan_length-2):])
    image = image[:, :-overscan_length] - bias_value
    if flip_x:
        image = np.flip(image, axis=1)
    if flip_y:
        image = np.flip(image, axis=0)
    if rotate:
        image = np.rot90(image)
    if add_rows:
        # Initialize a new image array with additional rows for padding
        new_image = np.zeros((len(image) + add_rows, image.shape[1]))

        # Copy the processed image into the new array, leaving the top rows as padding
        new_image[add_rows:, :] = image
        image = new_image

    return image


def base_reduction(data, master_bias, config_dict):
    """
    Perform basic image reduction by applying bias subtraction,
    gain correction, and calculating the error estimate.

    Args:
        data (np.ndarray): 2d numpy array, Raw input image to be reduced.
        master_bias (np.ndarray): 2d numpy array, Master bias frame to be subtracted from the image.
        config_dict: Dictionary of configuration parameters

    Returns:
        image (np.ndarray): 2d numpy array, Reduced image with bias subtracted and gain applied.
        error_estimate (np.ndarray): 2d numpy array, Error estimate for each pixel, including read noise and photon noise.
    """
    # TODO: update docstring to match input args, function signature

    # Preprocess the raw image (e.g., background subtraction, padding)
    image = prep_image(data, config_dict)

    # Subtract the master bias from the image
    image[:] -= master_bias

    # Apply gain correction to convert counts to electrons
    gain = config_dict['gain']
    image[:] *= gain

    # Calculate the error estimate (read noise + photon noise)
    read_noise = config_dict['read_noise']
    error_estimate = np.sqrt(read_noise**2 + np.where(image > 0., image, 0.))

    # Return the reduced image and the error estimate
    return image, error_estimate


def make_master_cal(filenames, config_dict):
    """
    Load a list of calibration FITS files, apply basic CCD preprocessing, and create a master calibration frame.

    Args:
        filenames (list(str)): list of filenames
        config_dict: Dictionary of configuration parameters
    Returns:
        master_cal (np.ndarray): median stacked calibration frame
        master_cal_time (float): average MJD of the frames that were stacked
    """
    # Extract from the files, re-oriented by prep_image()
    # Extract observation times (MJD) for frames in the current chunk
    frames, times = ([], [])
    for filename in filenames:
        frame, header = load_fits(filename)
        prepped_frame = prep_image(frame, config_dict)
        frames.append(prepped_frame)
        mjd = Time(header['DATE-OBS']).mjd
        times.append(mjd)

    # Compute median frame and the mean time for the current chunk
    master_cal      = np.nanmedian(frames, axis=0)  # maybe biweight() as an alternate method
    master_cal_time = np.mean(times)
    return master_cal, master_cal_time


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

def get_scattered_light(image, trace, fiber_profile_mask_size=11, nchunks=20):
    """
    Estimates the scattered light background in an image by masking out a region around the trace
    of each fiber and fitting splines in the row direction then column direction

    Args:
        image (np.ndarray): 2D numpy image array.
        trace (np.ndarray): 2D numpy array describing the trace position as a function of fiber
        fiber_profile_mask_size (int, optional): Vertical size (in pixels) of the mask region
            above and below the trace. Default is 11.
        nchunks (int, optional): Number of column chunks for Spline fitting
    Returns:
        scattered_light: Estimated background image with the same shape as `image`.
    """
    mask = make_mask_for_trace(image, trace, fiber_profile_mask_size=fiber_profile_mask_size)
    columns = np.arange(image.shape[1])
    row_fit = np.arange(image.shape[0]) / image.shape[0]

    masked_image = image * 1.
    masked_image[mask] = np.nan
    scattered_light = np.zeros_like(image)
    masked_image_chunks = np.array_split(masked_image, nchunks, axis=1)
    column_chunks = np.array_split(columns, nchunks)
    for masked_image_chunk, column_chunk in zip(masked_image_chunks, column_chunks):
        column_average = np.nanmedian(masked_image_chunk, axis=1)
        unmasked = np.isfinite(column_average)
        # Find break points to set spline knots
        diff = np.diff(unmasked.astype(int))
        # Start = 1 in diff → next index starts new unmasked chunk
        starts = np.where(diff == 1)[0] + 1
        # Stop = -1 in diff → current index ends unmasked chunk
        stops = np.where(diff == -1)[0] + 1
        # Handle edge cases:
        if unmasked[0]:
            starts = np.insert(starts, 0, 0)
        if unmasked[-1]:
            stops = np.append(stops, len(unmasked))
        # Get average normalized row value in each chunk
        row_chunk_means = []
        for start, stop in zip(starts, stops):
            indices = np.arange(start, stop)
            chunk_mean = row_fit[indices].mean()
            row_chunk_means.append(chunk_mean)

        # Fit a spline with knots at the average row chunk values
        spline = LSQUnivariateSpline(row_fit[unmasked], column_average[unmasked], t=row_chunk_means, k=3)
        scattered_light[:, column_chunk] = spline(row_fit)[:, np.newaxis]

    # Take initial 1D fit and fit a second spline in the column direction for smoothness
    column_fit = np.arange(image.shape[1]) / image.shape[1]
    number_fit_chunks = int(nchunks / 2)
    mean_column_chunks = [np.mean(column_fit_chunk) for column_fit_chunk in np.array_split(column_fit,
                                                                                           number_fit_chunks)]
    for row in np.arange(image.shape[0]):
        spline = LSQUnivariateSpline(column_fit, scattered_light[row], t=mean_column_chunks, k=3)
        scattered_light[row] = spline(column_fit)

    return scattered_light