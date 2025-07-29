import numpy as np
from astropy.io import fits
from astropy.stats import biweight_location as biweight
from astropy.time import Time

from scipy.interpolate import LSQUnivariateSpline
from antigen import config


def prep_image(image):
    """
    Measure the bias in the overscan region, trim the overscan and subtract the bias
    Perform any flips that are necessary

    Args:
        image (np.ndarray): FITS 2D image data array

    Returns:
        image (np.ndarray): Oriented fits 2D image data array, corrected for what amplifier it comes from.
    """
    # TODO: update docstring to match input args, function signature
    # TODO: the amplifier flip and over-scan should be read from a detector CONFIG
    # Should be read by config rather than hardcoded
    overscan_length = 32
    flipx = True
    flipy = False

    bias_value = biweight(image[:, -(overscan_length-2):])
    image = image[:, :-overscan_length] - bias_value
    if flipx:
        image = np.flip(image, axis=1)
    if flipy:
        image = np.flip(image, axis=0)

    return image


def base_reduction(data, masterbias, channel):
    """
    Perform basic image reduction by applying bias subtraction,
    gain correction, and calculating the error estimate.

    Args:
        data (np.ndarray): 2d numpy array, Raw input image to be reduced.
        masterbias (np.ndarray): 2d numpy array, Master bias frame to be subtracted from the image.
        channel (str): one of four channel char identifiers, e.g. 'g', 'b', 'r' or 'd'

    Returns:
        image (np.ndarray): 2d numpy array, Reduced image with bias subtracted and gain applied.
        error_estimate (np.ndarray): 2d numpy array, Error estimate for each pixel, including read noise and photon noise.
    """
    # TODO: update docstring to match input args, function signature

    # Preprocess the raw image (e.g., background subtraction, padding)
    image = prep_image(data)

    # Subtract the master bias from the image
    image[:] -= masterbias

    CHANNEL_DETECTOR, _ = config.get_channel_config_virus2()

    # Apply gain correction to convert counts to electrons
    gain = CHANNEL_DETECTOR[channel]['gain']
    image[:] *= gain

    # Calculate the error estimate (read noise + photon noise)
    rdnoise = CHANNEL_DETECTOR[channel]['rdnoise']
    error_estimate = np.sqrt(rdnoise**2 + np.where(image > 0., image, 0.))

    # Return the reduced image and the error estimate
    return image, error_estimate


def make_master_cal(filenames, channel):
    """
    Purpose: Load all files, slice array into 4 channels, select single channel slice, compute aggregate, return result

    Args:
        filenames (list(str)):
        channel (str): 'g', 'b'
    Returns:
    """
    # Extract from the files, re-oriented by prep_image()
    frames = [prep_image(fits.open(file)[0].data) for file in filenames]

    # Extract observation times (MJD) for frames in the current chunk
    times = [Time(fits.open(file)[0].header['DATE-OBS']).mjd for file in filenames]

    # Compute median frame and the mean time for the current chunk
    master_cal      = np.nanmedian(frames, axis=0)  # maybe biweight() as an alternate method
    master_cal_time = np.mean(times)
    return master_cal, master_cal_time


def make_mastercal_list(filenames, breakind, channel):
    """
    Creates a list of master calibration images and corresponding times
    by splitting the input list of filenames at given indices.

    Args:
        filenames (list(str)): List of FITS file paths containing calibration data.
        breakind (list(int)): List of indices to split the filenames into different chunks.
        channel (str): one of four channel char identifiers, e.g. 'g', 'b', 'r' or 'd'

    Returns
        masters (list(np.ndarray)): List of median calibration images (2D arrays) for each chunk.
        times (list(float)): List of mean observation times (MJD floats) corresponding to each chunk.
    """

    # Define break points for splitting the filenames into chunks
    breakind1 = np.hstack([0, breakind])  # Start indices for chunks
    breakind2 = np.hstack([breakind, len(filenames)+1])  # End indices for chunks

    masters = []  # List to store median calibration images
    times = []    # List to store mean observation times

    # Iterate over the chunks defined by breakind1 and breakind2
    for bk1, bk2 in zip(breakind1, breakind2):
        # Collect and preprocess frames within the current chunk
        frames = [prep_image(fits.open(f)[0].data, channel)
                  for cnt, f in enumerate(filenames)
                  if ((cnt > bk1) * (cnt < bk2))]  # Only include frames in the current chunk

        # Extract observation times (MJD) for frames in the current chunk
        t = [Time(fits.open(f)[0].header['DATE-OBS']).mjd
             for cnt, f in enumerate(filenames)
             if ((cnt > bk1) * (cnt < bk2))]


        # Append the median frame and the mean time for the current chunk
        masters.append(np.nanmedian(frames, axis=0))
        times.append(np.mean(t))

    return masters, times


def get_cal_index(mtime, time_list):
    """
    Finds the index of the closest calibration time to a given observation time, by
    finding the minimized absolute time difference between mtime and all time_list elements.

    Args:
        mtime (float): The observation time (MJD float) to compare against the calibration times.
        time_list (list(float)): A list of calibration times (MJD floats).

    Returns:
        index (int): The index of the closest calibration time in `time_list`.
    """
    return np.argmin(np.abs(mtime - np.array(time_list)))


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