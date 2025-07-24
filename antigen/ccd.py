import numpy as np
from astropy.io import fits
from astropy.stats import biweight_location as biweight
from astropy.time import Time

from antigen import config


def prep_image(image, channel):
    """
    Orient the images from blue to red (left to right). Fibers are oriented to match configuration files.
    Note: These flips are unique to the VIRUS/LRS2 amplifiers

    Args:
        image (np.ndarray): FITS 2D image data array
        channel (str): one of four channel char identifiers, e.g. 'g', 'b', 'r' or 'd'

    Returns:
        image (np.ndarray): Oriented fits 2D image data array, corrected for what amplifier it comes from.
    """
    # TODO: update docstring to match input args, function signature
    # TODO: the amplifier flip and over-scan should be read from a detector CONFIG

    overscan_length = 32
    bias_value = biweight(image[:, -(overscan_length-2):])
    image = image[:, :-overscan_length] - bias_value

    if channel == "b":
        image[:] = image[::-1, :] # flip-Y, due to how amplifier is reading the CCD
    if channel == "g":
        image[:] = image[:, ::-1] # flip-X
    if channel == "r":
        image[:] = image[::-1, :] # flip-Y
    if channel == "d":
        image[:] = image[:, ::-1] # flip-X
    # Overscan subtraction
    # TODO: comment line above implies something is missing?

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
    image = prep_image(data, channel)

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
    frames = [prep_image(fits.open(file)[0].data, channel) for file in filenames]

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
