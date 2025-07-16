import numpy as np
from astropy.io import fits
from astropy.stats import biweight_location as biweight
from astropy.time import Time

from antigen import config


def prep_image(image, channel):
    """
    Orient the images from blue to red (left to right)
    Fibers are oriented to match configuration files

    Parameters
    ----------
    image : 2d numpy array
        fits image
    amp : str
        Amplifier for the fits image
    ampname : str
        Amplifier name is the location of the amplifier

    Returns
    -------
    image : 2d numpy array
        Oriented fits image correcting for what amplifier it comes from
        These flips are unique to the VIRUS/LRS2 amplifiers
    """
    # TODO: update docstring to match input args, function signature

    overscan_length = 32
    bias_value = biweight(image[:, -(overscan_length-2):])
    image = image[:, :-overscan_length] - bias_value

    if channel == "b":
        image[:] = image[::-1, :]
    if channel == "g":
        image[:] = image[:, ::-1]
    if channel == "r":
        image[:] = image[::-1, :]
    if channel == "d":
        image[:] = image[:, ::-1]
    # Overscan subtraction


    return image


def base_reduction(data, masterbias, channel):
    """
    Perform basic image reduction by applying bias subtraction, gain correction,
    and calculating the error estimate.

    Parameters
    ----------
    data : 2d numpy array
        Raw input image to be reduced.
    masterbias : 2d numpy array
        Master bias frame to be subtracted from the image.

    Returns
    -------
    image : 2d numpy array
        Reduced image with bias subtracted and gain applied.
    E : 2d numpy array
        Error estimate for each pixel, including read noise and photon noise.
    """
    # TODO: update docstring to match input args, function signature

    # Preprocess the raw image (e.g., background subtraction, padding)
    image = prep_image(data, channel)

    # Subtract the master bias from the image
    image[:] -= masterbias

    # Apply gain correction to convert counts to electrons
    gain = config.CONFIG_CHANNEL_DETECTOR[channel]['gain']
    image[:] *= gain

    # Calculate the error estimate (read noise + photon noise)
    rdnoise = config.CONFIG_CHANNEL_DETECTOR[channel]['rdnoise']
    E = np.sqrt(rdnoise**2 + np.where(image > 0., image, 0.))

    # Return the reduced image and the error estimate
    return image, E


def make_mastercal_list(filenames, breakind, channel):
    """
    Creates a list of master calibration images and corresponding times
    by splitting the input list of filenames at given indices.

    Parameters
    ----------
    filenames : list of str
        List of FITS file paths containing calibration data.
    breakind : list of int
        List of indices to split the filenames into different chunks.

    Returns
    -------
    masters : list of 2D numpy arrays
        List of median calibration images for each chunk.
    times : list of float
        List of mean observation times (MJD) corresponding to each chunk.
    """
    # TODO: update docstring to match function signature, input args

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
    Finds the index of the closest calibration time to a given observation time.

    Parameters
    ----------
    mtime : float
        The observation time (MJD) to compare against the calibration times.
    time_list : list of float
        A list of calibration times (MJD).

    Returns
    -------
    int
        The index of the closest calibration time in `time_list`.
    """

    # Find the index of the closest time in time_list to mtime by minimizing the absolute time difference
    return np.argmin(np.abs(mtime - np.array(time_list)))
