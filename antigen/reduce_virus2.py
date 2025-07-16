#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import warnings
import sys
import traceback

from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.stats import mad_std
from astropy.stats import biweight_location as biweight
from astropy.table import Table
from astropy.time import Time

import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import percentile_filter
from sklearn.decomposition import PCA

from . import config
from . import io
from . import fiber
from .input_utils import setup_logging
from .plot import plot_wavelength, plot_trace

# Turn off annoying warnings (even though some deserve attention)
warnings.filterwarnings("ignore")


# TODO: mixed use of style guides, Google vs Numpy, especially in docstrings; pick one

def get_script_path():
    """
    Get script path, aka, where does Antigen live?
    """
    return os.path.dirname(os.path.realpath(sys.argv[0]))


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


def subtract_sky(spectra, good):
    """
    Subtract the sky background from spectra by identifying sky fibers 
    and performing a biweight calculation.

    Parameters
    ----------
    spectra : 2d numpy array
        The input spectra data where rows represent fibers and columns represent wavelengths.
    good : 1d numpy array of bools
        Boolean mask indicating which fibers are good (non-sky).

    Returns
    -------
    2d numpy array
        Spectra with the sky background subtracted for each fiber.
    """
    
    # Get the number of fibers and number of wavelength bins
    nfibs, N = spectra.shape

    # Define range for biweight calculation (middle third of the data)
    n1 = int(1. / 3. * N)
    n2 = int(2. / 3. * N)

    # Calculate the biweight of spectra over the middle third of each fiber's data
    y = biweight(spectra[:, n1:n2], axis=1, ignore_nan=True)

    # Identify sky pixels based on the biweighted data and apply a mask
    mask, cont = fiber.identify_sky_pixels(y[good], size=15)

    # Create a mask for fibers that are not good and are sky fibers
    m1 = ~good
    m1[good] = mask
    skyfibers = ~m1

    # Compute the biweighted sky spectrum based on sky fibers
    init_sky = biweight(spectra[skyfibers], axis=0, ignore_nan=True)

    # Subtract the sky spectrum from the original spectra
    return spectra - init_sky[np.newaxis, :]


def get_pca_sky_residuals(data, ncomponents=5):
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

    # Initialize the residual map with zeros
    res = data * 0.

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



def get_arc_pca(arcskysub, good, mask, components=15):
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


def get_continuum(skysub, masksky, nbins=50):
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

    # Initialize the output array for the continuum with zeros
    bigcont = skysub * 0.

    # Loop over each row (spectrum) in the sky-subtracted data
    for j in np.arange(skysub.shape[0]):
        # Copy the current spectrum and mask the sky pixels
        y = skysub[j] * 1.
        y[masksky] = np.nan

        # Divide the spectrum into bins and calculate the mean of each bin
        x = np.array([np.mean(chunk) for chunk in np.array_split(np.arange(len(y)), nbins)])

        # Calculate the biweighted median for each bin
        z = np.array([biweight(chunk, ignore_nan=True) for chunk in np.array_split(y, nbins)])

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

    
def reduce(fn, biastime_list, masterbias_list, masterflt_list, flttime_list,
           trace_list, wave_time, wave_list, ftf_list, channel, 
           pca=None, outfolder=None):
    """
    Reduce the raw data by performing a series of processing steps, 
    including bias subtraction, flat-fielding, sky subtraction, 
    and PCA-based residuals analysis.

    Parameters
    ----------
    fn : str
        The filename of the FITS file containing the data.
    biastime_list : list
        List of times associated with bias frames.
    masterbias_list : list
        List of master bias frames for bias correction.
    flttime_list : list
        List of times associated with flat field frames.
    trace_list : list
        List of fiber trace data.
    wave_time : list
        List of times associated with wavelength calibration.
    wave_list : list
        List of wavelength calibration data.
    ftf_list : list
        List of flat-field corrections.
    pca : PCA object, optional
        A pre-fitted PCA model for residual map analysis. Default is None.

    Returns
    -------
    pca : PCA object
        The fitted PCA model, returned if `pca` is None.
    continuum : 1d numpy array
        The computed continuum for the spectrum.
    """
    # TODO: update docstring to match function signature, input args

    # Open the FITS file and extract the observation time
    f = fits.open(fn)
    t = Time(f[0].header['DATE-OBS'] + 'T' + f[0].header['UT'])
    mtime = t.mjd

    # Select appropriate master bias frame based on observation time
    masterbias = get_element_with_closest_time(masterbias_list, biastime_list, mtime)

    # Perform basic image reduction (bias subtraction, gain adjustment)
    image, E = base_reduction(f[0].data, masterbias, channel)

    # Get the fiber trace and selection mask for the current observation
    trace, good = get_element_with_closest_time(trace_list, flttime_list, mtime)

    # Extract spectra from the image using the trace data
    spec = fiber.get_spectra(image, trace)

    # Calculate the spectrum error using the flat-field and error image
    specerr = fiber.get_spectra_error(E, trace)

    # Compute the chi-square of the spectrum to identify bad pixels
    masterflt = get_element_with_closest_time(masterflt_list, flttime_list, mtime)
    chi2 = fiber.get_spectra_chi2(masterflt - masterbias, image, E, trace)
    badpix = chi2 > 20.  # Pixels with chi2 > 20 are considered bad
    specerr[badpix] = np.nan
    spec[badpix] = np.nan

    # Retrieve the wavelength calibration data for the current observation
    wavelength = get_element_with_closest_time(wave_list, wave_time, mtime)

    # Rectify the spectrum and error based on the wavelength
    def_wave = config.CONFIG_CHANNEL_DEF_WAVE[channel]
    specrect, errrect = fiber.rectify(spec, specerr, wavelength, def_wave)

    # Apply flat-field correction
    ftf = get_element_with_closest_time(ftf_list, flttime_list, mtime)
    specrect[:] /= (ftf * f[0].header['EXPTIME'])
    errrect[:] /= (ftf * f[0].header['EXPTIME'])

    # Generate a sky mask and the continuum for sky subtraction
    skymask, cont = get_skymask(biweight(specrect, axis=0, ignore_nan=True), size=25)

    # Subtract the sky from the spectrum
    skysubrect = subtract_sky(specrect, good)

    # If PCA is not provided, compute it from the sky-subtracted data
    if pca is None:
        pca = get_arc_pca(skysubrect, good, skymask, components=config.CONFIG_PCA_COMPONENTS)
        return pca

    # Adjust the sky mask and compute the continuum
    skymask[1:] += skymask[:-1]
    skymask[:-1] += skymask[1:]
    cont1 = get_continuum(skysubrect, skymask, nbins=50)

    # Compute the residuals by subtracting the continuum
    Z = skysubrect - cont1
    res = get_residual_map(Z, pca)

    # Mask out residuals where sky mask is not valid
    res[:, ~skymask] = 0.0

    # Write the final reduced data to a new FITS file
    io.write_fits(skysubrect - res, skysubrect, specrect, errrect, f[0].header, channel, outfolder)

    # Return the biweighted spectrum and continuum
    return biweight(specrect, axis=0,ignore_nan=True), cont


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


def get_element_with_closest_time(element_list, time_list, target_time):
    """
    Purpose: Given a list of elements and a list of MJD times for those elements,
    find the element that has a time closest to the target_time

    Args:
        element_list (list): list of elements corresponding to the times in time_list
        time_list (list(numeric)): list of MJD times corresponding to the times for the elements in element_list
        target_time (numeric): MJD time

    Returns:
        closest_element: element for closest_time from time_list compared to target_time
    """
    # Find index of time closest to target_time
    index = np.argmin(np.abs(np.array(time_list) - target_time))

    # Use the index directly in the original file list
    closest_element = element_list[index]
    return closest_element


def get_filenames(gnames, typelist, names):
    """
    Finds filenames that match a list of keywords by checking if any of the keywords 
    appear in the associated types.

    Parameters
    ----------
    gnames : list of str
        List of filenames to check.
    typelist : list of str
        List of types or descriptions corresponding to the filenames.
    names : list of str
        List of keywords to search for within the types/descriptions.

    Returns
    -------
    np.ndarray
        Array of filenames from `gnames` that match any of the keywords in `names`.
    """
    # TODO: fix type, replace numpy.ndarray with lists when handling filenames; no vector operations are used
    matches = []  # List to store matched filenames
    # Iterate through each filename and its corresponding type
    for gn, tp in zip(gnames, typelist):
        # Check if any keyword appears in the type (case-insensitive)
        for name in names:
            if name.lower() in str(tp).lower():
                matches.append(gn)  # Append matching filename to the list
    return np.array(matches)  # Return matched filenames as a numpy array


def process(infolder, outfolder, obs_date, obs_name, reduce_all,
            bias_label, arc_label, dark_label, flat_label, twilight_flat_label):
    """
    Purpose: data reduction pipeline to process VIRUS2 observation files

    Args:
        infolder (str): Root path where reduction input file tree is located
        outfolder (str): Path where reduction output files will be written
        obs_date (str): Observation calendar date string formatted as YYYYMMDD
        obs_name (str): Observation object/target name, e.g. from FITS header card
        reduce_all (bool): Reduce all files found under infolder file tree
        bias_label (str): string label from FITS file header for bias frames
        arc_label (str): string label from FITS file header for arc frames
        dark_label (str): string label from FITS file header for dark frames
        flat_label (str): string label from FITS file header for flat frames
        twilight_flat_label (str): string label from FITS file header for twilight frames

    Returns:
        None
    """
    # TODO: dark_label is unused, in current and previous versions of this module

    log = setup_logging('virus2_reductions')

    os.makedirs(outfolder, exist_ok=True)

    ROOT_DATA_PATH = infolder
    metadata_records = io.parse_fits_file_tree(ROOT_DATA_PATH, date=obs_date, verbose=True)
    unit_list = [record['spec_id'] for record in metadata_records]
    units = list(set(unit_list))

    # =============================================================================
    # Get Line List
    # =============================================================================

    for unit in units:

        # =========================================
        # Get def_wave and line_list for given channel
        # TODO: implemented the other three channels, using green as a template
        # =========================================
        channel = unit[-1].lower()
        line_list = None
        limit = None
        xref = None
        if channel == 'b':
            def_wave = config.CONFIG_CHANNEL_DEF_WAVE[channel]
            continue
        if channel == 'g':
            def_wave = config.CONFIG_CHANNEL_DEF_WAVE[channel]
            line_list_filepath = config.get_config_filepath('line_list', 'virus2_green.txt')
            line_list = Table.read(line_list_filepath, format="ascii")
            limit = config.CONFIG_CHANNEL_DETECTOR[channel]['limit']
        if channel == 'r':
            def_wave = config.CONFIG_CHANNEL_DEF_WAVE[channel]
            continue
        if channel == 'd':
            def_wave = config.CONFIG_CHANNEL_DEF_WAVE[channel]
            continue

        lines = np.array(line_list['wavelength'])
        xref = np.array(line_list['col'])
        use_kernel = True

        # =============================================================================
        # Make a list of the objects for each filename, ignore those without 'OBJECT'
        # in the header
        # =============================================================================
        unit_filenames = [record['filename'] for record in metadata_records if record['spec_id'] == unit]

        typelist = []
        valid_files = []
        timelist = []

        for fn in unit_filenames:
            fn_meta_data = io.parse_fits_file_name(fn)
            try:
                obj_frame_type = fn_meta_data['frame_type']  # fits.open(f)[0].header['OBJECT']
                obs_date_str = fn_meta_data['utc_str_date']  # fits.open(f)[0].header['DATE-OBS']
            except:
                continue

            # Only store files for which the filename meta-data parsing is valid
            typelist.append(obj_frame_type)
            valid_files.append(fn)
            timelist.append(Time(obs_date_str))

        # =============================================================================
        # Get/sort subsets of filenames that are bias filenames, domeflat filenames, and arc lamp filenames
        # =============================================================================

        log.info('Sorting FITS Files by frame type')

        bias_filenames    = io.get_matching_filenames(valid_files, typelist, [bias_label])
        twiflt_filenames  = io.get_matching_filenames(valid_files, typelist, [twilight_flat_label])
        domeflt_filenames = io.get_matching_filenames(valid_files, typelist, [flat_label])
        arc_filenames     = io.get_matching_filenames(valid_files, typelist, [arc_label])

        if reduce_all:
            non_sci_files = set(bias_filenames) | set(twiflt_filenames) | set(domeflt_filenames) | set(arc_filenames)
            sci_filenames = [fn for fn in valid_files if fn not in non_sci_files]
        else:
            if obs_name is not None:
                sci_filenames =  io.get_matching_filenames(valid_files, typelist, [obs_name])
            else:
                sci_filenames = []

        if len(twiflt_filenames) > 0:
            flt_filenames = twiflt_filenames
        else:
            flt_filenames = domeflt_filenames

        # =============================================================================
        # Make a master bias, master dome flat, and master arc for the first set of OBS
        # Use the filename obs_id numbers for grouping/splitting contiguous blocks of observations files
        # =============================================================================

        log.info('Making master bias frames')
        bias_break_inds = io.get_file_block_break_indices(bias_filenames)
        masterbias_list, biastime_list = make_mastercal_list(bias_filenames, bias_break_inds, channel)

        log.info('Making master flat frames')
        flt_break_inds = io.get_file_block_break_indices(flt_filenames)
        masterflt_list, flttime_list = make_mastercal_list(flt_filenames, flt_break_inds, channel)

        log.info('Making master arc frames')
        arc_break_inds = io.get_file_block_break_indices(arc_filenames)
        masterarc_list, arctime_list = make_mastercal_list(arc_filenames, arc_break_inds, channel)

        # =============================================================================
        # Load reference fiber locations from a predefined file
        # =============================================================================

        ifu_cen_filepath = config.get_config_filepath('ifucen', 'IFUcen_VIRUS2_D3G.txt')
        ifu_cen_file_data = Table.read(ifu_cen_filepath, format="ascii")
        ref = ifu_cen_file_data
        ref.reverse()

        # =============================================================================
        # Get trace from the dome flat
        # =============================================================================
        trace_list = []
        fltspec = []
        log.info('Getting trace for each master flat')
        for masterflt, mtime in zip(masterflt_list, flttime_list):
            masterbias = get_element_with_closest_time(masterbias_list, biastime_list, mtime)
            trace, good, Tchunk, xchunk = fiber.get_trace(masterflt - masterbias, ref)
            plot_trace(trace, Tchunk, xchunk, outfolder=outfolder)
            trace_list.append([trace, good])
            domeflat_spec = fiber.get_spectra(masterflt - masterbias, trace)
            domeflat_error = 0. * domeflat_spec
            fltspec.append([domeflat_spec, domeflat_error])

        # =============================================================================
        # Get wavelength from arc lamps
        # =============================================================================
        wave_list = []
        wave_time = []
        bk1 = np.hstack([0, arc_break_inds+1])
        log.info('Getting wavelength for each master arc')

        for masterarc, mtime, bk in zip(masterarc_list, arctime_list, bk1):
            masterbias = get_element_with_closest_time(masterbias_list, biastime_list, mtime)
            trace, good = get_element_with_closest_time(trace_list, flttime_list, mtime)
            lamp_spec = fiber.get_spectra(masterarc - masterbias, trace)
            # TODO: move definition of test.fits to CONFIG param
            fits.PrimaryHDU(lamp_spec).writeto('test.fits',overwrite=True)
            try:
                wavelength, res, X, W = fiber.get_wavelength(lamp_spec, trace, good,
                                                       xref, lines, limit=limit,
                                                       use_kernel=use_kernel)

                # Plot wavelength solution for inspection
                plot_wavelength(lines, W, wavelength, outfolder)
            except:
                log.warning('Could not get wavelength solution for masterarc')
                log.warning('First file of failed masterarc included: %s' %
                            (arc_filenames[bk]))
                traceback.print_exc()
                continue
            wave_list.append(wavelength)
            wave_time.append(mtime)

        # =============================================================================
        # Rectify domeflat spectra and get fiber to fiber
        # =============================================================================
        ftf_list = []
        log.info('Getting fiber to fiber for each master domeFlat')
        for fltsp, mtime in zip(fltspec, flttime_list):
            wavelength = get_element_with_closest_time(wave_list, wave_time, mtime)
            domeflat_spec, domeflat_error = fltsp
            domeflat_rect, domeflat_error_rect = fiber.rectify(domeflat_spec, domeflat_error,
                                                               wavelength, def_wave)
            ftf, ftf_smooth = fiber.get_fiber_to_fiber(domeflat_rect)
            ftf_list.append(ftf)


        pca = reduce(arc_filenames[0], biastime_list, masterbias_list, masterflt_list, flttime_list,
                     trace_list, wave_time, wave_list, ftf_list, channel, pca=None, outfolder=outfolder)
        for fn in sci_filenames:
            log.info('Reducing: %s' % fn)
            sky, cont = reduce(fn, biastime_list, masterbias_list, masterflt_list, flttime_list,
                               trace_list, wave_time, wave_list, ftf_list,
                               channel, pca=pca, outfolder=outfolder)
    return None
