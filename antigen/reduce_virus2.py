import os
import warnings
import traceback

from astropy.io import fits
from astropy.stats import biweight_location as biweight
from astropy.table import Table
from astropy.time import Time

import numpy as np

from antigen import ccd
from antigen import config
from antigen import fiber
from antigen import io
from antigen import plot
from antigen import sky
from antigen import utils

# Turn off annoying warnings (even though some deserve attention)
warnings.filterwarnings("ignore")


def reduce(fn, biastime_list, masterbias_list, masterflt_list, flttime_list,
           trace_list, wave_time, wave_list, ftf_list, channel, 
           pca=None, outfolder=None):
    """
    Purpose: Reduce the raw data by performing a series of processing steps,
    including bias subtraction, flat-fielding, sky subtraction, 
    and PCA-based residuals analysis.

    Args:
        fn (str): The filename of the FITS file containing the data.
        biastime_list (list): List of times associated with bias frames.
        masterbias_list (list): List of master bias frames for bias correction.
        flttime_list (list): List of times associated with flat field frames.
        trace_list (list): List of fiber trace data.
        wave_time (list): List of times associated with wavelength calibration.
        wave_list (list): List of wavelength calibration data.
        ftf_list (list): List of flat-field corrections.
        pca (sklearn.decomposition.PCA), Default is None, optional, A pre-fitted PCA model for residual map analysis.

    Returns:
        pca (sklearn.decomposition.PCA): The fitted PCA model, returned if `pca` is None.  # TODO: this original comment does NOT match the data type returned by biweight()
            biweighted_spectrum (np.ndarray): This is what is actually returned.  # TODO: needs review!
        continuum (np.ndarray): 1d numpy array, The computed continuum for the spectrum.
    """
    # TODO: update docstring to match function signature, input args

    # Open the FITS file and extract the observation time
    f = fits.open(fn)
    t = Time(f[0].header['DATE-OBS'] + 'T' + f[0].header['UT'])
    mtime = t.mjd

    # Select appropriate master bias frame based on observation time
    masterbias = get_element_with_closest_time(masterbias_list, biastime_list, mtime)

    # Perform basic image reduction (bias subtraction, gain adjustment)
    image, E = ccd.base_reduction(f[0].data, masterbias, channel)

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

    # TODO: get clarification on the `cont` vs `cont1` handling here. Renamed `cont` to `continuum` to make it more obvious which is returned.

    # Generate a sky mask and the continuum for sky subtraction
    skymask, continuum = sky.get_skymask(biweight(specrect, axis=0, ignore_nan=True), size=25)

    # Subtract the sky from the spectrum
    skysubrect = sky.subtract_sky(specrect, good)

    # If PCA is not provided, compute it from the sky-subtracted data
    if pca is None:
        pca = sky.get_arc_pca(skysubrect, good, skymask, components=config.CONFIG_PCA_COMPONENTS)
        return pca

    # Adjust the sky mask and compute the continuum
    skymask[1:] += skymask[:-1]
    skymask[:-1] += skymask[1:]
    cont1 = sky.get_continuum(skysubrect, skymask, nbins=50)

    # Compute the residuals by subtracting the continuum
    Z = skysubrect - cont1
    res = sky.get_residual_map(Z, pca)

    # Mask out residuals where sky mask is not valid
    res[:, ~skymask] = 0.0

    # Write the final reduced data to a new FITS file
    io.write_fits(skysubrect - res, skysubrect, specrect, errrect, f[0].header, channel, outfolder)

    # Return the biweighted spectrum and continuum
    biweighted_spectrum = biweight(specrect, axis=0, ignore_nan=True)
    return biweighted_spectrum, continuum


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

    log = utils.setup_logging('virus2_reductions')

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
        masterbias_list, biastime_list = ccd.make_mastercal_list(bias_filenames, bias_break_inds, channel)

        log.info('Making master flat frames')
        flt_break_inds = io.get_file_block_break_indices(flt_filenames)
        masterflt_list, flttime_list = ccd.make_mastercal_list(flt_filenames, flt_break_inds, channel)

        log.info('Making master arc frames')
        arc_break_inds = io.get_file_block_break_indices(arc_filenames)
        masterarc_list, arctime_list = ccd.make_mastercal_list(arc_filenames, arc_break_inds, channel)

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
            plot.plot_trace(trace, Tchunk, xchunk, outfolder=outfolder)
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
                plot.plot_wavelength(lines, W, wavelength, outfolder)
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
