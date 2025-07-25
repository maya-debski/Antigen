import datetime
import os
import warnings

from astropy.io import fits
from astropy.stats import biweight_location as biweight
from astropy.table import Table

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


def reduce(data_filename, master_bias, master_flat, trace, good_fiber_mask, wavelength_cal, ftf_correction,
           channel, pca=None, pca_only=False, outfolder=None, debug=False):
    """
    Purpose: Reduce the raw data by performing a series of processing steps,
    including bias subtraction, flat-fielding, sky subtraction,
    and PCA-based residuals analysis.

    Args:
        data_filename (str): The filename of the FITS file containing the data to be reduced.
        master_bias (): master bias frames for bias correction.
        master_flat (): master flat frame
        trace (): fiber trace data for current observation
        good_fiber_mask (array(bool)): boolean selecttion mask of good (isfinite) fibers for current observation
        wavelength_cal (): wavelength calibration data for the current observation
        ftf_correction (list): flat-field (fiber_to_fiber) correction array.
        pca (sklearn.decomposition.PCA), Default is None, optional, A pre-fitted PCA model for residual map analysis.
        pca_only (bool): if True, return immediately after computing the PCA model fit object
        outfolder (str): file output path to write FITS files to
        debug (bool): default=False, if True, save PNG plot files for intermediate data artifacts

    Returns:
        pca (sklearn.decomposition.PCA): The fitted PCA model
        biweighted_spectrum (np.ndarray): bi-weighted spectrum, statistic for determining the
                                          central location of the specrect distribution
        continuum (np.ndarray): 1d numpy array, The computed continuum for the spectrum.
    """

    CONFIG_DETECTOR, CONFIG_DEF_WAVE = config.get_channel_config_virus2()

    obs_data, obs_header = io.load_fits(data_filename)

    # Perform basic image reduction (bias subtraction, gain adjustment)
    image, E = ccd.base_reduction(obs_data, master_bias, channel)

    # Extract spectra from the image using the trace data
    spec = fiber.get_spectra(image, trace)

    # Calculate the spectrum error using the flat-field and error image
    specerr = fiber.get_spectra_error(E, trace)

    # Compute the chi-square of the spectrum to identify bad pixels
    chi2 = fiber.get_spectra_chi2(master_flat - master_bias, image, E, trace)
    badpix = chi2 > 20.  # Pixels with chi2 > 20 are considered bad
    specerr[badpix] = np.nan
    spec[badpix] = np.nan

    # Rectify the spectrum and error based on the wavelength
    def_wave = CONFIG_DEF_WAVE[channel]
    specrect, errrect = fiber.rectify(spec, specerr, wavelength_cal, def_wave)

    # Apply flat-field correction
    obs_exp_time = obs_header['EXPTIME']
    specrect[:] /= (ftf_correction * obs_exp_time)
    errrect[:] /= (ftf_correction * obs_exp_time)

    # Generate a sky mask and the continuum for sky subtraction
    skymask, continuum = sky.get_skymask(biweight(specrect, axis=0, ignore_nan=True), size=25)

    # Subtract the sky from the spectrum
    skysubrect = sky.subtract_sky(specrect, good_fiber_mask)

    # If PCA is not provided, compute it from the sky-subtracted data
    if pca is None:
        biweighted_spectrum = None
        pca = sky.get_arc_pca(skysubrect, good_fiber_mask, skymask, components=config.VIRUS2_PCA_COMPONENTS)
        if pca_only:
            return pca, biweighted_spectrum, continuum, None

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
    skysubrect_adv = skysubrect - res
    output_fits_filename = io.write_fits(skysubrect_adv, skysubrect, specrect, errrect, obs_header, channel, outfolder)

    if debug:
        plot_title = 'skysubrect_adv'
        save_filename = os.path.abspath(os.path.join(outfolder, f'debug_{plot_title}.png'))
        plot.plot_frame(skysubrect_adv, save_file=save_filename, title=plot_title)

        plot_title = 'skysubrect'
        save_filename = os.path.abspath(os.path.join(outfolder, f'debug_{plot_title}.png'))
        plot.plot_frame(skysubrect, save_file=save_filename, title=plot_title)

        plot_title = 'specrect'
        save_filename = os.path.abspath(os.path.join(outfolder, f'debug_{plot_title}.png'))
        plot.plot_frame(specrect, save_file=save_filename, title=plot_title)

        plot_title = 'errrect'
        save_filename = os.path.abspath(os.path.join(outfolder, f'debug_{plot_title}.png'))
        plot.plot_frame(errrect, save_file=save_filename, title=plot_title)

    # Return the biweighted spectrum and continuum
    biweighted_spectrum = biweight(specrect, axis=0, ignore_nan=True)
    return pca, biweighted_spectrum, continuum, output_fits_filename


def build_manifest_records(infolder, obs_date, obs_name, reduce_all, time_radius,
                           bias_label, arc_label, dark_label, flat_label, twilight_flat_label):
    """
    Purpose: Search FITS file tree and generate a list of manifest records,
    one for each "obs id break block" in the files

    Args:
        infolder (str): Root path where reduction input file tree is located
        obs_date (str): Observation calendar date string formatted as YYYYMMDD
        obs_name (str): Observation object/target name, e.g. from FITS header card
        reduce_all (bool): Reduce all files found under infolder file tree
        time_radius (float) will group calibration files with times that fall within this distance from a given science file time
        bias_label (str): string label from FITS file header for bias frames
        arc_label (str): string label from FITS file header for arc frames
        dark_label (str): string label from FITS file header for dark frames
        flat_label (str): string label from FITS file header for flat frames
        twilight_flat_label (str): string label from FITS file header for twilight frames

    Returns:
        manifest_records (list(dict)): list of dicts/records, where each dict has the following keys:
            record['obs_time']: Obs time as MJD float
            record[observation_files']: filename of obs file to be reduced
            record['calibration_files']['bias']: list of bias FITS file names to use for calibration
            record['calibration_files']['flat']: list of flat FITS file names to use for calibration
            record['calibration_files']['arc']: list of arc FITS file names to use for calibration
    """

    ROOT_DATA_PATH = os.path.abspath(infolder)
    if not os.path.isdir(ROOT_DATA_PATH):
        raise NotADirectoryError(f'ERROR: user-specified input path does not exist: {ROOT_DATA_PATH}')

    # time filtering: Get FITS files in file-tree that match the obs_date
    metadata_records = io.parse_fits_file_tree(ROOT_DATA_PATH, date=obs_date, verbose=True)
    units_found = [record['spec_id'] for record in metadata_records]
    unique_units_found = list(set(units_found))

    manifest_records = []
    for unit in unique_units_found:
        # find subset of file records with matching spec_id, e.g. 'D3G' and get their filenames and frame types
        unit_records = [record for record in metadata_records if record['spec_id'] == unit]

        # Group subsets of filenames into lists based on the frame types
        unit_frame_types = [record['frame_type'] for record in unit_records]
        unit_filenames = [record['filename'] for record in unit_records]

        bias_filenames    = io.get_matching_filenames(unit_filenames, unit_frame_types, [bias_label])
        twiflt_filenames  = io.get_matching_filenames(unit_filenames, unit_frame_types, [twilight_flat_label])
        domeflt_filenames = io.get_matching_filenames(unit_filenames, unit_frame_types, [flat_label])
        arc_filenames     = io.get_matching_filenames(unit_filenames, unit_frame_types, [arc_label])
        dark_filenames    = io.get_matching_filenames(unit_filenames, unit_frame_types, [dark_label])

        if obs_name is not None:
            sci_filenames = io.get_matching_filenames(unit_filenames, unit_frame_types, [obs_name])
        else:
            sci_filenames = []

        if len(sci_filenames) == 0 and not reduce_all:
            warnings.warn(f'unit={unit}, found ZERO matching files for obs_name={obs_name}. Continuing to next unit...')
            continue

        if reduce_all:
            calibration_files = set(bias_filenames) | set(twiflt_filenames) | set(domeflt_filenames) | set(arc_filenames) | set(dark_filenames)
            sci_filenames = [name for name in unit_filenames if name not in calibration_files]

        if len(twiflt_filenames) > 0:
            flt_filenames = twiflt_filenames
        else:
            flt_filenames = domeflt_filenames

        # =============================================================================
        # Validate number of files found before attempting to use diffs on file obs_ids
        # Use exceptions to exit process if needed frame-types file counts were not found
        # =============================================================================
        fail_bias = False
        fail_flat = False
        fail_arc = False
        minimum_file_count_for_break = 1

        bias_minimum_count = minimum_file_count_for_break
        num_bias_files = len(bias_filenames)
        if num_bias_files < bias_minimum_count:
            fail_bias = True
            warnings.warn(f'WARNING: unit={unit}, Searched {ROOT_DATA_PATH}, found BIAS label = {bias_label}, '
                          f'found {num_bias_files}, needed >= {bias_minimum_count}')

        flat_minimum_count = minimum_file_count_for_break
        num_flt_files = len(flt_filenames)
        if num_flt_files < flat_minimum_count:
            fail_flat = True
            warnings.warn(f'WARNING: unit={unit}, Searched {ROOT_DATA_PATH}, FLAT label = {flat_label}, '
                          f'found {num_flt_files}, needed >= {flat_minimum_count}')

        arc_minimum_count = minimum_file_count_for_break
        num_arc_files = len(arc_filenames)
        if num_arc_files < arc_minimum_count:
            fail_arc = False
            warnings.warn(f'WARNING: unit={unit}, Searched {ROOT_DATA_PATH}, ARC label = {arc_label}, '
                          f'found {num_arc_files}, needed >= {arc_minimum_count}')

        if fail_bias or fail_flat or fail_arc:
            warnings.warn(f'WARNING: did not find enough calibration files to process unit={unit}. Continuing to next unit...')
            continue
        # =============================================================================
        # Use the filename obs_id numbers for grouping/splitting contiguous blocks of observations files
        # =============================================================================
        bias_times = [io.get_fits_file_time(name) for name in bias_filenames]
        flt_times  = [io.get_fits_file_time(name) for name in flt_filenames]
        arc_times  = [io.get_fits_file_time(name) for name in arc_filenames]

        # Generate manifest file for each science file
        print(f'Found {len(sci_filenames)} science files. '
              f'Attempting to find calibrations files withing time_radius of each ...')
        for sci_file in sci_filenames:
            fail_bias = False
            fail_flt = False
            fail_arc = False
            # For each science file, find the calibration files within a time radius
            obs_time_mjd = io.get_fits_file_time(sci_file)
            time_center = obs_time_mjd
            bias_files = io.get_elements_within_time_radius(bias_filenames, bias_times, time_center, time_radius)
            flt_files = io.get_elements_within_time_radius(flt_filenames, flt_times, time_center, time_radius)
            arc_files = io.get_elements_within_time_radius(arc_filenames, arc_times, time_center, time_radius)

            if len(bias_files) == 0:
                fail_bias = True
            if len(flt_files) == 0:
                fail_flt = True
            if len(arc_files) == 0:
                fail_arc = True

            if fail_bias or fail_arc or fail_flt:
                print(f'FAIL: found ZERO calibration files for sci_file={sci_file}: '
                      f'time_center={time_center}, time_radius={time_radius}')
                continue

            num_cals = len(bias_files) + len(flt_files) + len(arc_files)
            print(f'PASS: found {num_cals} calibration files for sci_file={sci_file}: '
                  f'time_center={time_center}, time_radius={time_radius}')

            record = dict()
            now_string = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
            record['reduction_name'] = f'antigen_manifest_{now_string}'
            record['unit_date'] = 'unknown'
            record['unit_instrument'] = 'VIRUS2'
            record['unit_id'] = unit
            record['obs_date'] = obs_date
            record['obs_name'] = obs_name
            record['in_folder'] = './'
            record['observation_files'] = sci_filenames
            record['calibration_files'] = dict()
            record['calibration_files']['bias'] = bias_filenames
            record['calibration_files']['flat'] = flt_filenames
            record['calibration_files']['arc'] = arc_filenames
            manifest_records.append(record)

    return manifest_records


def process_unit(manifest_record, output_path, debug=False):
    """
    Purpose: data reduction pipeline to process VIRUS2 observation files

    Args:
        manifest_record (dict): dict returned by yaml loading full-path filename to manifest.yaml containing lists of calibration files, etc
        output_path (str): Path where reduction output files will be written
        debug (bool): default=False, if True, save PNG plot files for intermediate data artifacts

    Returns:
        reduction_filename (str): full-path filename of FITS file written herein, containing obs file data reduction
    """
    log = utils.setup_logging('virus2_reductions')

    os.makedirs(output_path, exist_ok=True)

    CONFIG_DETECTOR, CONFIG_DEF_WAVE = config.get_channel_config_virus2()

    # =========================================
    # Get def_wave and line_list for given channel
    # TODO: implemented the other three channels, using green as a template
    # =========================================
    unit_id = manifest_record['unit_id']
    channel = unit_id[-1].lower()  # for VIRUS2, the unit is 'D3G', so 'D3G'[-1] returns 'G'

    line_list = None
    limit = None
    def_wave = None

    if channel == 'b':
        def_wave = CONFIG_DEF_WAVE[channel]
    if channel == 'g':
        def_wave = CONFIG_DEF_WAVE[channel]
        line_list_filepath = config.get_config_filepath('virus2', f'{unit_id}', 'virus2_lines_green.txt')
        line_list = Table.read(line_list_filepath, format="ascii")
        limit = CONFIG_DETECTOR[channel]['limit']
    if channel == 'r':
        def_wave = CONFIG_DEF_WAVE[channel]
    if channel == 'd':
        def_wave = CONFIG_DEF_WAVE[channel]

    lines = np.array(line_list['wavelength'])
    xref = np.array(line_list['col'])
    use_kernel = True

    # =============================================================================
    # Make a master bias, master dome flat, and master arc for the first set of OBS
    # Use the filename obs_id numbers for grouping/splitting contiguous blocks of observations files
    # =============================================================================

    bias_filenames = manifest_record['calibration_files']['bias']
    flat_filenames = manifest_record['calibration_files']['flat']
    arc_filenames = manifest_record['calibration_files']['arc']

    log.info('Making master bias frames')
    master_bias_data, master_bias_time = ccd.make_master_cal(bias_filenames, channel)

    log.info('Making master flat frames')
    master_flat_data, master_flat_time = ccd.make_master_cal(flat_filenames, channel)

    log.info('Making master arc frames')
    master_arc_data, master_arc_time = ccd.make_master_cal(arc_filenames, channel)

    # =============================================================================
    # Load reference fiber locations from a predefined file
    # =============================================================================

    ifu_cen_filepath = config.get_config_filepath('virus2', f'{unit_id}', f'virus2_ifucen_{unit_id}.txt')
    ifu_cen_file_data = Table.read(ifu_cen_filepath, format="ascii")
    ref = ifu_cen_file_data
    ref.reverse()

    # =============================================================================
    # Get trace from the dome flat
    # =============================================================================
    log.info('Getting trace for each master flat')

    trace, good_fiber_mask, Tchunk, xchunk = fiber.get_trace(master_flat_data - master_bias_data, ref)
    _, _ = plot.plot_trace(trace, Tchunk, xchunk, outfolder=output_path)

    domeflat_spec = fiber.get_spectra(master_flat_data - master_bias_data, trace)
    domeflat_error = 0. * domeflat_spec

    # =============================================================================
    # Get wavelength from arc lamps
    # =============================================================================
    log.info('Getting wavelength for each master arc')

    lamp_spec = fiber.get_spectra(master_arc_data - master_bias_data, trace)

    # save lamp spec data to FITS and PNG
    if debug:
        lamp_spec_test_fits_filename = os.path.abspath(os.path.join(output_path, 'debug_lamp_spec.fits'))
        fits.PrimaryHDU(lamp_spec).writeto(lamp_spec_test_fits_filename, overwrite=True)

        lamp_spec_test_plot_filename = os.path.abspath(os.path.join(output_path, 'debug_lamp_spec.png'))
        plot.plot_frame(lamp_spec, save_file=lamp_spec_test_plot_filename, title='Lamp Spec')

    try:
        wavelength, res, X, W = fiber.get_wavelength(lamp_spec, trace, good_fiber_mask,
                                                     xref, lines, limit=limit,
                                                     use_kernel=use_kernel)

        # Plot wavelength solution for inspection
        plot.plot_wavelength(lines, W, wavelength, output_path)
    except:
        error_message = 'Could not get wavelength solution for arc_filenames from manifest'
        log.warning(error_message)
        raise RuntimeError(error_message)

    # =============================================================================
    # Rectify domeflat spectra and get fiber to fiber
    # =============================================================================
    log.info('Getting fiber to fiber for each master domeFlat')

    domeflat_rect, domeflat_error_rect = fiber.rectify(domeflat_spec, domeflat_error,
                                                       wavelength, def_wave)
    ftf, ftf_smooth = fiber.get_fiber_to_fiber(domeflat_rect)

    # =============================================================================
    # Reduce!
    # =============================================================================
    arc_filename = arc_filenames[0]
    log.info(f'Reducing Arc Frame to generate PCS model: {arc_filename}')

    pca, _, _, _ = reduce(arc_filename, master_bias_data, master_flat_data,
                       trace, good_fiber_mask, wavelength, ftf, channel,
                       pca=None, pca_only=True, outfolder=output_path, debug=debug)

    science_file = manifest_record['observation_files'][0]  # TODO: process ALL files from the obs list, not just [0]
    log.info(f'Reducing Science Frame: {science_file}')

    _, sky, cont, reduction_filename = reduce(science_file, master_bias_data, master_flat_data,
                                              trace, good_fiber_mask, wavelength, ftf, channel,
                                              pca=pca, outfolder=output_path, debug=debug)

    return reduction_filename


def process(infolder, outfolder, obs_date, obs_name, reduce_all, time_radius,
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
    log = utils.setup_logging('virus2_reductions')

    manifest_records = build_manifest_records(infolder, obs_date, obs_name, reduce_all, time_radius,
                                              bias_label, arc_label, dark_label, flat_label, twilight_flat_label)

    os.makedirs(outfolder, exist_ok=True)

    for record in manifest_records:
        unit_id = record['unit_id']
        log.info(f'Processing reduction for unit = {unit_id}:')
        try:
            output_fits_filename = process_unit(record, outfolder)
            log.error(f'Processing reduction for unit = {unit_id}: PASS: wrote reduction to FITS file {output_fits_filename}')
        except Exception as error:
            log.error(f'Processing reduction for unit = {unit_id}: FAILED: {error}')

    return None