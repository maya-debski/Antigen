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


def reduce_science(data_filename, master_bias, master_flat, trace, good_fiber_mask, wavelength_cal, ftf_correction,
                   channel, pca=None, pca_only=False, outfolder=None):
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

    # Return the biweighted spectrum and continuum
    biweighted_spectrum = biweight(specrect, axis=0, ignore_nan=True)
    return pca, biweighted_spectrum, continuum, output_fits_filename


def process_calibration(manifest_record, output_path):
    """
    Purpose: Process calibration files needed for data reduction for VIRUS2 observation files

    Args:
        manifest_record (dict): dict returned by yaml loading full-path filename to manifest.yaml
                                containing lists of calibration files, etc
        output_path (str): Path where reduction output files will be written

    Returns:
        master_bias (arraylike): master bias frames for bias correction.
        master_flat (arraylike): master flat frame
        master_arc (arraylike): master arc frame
        trace (arraylike): fiber trace data for current observation
        good_fiber_mask (array(bool)): boolean selection mask of good (isfinite) fibers for current observation
        wavelength (arraylike): wavelength calibration data for the current observation
        ftf (list): flat-field (fiber_to_fiber) correction array.
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
    lamp_spec_test_fits_filename = os.path.abspath(os.path.join(output_path, 'lamp_spec.fits'))
    fits.PrimaryHDU(lamp_spec).writeto(lamp_spec_test_fits_filename, overwrite=True)
    lamp_spec_test_plot_filename = os.path.abspath(os.path.join(output_path, 'lamp_spec.png'))
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

    return master_bias_data, master_flat_data, master_arc_data, trace, good_fiber_mask, wavelength, ftf


def reduction_pipeline(dataset_manifest, output_path):
    """
    Purpose: data reduction pipeline to process VIRUS2 observation files

    Args:
        dataset_manifest (dict): dataset manifest dictionary returned from e.g. dataset.find_datasets()
        output_path (str): output file path to which this method will write a reduced FITS file
    Returns:
        reduction_filename (str): full-path filename of FITS file written herein, containing obs file data reduction
    """

    log = utils.setup_logging('virus2_reductions')

    # TODO: replace tuple unpack with a more intentional data structure, e.g. dict, namedtuple, data-class, etc
    log.info(f'Processing calibration for reduction.')
    calibration_tuple = process_calibration(dataset_manifest, output_path)
    (master_bias_data,
     master_flat_data,
     master_arc_data,
     trace, good_fiber_mask,
     wavelength, ftf) = calibration_tuple

    arc_file = dataset_manifest['calibration_files']['arc'][0]
    channel = dataset_manifest['unit_id'][-1].lower()
    log.info(f'Reducing Arc Frame to generate PCA model: arc_file={arc_file}')
    pca, _, _, _ = reduce_science(arc_file, master_bias_data, master_flat_data,
                                  trace, good_fiber_mask, wavelength, ftf, channel,
                                  pca=None, pca_only=True, outfolder=output_path)

    science_file = dataset_manifest['observation_files'][0]
    log.info(f'Reducing Science Frame: science_file={science_file}')

    _, sky, cont, reduction_filename = reduce_science(science_file, master_bias_data, master_flat_data,
                                                      trace, good_fiber_mask, wavelength, ftf, channel,
                                                      pca=pca, outfolder=output_path)
    return reduction_filename
