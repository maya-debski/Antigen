import os
from pathlib import Path
import logging
import warnings

from astropy.io import fits
from astropy.stats import biweight_location as biweight
import numpy as np

from antigen import ccd
from antigen import config
from antigen import fiber
from antigen import io
from antigen import plot
from antigen import sky
from antigen import spectra
from antigen import trace
from antigen import wavelength

# Turn off annoying warnings (even though some deserve attention)
warnings.filterwarnings("ignore")
logger = logging.getLogger('antigen.reduce')


def reduce_science(data_filename, master_bias, master_flat, trace_array, good_fiber_mask,
                   wavelength_cal, ftf_correction,
                   config_dict, pca=None, pca_only=False, outfolder=None):
    """
    Purpose: Reduce the raw data by performing a series of processing steps,
    including bias subtraction, flat-fielding, sky subtraction,
    and PCA-based residuals analysis.

    Args:
        data_filename (str): The filename of the FITS file containing the data to be reduced.
        master_bias (): master bias frames for bias correction.
        master_flat (): master flat frame
        trace_array (): fiber trace data for current observation
        good_fiber_mask (array(bool)): boolean selecttion mask of good (isfinite) fibers for current observation
        wavelength_cal (): wavelength calibration data for the current observation
        ftf_correction (list): flat-field (fiber_to_fiber) correction array.
        config_dict (dict): configuration dictionary for fits file
        pca (sklearn.decomposition.PCA), Default is None, optional, A pre-fitted PCA model for residual map analysis.
        pca_only (bool): if True, return immediately after computing the PCA model fit object
        outfolder (str): file output path to write FITS files to

    Returns:
        pca (sklearn.decomposition.PCA): The fitted PCA model
        biweighted_spectrum (np.ndarray): bi-weighted spectrum, statistic for determining the
                                          central location of the specrect distribution
        continuum (np.ndarray): 1d numpy array, The computed continuum for the spectrum.
    """


    obs_data, obs_header = io.load_fits(data_filename)

    # Perform basic image reduction (bias subtraction, gain adjustment)
    image, E = ccd.base_reduction(obs_data, master_bias, config_dict)

    # Scattered light correction
    scattered_light = ccd.get_scattered_light(image, trace_array)
    image -= scattered_light

    # Extract spectra from the image using the trace data
    spec = spectra.get_spectra(image, trace_array)

    # Calculate the spectrum error using the flat-field and error image
    specerr = spectra.get_spectra_error(E, trace_array)

    # Compute the chi-square of the spectrum to identify bad pixels
    chi2 = spectra.get_spectra_chi2(master_flat - master_bias, image, E, trace_array)
    badpix = chi2 > 20.  # Pixels with chi2 > 20 are considered bad
    specerr[badpix] = np.nan
    spec[badpix] = np.nan

    # Rectify the spectrum and error based on the wavelength
    def_wave = wavelength.get_rectified_wavelength(config_dict)

    specrect, errrect = spectra.rectify(spec, specerr, wavelength_cal, def_wave)

    # Apply flat-field correction
    obs_exp_time = obs_header['EXPTIME']
    specrect = spectra.convert_spectral_units(specrect, def_wave, obs_exp_time)
    errrect = spectra.convert_spectral_units(errrect, def_wave, obs_exp_time)

    specrect[:] /= (ftf_correction)
    errrect[:] /= (ftf_correction)

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
    output_fits_filename = io.write_fits(skysubrect_adv, skysubrect, specrect, errrect, obs_header, config_dict,
                                         outfolder)

    # Return the biweighted spectrum and continuum
    biweighted_spectrum = biweight(specrect, axis=0, ignore_nan=True)
    return pca, biweighted_spectrum, continuum, output_fits_filename


def process_calibration(manifest_record, output_path, config_dict):
    """
    Purpose: Process calibration files needed for data reduction for VIRUS2 observation files

    Args:
        manifest_record (dict): dict returned by yaml loading full-path filename to manifest.yaml
                                containing lists of calibration files, etc
        output_path (str): Path where reduction output files will be written
        config_dict (dict): Dictionary of configuration parameters

    Returns:
        master_bias (arraylike): master bias frames for bias correction.
        master_flat (arraylike): master flat frame
        master_arc (arraylike): master arc frame
        trace_array (arraylike): fiber trace data for current observation
        good_fiber_mask (array(bool)): boolean selection mask of good (isfinite) fibers for current observation
        wavelength (arraylike): wavelength calibration data for the current observation
        ftf (list): flat-field (fiber_to_fiber) correction array.
    """

    os.makedirs(output_path, exist_ok=True)

    def_wave = np.linspace(config_dict['start_wavelength'], config_dict['end_wavelength'],
                           config_dict['detector_dimensions']['X'])

    lines = config_dict['wavelength']
    xref = config_dict['column']
    peak_threshold = config_dict['arc_flux_limit']
    fiber_ref = config_dict['reference_fiber_index']
    use_kernel = True
    #TODO: Decide how to do the VIRUS2 ifucen file
    trace_rows = config_dict['trace_row']
    exclude_fiber = config_dict['exclude_fiber']

    # =============================================================================
    # Make a master bias, master dome flat, and master arc for the first set of OBS
    # Use the filename obs_id numbers for grouping/splitting contiguous blocks of observations files
    # =============================================================================

    bias_filenames = manifest_record['calibration_files']['bias']
    flat_filenames = manifest_record['calibration_files']['flat']
    arc_filenames = manifest_record['calibration_files']['arc']

    logger.info('Making master bias frames')
    master_bias_data, master_bias_time = ccd.make_master_cal(bias_filenames, config_dict)
    fits.PrimaryHDU(master_bias_data).writeto(Path(output_path) / 'master_bias.fits', overwrite=True)

    logger.info('Making master flat frames')
    master_flat_data, master_flat_time = ccd.make_master_cal(flat_filenames, config_dict)
    fits.PrimaryHDU(master_flat_data).writeto(Path(output_path) / 'master_flat.fits', overwrite=True)

    logger.info('Making master arc frames')
    master_arc_data, master_arc_time = ccd.make_master_cal(arc_filenames, config_dict)
    fits.PrimaryHDU(master_arc_data).writeto(Path(output_path) / 'master_arc.fits', overwrite=True)

    # =============================================================================
    # Get trace from the dome flat
    # =============================================================================
    logger.info('Getting trace for each master flat')
    trace_array, good_fiber_mask, raw_trace_matrix, x_chunk_centers = trace.get_trace(master_flat_data - master_bias_data,
                                                                                trace_rows, exclude_fiber)
    _, _ = plot.plot_trace(trace_array, raw_trace_matrix, x_chunk_centers,
                           fiber_indices=config_dict['sample_fiber_indices'],
                           outfolder=output_path)

    domeflat_spec = spectra.get_spectra(master_flat_data - master_bias_data, trace_array)
    domeflat_error = 0. * domeflat_spec

    # =============================================================================
    # Get wavelength from arc lamps
    # =============================================================================
    logger.info('Getting wavelength for each master arc')

    lamp_spec = spectra.get_spectra(master_arc_data - master_bias_data, trace_array)

    # save lamp spec data to FITS and PNG
    lamp_spec_test_fits_filename = os.path.abspath(os.path.join(output_path, 'lamp_spec.fits'))
    fits.PrimaryHDU(lamp_spec).writeto(lamp_spec_test_fits_filename, overwrite=True)
    lamp_spec_test_plot_filename = os.path.abspath(os.path.join(output_path, 'lamp_spec.png'))
    plot.plot_frame(lamp_spec, save_file=lamp_spec_test_plot_filename, title='Lamp Spec')

    #TODO: Need to test if binned in a better way.
    if config_dict['detector_dimensions']['X'] < config_dict['detector_dimensions']['Y']:
        xref = xref / 2

    # TODO: peak_threshold is now a function of the noise in the arc spectrum
    wavelength_array, res, X, W = wavelength.get_wavelength(lamp_spec, trace_array, good_fiber_mask,
                                                 xref, lines, peak_threshold=peak_threshold,
                                                 reference_fiber_index=fiber_ref,
                                                 use_kernel=use_kernel)

    # Plot wavelength solution for inspection
    plot.plot_wavelength(lines, W, wavelength_array, output_path)
    #except:
    #    error_message = 'Could not get wavelength solution for arc_filenames from manifest'
    #    raise RuntimeError(error_message)

    # =============================================================================
    # Rectify domeflat spectra and get fiber to fiber
    # =============================================================================
    logger.info('Getting fiber to fiber for each master domeFlat')

    domeflat_rect, domeflat_error_rect = spectra.rectify(domeflat_spec, domeflat_error,
                                                       wavelength_array, def_wave)
    ftf, ftf_smooth = fiber.get_fiber_to_fiber(domeflat_rect)

    return master_bias_data, master_flat_data, master_arc_data, trace_array, good_fiber_mask, wavelength_array, ftf


def reduction_pipeline(dataset_manifest, output_path):
    """
    Purpose: data reduction pipeline to process VIRUS2 observation files

    Args:
        dataset_manifest (dict): dataset manifest dictionary returned from e.g. dataset.find_datasets()
        output_path (str): output file path to which this method will write a reduced FITS file

    Returns:
        reduction_filename (str): full-path filename of FITS file written herein, containing obs file data reduction
    """
    # TODO: config file naming conventions are lower case instrument and upper case element.
    config_dict = config.build_config_for_element(dataset_manifest['unit_instrument'].lower(),
                                                  dataset_manifest['unit_id'].upper())

    # TODO: replace tuple unpack with a more intentional data structure, e.g. dict, namedtuple, data-class, etc
    logger.info(f'Processing calibration for reduction.')
    calibration_tuple = process_calibration(dataset_manifest, output_path, config_dict)
    (master_bias_data,
     master_flat_data,
     master_arc_data,
     trace_array, good_fiber_mask,
     wavelength_array, ftf) = calibration_tuple


    arc_file = dataset_manifest['calibration_files']['arc'][0]
    logger.info(f'Reducing Arc Frame to generate PCA model: arc_file={arc_file}')
    pca, _, _, _ = reduce_science(arc_file, master_bias_data, master_flat_data,
                                  trace_array, good_fiber_mask, wavelength_array, ftf, config_dict,
                                  pca=None, pca_only=True, outfolder=output_path)

    for science_file in dataset_manifest['observation_files']:
        logger.info(f'Reducing Science Frame: science_file={science_file}')
        _, sky, cont, reduction_filename = reduce_science(science_file, master_bias_data, master_flat_data,
                                                          trace_array, good_fiber_mask, wavelength_array, ftf, config_dict,
                                                          pca=pca, outfolder=output_path)
        logger.info(f'Wrote reduction to FITS file {reduction_filename}')

    return reduction_filename
