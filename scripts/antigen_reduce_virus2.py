#!/usr/bin/env python

import argparse
import os
import sys

from antigen.cli import add_calibration_args, add_common_args
from antigen.datasets import find_datasets
from antigen.reduce_virus2 import reduction_pipeline
from antigen.utils import setup_logging

DESCRIPTION = r"""
Purpose:
    This script runs the full reduction pipeline for VIRUS2 instrument datasets
    for a given night of observations. It uses manifest files to locate input
    science and calibration data, applies the standard reduction steps, and
    writes reduced science-ready FITS files to the specified output directory.

What it does:
    - Searches the input folder for datasets matching the given observation date,
      observation name, and calibration frame labels.
    - Groups raw files into logical datasets for reduction.
    - Runs the VIRUS2 reduction pipeline on each dataset, which includes:
        * Bias correction
        * Dark subtraction
        * Flat fielding
        * Wavelength calibration
        * Science frame spectral extraction
    - Writes the reduced, calibrated science frames to the output folder.

Inputs:
    - Input folder: Directory containing raw VIRUS2 FITS files.
    - Output folder: Destination for reduced FITS files.
    - Observation date: Date of the observation run (YYYYMMDD format).
    - Optional: Observation name, time radius, and custom labels for calibration frames.
    - Other common flags: e.g., --reduce-all to include all matching datasets.

Outputs:
    - One reduced, science-ready FITS file per dataset found.
    - Files are saved in the output folder with names indicating unit ID and processing details.

Example:
    Reduce all VIRUS2 data for the night of 20250801:

    $ antigen_reduce_virus2.py -i /data/virus2/raw/ -o /data/virus2/reduced/ -c 20250801

Notes:
    - The script logs all processing steps and reports any datasets that failed to reduce.
    - Reduction uses calibration frames as specified by bias, dark, arc, flat, and twilight labels.
    - This script assumes a valid configuration of the VIRUS2 reduction pipeline.
"""



def get_args():
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_args(parser)
    add_calibration_args(parser)
    return parser.parse_args()


def main():
    args = get_args()

    logger = setup_logging('antigen', verbose=args.verbose)
    logger.info(f'Starting application...')

    os.makedirs(args.out_folder, exist_ok=True)

    dataset_manifests = find_datasets(args.in_folder, args.obs_date, args.obs_name,
                                      args.reduce_all, args.time_radius,
                                      args.bias_label, args.arc_label, args.dark_label,
                                      args.flat_label, args.twilight_flat_label)

    logger.info(f'Found {len(dataset_manifests)} datasets to reduce.')

    save_files = []
    for manifest in dataset_manifests:
        unit_id = manifest['unit_id']
        sci_file = manifest['observation_files'][0]
        logger.info(f'Processing reduction for unit={unit_id}, science_file={sci_file}')
        try:
            output_fits_filename = reduction_pipeline(manifest, args.out_folder)
            save_files.append(output_fits_filename)
            logger.info(f'Processing reduction for unit = {unit_id}: PASS: wrote reduction to FITS file {output_fits_filename}')
        except Exception as error:
            logger.error(f'Processing reduction for unit = {unit_id}: FAILED: {error}')

    logger.info(f'Application completed: Completed reduction and FITS save for {len(save_files)} out of {len(dataset_manifests)} datasets found.')

    return None


if __name__ == '__main__':
    sys.exit(main())
