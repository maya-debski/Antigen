#!/usr/bin/env python

import argparse
from pathlib import Path
import sys

from antigen.cli import add_calibration_args, add_common_args
from antigen.datasets import find_datasets
from antigen.manifest import save_manifest
from antigen.reduce import reduction_pipeline
from antigen.utils import setup_logging

DESCRIPTION = r"""
Purpose:
    This script runs the full reduction pipeline for VIRUS2/GCMS instrument datasets
    for a given night of observations. It uses manifest files to locate input
    science and calibration data, applies the standard reduction steps, and
    writes reduced science-ready FITS files to the specified output directory.

What it does:
    - Searches the input folder for datasets matching the given observation date,
      observation name, and calibration frame labels.
    - Groups raw files into logical datasets for reduction.
    - Runs the VIRUS2 reduction pipeline on each dataset, which includes:
        * Bias correction
        * Trace calibration
        * Flat fielding
        * Wavelength calibration
        * Science frame spectral extraction
    - Writes the reduced, calibrated science frames to the output folder.

Inputs:
    - Input folder: Directory containing raw spectrograph FITS files
    - Output folder: Destination for reduced FITS files.
    - Observation date: Date of the observation run (YYYYMMDD format).
    - Optional: Observation name, time radius, and custom labels for calibration frames.
    - Other common flags: e.g., --reduce-all to include all matching datasets.

Outputs:
    - One reduced, science-ready FITS file per dataset found.
    - Files are saved in the output folder with names indicating unit ID and processing details.

Example:
    Reduce all VIRUS2 data for the night of 20250801:

    $ antigen_reduce.py -i ~/data -o ~/data/VIRUS/reduced -c 20250801 -m VIRUS2

Notes:
    - The script logs all processing steps and reports any datasets that failed to reduce.
    - Reduction uses calibration frames as specified by bias, dark, arc, flat, and twilight labels.
    - This script assumes a valid configuration of the VIRUS2 reduction pipeline.
"""



def get_args():
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_args(parser)
    add_calibration_args(parser)

    parser.add_argument(
        '-m', '--instrument',
        type=str,
        default='GCMS',
        help='Name of the instrument used for the observation '
             '(default: %(default)s). Example: GCMS.'
    )

    return parser.parse_args()


def main():
    args = get_args()

    logger = setup_logging('antigen', verbose=args.verbose)
    logger.info(f'Starting application...')

    save_path = Path(args.out_folder).expanduser().resolve()
    save_path.mkdir(parents=True, exist_ok=True)

    dataset_manifests = find_datasets(args.in_folder, args.obs_date, args.obs_name, args.time_radius,
                                      args.bias_label, args.arc_label, args.dark_label,
                                      args.flat_label, args.twilight_flat_label, instrument=args.instrument)

    logger.info(f'Found {len(dataset_manifests)} datasets to reduce.')

    save_files = []
    for nr, manifest in enumerate(dataset_manifests):
        reduction_name = manifest['reduction_name']
        manifest_filename = f'manifest_{args.obs_date}_{args.obs_name}_record{nr}.yml'
        save_filepath = save_path / manifest_filename
        save_manifest(manifest, str(save_filepath))
        element = manifest['unit_id']
        instrument = manifest['unit_instrument']
        logger.info(f'Processing {reduction_name} for unit={element}')
        output_fits_filename = reduction_pipeline(manifest, args.out_folder)
        save_files.append(output_fits_filename)

    logger.info(f'Application completed: Completed reduction and FITS save for {len(save_files)} out of {len(dataset_manifests)} datasets found.')

    return None


if __name__ == '__main__':
    sys.exit(main())
