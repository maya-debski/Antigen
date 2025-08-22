#!/usr/bin/env python

import argparse
from pathlib import Path
import sys

from antigen.cli import add_calibration_args, add_common_args
from antigen.datasets import build_dataset_from_reduced_files
from antigen.manifest import save_manifest
from antigen.utils import setup_logging
from antigen.calibrate import build_response

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

    $ antigen_reduce_gcms.py -i /data/virus2/raw/ -o /data/virus2/reduced/ -c 20250801

Notes:
    - The script logs all processing steps and reports any datasets that failed to reduce.
    - Reduction uses calibration frames as specified by bias, dark, arc, flat, and twilight labels.
    - This script assumes a valid configuration of the VIRUS2 reduction pipeline.
"""



def get_args():
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-r', '--reduced_dir', required=True,
                        help='Path to directory with reduced standard star frames.')
    parser.add_argument('-o', '--output_folder', default=None,
                        help='Path to output folder for response function.')
    parser.add_argument('-s', '--standard_name', required=True,
                        help='Name for standard star files in reduced-dir (ex: Feige).')

    # PSF/extraction settings
    parser.add_argument('-e', '--extraction_radius', type=float, default=10.0,
                        help='Extraction radius in arcsec.')
    parser.add_argument('-p', '--pixel_size', type=float, default=1.0,
                        help='Cube pixel size in arcsec.')
    parser.add_argument('-v', '--verbose', action='store_true',
                   help='if True, print more process details and logger.info to console')
    return parser.parse_args()


def main():
    args = get_args()

    logger = setup_logging('antigen', verbose=args.verbose)
    logger.info(f'Starting application...')

    if args.output_folder is None:
        args.output_folder = args.reduced_dir
    save_path = Path(args.output_folder).expanduser().resolve()
    save_path.mkdir(parents=True, exist_ok=True)

    # Build response manifests
    response_manifest = build_dataset_from_reduced_files(args.reduced_dir)

    for nr, manifest in enumerate(response_manifest):
        reduction_name = manifest['reduction_name']
        element = manifest['unit_id']
        instrument = manifest['unit_instrument']
        target_name = manifest['target']
        logger.info(f'Manifest for {target_name} with instrument={instrument} and unit={element}')
        if args.standard_name.lower() in target_name.lower():
            manifest_filename = f'manifest_response_{target_name}.yml'
            save_filepath = save_path / manifest_filename
            save_manifest(manifest, str(save_filepath))
            logger.info(f'Processing response for {target_name} with instrument={instrument} and unit={element}')
            output_response_filename = build_response(manifest, args.output_folder, args.standard_name,
                                                      args.pixel_size)
            logger.info(f'Saved response to {output_response_filename}')

    return None


if __name__ == '__main__':
    sys.exit(main())
