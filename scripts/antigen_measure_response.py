#!/usr/bin/env python

import argparse
from pathlib import Path
import sys

from antigen.datasets import build_dataset_from_reduced_files
from antigen.manifest import save_manifest
from antigen.utils import setup_logging
from antigen.calibrate import build_response

DESCRIPTION = r"""
Purpose:
    This script builds an instrument response function from reduced
    standard star observations for the VIRUS2/GCMS instrument.

What it does:
    - Scans a directory of reduced datasets for standard star frames.
    - Generates a dataset manifest for each candidate standard star.
    - Selects the datasets matching the requested standard star name.
    - For each match:
        * Builds an instrument response function by comparing the observed
          spectrum with a reference CALSPEC spectrum.
        * Saves the response manifest and response products to the output folder.

Inputs:
    - Reduced directory: Path containing reduced standard star datasets
      (already processed with bias, flats, arcs, etc.).
    - Output folder: Destination for response files (defaults to reduced directory).
    - Standard star name: Substring to identify the target (e.g., "Feige").
    - PSF and cube-building settings such as extraction radius and pixel size.

Outputs:
    - A YAML manifest file for each standard star response dataset.
    - Instrument response function FITS and diagnostic plots, written to the output folder.

Example:
    Build a response function from reduced Feige standard star frames:

    $ antigen_build_response.py -r /data/VIRUS2/reduced/ -o /data/VIRUS2/response/ -s Feige

Notes:
    - This script assumes that the input data have already been reduced
      (use the main reduction pipeline first).
    - Logs will report any datasets that do not match the requested standard.
    - Response functions are required for flux calibration of science targets.
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
            output_response_filename = build_response(manifest, args)
            logger.info(f'Saved response to {output_response_filename}')

    return None


if __name__ == '__main__':
    sys.exit(main())
