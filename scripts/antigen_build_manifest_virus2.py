#!/usr/bin/env python

import argparse
import os
import sys

from antigen.cli import add_calibration_args, add_common_args
from antigen.datasets import find_datasets
from antigen.manifest import save_manifest
from antigen.utils import setup_logging

DESCRIPTION = r"""
Purpose:
    This script builds a manifest file for the VIRUS2 instrument for a given night of observations.
    The manifest serves as an index of all relevant datasets found in the specified input folder,
    grouped and labeled according to observation date, target name, and calibration frames.
    It is intended to assist in organizing and tracking raw data before further reduction.

What it does:
    - Searches a given input directory for raw VIRUS2 datasets.
    - Groups exposures by science target and calibration type (bias, dark, arc, flat, twilight flat).
    - Writes one or more manifest files in YAML format to the specified output directory.
    - Each manifest file includes paths and metadata for the matched datasets.

Inputs:
    - Input folder: Root directory containing raw FITS files.
    - Output folder: Destination for generated manifest YAML files.
    - Observation date: Required to identify datasets from a specific night.
    - Optional: Observation name, time radius, and custom labels for calibration frames.

Outputs:
    - One or more YAML manifest files named like:
      manifest_<obs_date>_<obs_name>_record<N>.yml

Examples:
    Build a manifest for the night of 20240801 for object M57:

    $ antigen_build_manifest_virus2.py -i /data/virus2/raw/ -o /data/virus2/manifests/ -c 20240801 -n M57

Notes:
    - You can use the --reduce-all flag to include all datasets regardless of name match.
    - The manifest files are used as input for downstream reduction pipelines.
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

    dataset_manifests = find_datasets(args.in_folder, args.obs_date,
                                      args.obs_name, args.reduce_all, args.time_radius,
                                      args.bias_label, args.arc_label, args.dark_label,
                                      args.flat_label, args.twilight_flat_label
                                      )

    for nr, record in enumerate(dataset_manifests):
        manifest_filename = f'manifest_{args.obs_date}_{args.obs_name}_record{nr}.yml'
        save_path = os.path.abspath(args.out_folder)
        os.makedirs(save_path, exist_ok=True)
        save_filepath = os.path.join(save_path, manifest_filename)
        save_manifest(record, save_filepath)
        logger.info(f'Application completed: Wrote manifest file: {save_filepath}')


    return None


if __name__ == '__main__':
    sys.exit(main())
