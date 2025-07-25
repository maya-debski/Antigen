#!/usr/bin/env python

from argparse import ArgumentParser
import datetime
import os
import sys

from antigen import manifest
from antigen import reduce_virus2
from antigen import utils


def get_args():
    """
    Purpose: Input arg handling for VIRUS2 data reduction CLI tools, e.g. antigen_reduce_virus2.py
    """

    defaults = {
        'manifest_filename' : os.path.abspath(os.curdir),
        'out_folder': datetime.datetime.now().strftime('antigen_reduce_virus2_%Y%m%d_%H%M%S'),
        'verbose': False,
        'debug': False,
    }

    helps = {
        'manifest_filename' : 'Full path filename of YAML manifest containing params and lists of calibration files, (default: %(default)s)',
        'out_folder' : 'Path where reduction output files will be written, (default: %(default)s)',
        'verbose': 'if True, validate manifest file and print process details to console, (default: %(default)s)',
        'debug': 'if True, generate and save PNG plot files of intermediate data artifacts, (default: %(default)s)',
    }

    parser = ArgumentParser(add_help=True)

    parser.add_argument('-m', '--manifest_filename', type=str, help=helps['manifest_filename'], default=defaults['manifest_filename'])
    parser.add_argument('-o', '--out_folder', type=str, help=helps['out_folder'], default=defaults['out_folder'])
    parser.add_argument('-v', '--verbose', action='store_true', help=helps['verbose'], default=defaults['verbose'])
    parser.add_argument('-d', '--debug', action='store_true', help=helps['debug'], default=defaults['debug'])

    args = parser.parse_args()

    return args


def main():


    args = get_args()

    logger = utils.setup_logging('antigen', verbose=args.verbose, debug=args.debug)
    logger.info(f'Starting application...')

    if args.out_folder:
        if not os.path.isdir(args.out_folder):
            os.makedirs(args.out_folder)

    logger.info(f'Reading dataset manifest for reduction: {args.manifest_filename}')
    manifest_record = manifest.read_manifest(args.manifest_filename, validate=args.verbose, verbose=args.verbose)
    manifest_savefile = os.path.join(os.path.abspath(args.out_folder), 'reduction_manifest.yml')
    manifest.save_manifest(manifest_record, manifest_savefile)

    obs_filenames = manifest_record[manifest.MANIFEST_OBSERVATION_LIST_KEY]
    for file in obs_filenames:
        logger.info(f'Processing reduction for FITS obs file: {file.name}')
        output_fits_filename = reduce_virus2.reduction_pipeline(manifest_record, output_path=args.out_folder)

        logger.info(f'Wrote reduction FITS file to {output_fits_filename}')

    return None


if __name__ == '__main__':
    sys.exit(main())
