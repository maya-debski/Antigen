#!/usr/bin/env python

from argparse import ArgumentParser
import datetime
import os
import sys

from antigen.config import get_config_filepath
from antigen.manifest import read_manifest, save_manifest
from antigen.utils import setup_logging


def get_args():
    """
    Purpose: Input arg handling for VIRUS2 data reduction CLI tools, e.g. antigen_reduce_virus2.py
    """
    defaults = {
        'filename': None,
        'validate': False,
        'example': False,
    }

    helps = {
        'filename': 'Name of yaml file containing reduction input file and parameter manifest, (default: %(default)s)',
        'validate': 'if True, validate manifest file structure, test existence of listed files, (default: %(default)s)',
        'example': 'if True, SKIP --filename, and instead print an EXAMPLE manifest file, (default: %(default)s)'
    }

    parser = ArgumentParser(add_help=True)
    parser.add_argument('-f', '--filename', type=str, help=helps['filename'], default=defaults['filename'])
    parser.add_argument('-v', '--validate', action='store_true', help=helps['validate'], default=defaults['validate'])
    parser.add_argument('-e', '--example', action='store_true', help=helps['example'], default=defaults['example'])
    args = parser.parse_args()
    return args


def print_example_manifest():
    example_file = get_config_filepath('virus2', 'virus2_manifest_template.yml')
    print('*** EXAMPLE MANIFEST, INPUT --filename IGNORED ***')
    example_manifest_text = example_file.read_text()
    print(example_manifest_text)
    print('*** EXAMPLE MANIFEST, INPUT --filename IGNORED ***')
    return example_manifest_text


def main():
    args = get_args()

    logger = setup_logging('antigen', debug=args.validate, verbose=True)
    logger.info(f'Starting application...')

    if args.example:
        print_example_manifest()
        sys.exit(0)
    else:
        manifest = read_manifest(filename=args.filename, validate=args.validate, verbose=True)

        # save manifest dict back to disk to test/validate
        manifest_save_name_stem = datetime.datetime.now().strftime('antigen_read_manifest_save_%Y%m%d_%H%M%S')
        manifest_savefile = os.path.join(os.path.abspath(os.curdir), f'{manifest_save_name_stem}.yml')
        save_manifest(manifest, manifest_savefile)

    logger.info(f'Completed application.')

    return None


if __name__ == '__main__':
    sys.exit(main())
