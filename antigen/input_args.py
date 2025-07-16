#!/usr/bin/env python

import argparse as ap
import datetime
import os


def get_args():

    defaults = {
        'infolder' : os.curdir,
        'outfolder': os.curdir,
        'date': datetime.datetime.now().strftime('%Y%m%d'),
        'name': None,
        'reduce_all': False,
        'bias_label': 'bias',
        'arc_label': 'arc',
        'dark_label': 'dark',
        'flat_label': 'flat',
        'twilight_flat_label': 'twi',
    }

    helps = {
        'infolder' : 'Root path where reduction input file tree is located, (default: %(default)s)',
        'outfolder' : 'Path where reduction output files will be written, (default: %(default)s)',
        'date': 'Observation calendar date string formatted as YYYYMMDD, ex: 20250613, (default: %(default)s)',
        'name': 'Observation object/target name, e.g. from FITS header card, (default: %(default)s)',
        'reduce_all': 'Reduce all files found under infolder file tree, (default: %(default)s)',
        'bias_label': 'The object name from the FITS header card for bias files, (default: %(default)s)',
        'arc_label': 'The object name from the FITS header card for arc files, (default: %(default)s)',
        'dark_label': 'The object name from the FITS header card for dark files, (default: %(default)s)',
        'flat_label': 'The object name from the FITS header card for flat files, (default: %(default)s)',
        'twilight_flat_label': 'The object name from the FITS header card for twilight flat files, (default: %(default)s)',
    }

    parser = ap.ArgumentParser(add_help=True)

    parser.add_argument('-i', '--infolder', type=str, help=helps['infolder'], default=defaults['infolder'])
    parser.add_argument('-o', '--outfolder', type=str, help=helps['outfolder'], default=defaults['outfolder'])
    parser.add_argument('-c', '--date', type=str, help=helps['date'], default=defaults['date'])
    parser.add_argument('-n', '--name', type=str, help=helps['name'], default=defaults['name'])
    parser.add_argument('-r', '--reduce_all', action='store_true', help=helps['reduce_all'], default=defaults['reduce_all'])
    parser.add_argument('-b', '--bias', type=str, help=helps['bias_label'], default=defaults['bias_label'])
    parser.add_argument('-a', '--arc', type=str, help=helps['arc_label'], default=defaults['arc_label'])
    parser.add_argument('-d', '--dark', type=str, help=helps['dark_label'], default=defaults['dark_label'])
    parser.add_argument('-f', '--flat', type=str, help=helps['flat_label'], default=defaults['flat_label'])
    parser.add_argument("-t", "--twilight_flat_label", type=str, help=helps['twilight_flat_label'], default=defaults['twilight_flat_label'])
    argv = None
    args = parser.parse_args(args=argv)

    return args

