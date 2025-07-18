from argparse import ArgumentParser
import datetime
import os
import sys
from datetime import datetime as dt


def get_args():
    """
    Purpose: Input arg handling for VIRUS2 data reduction CLI tools, e.g. antigen_reduce_virus2.py
    """

    defaults = {
        'in_folder' : os.curdir,
        'out_folder': datetime.datetime.now().strftime('antigen_reduce_virus2_%Y%m%d_%H%M%S'),
        'obs_date': datetime.datetime.now().strftime('%Y%m%d'),
        'obs_name': None,
        'reduce_all': False,
        'bias_label': 'bias',
        'arc_label': 'arc',
        'dark_label': 'dark',
        'flat_label': 'flat',
        'twilight_flat_label': 'twi',
    }

    helps = {
        'in_folder' : 'Root path where reduction input file tree is located, (default: %(default)s)',
        'out_folder' : 'Path where reduction output files will be written, (default: %(default)s)',
        'obs_date': 'Observation calendar date string formatted as YYYYMMDD, ex: 20250613, (default: %(default)s)',
        'obs_name': 'Observation object/target name, e.g. from FITS header card, (default: %(default)s)',
        'reduce_all': 'Reduce all files found under infolder file tree, (default: %(default)s)',
        'bias_label': 'The object name from the FITS header card for bias files, (default: %(default)s)',
        'arc_label': 'The object name from the FITS header card for arc files, (default: %(default)s)',
        'dark_label': 'The object name from the FITS header card for dark files, (default: %(default)s)',
        'flat_label': 'The object name from the FITS header card for flat files, (default: %(default)s)',
        'twilight_flat_label': 'The object name from the FITS header card for twilight flat files, (default: %(default)s)',
    }

    parser = ArgumentParser(add_help=True)

    parser.add_argument('-i', '--in_folder', type=str, help=helps['in_folder'], default=defaults['in_folder'])
    parser.add_argument('-o', '--out_folder', type=str, help=helps['out_folder'], default=defaults['out_folder'])
    parser.add_argument('-c', '--obs_date', type=str, help=helps['obs_date'], default=defaults['obs_date'])
    parser.add_argument('-n', '--obs_name', type=str, help=helps['obs_name'], default=defaults['obs_name'])
    parser.add_argument('-r', '--reduce_all', action='store_true', help=helps['reduce_all'], default=defaults['reduce_all'])
    parser.add_argument('-b', '--bias_label', type=str, help=helps['bias_label'], default=defaults['bias_label'])
    parser.add_argument('-a', '--arc_label', type=str, help=helps['arc_label'], default=defaults['arc_label'])
    parser.add_argument('-d', '--dark_label', type=str, help=helps['dark_label'], default=defaults['dark_label'])
    parser.add_argument('-f', '--flat_label', type=str, help=helps['flat_label'], default=defaults['flat_label'])
    parser.add_argument('-t', '--twilight_flat_label', type=str, help=helps['twilight_flat_label'], default=defaults['twilight_flat_label'])
    argv = None
    args = parser.parse_args(args=argv)

    return args


def setup_parser():
    ''' BRIEF DESCRIPTION '''
    # TODO: update docstring
    # TODO: specify what CLI tool or script this is used for?
    parser = ArgumentParser(add_help=True)

    parser.add_argument("-sd", "--start_date",
                        help='''Start Date, e.g., 20170321, YYYYMMDD''',
                        type=str, default=None)

    parser.add_argument("-ed", "--end_date",
                        help='''Start Date, e.g., 20170326, YYYYMMDD''',
                        type=str, default=None)

    parser.add_argument("-dl", "--date_length",
                        help='''Days after/before start/end date, e.g., 10''',
                        type=int, default=None)

    parser.add_argument("-r", "--rootdir",
                        help='''Root Directory for Date''',
                        type=str, default='/work/03946/hetdex/maverick')

    parser.add_argument("-in", "--instrument",
                        help='''Instrument, e.g., virus''',
                        type=str, default='virus')

    return parser


def setup_basic_parser():
    ''' BRIEF DESCRIPTION '''
    # TODO: update docstring
    # TODO: specify what CLI tool or script this is used for?
    parser = ArgumentParser(add_help=True)

    parser.add_argument("-d", "--date",
                        help='''Date, e.g., 20170321, YYYYMMDD''',
                        type=str, default=None)

    parser.add_argument("-o", "--observation",
                        help='''Observation number, "00000007" or "7"''',
                        type=str, default=None)

    parser.add_argument("-e", "--exposure_number",
                        help='''Exposure number, 10''',
                        type=int, default=None)

    parser.add_argument("-r", "--rootdir",
                        help='''Root Directory for Reductions''',
                        type=str, default='/work/03946/hetdex/maverick')

    parser.add_argument("-in", "--instrument",
                        help='''Instrument, e.g., lrs2''',
                        type=str, default='lrs2')

    parser.add_argument("-i", "--ifuslot",
                        help='''Ifuslot, e.g., 066''',
                        type=str, default='066')

    parser.add_argument("-s", "--side",
                        help='''Instrument Side, e.g., L''',
                        type=str, default='L')

    return parser


def set_daterange(args):
    # TODO: add docstring
    # TODO: specify what CLI tool or script this is used for? location here is based on input/output of Argparser object
    dateatt = ['start_date', 'end_date']
    if args.date_length is None:
        if args.start_date is None:
            args.log.error('You must include two of the following: '
                           '"start_date", "end_date", or "date_length"')
            sys.exit(1)
        if args.end_date is None:
            args.log.error('You must include two of the following: '
                           '"start_date", "end_date", or "date_length"')
            sys.exit(1)
        dates = {}
        for da in dateatt:
            dates[da] = dt(int(getattr(args, da)[:4]),
                           int(getattr(args, da)[4:6]),
                           int(getattr(args, da)[6:]))

        args.daterange = [datetime.date.fromordinal(i)
                          for i in range(dates[dateatt[0]].toordinal(),
                                         dates[dateatt[1]].toordinal())]
    else:
        if args.start_date is not None and args.end_date is not None:
            args.log.warning('Using "start_date" and "date_length", '
                             'however, you specified "end_date" as well '
                             'which will not be used.')
            args.end_date = None
        if args.start_date is not None:
            base = dt(int(args.start_date[:4]),
                      int(args.start_date[4:6]),
                      int(args.start_date[6:]))
            args.daterange = [base + datetime.timedelta(days=x)
                              for x in range(0, args.date_length)]

        if args.end_date is not None:
            base = dt(int(args.end_date[:4]),
                      int(args.end_date[4:6]),
                      int(args.end_date[6:]))
            args.daterange = [base - datetime.timedelta(days=x)
                              for x in range(0, args.date_length)]

    return args
