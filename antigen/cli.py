from argparse import ArgumentParser
import datetime
import os

DEFAULTS = {
    'in_folder': os.path.abspath(os.curdir),
    'out_folder': datetime.datetime.now().strftime('antigen_reduce_virus2_%Y%m%d_%H%M%S'),
    'obs_date': datetime.datetime.now().strftime('%Y%m%d'),
    'obs_name': None,
    'reduce_all': False,
    'time_radius': 2.0,
    'bias_label': 'bias',
    'arc_label': 'arc',
    'dark_label': 'dark',
    'flat_label': 'flat',
    'twilight_flat_label': 'twi',
    'verbose': False,
}

TYPES = {
    'in_folder': str,
    'out_folder': str,
    'obs_date': str,
    'obs_name': str,
    'reduce_all': bool,
    'time_radius': float,
    'bias_label': str,
    'arc_label': str,
    'dark_label': str,
    'flat_label': str,
    'twilight_flat_label': str,
    'verbose': bool,
}

HELPS = {
    'in_folder': 'Root path where reduction input file tree is located, (default: %(default)s)',
    'out_folder': 'Path where reduction output files will be written, (default: %(default)s)',
    'obs_date': 'Observation calendar date string formatted as YYYYMMDD, ex: 20250613, (default: %(default)s)',
    'obs_name': 'Observation object/target name, e.g. from FITS header card, (default: %(default)s)',
    'reduce_all': 'Reduce all files found under infolder file tree, (default: %(default)s)',
    'time_radius': 'All calibration files within this MJD radius of a science file will be added to its manifest, (default: %(default)s)',
    'bias_label': 'The object name from the FITS header card for bias files, (default: %(default)s)',
    'arc_label': 'The object name from the FITS header card for arc files, (default: %(default)s)',
    'dark_label': 'The object name from the FITS header card for dark files, (default: %(default)s)',
    'flat_label': 'The object name from the FITS header card for flat files, (default: %(default)s)',
    'twilight_flat_label': 'The object name from the FITS header card for twilight flat files, (default: %(default)s)',
    'verbose': 'if True, print more process details and logger.info to console, (default: %(default)s)',
}

ABBREVS = {
    'in_folder': 'i',
    'out_folder': 'o',
    'obs_date': 'c',
    'obs_name': 'n',
    'reduce_all': 'r',
    'time_radius': 'w',
    'bias_label': 'b',
    'arc_label': 'a',
    'dark_label': 'd',
    'flat_label': 'f',
    'twilight_flat_label': 't',
    'verbose': 'v',
}

def add_args_from_dicts(parser, arg_names):
    """Add arguments to an ArgumentParser based on provided defaults and help strings.

    Args:
        parser (ArgumentParser): The parser to which arguments will be added.
        arg_names (list[str]): A list of argument names (keys in `defaults` and `helps`).

    Returns:
        parser (ArgumentParser): The updated parser.
    """
    for name in arg_names:
        default = DEFAULTS[name]
        help_msg = HELPS[name]
        abbrev = ABBREVS[name]
        argtype = TYPES[name]

        if isinstance(argtype, bool):
            parser.add_argument(f'-{abbrev}', f'--{name}', action='store_true',  default=default,  help=help_msg)
        else:
            parser.add_argument(f'-{abbrev}', f'--{name}', type=argtype, default=default, help=help_msg)

    return parser


def add_common_args(parser):
    """Add common arguments to the parser.

    Args:
        parser (ArgumentParser): The parser to extend.

    Returns:
        parser (ArgumentParser): The updated parser with common args.
    """
    common_args = ['in_folder', 'out_folder', 'obs_date', 'obs_name', 'reduce_all', 'time_radius', 'verbose']
    return add_args_from_dicts(parser)


def add_calibration_args(parser):
    """Add common arguments to the parser.

    Args:
        parser (ArgumentParser): The parser to extend.

    Returns:
        parser (ArgumentParser): The updated parser with common args.
    """
    common_args = ['bias_label', 'arc_label', 'dark_label', 'flat_label', 'twilight_flat_label']
    return add_args_from_dicts(parser, common_args)

