#!/usr/bin/env python

import argparse
import os
import glob
from pathlib import Path
import shutil

from antigen.cli import add_calibration_args, add_common_args
from antigen.io import load_fits_header
from antigen.utils import setup_logging

DESCRIPTION = r"""
Purpose: Create a folder structure for a night of GCMS data from the following pattern:
IN_FOLDER/
    *0001.fits
    *0002.fits
    *0003.fits
    *0004.fits

to this pattern:
OUT_FOLDER/GCMS/<obsdate>/<obsid>/<specid>/
    GCMS_<obsdate>_<obsid>_<frametype>_<specid>_exp<exposureindex>_<utctime>_<objectname>.fits

Examples:
    antigen_prepare_gcms_dataset.py -i <input_folder> -o <output_folder> -c <obs_date>
"""


def get_args():
    parser = argparse.ArgumentParser(description=DESCRIPTION, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_common_args(parser)
    add_calibration_args(parser)
    return parser.parse_args()

def main():
    args = get_args()
    logger = setup_logging('antigen', verbose=args.verbose)

    filenames = sorted(glob.glob(os.path.join(args.in_folder, '*.fits')))

    # TODO: Change this to be read from the arguments later
    instrument = 'GCMS'
    config_element = 'VP1B'

    required_keywords = ['OBJECT', 'DATE-OBS', 'UT']
    bias_labels = ['zero', 'bias', args.bias_label]
    flat_labels = ['ldls', 'flat', 'flt', 'dome', 'twilight', args.flat_label]
    lamp_labels = ['arc', 'comp', 'HgNe', 'FeAr', args.arc_label]
    dark_labels = ['dark', args.dark_label]

    exposure_counter = 1 # Starts at 1 and moves up if there is an object repeat
    observation_counter = 0 # Starts at zero and goes to 1 for first file
    previous_object = ''
    for filename in filenames:
        header, is_header_valid = load_fits_header(filename)

        object_name = header['OBJECT']
        date = ''.join(header['DATE-OBS'].split('-'))
        ut = ''.join(header['UT'].split(':'))

        if object_name == previous_object:
            exposure_counter += 1
        else:
            exposure_counter = 1
            observation_counter += 1
            previous_object = object_name

        obs_string = f"{observation_counter:07d}"
        exp_string = f"exp{exposure_counter:02d}"

        frametype = 'sci'
        for labels, framelabel in zip([bias_labels, flat_labels, lamp_labels, dark_labels],
                                      ['bias', 'flat', 'arc', 'dark']):
            for label in labels:
                if label.lower() in object_name.lower():
                    frametype = framelabel

        clean_object_name = ''.join(object_name.split())
        output_dir = Path(args.out_folder) / instrument / args.obs_date / obs_string / config_element

        output_filename = (
            f"{instrument}_{args.obs_date}_{obs_string}_"
            f"{frametype}_{config_element}_{exp_string}_{date}T{ut}_{clean_object_name}.fits"
        )

        output_path = output_dir / output_filename

        dst_file = Path(output_path)
        dst_file.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(filename, dst_file)
        logger.info(f"Copied {filename} -> {dst_file}")

if __name__ == "__main__":
    main()
