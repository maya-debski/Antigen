#!/usr/bin/env python

import os
import glob
import shutil
from pathlib import Path
from astropy.io import fits
from antigen.cli import get_args
from antigen.utils import setup_logging

def main():
    args = get_args()
    logger = setup_logging('antigen', verbose=args.verbose)

    filenames = sorted(glob.glob(os.path.join(args.input_folder, '*.fits')))

    # TODO: Change this to be read from the arguments later
    instrument = 'GCMS'
    setup = 'VP1B'
    
    required_keywords = ['OBJECT', 'DATE-OBS', 'UT']
    bias_labels = ['zero', 'bias', args.bias_label]
    flat_labels = ['ldls', 'flat', 'flt', 'dome', 'twilight', args.flat_label]
    lamp_labels = ['arc', 'comp', 'HgNe', 'FeAr', args.arc_label]
    dark_labels = ['dark', args.dark_label]

    exposure_counter = 1
    observation_counter = 1
    previous_object = ''

    for filename in filenames:
        with fits.open(filename) as hdu:
            header = hdu[0].header

            missing = [key for key in required_keywords if key not in header]
            if missing:
                logger.warning(f"{filename} missing required keywords: {', '.join(missing)}")
                continue

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
            output_dir = Path(args.out_folder) / instrument / args.obsdate / obs_string / args.setup

            filename = (
                f"{instrument}_{args.obsdate}_{obs_string}_"
                f"{frametype}_{setup}_{exp_string}_{date}T{ut}_{clean_object_name}.fits"
            )

            output_path = output_dir / filename

            dst_file = Path(output_path)
            dst_file.parent.mkdir(parents=True, exist_ok=True)

            shutil.copy2(filename, dst_file)
            logger.info(f"Copied {filename} -> {dst_file}")

if __name__ == "__main__":
    main()
