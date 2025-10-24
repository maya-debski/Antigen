#!/usr/bin/env python

import argparse
from pathlib import Path
import re
import shutil

from antigen.cli import add_calibration_args, add_common_args
from antigen.io import load_fits_header, get_fits_files_in_path
from antigen.utils import setup_logging
from astropy.io import fits

DESCRIPTION = r"""
Purpose:
    This script organizes raw GCMS (George and Cynthia Mitchell Spectrograph) data for a given night of observations.
    It scans a folder of raw FITS files, inspects their headers, determines frame types,
    and copies them into a structured output directory with standardized filenames.
    This ensures consistency and traceability for later data reduction steps.

What it does:
    - Scans the specified input folder for all FITS files matching *.fits.
    - Reads each file’s FITS header to extract the target name, observation date, and UT time.
    - Determines the frame type (science, bias, flat, arc, or dark) using keywords in the OBJECT header.
    - Creates an organized output folder hierarchy:
      OUT_FOLDER/<instrument>/<obsdate>/<obsid>/<element>/
    - Renames and copies each file into its target folder using a clear, standardized naming convention:
      <instrument>_<obsdate>_<obsid>_<frametype>_<element>_exp<exposureindex>_<utctime>_<objectname>.fits

Inputs:
    - Input folder: Directory containing raw GCMS FITS files.
    - Output folder: Destination for organized data.
    - Observation date: Date of the observation run (YYYYMMDD format).
    - Instrument name: Name of the instrument (default: GCMS).
    - Spectrograph configuration element: e.g., VP1B or VP1R (default: VP1B).
    - Optional: Custom labels for calibration frames (bias, flat, arc, dark).

Outputs:
    - An organized folder tree under the output folder.
    - Renamed FITS files with descriptive filenames that include:
      instrument, date, observation ID, frame type, configuration element, exposure index, UT time, and object name.

Example:
    Prepare GCMS data for the night of 20240609 taken with the VP1R configuration:

    $ antigen_prepare_dataset.py -i ~/Downloads/VIRUS-P_Data/20240609 -c 20240609 -m GCMS -l VP1R -o ~/data -v

Notes:
    - Frame types are inferred by matching keywords in the OBJECT header card.
    - The observation ID increments each time a new object is found in the sequence.
    - The exposure index counts repeated exposures of the same target.
    - This script is typically run before building a manifest or starting data reduction.
    - Supported setups for GCMS include VP1B and VP1R.
"""

def get_args():
    parser = argparse.ArgumentParser(
        description=DESCRIPTION,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    add_common_args(parser)
    add_calibration_args(parser)

    # Override the generic in_folder help text from add_common_args for this script's context
    # We only change the help string, leaving defaults and behavior unchanged.
    for action in parser._actions:
        if '--in_folder' in action.option_strings:
            action.help = 'Root path where raw data folder is located, (default: %(default)s)'
            break

    for action in parser._actions:
        if '--out_folder' in action.option_strings:
            action.help = 'Root path modified raw data folder will go, (default: %(default)s)'
            break

    parser.add_argument(
        '-m', '--instrument',
        type=str,
        default='GCMS',
        help='Name of the instrument used for the observation '
             '(default: %(default)s). Example: GCMS.'
    )

    parser.add_argument(
        '-j', '--element',
        type=str,
        default='VP1B',
        help='Spectrograph setup or configuration element, such as VP1B or VP1R '
             '(default: %(default)s).'
    )

    return parser.parse_args()


def detect_and_fix_virus_w_header(header, filename, logger):
    """
    Detect virus-w style DATE-OBS and split into DATE-OBS and UT if needed.

    Args:
        header (fits.Header): FITS header to check and potentially modify
        filename (str): Filename for logging purposes
        logger: Logger instance

    Returns:
        bool: True if header was modified, False otherwise
    """
    if 'DATE-OBS' not in header:
        return False

    date_obs = header['DATE-OBS']

    # Check if DATE-OBS contains time (virus-w format): YYYYMMDDTHH:MM:SS.sss or YYYY-MM-DDTHH:MM:SS.sss
    time_pattern = r'(\d{4}[-]?\d{2}[-]?\d{2})T(\d{2}:\d{2}:\d{2}(?:\.\d+)?)'
    match = re.match(time_pattern, str(date_obs))

    if match:
        date_part = match.group(1)
        time_part = match.group(2)

        # Update header
        header['DATE-OBS'] = date_part
        header['UT'] = time_part

        logger.info(
            f"Fixed virus-w header format in {filename}: split '{date_obs}' into DATE-OBS='{date_part}' and UT='{time_part}'")
        return header, True

    return header, False


def main():
    args = get_args()
    logger = setup_logging('antigen', verbose=args.verbose)

    filenames = get_fits_files_in_path(args.in_folder, pattern='*.fits')

    instrument = args.instrument
    config_element = args.element

    bias_labels = ['zero', 'bias', args.bias_label]
    flat_labels = ['ldls', 'flat', 'flt', 'dome', 'twilight', args.flat_label]
    lamp_labels = ['arc', 'comp', 'HgNe', 'FeAr', args.arc_label]
    dark_labels = ['dark', args.dark_label]

    exposure_counter = 1 # Starts at 1 and moves up if there is an object repeat
    observation_counter = 0 # Starts at zero and goes to 1 for first file
    previous_object = ''
    for filename in filenames:
        header, is_header_valid = load_fits_header(filename)

        if not is_header_valid:
            logger.warning(f"Invalid header for {filename}")
            continue

        header, header_was_modified = detect_and_fix_virus_w_header(header, filename, logger)
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

        try:
            dst_file = Path(output_path)
            dst_file.parent.mkdir(parents=True, exist_ok=True)
            
            # If header was modified, write the corrected version
            if header_was_modified:
                # Read the full FITS file to preserve data
                with fits.open(filename) as hdul:
                    # Update the primary HDU header with corrected values
                    hdul[0].header['DATE-OBS'] = header['DATE-OBS']
                    hdul[0].header['UT'] = header['UT']
                    # Write to destination with corrected header
                    hdul.writeto(dst_file, overwrite=True)
                logger.info(f"Copied and corrected header: {filename} -> {dst_file}")
            else:
                # Just copy the file as-is
                shutil.copy2(filename, dst_file)
                logger.info(f"Copied {filename} -> {dst_file}")
                
        except (OSError, IOError, PermissionError) as e:
            logger.error(f"Failed to copy {filename}: {e}")
            continue
        except Exception as e:
            logger.error(f"Unexpected error copying {filename}: {e}")
            continue

if __name__ == "__main__":
    main()