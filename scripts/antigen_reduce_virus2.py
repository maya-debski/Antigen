#!/usr/bin/env python
import os
import sys

from antigen import utils
from antigen.cli import get_args
from antigen.datasets import find_datasets
from antigen.reduce_virus2 import process_unit


def main():
    args = get_args()

    log = utils.setup_logging('virus2_reductions')
    os.makedirs(args.out_folder, exist_ok=True)

    dataset_manifests = find_datasets(args.in_folder, args.obs_date, args.obs_name,
                                      args.reduce_all, args.time_radius,
                                      args.bias_label, args.arc_label, args.dark_label,
                                      args.flat_label, args.twilight_flat_label)

    for record in dataset_manifests:
        unit_id = record['unit_id']
        log.info(f'Processing reduction for unit = {unit_id}:')
        try:
            output_fits_filename = process_unit(record, args.out_folder)
            log.info(f'Processing reduction for unit = {unit_id}: PASS: wrote reduction to FITS file {output_fits_filename}')
        except Exception as error:
            log.error(f'Processing reduction for unit = {unit_id}: FAILED: {error}')

    return None


if __name__ == '__main__':
    sys.exit(main())
