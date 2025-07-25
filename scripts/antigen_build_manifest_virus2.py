#!/usr/bin/env python

import os
import sys
from antigen.cli import get_args
from antigen.process import build_manifest_records
from antigen.manifest import save_manifest


def main():

    args = get_args()

    manifest_records = build_manifest_records(args.in_folder, args.obs_date,
                                              args.obs_name, args.reduce_all, args.time_radius,
                                              args.bias_label, args.arc_label, args.dark_label,
                                              args.flat_label, args.twilight_flat_label
                                              )

    for nr, record in enumerate(manifest_records):
        manifest_filename = f'manifest_{args.obs_date}_{args.obs_name}_record{nr}.yml'
        save_path = os.path.abspath(args.out_folder)
        os.makedirs(save_path, exist_ok=True)
        save_filepath = os.path.join(save_path, manifest_filename)
        save_manifest(record, save_filepath)
        print(f'Wrote manifest file: {save_filepath}')


    return None


if __name__ == '__main__':
    sys.exit(main())
