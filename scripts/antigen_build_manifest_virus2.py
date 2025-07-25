#!/usr/bin/env python

import os
import sys
from antigen.cli import get_args
from antigen.datasets import find_datasets
from antigen.manifest import save_manifest


def main():

    args = get_args()

    dataset_manifests = find_datasets(args.in_folder, args.obs_date,
                                      args.obs_name, args.reduce_all, args.time_radius,
                                      args.bias_label, args.arc_label, args.dark_label,
                                      args.flat_label, args.twilight_flat_label
                                      )

    for nr, record in enumerate(dataset_manifests):
        manifest_filename = f'manifest_{args.obs_date}_{args.obs_name}_record{nr}.yml'
        save_path = os.path.abspath(args.out_folder)
        os.makedirs(save_path, exist_ok=True)
        save_filepath = os.path.join(save_path, manifest_filename)
        save_manifest(record, save_filepath)
        print(f'Wrote manifest file: {save_filepath}')


    return None


if __name__ == '__main__':
    sys.exit(main())
