#!/usr/bin/env python

import sys
from antigen.cli import get_args
from antigen.process import process


def main():

    args = get_args()

    process(args.in_folder, args.out_folder,
            args.obs_date, args.obs_name, args.reduce_all, args.time_radius,
            args.bias_label, args.arc_label, args.dark_label, args.flat_label, args.twilight_flat_label
            )

    return None


if __name__ == '__main__':
    sys.exit(main())
