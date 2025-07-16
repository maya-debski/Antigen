#!/usr/bin/env python

import sys
from antigen.cli import get_args
from antigen.reduce_virus2 import process


def main():

    args = get_args()

    process(args.infolder, args.outfolder,
            args.date, args.target_name, args.reduce_all,
            args.bias_label, args.arc_label, args.dark_label, args.flat_label, args.twilight_flat_label
            )

    return None


if __name__ == '__main__':
    sys.exit(main())
