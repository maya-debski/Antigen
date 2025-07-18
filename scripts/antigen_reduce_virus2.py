#!/usr/bin/env python

import sys
from antigen.cli import get_args
from antigen.reduce_virus2 import process


def main():

    args = get_args()

    process(args.infolder, args.outfolder,
            args.date, args.name, args.reduce_all,
            args.bias, args.arc, args.dark, args.flat, args.twilight
            )

    return None


if __name__ == '__main__':
    sys.exit(main())
