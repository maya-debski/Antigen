#!/usr/bin/env python
"""
Purpose: View a FITS file header and ...
         - optionally save header to TXT file
         - optionally load, plot, and display data
         - optionally save data plot to PNG file

Note: Why? Because (1) ds9 is slow, and (2) ds9 is a GUI and cannot be easily scripted.
"""
from argparse import ArgumentParser, RawTextHelpFormatter
import json
import os

# Set the rendering graphics back-end to ensure any plots will be displayed on all systems
import platform
import matplotlib
if platform.system() == 'Darwin':
    matplotlib.use('MacOSX')
elif platform.system() == 'Linux':
    matplotlib.use('TkAgg')
elif platform.system() == 'Windows':
    matplotlib.use('TkAgg')

from antigen.io import read_fits, write_fits_header_txt
from antigen.plot import plot_frame, get_array_data_summary


def get_input_defaults():
    default_dict = dict(file=None,
                        data=False,
                        plot=False,
                        save=False,
                        outpath=".",
                        verbose=False)
    return default_dict


def get_inputs():
    parser_desc = """Purpose: View a FITS file header and ...
         - optionally save header to TXT file
         - optionally load data, display data statistics, plot data
         - optionally save data stats and data plot to files
    """
    parser = ArgumentParser(description=parser_desc, formatter_class=RawTextHelpFormatter)

    arg_defaults = get_input_defaults()

    arg_helps = dict(file='FITS file name to read (only loads header, unless data==True), (default: %(default)s)',
                     data='Load FITS data, display summary statistics (if plot==True, display) (default: %(default)s)',
                     plot='Plot FITS data, display (if save==True, save to PNG), (default: %(default)s)',
                     save='Save header (if data==True, also save stats,plot) to outpath, (default: %(default)s)',
                     outpath='Output path to write files (header, data stats, data plot) (default: %(default)s)',
                     verbose='if verbose==True, display more process detail to console, (default: %(default)s)')

    parser.add_argument("-f", "--file", default=arg_defaults['file'], help=arg_helps['file'])
    parser.add_argument("-d", "--data", action='store_true', default=arg_defaults['data'], help=arg_helps['data'])
    parser.add_argument("-p", "--plot", action='store_true', default=arg_defaults['plot'], help=arg_helps['plot'])
    parser.add_argument("-s", "--save", action='store_true', default=arg_defaults['save'], help=arg_helps['save'])
    parser.add_argument("-o", "--outpath", default=arg_defaults['outpath'], help=arg_helps['outpath'])
    parser.add_argument("-v", "--verbose", action='store_true', default=arg_defaults['verbose'], help=arg_helps['verbose'])

    args = parser.parse_args()

    return args


def main():
    args = get_inputs()

    if args.file and not os.path.isfile(args.file):
        error_msg = f'Specified FITS file does not exist: {args.file}'
        raise FileNotFoundError(error_msg)
    if args.verbose:
        print(f'FITS file name = {args.file}')

    if args.save:
        if os.path.isdir(args.outpath):
            error_msg = f'ERROR: directory already exists: too risky to write to it: {args.outpath}'
            raise IsADirectoryError(error_msg)
        else:
            os.makedirs(args.outpath)

    if args.data or args.plot:
        fits_data, fits_header = read_fits(file_name=args.file, read_data=True)
        fits_data_stats = get_array_data_summary(fits_data)
    else:
        fits_data, fits_header = read_fits(args.file, read_data=False)
        fits_data_stats = None

    if fits_header:
        if args.verbose:
            print(fits_header.tostring(sep='\n', endcard=False))
        if args.save:
            save_file = os.path.join(args.outpath, 'fits_header.txt')
            _ = write_fits_header_txt(fits_header, save_file)


    if fits_data_stats:
        if args.verbose:
            print(f'\nFITS data stats:')
            for key, val in fits_data_stats.items():
                print(f'{key:10} = {val}')
        if args.save:
            save_file = os.path.join(args.outpath, 'fits_data_stats.json')
            if save_file:
                with open(save_file, 'w') as fob:
                    json.dump(fits_data_stats, fob, indent=4, sort_keys=True)

    if args.data and args.plot:
        import matplotlib.pyplot as plt  # needed only for plt.show()

        naxis = fits_header.get('NAXIS', None)    # e.g. 2 or 3 / number of data axes
        naxis1 = fits_header.get('NAXIS1', None)  # e.g. 2096 / length of data axis 1 / number of x-bins
        naxis2 = fits_header.get('NAXIS2', None)  # e.g. 2064 / length of data axis 2 / number of y-bins
        naxis3 = fits_header.get('NAXIS3', None)  # e.g. 4 / length of data axis 3 / number of channels

        if naxis and naxis == 2:
            plot_fig, plot_ax = plot_frame(fits_data)
            plt.show()
        else:
            raise Exception(f'Unexpected number of axes, dimensions: naxis={naxis}, naxis1={naxis1}, naxis2={naxis2}, naxis3={naxis3}')

        if args.save:
            save_file = os.path.join(args.outpath, 'fits_data_plot.png')
            if os.path.isfile(save_file):
                raise FileExistsError('File name already exists. Will not over-write')
            else:
                plot_fig.savefig(save_file)

    if args.plot and not args.data:
        print(f'WARNING: You asked to plot the data, but also asked to NOT read the data. Please see --help')



if __name__ == "__main__":
    main()


