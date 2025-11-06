# Antigen

Data Reduction Pipeline for GCMS (VIRUS-P), VIRUS-W, and VIRUS-2

## Overview

Antigen is designed to reduce data from the George and Cynthia Mitchell Spectrograph (GCMS), VIRUS-W, and the forthcoming VIRUS-2 instruments on the 2.7m Harlan J. Smith Telescope at McDonald Observatory. Antigen automates CCD reduction, wavelength and flux calibration, and the generation of science-ready datacubes and outputs FITS files and diagnostic plots at each step of the process (base CCD reduction, response measurement (optional), and cube creation).

## Contributing

- for collaborative development, see git workflow, see [git-workflow.md](./git-workflow.md)
- for Issues/Requests, see the data reduction template presented by GitHub when creating a new issue.

## Installation

- for installing conda, see [install_conda.sh](.docs/installation/install_conda.sh)
- for installing the `antigen` python package, see [INSTALL.md](./docs/installation/INSTALL.md)

## Usage

### Dataset Preparation

Antigen requires a particular folder and file tree structure. Given a folder of a night of raw data, the script [antigen_prepare_dataset.py](./scripts/antigen_prepare_dataset.py) will:
- Scan the specified input folder for all FITS files matching *.fits.
- Read each file’s FITS header to extract the target name, observation date, and UT time.
- Determine the frame type (science, bias, flat, arc, or dark) using keywords in the OBJECT header.
- Create an organized output folder hierarchy:
  `OUT_FOLDER/<instrument>/<obsdate>/<obsid>/<element>/`
- Rename and copy each file into its target folder using a clear, standardized naming convention:
  `<instrument>_<obsdate>_<obsid>_<frametype>_<element>_exp<exposureindex>_<utctime>_<objectname>.fits`

An example call might be:
```bash
$ conda activate env_antigen
$ antigen_prepare_dataset.py -i ~/Downloads/VIRUS-P_Data/20240609 -o ~/data -c 20240609 -m GCMS -j VP1R -v
```
The arguments for this script are:
```plaintext
  -h, --help            show this help message and exit
  -i, --in_folder IN_FOLDER
                        Root path where raw data folder is located, (default: /Users/mhd674/Code/github-com/Antigen)
  -o, --out_folder OUT_FOLDER
                        Root path modified raw data folder will go, (default: antigen_reduce_virus2_20251106_124228)
  -c, --obs_date OBS_DATE
                        Observation calendar date string formatted as YYYYMMDD, ex: 20250613, (default: 20251106)
  -n, --obs_name OBS_NAME
                        Observation object/target name, e.g. from FITS header card, (default: None)
  -r, --reduce_all      Reduce all files found under in_folder file tree, (default: False)
  -w, --time_radius TIME_RADIUS
                        All calibration files within this MJD radius of a science file will be added to its manifest, (default: 2.0)
  -v, --verbose         if True, print more process details and logger.info to console, (default: False)
  -d, --debug           if True, print more process details and logger.debug to console, (default: False)
  -b, --bias_label BIAS_LABEL
                        The object name from the FITS header card for bias files, (default: bias)
  -a, --arc_label ARC_LABEL
                        The object name from the FITS header card for arc files, (default: arc)
  -l, --dark_label DARK_LABEL
                        The object name from the FITS header card for dark files, (default: dark)
  -f, --flat_label FLAT_LABEL
                        The object name from the FITS header card for flat files, (default: flat)
  -t, --twilight_flat_label TWILIGHT_FLAT_LABEL
                        The object name from the FITS header card for twilight flat files, (default: twi)
  -m, --instrument INSTRUMENT
                        Name of the instrument used for the observation (default: GCMS). Example: GCMS.
  -j, --element ELEMENT
                        Spectrograph setup or configuration element, such as VP1B or VP1R (default: VP1B).
  ```

See [prepare_dataset.md](./docs/scripts/prepare_dataset.md) for more information about `antigen_prepare_dataset.py

### Base CCD Reduction

```plaintext
options:
  -h, --help            show this help message and exit
  -i, --in_folder IN_FOLDER
                        Root path where reduction input file tree is located, (default: /Users/mhd674/Code/github-com/Antigen)
  -o, --out_folder OUT_FOLDER
                        Path where reduction output files will be written, (default: antigen_reduce_virus2_20251106_124309)
  -c, --obs_date OBS_DATE
                        Observation calendar date string formatted as YYYYMMDD, ex: 20250613, (default: 20251106)
  -n, --obs_name OBS_NAME
                        Observation object/target name, e.g. from FITS header card, (default: None)
  -r, --reduce_all      Reduce all files found under in_folder file tree, (default: False)
  -w, --time_radius TIME_RADIUS
                        All calibration files within this MJD radius of a science file will be added to its manifest, (default: 2.0)
  -v, --verbose         if True, print more process details and logger.info to console, (default: False)
  -d, --debug           if True, print more process details and logger.debug to console, (default: False)
  -b, --bias_label BIAS_LABEL
                        The object name from the FITS header card for bias files, (default: bias)
  -a, --arc_label ARC_LABEL
                        The object name from the FITS header card for arc files, (default: arc)
  -l, --dark_label DARK_LABEL
                        The object name from the FITS header card for dark files, (default: dark)
  -f, --flat_label FLAT_LABEL
                        The object name from the FITS header card for flat files, (default: flat)
  -t, --twilight_flat_label TWILIGHT_FLAT_LABEL
                        The object name from the FITS header card for twilight flat files, (default: twi)
  -m, --instrument INSTRUMENT
                        Name of the instrument used for the observation (default: GCMS). Example: GCMS.
  -j, --good_arc_residual_limit GOOD_ARC_RESIDUAL_LIMIT
                        Residual limit for an arc line to be considered well-fit(default: 0.2). Example: 0.2.
  -g, --binned          Data is binned in the x-direction?
```

### Measuring Response

```plaintext
options:
  -h, --help            show this help message and exit
  -r, --reduced_dir REDUCED_DIR
                        Path to directory with reduced standard star frames.
  -o, --output_folder OUTPUT_FOLDER
                        Path to output folder for response function.
  -s, --standard_name STANDARD_NAME
                        Name for standard star files in reduced-dir (ex: Feige).
  -e, --extraction_radius EXTRACTION_RADIUS
                        Extraction radius in arcsec.
  -p, --pixel_size PIXEL_SIZE
                        Cube pixel size in arcsec.
  -v, --verbose         if True, print more process details and logger.info to console
  -d, --debug           if True, print more process details and logger.debug to console
  -g, --binned          Data is binned in the x-direction?
```

### Cube Creation

```plaintext
options:
  -h, --help            show this help message and exit
  -r, --reduced_dir REDUCED_DIR
                        Path to directory with reduced standard star frames.
  -o, --output_folder OUTPUT_FOLDER
                        Path to output folder for response function.
  -s, --object_name OBJECT_NAME
                        Name for standard star files in reduced-dir (ex: Feige).
  -e, --extraction_radius EXTRACTION_RADIUS
                        Extraction radius in arcsec.
  -p, --pixel_size PIXEL_SIZE
                        Cube pixel size in arcsec.
  -v, --verbose         if True, print more process details and logger.info to console
  -d, --debug           if True, print more process details and logger.debug to console
  -i, --interpolation_method {linear,nearest,cubic,rbf,gdw}
                        Interpolation method for cube construction (default: gdw)
  -k, --k_neighbors K_NEIGHBORS
                        Number of neighbors for GDW interpolation (default: 11)
  --sigma_gdw SIGMA_GDW
                        Standard deviation for GDW interpolation in arcseconds (default: 2.5)
  --rbf_function {multiquadric,inverse,gaussian,linear,cubic,quintic,thin_plate}
                        Radial basis function type for RBF interpolation (default: multiquadric)
  -g, --binned          Data is binned in the x-direction?
```
    ```
