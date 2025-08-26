# Antigen

Data Reduction Pipelines for GCMS (VIRUS-P), VIRUS-W, and VIRUS-2

## Overview

- Antigen is designed to reduce data from the George and Cynthia Mitchell Spectrograph (GCMS) and VIRUS-W on the 2.7m Harlan J. Smith Telescope at McDonald Observatory. 
- In its current state, Antigen outputs a fits file for each science exposure in a night of data taken with GCMS that is the fiber-extracted, wavelength-calibrated, but "raw" spectra.  

## Contributing

- for collaborative development, see git workflow, see [./git-workflow.md](./git-workflow.md)
- for Issues/Requests, see the data reduction template presented by GitHub when creating a new issue.

## Installation

- for installing conda, see [./install_conda.sh](./install_conda.sh)
- for installing the `antigen` python package, see [./INSTALL.md](INSTALL.md)

## Usage

### GCMS Usage
- for GCMS, you will want to copy the data to a folder structure to help the antigen reduction
aAn example call might be:

    ```bash
    $ conda activate env_antigen
    $ antigen_prepare_gcms_dataset.py -i ~/Downloads/VIRUS-P_Data/20240609 -c 20240609 -m GCMS -l VP1R -o ~/data -v
    ```
- then we can begin the basic reduction
- an example call might be:

    ```bash
    $ antigen_reduce.py -i ~/data -c 20240609 -m GCMS -v
    ```

- The arguments for the script are:

    ```plaintext
    usage: antigen_reduce.py [-h] [-i IN_FOLDER] [-o OUT_FOLDER] [-c OBS_DATE] [-n OBS_NAME] [-r] [-w TIME_RADIUS] [-v] [-b BIAS_LABEL] [-a ARC_LABEL] [-d DARK_LABEL] [-f FLAT_LABEL] [-t TWILIGHT_FLAT_LABEL]
                             [-m INSTRUMENT]
  
    Purpose:
        This script runs the full reduction pipeline for VIRUS2/GCMS instrument datasets
        for a given night of observations. It uses manifest files to locate input
        science and calibration data, applies the standard reduction steps, and
        writes reduced science-ready FITS files to the specified output directory.
  
    What it does:
        - Searches the input folder for datasets matching the given observation date,
          observation name, and calibration frame labels.
        - Groups raw files into logical datasets for reduction.
        - Runs the VIRUS2 reduction pipeline on each dataset, which includes:
            * Bias correction
            * Trace calibration
            * Flat fielding
            * Wavelength calibration
            * Science frame spectral extraction
        - Writes the reduced, calibrated science frames to the output folder.
  
    Inputs:
        - Input folder: Directory containing raw spectrograph FITS files
        - Output folder: Destination for reduced FITS files.
        - Observation date: Date of the observation run (YYYYMMDD format).
        - Optional: Observation name, time radius, and custom labels for calibration frames.
        - Other common flags: e.g., --reduce-all to include all matching datasets.
  
    Outputs:
        - One reduced, science-ready FITS file per dataset found.
        - Files are saved in the output folder with names indicating unit ID and processing details.
  
    Example:
        Reduce all VIRUS2 data for the night of 20250801:
  
        $ antigen_reduce.py -i ~/data -o ~/data/VIRUS/reduced -c 20250801 -m VIRUS2
  
    Notes:
        - The script logs all processing steps and reports any datasets that failed to reduce.
        - Reduction uses calibration frames as specified by bias, dark, arc, flat, and twilight labels.
        - This script assumes a valid configuration of the VIRUS2 reduction pipeline.
  
    options:
      -h, --help            show this help message and exit
      -i IN_FOLDER, --in_folder IN_FOLDER
                            Root path where reduction input file tree is located, (default: /Users/grz85/work/code/Antigen)
      -o OUT_FOLDER, --out_folder OUT_FOLDER
                            Path where reduction output files will be written, (default: antigen_reduce_virus2_20250826_095917)
      -c OBS_DATE, --obs_date OBS_DATE
                            Observation calendar date string formatted as YYYYMMDD, ex: 20250613, (default: 20250826)
      -n OBS_NAME, --obs_name OBS_NAME
                            Observation object/target name, e.g. from FITS header card, (default: None)
      -r, --reduce_all      Reduce all files found under in_folder file tree, (default: False)
      -w TIME_RADIUS, --time_radius TIME_RADIUS
                            All calibration files within this MJD radius of a science file will be added to its manifest, (default: 2.0)
      -v, --verbose         if True, print more process details and logger.info to console, (default: False)
      -b BIAS_LABEL, --bias_label BIAS_LABEL
                            The object name from the FITS header card for bias files, (default: bias)
      -a ARC_LABEL, --arc_label ARC_LABEL
                            The object name from the FITS header card for arc files, (default: arc)
      -d DARK_LABEL, --dark_label DARK_LABEL
                            The object name from the FITS header card for dark files, (default: dark)
      -f FLAT_LABEL, --flat_label FLAT_LABEL
                            The object name from the FITS header card for flat files, (default: flat)
      -t TWILIGHT_FLAT_LABEL, --twilight_flat_label TWILIGHT_FLAT_LABEL
                            The object name from the FITS header card for twilight flat files, (default: twi)
      -m INSTRUMENT, --instrument INSTRUMENT
                            Name of the instrument used for the observation (default: GCMS). Example: GCMS.
      ```


### VIRUS2 Usage

- Example 1: How to build and read a manifest

    ```bash
    $ conda activate env_antigen
    
    $ antigen_build_manifest_virus2.py \
        -c 20250619  \
        -f flatp  \
        --obs_name standard  \
        -w 2  \
        --verbose
    
    $ antigen_read_manifest_virus2.py -f antigen/config_files/virus2/virus2_manifest_template.yml --validate
    ```


- Example 2: How to run a reduction using an existing manifest (recommended!)

    ```bash
    $ conda activate env_antigen
    $ antigen_reduce_manifest.py -m antigen/config_files/virus2/virus2_manifest_template.yml -v
    ```

- Example 3: How to build a manifest and reduce it, without pre-existing manifest

    ```bash
    $ conda activate env_antigen
    $ dir_containing_VIRUS2=/home/user/my_data
    $ antigen_reduce.py -i $dir_containing_VIRUS2 -c 20240609 -m VIRUS2 -v
    ```
