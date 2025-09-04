# Usage

## GCMS

GCMS data for a single night of observations is output in a folder containing all files, including calibrations, standard stars, and science frames. In order to reduce the data for a single night:

- Copy the data to a folder structure to help the antigen reduction using <tt> antigen\_prepare\_gcms\_dataset.py</tt>. An example call might be:
   ```bash
   $ conda activate env_antigen
   $ antigen_prepare_gcms_dataset.py -i ~/Downloads/VIRUS-P_Data/20240609 -c 20240609 -m GCMS -l VP1R -o ~/data -v
   ```
- Run the basic reduction using <tt> antigen\_reduce.py</tt>. An example call might be:
   ```bash
   $ antigen_reduce.py -i ~/data -c 20240609 -m GCMS -v
   ```
- The arguments for <tt> antigen_reduce.py</tt> are:
  ```bash
  usage: antigen_reduce.py [-h] [-i IN_FOLDER] [-o OUT_FOLDER] [-c OBS_DATE] [-n OBS_NAME] [-r] [-w TIME_RADIUS] [-v] [-b BIAS_LABEL] [-a ARC_LABEL] [-d DARK_LABEL] [-f FLAT_LABEL] [-t TWILIGHT_FLAT_LABEL]
                         [-m INSTRUMENT]

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
    -r, --reduce_all    Reduce all files found under in_folder file tree, (default: False)
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

You can also run antigen on GCMS data by building a manifest file from your restructured data folder through these steps:

- Build the manifest using <tt> antigen\_build\_manifest.py</tt>. An example call might be:
   ```bash
   $ antigen_build_manifest.py -i ~/data -o ~/data/manifests -c 20240609
   ```
- Run a basic reduction using <tt> antigen\_reduce\_manifest.py</tt>. An example call might be:
   ```bash
   $ antigen_reduce_manifest.py -f manifest.yml -o ~/data/reduced -v
   ```

## VIRUS2
