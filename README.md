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

- An example call might be:

    ```bash
    $ conda activate env_antigen
    $ python ./scripts/reduce_virusp.py ~/Downloads/VIRUS-P_Data/20240606 ~/Downloads/VIRUS-P_Data/20240606/reduced -bn -ra -bl "zero" -tfl "Twilight flat" -al "comp" -b
    ```

- The arguments for the script are:

    ```plaintext
    usage: reduce_virusp.py [-h] [-n NAME] [-b] [-bn] [-ra] [-bl BIAS_LABEL] [-al ARC_LABEL] [-dfl DOME_FLAT_LABEL] [-tfl TWILIGHT_FLAT_LABEL] folder outfolder
    
    positional arguments:
    
      folder                Input folder
      
      outfolder             name of the output file
    
    options:
    
      -h, --help            show this help message and exit
      
      -n NAME, --name NAME  Name of the science target
      
      -b, --blue            blue Side?
      
      -bn, --binned         Binned?
      
      -ra, --reduce_all     Reduce all files in folder
      
      -bl BIAS_LABEL, --bias_label BIAS_LABEL
                            The objet name for bias files
                            
      -al ARC_LABEL, --arc_label ARC_LABEL
                            The objet name for arc files
                            
      -dfl DOME_FLAT_LABEL, --dome_flat_label DOME_FLAT_LABEL
                            The objet name for dome flat files
                            
      -tfl TWILIGHT_FLAT_LABEL, --twilight_flat_label TWILIGHT_FLAT_LABEL
                            The objet name for twilight flat files
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
    
    $ antigen_read_manifest_virus2.py \
        -f ./config_files/virus2/virus2_manifest_template.yml \
        --validate
    ```


- Example 2: How to run a reduction using an existing manifest (recommended!)

    ```bash
    $ conda activate env_antigen
    $ antigen_reduce_virus2_manifest.py \ 
             --manifest_file ./config_files/virus2/virus2_manifest_template.yml \
             --verbose \
             --debug
    ```

- Example 3: How to build a manifest and reduce it, without pre-existing manifest (worst case, avoid if possible)

    ```bash
    $ conda activate env_antigen
    $ dir_containing_VIRUS2=/home/user/my_data
    $ antigen_reduce_virus2.py \ 
            -i $dir_containing_VIRUS2 \
            -c 20250619 \ 
            -f flatp \ 
            --obs_name standard \ 
            -w 2 \ 
            --verbose
    ```
