# Antigen
## Reduction Pipelines for GCMS (VIRUS-P) and VIRUS-W

Antigen is designed to reduce data from the George and Cynthia Mitchell Spectrograph (GCMS) and VIRUS-W on the 2.7m Harlan J. Smith Telescope at McDonald Observatory. In its current state, Antigen  outputs a fits file for each science exposure in a night of data taken with GCMS that is the fiber-extracted, wavelength-calibrated, but "raw" spectra.  

An example call might be:

python reduce_virusp.py ~/Downloads/VIRUS-P_Data/20240606 ~/Downloads/VIRUS-P_Data/20240606/reduced -bn -ra -bl "zero" -tfl "Twilight flat" -al "comp" -b

The arguments for the script are:

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
