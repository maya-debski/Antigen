# Configuration and Packaged Resources

Antigen bundles instrument configuration resources with the Python package under antigen/config_files. At runtime, code locates these files via importlib.resources, so you normally do not need to manage absolute paths.

What’s inside config_files repository:

The config_files repository is organized into sub-folders for each instrument (GCMS, VIRUS2, VIRUSW), with additional sub-folders for each instrument element (i.e., VP1B, etc.).
Within each sub-folder for each instrument element, there is:

1) Instrument element config YAML files (used to define element-specific information, i.e., wavelength coverage, detector dimensions, etc.)
- gcms/VP1B/gcms_config_VP1B.yml
- gcms/VP1R/gcms_config_VP1R.yml
- virus2/D3G/virus2_config_D3G.yml
- virusw/HRES/virusw_config_HRES.yml
- virusw/LRES/virusw_config_LRES.yml

Where used: **FILL IN**

2) IFU fiber-center mappings (used for spatial geometry)
- gcms/VP1B/gcms_ifucen_VP1B.txt
- gcms/VP1R/gcms_ifucen_VP1R.txt
- virus2/D3G/virus2_ifucen_D3G.txt
- virusw/HRES/virusw_ifucen_HRES.txt
- virusw/LRES/virusw_ifucen_LRES.txt

Where used: **FILL IN**

3) Arc line lists (used for wavelength calibration)
- gcms/VP1B/gcms_lines_VP1B.txt
- gcms/VP1R/gcms_lines_VP1R.txt
- virus2/D3G/virus2_lines_D3G.txt
- virusw/HRES/virusw_lines_HRES.txt
- virusw/LRES/virusw_lines_LRES.txt

Where used: **FILL IN**

4) Throughput curves (used for ...)
- gcms/VP1B/gcms_throughput_VP1B.txt
- gcms/VP1R/gcms_throughput_VP1R.txt
- virus2/D3G/virus2_throughput_D3G.txt
- virusw/HRES/virusw_throughput_HRES.txt
- virusw/LRES/virusw_throughput_LRES.txt

Where used: **FILL IN**

Additionally, within each instrument sub-folder, there is:

1) Dither files (used for ...)
- gcms/gcms_dither_3pt.lis
- gcms/gcms_dither_6pt.lis
- virusw/virusw_dither_3pt.lis
- virusw/virusw_dither_6pt.lis

Where used: **FILL IN**

Also in the config_files repository:

1) McDonald Observatory extinction table (used for ...)
- extinction/mcdonald_extinction.dat

Where used: **FILL IN**

2) Recipes (used for defining inputs, workflow, and outputs for a given antigen script)
- recipes/base_reduction.yml
- recipes/make_cubes.yml
- recipes/measure_response.yml

Where used: antigen_*.py scripts with same name (antigen_base_reduction.py, antigen_make_cubes.py, antigen_measure_response.py)

3) Schema (used for defining validation rules for recipe inputs)
- schema/base_reduction_schema.yml
- schema/make_cubes_schema.yml
- schema/measure_response_schema.yml

Where used: antigen_*.py scripts with same name

4) Operations (used to define implementation details for each step of a recipe workflow)
- schema/base_reduction_operations.yml
- schema/make_cubes_operations.yml
- schema/measure_response_operations.yml

Where used: antigen_*.py scripts with same name

Overriding the packaged defaults
- Within Python: You can open and use your own files, but for the CLI quicklook, the simplest “override” is to place a file with the same name in the installed package data. This is easiest in an editable/development install. Example for replacing the VP1B arc lines:
  cp my_lines_VP1B.dat antigen/config_files/gcms/VP1B/gcms_lines_VP1B.txt
  (Reinstall or rebuild the package as needed.)

File formats (summary)
- Instrument element config files (\*\_config\_\*.yml): YAML with keys [instrument_element, gain, read_noise, ...]; **HOW IT IS READ**
- Mapping files (\*\_ifucen\_\*.txt): ASCII with columns [fiber_id, head_id, ifu_x, ifu_y, trace_row, exclude_fiber]; **HOW IT IS READ**
- Arc line lists (\*\_lines\_\*.txt): ASCII with columns [wavelength, column]; **HOW IT IS READ**
- Throughput curves (\*\_throughput\_\*.txt): ASCII with columns [rectified_wavelength, throughput]; **HOW IT IS READ**
- Dither files (\*\_dither\_\*pt.lis): Whitespace-separated columns describing change in IFU x and y center positions for each dither; **HOW IT IS READ**
- Extinction table (mcdonald_extinction.dat): ASCII with columns [wavelength, mag_per_airmass]; **HOW IT IS READ**
- Recipes: YAML with keys [description, inputs, workflow, outputs]; **HOW IT IS READ**
- Schema: YAML with keys [cli, dataset, config];**HOW IT IS READ**
- Operations: YAML with keys [operations]; **HOW IT IS READ**


See also
- [CLI Usage](../user-guide/cli.md)
- [API Reference for utils.get_config_file (**this needs to be changed for Antigen**) and related helpers](../api/index.md)
- [Recipe Building](../recipes/Recipe_Building.md)