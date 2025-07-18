"""
Purpose: Module for handling reduction input manifest files
Uses a yaml file to store information likely only known by observer that are needed to process reduction
Provides reader, structure validate, file validation

Example MANIFEST yaml File Format:

    reduction_name: test_01
    in_folder: ./
    out_folder: reduc

    parameters:
      instrument: VIRUS2
      obs_date: 20250619
      obs_name: M81

    input_files:
      bias:
        - VIRUS2/20250619/0000029/D3G/VIRUS2_20250619_0000029_bias_D3G_exp01_20250620T021157.0_zero.fits
        - VIRUS2/20250619/0000029/D3G/VIRUS2_20250619_0000029_bias_D3G_exp03_20250620T021655.1_zero.fits
      flat:
        - VIRUS2/20250619/0000030/D3G/VIRUS2_20250619_0000030_flatp_D3G_exp01_20250620T022204.3_twilight.fits
      arc:
        - VIRUS2/20250619/0000025/D3G/VIRUS2_20250619_0000025_arc_D3G_exp01_20250620T010914.6_fear.fits
        - VIRUS2/20250619/0000025/D3G/VIRUS2_20250619_0000025_arc_D3G_exp02_20250620T011523.7_fear.fits
"""

from pathlib import Path
import yaml


def load_manifest(filename):
    """
    Purpose: Load a YAML config and resolve input file paths.
    Args:
        filename (str): Full-path filename of manifest file
    Returns:
        manifest (dict): data structure loaded from manifest YAML file
    """
    with open(filename, 'r') as f:
        manifest = yaml.safe_load(f)
    return manifest


def normalize_manifest(manifest):
    """
    Purpose: Replace relative file paths with fully resolved Path objects.
    Args:
        manifest (dict): data structure loaded from manifest YAML file
    Returns:
        manifest (dict): same structure, but with manifest['input_files'] paths updated to prepend 'input_root'
    """
    input_files = manifest['input_files']
    for category, relative_paths in input_files.items():
        input_files[category] = [Path(manifest['in_folder']) / Path(file) for file in relative_paths]
    manifest['input_files'] = input_files
    return manifest


def validate_manifest(manifest):
    """
    Purpose: Validate presence of required keys and check that all input files exist and are files.

    Args:
        manifest (dict): dict returned from yaml.safe_load(manifest_filename)
    Returns:
        pass (bool): returned True only if all tests passed, else an exception is raised
    Raises:
        ValueError
        TypeError
    """
    required_top_level = ['reduction_name', 'in_folder', 'out_folder', 'parameters', 'input_files']
    required_parameters = ['instrument', 'obs_date', 'obs_name']
    required_input_file_categories = ['bias', 'flat', 'arc']

    # Top-level key check
    for key in required_top_level:
        if key not in manifest:
            raise ValueError(f"Missing top-level config key: '{key}'")

    # File categories
    input_files = manifest['input_files']
    for category in required_input_file_categories:
        if category not in input_files:
            raise ValueError(f"Missing file category in input_files: '{category}'")

    # Ensure input files are a list of valid files
    for category, file_list in input_files.items():
        if not isinstance(file_list, list):
            raise TypeError(f"input_files[{category}] must be a list")
        for file_path in file_list:
            if not isinstance(file_path, Path):
                raise TypeError(f"input_files[{category}] contains non-Path item: {file_path}")
            if not file_path.exists():
                raise FileNotFoundError(f"File does not exist: {file_path}")
            if not file_path.is_file():
                raise ValueError(f"Path is not a file: {file_path}")

    # Parameter key check
    for param in required_parameters:
        if param not in manifest['parameters']:
            raise ValueError(f"Missing parameter: '{param}'")

    return True


def print_manifest(manifest):
    """
    Purpose: Pretty-print summary of the loaded and normalized config.
    Args:
        manifest (dict): data structure loaded from manifest YAML file
    Returns:
        None
    """
    print('\n')
    print('Reduction Manifest:')
    print(f"- Reduction Name : {manifest.get('reduction_name', '<unnamed>')}")
    print(f"- In Folder  : {manifest.get('in_folder', '<unnamed>')}")
    print(f"- Out Folder : {manifest.get('out_folder', '<unnamed>')}")

    params = manifest.get('parameters', {})
    print(f"- Parameters:")
    print(f"  - Instrument Name  : {params.get('instrument', 'N/A')}")
    print(f"  - Observation Name : {params.get('obs_name', 'N/A')}")
    print(f"  - Observation Date : {params.get('obs_date', 'N/A')}")

    input_files = manifest.get('input_files', {})
    print("- Input Files:")
    for category, files in input_files.items():
        print(f"  - {category} (file count = {len(files)}):")
        for file in files:
            print(f"    - {file}")
    print()
    return None


def read_manifest(filename, validate=False, verbose=False):
    """
    Purpose: Read, Normalize, and validate manifest file

    Args:
        filename (str): Full-path filename of manifest file
        validate (bool): if True, test manifest structure and existence and all files listed therein
        verbose (bool): if True, pretty-print manifest structure to console
    Returns:
        manifest_normalized (dict): data structure loaded from manifest YAML file
    """
    manifest_raw = load_manifest(filename)
    manifest_normalized = normalize_manifest(manifest_raw)
    if verbose:
        print_manifest(manifest_normalized)
    if validate:
        valid = validate_manifest(manifest_normalized)
        if verbose:
            result = 'PASS' if valid else 'FAIL'
            result_message = f'{result}: manifest validation: valid={valid}, filename={filename}'
            print(result_message + '\n')
    return manifest_normalized

