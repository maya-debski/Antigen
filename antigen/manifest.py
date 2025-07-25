"""
Purpose: Module for handling reduction input manifest files
Uses a yaml file to store information likely only known by observer that are needed to process reduction
Provides reader, structure validate, file validation

Example MANIFEST yaml File Format:

reduction_name: antigen_test_02
in_folder: ~/Antigen/tests/test_02/
out_folder: reduction
unit_instrument: VIRUS2
unit_id: D3G
unit_date: 20250618
obs_date: 20250619
obs_name: M57
observation_files:
  - VIRUS2/20250619/0000099/D3G/VIRUS2_20250620_0000017_standard_D3G_exp01_20250621T040909.6_HR4963.fits
calibration_files:
  bias:
    - VIRUS2/20250619/0000029/D3G/VIRUS2_20250619_0000029_bias_D3G_exp01_20250620T021157.0_zero.fits
    - VIRUS2/20250619/0000029/D3G/VIRUS2_20250619_0000029_bias_D3G_exp02_20250620T021608.9_zero.fits
    - VIRUS2/20250619/0000029/D3G/VIRUS2_20250619_0000029_bias_D3G_exp03_20250620T021655.1_zero.fits
    - VIRUS2/20250619/0000029/D3G/VIRUS2_20250619_0000029_bias_D3G_exp04_20250620T021744.5_zero.fits
    - VIRUS2/20250619/0000029/D3G/VIRUS2_20250619_0000029_bias_D3G_exp05_20250620T021833.9_zero.fits
  flat:
    - VIRUS2/20250619/0000030/D3G/VIRUS2_20250619_0000030_flatp_D3G_exp01_20250620T022204.3_twilight.fits
  arc:
    - VIRUS2/20250619/0000025/D3G/VIRUS2_20250619_0000025_arc_D3G_exp01_20250620T010914.6_fear.fits
    - VIRUS2/20250619/0000025/D3G/VIRUS2_20250619_0000025_arc_D3G_exp02_20250620T011523.7_fear.fits
    - VIRUS2/20250619/0000025/D3G/VIRUS2_20250619_0000025_arc_D3G_exp03_20250620T012133.2_fear.fits
"""

from collections import OrderedDict
from pathlib import Path
import yaml


MANIFEST_PARAMETER_KEYS = ('reduction_name', 'in_folder', 'out_folder',
                           'unit_instrument', 'unit_id', 'unit_date',
                           'obs_date', 'obs_name')
MANIFEST_OBSERVATION_LIST_KEY = 'observation_files'
MANIFEST_CALIBRATION_LIST_KEY = 'calibration_files'
MANIFEST_CALIBRATION_TYPE_KEYS = ('bias', 'flat', 'arc')


def load_manifest(filename):
    """
    Purpose: Load a YAML config and resolve input file paths.
    Args:
        filename (str): Full-path filename of manifest file
    Returns:
        manifest (dict): data structure loaded from manifest YAML file
    """
    with open(filename, 'r') as f:
        manifest = ordered_yaml_loader(f)
    return manifest


def ordered_yaml_loader(stream):
    class OrderedLoader(yaml.SafeLoader):
        pass

    def construct_mapping(loader, node):
        loader.flatten_mapping(node)
        return OrderedDict(loader.construct_pairs(node))

    OrderedLoader.add_constructor(
        yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
        construct_mapping)

    return yaml.load(stream, OrderedLoader)


def ordered_yaml_dumper(data, stream, **kwargs):
    class OrderedDumper(yaml.SafeDumper):
        pass

    def represent_ordered_dict(dumper, data):
        return dumper.represent_mapping(
            yaml.resolver.BaseResolver.DEFAULT_MAPPING_TAG,
            data.items())

    OrderedDumper.add_representer(OrderedDict, represent_ordered_dict)

    # Default: block-style output like safe_dump
    yaml.dump(data, stream, OrderedDumper, sort_keys=False, default_flow_style=False, **kwargs)


def stringify_paths(obj):
    if isinstance(obj, dict):
        return {k: stringify_paths(v) for k, v in obj.items()}
    elif isinstance(obj, list):
        return [stringify_paths(i) for i in obj]
    elif isinstance(obj, Path):
        return str(obj)
    else:
        return obj


def save_manifest(manifest, filename):
    """
    Purpose: Save a manifest dict to a YAML file
    Args:
        manifest (dict): data structure loaded from manifest YAML file
        filename (str): Full-path filename of manifest file
    Returns:
        None
    """
    path = Path(filename).expanduser()
    with open(path, 'w') as fob:
        ordered_yaml_dumper(stringify_paths(manifest), fob)
        # yaml.safe_dump(stringify_paths(manifest), fob, default_flow_style=False)


def normalize_manifest(manifest):
    """
    Purpose: Replace file strings and relative file paths with fully resolved file paths as pathlib.Path objects.

    Args:
        manifest (dict): data structure loaded from manifest YAML file
    Returns:
        manifest (dict): same structure, but with all file paths updated to prepend 'in_folder'
    """
    manifest['in_folder'] = Path(manifest['in_folder']).expanduser().resolve(strict=False)

    calibration_files = manifest[MANIFEST_CALIBRATION_LIST_KEY]
    for category, relative_paths in calibration_files.items():
        calibration_files[category] = [manifest['in_folder'] / Path(file) for file in relative_paths]
    manifest[MANIFEST_CALIBRATION_LIST_KEY] = calibration_files

    obs_files_relative_paths = manifest[MANIFEST_OBSERVATION_LIST_KEY]
    observation_files = [manifest['in_folder'] / Path(file) for file in obs_files_relative_paths]
    manifest[MANIFEST_OBSERVATION_LIST_KEY] = observation_files
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

    # Top-level key check
    for key in MANIFEST_PARAMETER_KEYS:
        if key not in manifest:
            raise ValueError(f"Missing top-level config key: '{key}'")

    # Calibration file categories
    for category in MANIFEST_CALIBRATION_TYPE_KEYS:
        if category not in manifest[MANIFEST_CALIBRATION_LIST_KEY]:
            raise ValueError(f"Missing file category in calibration_files: '{category}'")

    # Ensure in_folder directory exists:
    in_folder_path = manifest['in_folder']
    if not isinstance(in_folder_path, Path):
        raise TypeError(f"manifest['in_folder'] contains non-Path item: {in_folder_path}")
    if not in_folder_path.exists():
        raise FileNotFoundError(f"Path does not exist: {in_folder_path}")
    if not in_folder_path.is_dir():
        raise ValueError(f"Path is not a directory: {in_folder_path}")

    # Ensure observation files are a list of valid files
    for file_path in manifest[MANIFEST_OBSERVATION_LIST_KEY]:
        if not isinstance(file_path, Path):
            raise TypeError(f"manifest[{MANIFEST_OBSERVATION_LIST_KEY}] contains non-Path item: {file_path}")
        if not file_path.exists():
            raise FileNotFoundError(f"File does not exist: {file_path}")
        if not file_path.is_file():
            raise ValueError(f"Path is not a file: {file_path}")

    # Ensure calibration files are a list of valid files
    for category, file_list in manifest[MANIFEST_CALIBRATION_LIST_KEY].items():
        if not isinstance(file_list, list):
            raise TypeError(f"input_files[{category}] must be a list")
        for file_path in file_list:
            if not isinstance(file_path, Path):
                raise TypeError(f"input_files[{category}] contains non-Path item: {file_path}")
            if not file_path.exists():
                raise FileNotFoundError(f"File does not exist: {file_path}")
            if not file_path.is_file():
                raise ValueError(f"Path is not a file: {file_path}")

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
    for key in MANIFEST_PARAMETER_KEYS:
        print(f"{key} : {manifest.get(key, '<unnamed>')}")

    observation_files = manifest.get(MANIFEST_OBSERVATION_LIST_KEY, {})
    print("observation_files:")
    for file in observation_files:
        print(f"  - {file}")

    calibration_file_dict = manifest.get(MANIFEST_CALIBRATION_LIST_KEY, {})
    print("calibration_files:")
    for category, files in calibration_file_dict.items():
        print(f"  - {category}: (file count = {len(files)}):")
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

