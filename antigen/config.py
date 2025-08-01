from importlib.resources import files
from pathlib import Path
import logging
import yaml

import numpy as np
from astropy.table import Table

VIRUS2_PCA_COMPONENTS = 15
VIRUS2_FIBER_REF = 130
VIRUS2_FIBER_RADIUS = 2.483 / 2.

DEFAULT_PCA_COMPONENTS_SKY = 5
DEFAULT_PCA_COMPONENTS_ARC = 15
DEFAULT_SKY_CONTINUUM_BINS = 50


def get_channel_config_virus2():
    CHANNEL_DETECTOR = {'g': {'gain': 2.017, 'rdnoise': 3.09, 'limit': 5},
                        'b': {'gain': None, 'rdnoise': None, 'limit': None},
                        'r': {'gain': None, 'rdnoise': None, 'limit': None},
                        'd': {'gain': None, 'rdnoise': None, 'limit': None}}

    CHANNEL_DEF_WAVE = {'g': np.linspace(4610., 5925., 2064),
                        'b': np.linspace(3700., 4630., 2064),
                        'r': np.linspace(5900., 7610., 2064),
                        'd': np.linspace(7590., 9300., 2064)}

    return CHANNEL_DETECTOR, CHANNEL_DEF_WAVE

def get_base_config_path():
    """
    Purpose: Returns the base path to the config_files directory.

    Returns:
        base_path: Path object pointing to the root of config_files/.
    """
    base_path = files(__package__) / "config_files"
    return base_path


def get_config_filepath(*path_parts):
    """
    Purpose: Returns the filename path to a specific config file under config_files/.

    Args:
        *path_parts (str): One or more subpaths under config_files/.
                           For example: get_config_file("lines", "virus2_green.txt")

    Returns:
        config_file_path (Path): A Path object pointing to the config file.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    installed_config_path_base = get_base_config_path()
    config_file_path = installed_config_path_base.joinpath(*path_parts)
    if not config_file_path.exists():
        raise FileNotFoundError(f"Config file does not exist: {config_file_path}")
    if not config_file_path.is_file():
        raise FileNotFoundError(f"Config file is not a file: {config_file_path}")
    return config_file_path


def read_config_file(*path_parts):
    """
    Purpose: Reads and returns the contents of a config file_path as a string.

    Args:
        *path_parts (str): *path_parts (str): One or more subpaths under config_files/,
                           for example pass in path_parts=('lines', 'virus2_green.txt')

    Returns:
        str: The contents of the config file_path.
    """
    full_file_path = get_config_filepath(*path_parts)
    with full_file_path.open("r", encoding="utf-8") as f:
        file_string_contents = f.read()
    return file_string_contents


def list_available_config_files(subdir=""):
    """
    Purpose: Get a list of all available config files by searching a subdirectory under config_files/

    Args:
        subdir (str): The subdirectory to scan (e.g. 'instruments').
            Defaults to the top-level config_files directory.

    Returns:
        list[str]: A list of filenames in the subdirectory.
    """
    dir_path = get_base_config_path().joinpath(subdir)
    files_found = [entry.name for entry in dir_path.iterdir() if entry.is_file()]
    return files_found


logger = logging.getLogger('antigen.config')

CONFIG_YAML_PARAMETER_KEYS = ('instrument_element', 'gain', 'read_noise',
                              'fiber_radius', 'pca_components_arc', 'reference_fiber_index',
                              'sample_fiber_indices', 'start_wavelength', 'end_wavelength',
                              'detector_dimensions', 'arc_flux_limit', 'telescope_diameter',
                              'overscan_length', 'flip_x', 'flip_y', 'rotate')

CONFIG_IFUCEN_PARAMETER_KEYS = ('fiber_id', 'head_id', 'ifu_x', 'ifu_y',
                                'trace_row', 'exclude_fiber')

CONFIG_LINES_PARAMETER_KEYS = ('wavelength', 'column')

CONFIG_THROUGHPUT_PARAMETER_KEYS = ('rectified_wavelength', 'throughput')

def get_base_config_path_new(instrument, instrument_element):
    """
    Purpose: Returns the base path to the config_files directory.
    Args:
        instrument (str): The instrument used in observations (e.g., gcms, virus2, virusw)
        instrument_element (str): The instrument element used in observations (e.g., VP1R, VP1B, D3G)
    Returns:
        base_path: Path object pointing to the root of config_files/.
    """
    base_path = files(__package__) / "config_files" / f"{instrument}" / f"{instrument_element}"
    return base_path

def get_config_filepath_new(instrument, instrument_element, config_file):
    """
    Purpose: Returns the filename path to a specific config file under config_files/.

    Args:
        instrument (str): The instrument used in observations (e.g., gcms, virus2, virusw)
        instrument_element (str): The instrument element used in observations (e.g., VP1R, VP1B, D3G)
        config_file (str): The exact config file to open (e.g., ifucen, lines, throughput)

    Returns:
        config_file_path (Path): A Path object pointing to the config file.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    installed_config_path_base = get_base_config_path_new(instrument, instrument_element)
    config_file_path = installed_config_path_base.joinpath(f"{instrument}_{config_file}_{instrument_element}.txt")
    if not config_file_path.exists():
        raise FileNotFoundError(f"Config file does not exist: {config_file_path}")
    if not config_file_path.is_file():
        raise FileNotFoundError(f"Config file is not a file: {config_file_path}")
    return config_file_path

def load_config_yaml(instrument, instrument_element):
    """
    Purpose: Load a YAML config
    Args:
        instrument (str): The instrument used in observations (e.g., gcms, virus2, virusw)
        instrument_element (str): The instrument element used in observations (e.g., VP1R, VP1B, D3G)
    Returns:
        config_yaml (dict): data structure loaded from config YAML file
    """
    base_config_path = get_base_config_path_new(instrument, instrument_element)
    config_file_path = base_config_path.joinpath(f"{instrument}_config_{instrument_element}.yml")
    with open(config_file_path, 'r') as f:
        config_yaml = yaml.safe_load(f)
    return config_yaml

def load_config_file_new(instrument, instrument_element, config_file):
    """
    Purpose: Reads and returns the contents of a config file_path as a string.

    Args:
        instrument (str): The instrument used in observations (e.g., gcms, virus2, virusw)
        instrument_element (str): The instrument element used in observations (e.g., VP1R, VP1B, D3G)
        config_file (str): The exact config file to open (e.g., ifucen, lines, throughput)

    Returns:
        file_contents (astropy Table): The contents of the config file_path.
    """
    full_file_path = get_config_filepath_new(instrument, instrument_element, config_file)
    file_contents = Table.read(full_file_path, format = 'ascii')
    return file_contents

def validate_config_yaml(config_yaml):
    """
    Purpose: Validate presence of required keys

    Args:
        config_yaml (dict): dict returned from yaml.safe_load(config_filename.yml)
    Returns:
        pass (bool): returned True only if all tests passed, else an exception is raised
    Raises:
        ValueError
        TypeError
    """

    # Top-level key check
    for key in CONFIG_YAML_PARAMETER_KEYS:
        if key not in config_yaml:
            raise ValueError(f"Missing top-level config key: '{key}'")
    return True

def validate_config_lines(line_file_contents):
    """
    Purpose: Validate presence of required columns

    Args:
        line_file_contents (astropy Table): The contents of the line_file
    Returns:
        pass (bool): returned True only if all tests passed, else an exception is raised
    Raises:
        ValueError
    """

    column_names = line_file_contents.columns

    for key in CONFIG_LINES_PARAMETER_KEYS:
        if key not in column_names:
            raise ValueError(f"Missing top-level line list config key: '{key}'")

    return True


def validate_config_ifucen(ifucen_file_contents):
    """
    Purpose: Validate presence of required columns

    Args:
        ifucen_file_contents (astropy Table): The contents of the ifucen_file
    Returns:
        pass (bool): returned True only if all tests passed, else an exception is raised
    Raises:
        ValueError
    """

    column_names = ifucen_file_contents.columns

    for key in CONFIG_IFUCEN_PARAMETER_KEYS:
        if key not in column_names:
            raise ValueError(f"Missing top-level ifucen config key: '{key}'")

    return True

def validate_config_throughput(throughput_file_contents):
    """
    Purpose: Validate presence of required columns

    Args:
        throughput_file_contents (astropy Table): The contents of the throughput_file
    Returns:
        pass (bool): returned True only if all tests passed, else an exception is raised
    Raises:
        ValueError
    """

    column_names = throughput_file_contents.columns

    for key in CONFIG_THROUGHPUT_PARAMETER_KEYS:
        if key not in column_names:
            raise ValueError(f"Missing top-level throughput config key: '{key}'")

    return True

def read_config_yaml(instrument, instrument_element, validate=False):
    """
    Purpose: Read and validate config YAML file

    Args:
        instrument (str): The instrument used in observations (e.g., gcms, virus2, virusw)
        instrument_element (str): The instrument element used in observations (e.g., VP1R, VP1B, D3G)
        validate (bool): if True, test config YAML structure
    Returns:
        config_yaml (dict): data structure loaded from config YAML file
    """

    config_yaml = load_config_yaml(instrument, instrument_element)
    if validate:
        valid = validate_config_yaml(config_yaml)
        result = 'PASS' if valid else 'FAIL'
        logger.info(f'{result}: config YAML validation: valid={valid}')
    return config_yaml

def read_config_files(instrument, instrument_element, validate=False):
    """
    Purpose: Read and validate config files for lines, throughput, and ifucen

    Args:
        instrument (str): The instrument used in observations (e.g., gcms, virus2, virusw)
        instrument_element (str): The instrument element used in observations (e.g., VP1R, VP1B, D3G)
        validate (bool): if True, test config file structure
    Returns:
        config_lines (astropy Table): file contents loaded from line list file
        config_ifucen (astropy Table): ifucen file contents loaded from ifucen file
        config_throughput (astropy Table): throughput file contents loaded from throughput
    """

    config_lines = load_config_file_new(instrument, instrument_element, "lines")
    config_ifucen = load_config_file_new(instrument, instrument_element, "ifucen")
    config_throughput = load_config_file_new(instrument, instrument_element, "throughput")
    if validate:
        valid_lines = validate_config_lines(config_lines)
        result = 'PASS' if valid_lines else 'FAIL'
        logger.info(f'{result}: config line list validation: valid={valid_lines}')

        valid_ifucen = validate_config_ifucen(config_ifucen)
        result = 'PASS' if valid_ifucen else 'FAIL'
        logger.info(f'{result}: config ifucen file validation: valid={valid_ifucen}')

        valid_throughput = validate_config_throughput(config_throughput)
        result = 'PASS' if valid_throughput else 'FAIL'
        logger.info(f'{result}: config throughput file validation: valid={valid_throughput}')
    return config_lines, config_ifucen, config_throughput

def build_config_for_element(instrument, instrument_element,validate = False):
    """
    Purpose: Build configuration object for a given instrument and instrument element
    
    Args:
        instrument (str): The instrument used in observations (e.g., gcms, virus2, virusw)
        instrument_element (str): The instrument element used in observations (e.g., VP1R, VP1B, D3G)
        validate (bool): if True, test config file and YAML structures
    Returns:
        config_dict (dict): Python dictionary containing all configuration data from config files
    """

    config_yaml = read_config_yaml(instrument, instrument_element, validate)
    config_lines, config_ifucen, config_throughput = read_config_files(instrument, instrument_element, validate)

    config_dict = {'instrument_element': config_yaml['instrument_element'],
                   'gain': config_yaml['gain'],
                   'read_noise': config_yaml['read_noise'],
                   'fiber_radius': config_yaml['fiber_radius'],
                   'pca_components_arc': config_yaml['pca_components_arc'],
                   'reference_fiber_index': config_yaml['reference_fiber_index'],
                   'start_wavelength': config_yaml['start_wavelength'],
                   'end_wavelength': config_yaml['end_wavelength'],
                   'detector_dimensions': config_yaml['detector_dimensions'],
                   'arc_flux_limit': config_yaml['arc_flux_limit'],
                   'telescope_diameter': config_yaml['telescope_diameter'],
                   'overscan_length': config_yaml['overscan_length'],
                   'flip_x': config_yaml['flip_x'],
                   'flip_y': config_yaml['flip_y'],
                   'rotate': config_yaml['rotate'],
                   'sample_fiber_indices': config_yaml['sample_fiber_indices'],
                   'fiber_id': np.array(config_ifucen['fiber_id']),
                   'head_id': np.array(config_ifucen['head_id']),
                   'ifu_x': np.array(config_ifucen['ifu_x']),
                   'ifu_y': np.array(config_ifucen['ifu_y']),
                   'trace_row': np.array(config_ifucen['trace_row']),
                   'exclude_fiber': np.array(config_ifucen['exclude_fiber']),
                   'wavelength': np.array(config_lines['wavelength']),
                   'column': np.array(config_lines['column']),
                   'throughput': np.array(config_throughput['throughput'])}

    return config_dict





