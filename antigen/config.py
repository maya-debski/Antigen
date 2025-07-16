from importlib.resources import files
from pathlib import Path


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
        Path: A Path object pointing to the config file.

    Raises:
        FileNotFoundError: If the file does not exist.
    """
    config_filepath = get_base_config_path().joinpath(*path_parts)
    if not config_filepath.is_file():
        raise FileNotFoundError(f"Config file not found: {config_filepath}")
    return config_filepath


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
