import logging

logger = logging.getLogger('antigen.utils')

def setup_logging(log_name='antigen', debug=False, verbose=False):
    """
    Purpose: Set up a logger object with specific log level time and message string format
    Note: Use a StreamHandler to write to stdout and set the level to DEBUG if verbose is set from the command line

    Args:
        log_name (str): name of logger

    Returns:
        log (logging.Logger): Logger object
    """

    logger = logging.getLogger(log_name)

    if debug:
        level = logging.DEBUG
    elif verbose:
        level = logging.INFO
    else:
        level = logging.WARNING
    logger.setLevel(level)

    if not logger.handlers:
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            "[%(asctime)s] [%(name)s] %(levelname)s: %(message)s",
            datefmt="%H:%M:%S"
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger

def validate_inputs(dataset_manifest, args, config_dict):
    """Validate required inputs for building the instrument response function.

    This function checks that the dataset manifest, command-line arguments,
    and configuration dictionary contain the attributes needed for downstream
    processing. If any are missing or invalid, a ValueError is raised with a
    descriptive message. Validation progress is logged.

    Args:
        dataset_manifest (dict): Dictionary containing dataset information.
            Must include:
                - 'unit_instrument' (str)
                - 'unit_id' (str)
                - 'ndithers' (int)
                - 'dither_number' (int)
                - 'in_folder' (str)
                - 'reduced_files' (list[str])
        args (argparse.Namespace): Parsed command-line arguments.
            Must include:
                - extraction_radius (float)
                - pixel_size (float)
        config_dict (dict): Instrument configuration dictionary.
            Must include:
                - 'fiber_radius' (float)

    Raises:
        ValueError: If any required key or attribute is missing, None, or invalid.
    """
    # Dataset manifest checks
    required_manifest_keys = [
        "unit_instrument", "unit_id", "ndithers", "dither_number",
        "in_folder", "reduced_files"
    ]
    for key in required_manifest_keys:
        if key not in dataset_manifest:
            logger.error(f"Dataset manifest missing required key: '{key}'")
            raise ValueError(f"Missing required dataset_manifest key: '{key}'")
        logger.debug(f"Dataset manifest key '{key}' present: {dataset_manifest[key]}")

    # Args checks
    required_args = ["extraction_radius", "pixel_size"]
    for attr in required_args:
        if not hasattr(args, attr):
            logger.error(f"Args missing required attribute: '{attr}'")
            raise ValueError(f"Missing required argument: '{attr}'")
        value = getattr(args, attr)
        if value is None:
            logger.error(f"Argument '{attr}' is None")
            raise ValueError(f"Argument '{attr}' cannot be None")
        logger.debug(f"Argument '{attr}' validated: {value}")

    # Config dict checks
    required_config_keys = ["fiber_radius"]
    for key in required_config_keys:
        if key not in config_dict:
            logger.error(f"Config dict missing required key: '{key}'")
            raise ValueError(f"Missing required config_dict key: '{key}'")
        if config_dict[key] is None:
            logger.error(f"Config value '{key}' is None")
            raise ValueError(f"Config value '{key}' cannot be None")
        logger.debug(f"Config dict key '{key}' validated: {config_dict[key]}")

    # Type / value sanity checks
    if args.extraction_radius <= 0:
        logger.error("Invalid extraction_radius: must be positive")
        raise ValueError("extraction_radius must be positive")
    if args.pixel_size <= 0:
        logger.error("Invalid pixel_size: must be positive")
        raise ValueError("pixel_size must be positive")
    if config_dict["fiber_radius"] <= 0:
        logger.error("Invalid fiber_radius: must be positive")
        raise ValueError("fiber_radius must be positive")

    logger.info("All inputs successfully validated.")