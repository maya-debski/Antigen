import logging


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
