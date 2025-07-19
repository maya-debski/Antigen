import logging


def setup_logging(log_name='antigen.logger'):
    """
    Purpose: Set up a logger object with specific log level time and message string format
    Note: Use a StreamHandler to write to stdout and set the level to DEBUG if verbose is set from the command line

    Args:
        log_name (str): name of logger

    Returns:
        log (logging.Logger): Logger object
    """
    log = logging.getLogger(log_name)
    if not len(log.handlers):
        fmt = '[%(levelname)s - %(asctime)s] %(message)s'
        fmt = logging.Formatter(fmt)

        level = logging.INFO

        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        handler.setLevel(level)

        log = logging.getLogger('input_utils')
        log.setLevel(logging.DEBUG)
        log.addHandler(handler)
    return log
