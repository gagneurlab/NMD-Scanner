import logging
from .scanner import (
    annotate_nmd
)

def setup_logging(verbosity=0):
    """
    Configure logging for the nmd_scanner package.
    
    :param verbosity: 0 = WARNING, 1 = INFO, 2+ = DEBUG
    """
    if verbosity == 0:
        level = logging.WARNING
    elif verbosity == 1:
        level = logging.INFO
    else:
        level = logging.DEBUG
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    
    # Get the package logger
    logger = logging.getLogger('nmd_scanner')
    logger.setLevel(level)
    
    return logger