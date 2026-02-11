"""Logging configuration for SDRF Refiner."""

import logging
import sys


def setup_logging(verbose: int = 0) -> logging.Logger:
    """
    Set up logging configuration based on verbosity level.

    Args:
        verbose: Verbosity level (0=WARNING, 1=INFO, 2+=DEBUG)

    Returns:
        Configured logger instance
    """
    if verbose == 0:
        log_level = logging.WARNING
    elif verbose == 1:
        log_level = logging.INFO
    else:
        log_level = logging.DEBUG

    # Configure root logger
    logging.basicConfig(
        level=log_level,
        format="%(levelname)s: %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    # Get logger for this package
    logger = logging.getLogger("sdrf_refiner")
    logger.setLevel(log_level)

    return logger


def get_logger(name: str) -> logging.Logger:
    """Get a logger instance for a module."""
    return logging.getLogger(f"sdrf_refiner.{name}")
