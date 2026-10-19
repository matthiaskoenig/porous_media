"""Module for logging.

Using rich for output formating.
"""

import logging

from rich.logging import RichHandler

from porous_media.console import console


def get_logger(name: str, level: int = logging.INFO) -> logging.Logger:
    """Get a custom logger without adding duplicate console handlers."""
    logger = logging.getLogger(name)
    logger.setLevel(level)
    if any(isinstance(handler, RichHandler) for handler in logger.handlers):
        return logger
    formatter = logging.Formatter(
        fmt="%(message)s",
        datefmt="[%X]",
    )

    # handler = logging.StreamHandler()
    handler = RichHandler(
        markup=False, rich_tracebacks=True, show_time=False, console=console
    )
    handler.setFormatter(formatter)

    logger.addHandler(handler)
    return logger
