"""Pipeline logging setup."""

from __future__ import annotations

import logging
import sys
from pathlib import Path


def configure_logging(log_path: Path, verbose: bool = False) -> logging.Logger:
    logger = logging.getLogger("ulm_3d_vascular")
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG)
    logger.propagate = False

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s | %(levelname)-8s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S")
    )

    console_handler = logging.StreamHandler(sys.stdout)
    console_handler.setLevel(logging.DEBUG if verbose else logging.INFO)
    console_handler.setFormatter(logging.Formatter("%(levelname)-8s %(message)s"))

    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    logging.captureWarnings(True)
    return logger

