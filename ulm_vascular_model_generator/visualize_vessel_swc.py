#!/usr/bin/env python
"""
Visualize the generated SWC vessel tree as a cylinder mesh.

By default, the newest timestamped generator result is selected automatically.
Pass ``--result-folder PATH`` to pin a specific run. The exported mesh,
quantitative figure, and caption are written back to the selected folder.
Rendering settings come from ``configs/visualization_config.yaml``.
"""

from __future__ import annotations
from pathlib import Path
import sys
if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from ulm_vascular_model_generator.utils.visualization.cli import (
    main,
    newest_result_folder,
    parse_args,
)


# Keep this as None for automatic newest-result selection. Set a folder name
# only when the script should remain pinned to one historical result.
RESULT_FOLDER_NAME = None       # '20260729_020617_759618_volumetric_3d_seed_105'
RESULT_FOLDER = (
    None
    if RESULT_FOLDER_NAME is None
    else Path(__file__).resolve().parent
    / "vessel_swc_models"
    / RESULT_FOLDER_NAME
)

__all__ = [
    "RESULT_FOLDER",
    "RESULT_FOLDER_NAME",
    "main",
    "newest_result_folder",
    "parse_args",
]

if __name__ == "__main__":
    main(result_folder=RESULT_FOLDER)
