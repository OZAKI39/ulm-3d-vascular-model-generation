"""Create one reproducible output directory for each generator run."""

from datetime import datetime
from pathlib import Path
import shutil


def create_run_directory(output_root, seed, geometry_mode):
    """Create a timestamped run folder whose name states its geometry mode."""
    timestamp = datetime.now().astimezone().strftime("%Y%m%d_%H%M%S_%f")
    mode = str(geometry_mode).strip().lower()
    if mode not in {"planar_2d", "volumetric_3d"}:
        raise ValueError(
            "geometry_mode must be 'volumetric_3d' or 'planar_2d' when "
            "creating a generator result directory."
        )
    run_directory = Path(output_root) / (
        f"{timestamp}_{mode}_seed_{int(seed)}"
    )
    run_directory.mkdir(parents=True, exist_ok=False)
    return run_directory


def sync_swc_to_output_root(run_swc_path, output_root_swc_path):
    """Copy the run SWC to the stable latest-model path."""
    output_root_swc_path = Path(output_root_swc_path)
    shutil.copy2(run_swc_path, output_root_swc_path)
    return output_root_swc_path
