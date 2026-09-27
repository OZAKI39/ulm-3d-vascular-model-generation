from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

from utils.nne2.config import NNE2Config
from utils.nne2.segmentation import segment_vessels
from utils.nne2.stack_io import inspect_stack, list_stack_frames, load_stack


def _stack(tmp_path: Path) -> Path:
    root = tmp_path / "Stack-Test"
    root.mkdir()
    yy, xx = np.indices((48, 48))
    for index in range(1, 9):
        image = np.full((48, 48), 100, dtype=np.uint16)
        image[(xx - 24) ** 2 + (yy - 24) ** 2 <= 4**2] = 3500
        Image.fromarray(image).save(
            root / f"Stack-Test_Cycle00001_CurrentSettings_Ch2_{index:06d}.tif"
        )
    (root / "Stack-Test.xml").write_text(
        '<PVScan><Key key="micronsPerPixel_XAxis" value="1.0"/>'
        '<Key key="micronsPerPixel_YAxis" value="1.0"/></PVScan>',
        encoding="utf-8",
    )
    return root


def test_stack_loader_sorts_and_downsamples(tmp_path: Path) -> None:
    stack = _stack(tmp_path)
    assert len(list_stack_frames(stack)) == 8
    metadata = inspect_stack(
        stack, xy_spacing_um=1.0, z_spacing_um=2.0, target_xy_spacing_um=2.0
    )
    volume = load_stack(metadata, workers=2)

    assert volume.shape == (8, 24, 24)
    assert metadata.processed_spacing_xyz_um == (2.0, 2.0, 2.0)


def test_bright_tube_is_segmented(tmp_path: Path) -> None:
    stack = _stack(tmp_path)
    metadata = inspect_stack(
        stack, xy_spacing_um=1.0, z_spacing_um=2.0, target_xy_spacing_um=1.0
    )
    volume = load_stack(metadata, workers=2)
    config = NNE2Config(
        input_dir=tmp_path,
        output_root=tmp_path,
        foreground_quantile=0.90,
        min_component_voxels=8,
    )
    result = segment_vessels(volume, metadata.processed_spacing_xyz_um, config)

    assert np.any(result.mask_zyx[:, 24, 24])
    assert 0 < np.mean(result.mask_zyx) < 0.35
