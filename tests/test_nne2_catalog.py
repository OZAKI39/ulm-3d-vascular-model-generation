from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image
from scipy.io import savemat

from utils.nne2.catalog import load_nne2_catalog


def _write_minimal_nne2(root: Path) -> Path:
    refs = root / "hana_refs" / "hana_refs"
    stack = root / "hana_stk" / "hana_stk" / "Stack-A"
    maps = root / "maps" / "maps"
    refs.mkdir(parents=True)
    stack.mkdir(parents=True)
    maps.mkdir(parents=True)
    Image.fromarray(np.zeros((16, 16, 3), dtype=np.uint8)).save(
        refs / "Ref-A_ImageWindow1Ch1Ch2-SingleImageReference.tiff"
    )
    for index in range(1, 4):
        Image.fromarray(np.zeros((16, 16), dtype=np.uint16)).save(
            stack / f"Stack-A_Cycle00001_CurrentSettings_Ch2_{index:06d}.tif"
        )
    (stack / "Stack-A.xml").write_text(
        '<PVScan><Key key="micronsPerPixel_XAxis" value="1.0"/>'
        '<Key key="micronsPerPixel_YAxis" value="1.0"/></PVScan>',
        encoding="utf-8",
    )
    Image.fromarray(np.zeros((16, 16, 3), dtype=np.uint8)).save(maps / "010101.jpg")
    table = np.empty((3, 24), dtype=object)
    table[:] = ""
    headers = [f"field_{index}" for index in range(24)]
    table[0] = headers
    complete = [""] * 24
    complete[0] = "010101"
    complete[1] = 2
    complete[5] = 0
    complete[8] = 50
    complete[15] = "/old/path/Ref-A/"
    complete[16] = "/old/path/Stack-A"
    complete[18] = 2
    complete[20] = 1.0
    complete[21] = 1.0
    complete[22] = "/old/path/map.bmp"
    complete[23] = 3.0
    table[1] = complete
    incomplete = list(complete)
    incomplete[15] = ""
    table[2] = incomplete
    savemat(root / "vdb.mat", {"vdb_gnu": table})
    return root


def test_catalog_skips_missing_records_before_processing(tmp_path: Path) -> None:
    catalog = load_nne2_catalog(_write_minimal_nne2(tmp_path / "NNE2"))

    assert len(catalog.records) == 2
    assert len(catalog.complete_records) == 1
    assert catalog.complete_records[0].tree_key == "010101_tree_2"
    assert len(catalog.skipped_records) == 1
    assert "missing_reference_image" in catalog.skipped_records[0].skip_reasons


def test_tree_id_is_combined_with_subject_id(tmp_path: Path) -> None:
    catalog = load_nne2_catalog(_write_minimal_nne2(tmp_path / "NNE2"))

    assert catalog.complete_records[0].tree_key == "010101_tree_2"
