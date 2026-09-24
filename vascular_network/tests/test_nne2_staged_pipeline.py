from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from PIL import Image
from scipy.io import savemat

from utils.nne2.components import clean_components
from utils.nne2.config import NNE2Config
from utils.nne2.manifest import load_source_manifests, resolve_and_verify_artifacts
from utils.nne2.pipeline import run_nne2_pipeline


def _write_staged_nne2(root: Path) -> Path:
    refs = root / "hana_refs" / "hana_refs"
    stack = root / "hana_stk" / "hana_stk" / "Stack-A"
    maps = root / "maps" / "maps"
    refs.mkdir(parents=True)
    stack.mkdir(parents=True)
    maps.mkdir(parents=True)
    yy, xx = np.indices((40, 40))
    reference = np.zeros((40, 40, 3), dtype=np.uint8)
    reference[(xx - 20) ** 2 + (yy - 20) ** 2 <= 5**2, 1] = 230
    reference[19:22, 19:22, 0] = 255
    Image.fromarray(reference).save(
        refs / "Ref-A_ImageWindow1Ch1Ch2-SingleImageReference.tiff"
    )
    for index in range(1, 10):
        image = np.full((40, 40), 100, dtype=np.uint16)
        image[(xx - 20) ** 2 + (yy - 20) ** 2 <= 4**2] = 3500
        Image.fromarray(image).save(
            stack / f"Stack-A_Cycle00001_CurrentSettings_Ch2_{index:06d}.tif"
        )
    (stack / "Stack-A.xml").write_text(
        '<PVScan><Key key="micronsPerPixel_XAxis" value="1.0"/>'
        '<Key key="micronsPerPixel_YAxis" value="1.0"/></PVScan>',
        encoding="utf-8",
    )
    Image.fromarray(np.zeros((40, 40, 3), dtype=np.uint8)).save(maps / "010101.jpg")
    table = np.empty((2, 24), dtype=object)
    table[:] = ""
    table[0] = [f"field_{index}" for index in range(24)]
    row = [""] * 24
    row[0] = "010101"
    row[1] = 2
    row[5] = 0
    row[8] = 50
    row[15] = "/old/path/Ref-A/"
    row[16] = "/old/path/Stack-A"
    row[18] = 5
    row[20] = 1.0
    row[21] = 1.0
    row[22] = "/old/path/map.bmp"
    row[23] = 2.0
    table[1] = row
    savemat(root / "vdb.mat", {"vdb_gnu": table})
    return root


def test_component_decisions_account_for_kept_and_removed_voxels() -> None:
    candidate = np.zeros((8, 12, 12), dtype=bool)
    candidate[1:7, 5:7, 5:7] = True
    candidate[0, 0, 0] = True
    cleaned, removed, decisions = clean_components(
        candidate, (1.0, 1.0, 2.0), min_component_voxels=4
    )

    assert np.count_nonzero(cleaned) == 24
    assert np.count_nonzero(removed) == 1
    assert sum(item["voxel_count"] for item in decisions) == np.count_nonzero(candidate)
    assert {item["decision"] for item in decisions} == {"keep", "remove"}


def test_preprocess_manifest_can_drive_independent_step3(tmp_path: Path) -> None:
    input_dir = _write_staged_nne2(tmp_path / "NNE2")
    output_dir = tmp_path / "outputs"
    preprocess = run_nne2_pipeline(
        NNE2Config(
            input_dir=input_dir,
            output_root=output_dir,
            stage="preprocess",
            target_xy_spacing_um=1.0,
            foreground_quantile=0.90,
            min_component_voxels=8,
            visualizations_enabled=False,
        )
    )

    assert preprocess.status == "completed"
    manifests = load_source_manifests(preprocess.run_root)
    assert set(manifests) == {"Stack-A"}
    artifacts = resolve_and_verify_artifacts(preprocess.run_root, manifests["Stack-A"])
    assert artifacts["preprocess_arrays"].is_file()

    graph = run_nne2_pipeline(
        NNE2Config(
            input_dir=input_dir,
            output_root=output_dir,
            stage="hierarchical-graph",
            source_run=preprocess.run_root,
            target_xy_spacing_um=1.0,
            min_registration_score=-1.0,
            max_anchor_distance_um=100.0,
            visualizations_enabled=False,
            save_vtp=False,
        )
    )

    assert graph.status == "completed"
    assert graph.summary["successful_hierarchy_count"] == 1
    tree_json = next(graph.run_root.glob("trees/*/graphs/directed_hierarchy.json"))
    payload = json.loads(tree_json.read_text(encoding="utf-8"))
    assert payload["representation"]["directed"] is True
    assert payload["summary"]["root_branch_id"] >= 0


def test_manifest_detects_changed_source_artifact(tmp_path: Path) -> None:
    input_dir = _write_staged_nne2(tmp_path / "NNE2")
    preprocess = run_nne2_pipeline(
        NNE2Config(
            input_dir=input_dir,
            output_root=tmp_path / "outputs",
            stage="preprocess",
            target_xy_spacing_um=1.0,
            foreground_quantile=0.90,
            min_component_voxels=8,
            visualizations_enabled=False,
            write_nifti=False,
        )
    )
    manifest = load_source_manifests(preprocess.run_root)["Stack-A"]
    path = preprocess.run_root / manifest["artifacts"]["step1_report"]["relative_path"]
    path.write_text(path.read_text(encoding="utf-8") + " ", encoding="utf-8")

    try:
        resolve_and_verify_artifacts(preprocess.run_root, manifest)
    except (FileNotFoundError, ValueError) as exc:
        assert "mismatch" in str(exc)
    else:
        raise AssertionError("Expected changed Step 1 artifact to be rejected")


def test_all_stage_writes_preprocess_and_directed_outputs(tmp_path: Path) -> None:
    input_dir = _write_staged_nne2(tmp_path / "NNE2")
    run = run_nne2_pipeline(
        NNE2Config(
            input_dir=input_dir,
            output_root=tmp_path / "outputs",
            stage="all",
            target_xy_spacing_um=1.0,
            foreground_quantile=0.90,
            min_component_voxels=8,
            min_registration_score=-1.0,
            max_anchor_distance_um=100.0,
            visualizations_enabled=False,
            save_vtp=False,
        )
    )

    assert run.status == "completed"
    assert len(list(run.run_root.glob("stacks/*/reports/step1_step2_manifest.json"))) == 1
    assert len(list(run.run_root.glob("trees/*/graphs/directed_hierarchy.json"))) == 1
