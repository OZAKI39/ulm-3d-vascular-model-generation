"""Physical-unit adaptation for the existing mouse single-network workflow.

Only the selected analysis component is modeled. Other reference components,
image data and masks are retained; no alternative radius filter is used.
"""
from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import numpy as np

from utils.rodent_vasculature.catalog import build_catalog, select_records
from utils.rodent_vasculature.swc_analysis import select_analysis_swc
from utils.rodent_vasculature.swc_io import load_swc
from utils.swc_roi_yaml_config import load_swc_roi_yaml_config

from .pipeline import Options, process_file
from .swc_export import read_source


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def backup_file(source: Path, target: Path):
    """Create a byte-for-byte independent backup before fitting, never overwrite."""
    target.parent.mkdir(parents=True, exist_ok=True)
    with source.open("rb") as original, target.open("xb") as backup:
        shutil.copyfileobj(original, backup)
    shutil.copystat(source, target)
    digest = sha256(source)
    if sha256(target) != digest:
        raise ValueError(f"Backup differs from original: {source}")
    return {"source": str(source), "path": str(target), "sha256": digest, "byte_identical": True}


def write_rows(path: Path, rows: np.ndarray, header: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        np.savetxt(handle, rows, fmt=["%d", "%d", "%.17g", "%.17g", "%.17g", "%.17g", "%d"], header=header)


def replace_analysis_component(reference, selected_ids, smooth_rows, spacing_xyz_um):
    """Return a forest in legacy voxel-XYZ / micrometre-radius representation.

    Keep the untouched components' IDs and parents. Allocate fresh IDs above
    the original maximum to the complete modeled component, avoiding collisions.
    """
    rows = np.column_stack((reference.node_ids, reference.type_codes, reference.points_voxel_xyz,
                            reference.radius_raw_um, reference.parent_ids))
    keep = ~np.isin(reference.node_ids, selected_ids)
    modeled = smooth_rows.copy()
    modeled[:, 2:5] /= np.asarray(spacing_xyz_um)
    offset = int(reference.node_ids.max()) + 1
    mapping = {int(row[0]): offset + i for i, row in enumerate(modeled)}
    modeled[:, 0] = [mapping[int(n)] for n in smooth_rows[:, 0]]
    modeled[:, 6] = [-1 if p == -1 else mapping[int(p)] for p in smooth_rows[:, 6]]
    # Retain the component's original position in the forest's traversal order.
    insert = int(np.flatnonzero(~keep)[0])
    before = rows[:insert][keep[:insert]]
    after = rows[insert:][keep[insert:]]
    return np.vstack((before, modeled, after)), mapping


def prepare_mouse_dataset(config_path: Path, output_dir: Path, *, project_root: Path,
                          accept_native_merges: bool = True, auto_resample: bool = True,
                          model_mode: str = "network") -> dict:
    options = Options(surface=False, qc_plot=True, input_units="um", model_mode=model_mode,
                      accept_native_merges=accept_native_merges if model_mode == "network" else False,
                      auto_resample=auto_resample)
    options.validate()
    settings = load_swc_roi_yaml_config(config_path, project_root=project_root)
    config = settings.rodent
    if config.swc_units != "voxel" or config.cohort not in {"raw-analysis", "raw-total", "train"}:
        raise ValueError("Mouse preparation requires a voxel-coordinate mouse dataset and one explicit cohort")
    output_dir = output_dir.resolve()
    if output_dir.is_relative_to(config.input_dir.resolve()):
        raise ValueError("Derived dataset must be outside the original dataset directory")
    records = select_records(build_catalog(config.input_dir, config.cohort), sample_id=config.sample_id,
                             parent_group_id=config.parent_group_id, split=config.split, max_samples=config.max_samples)
    if not records:
        raise ValueError("No selected mouse SWC files")
    output_dir.mkdir(parents=True, exist_ok=False)
    report = {"status": "failed", "source_config": str(settings.source_path),
              "model_mode": model_mode,
              "source_config_sha256": sha256(settings.source_path), "output_dataset": str(output_dir),
              "scope": "Only each source configuration's selected analysis component is modeled; other components remain unchanged reference data",
              "input_coordinates": "voxel XYZ; radius already in micrometres",
              "model_coordinates": "XYZ and radius in micrometres; native density 0.4..0.6 points/um; native max_distance=6 is a dimensionless factor multiplied by mean radius internally",
              "output_coordinates": "XYZ divided by original voxel spacing; radius stays in micrometres",
              "spacing_xyz_um": list(config.spacing_xyz_um), "samples": []}
    report["config_backup"] = backup_file(settings.source_path, output_dir / "backups/configs" / settings.source_path.name)
    for record in records:
        item = {"source_swc": str(record.swc_path), "sample_id": record.sample_id, "status": "failed",
                "source_sha256": sha256(record.swc_path)}
        report["samples"].append(item)
        try:
            item["original_backup"] = backup_file(record.swc_path, output_dir / "backups" / record.swc_path.relative_to(config.input_dir))
            reference = load_swc(record.swc_path, spacing_xyz_um=config.spacing_xyz_um,
                                 volume_shape_zyx=None, swc_units="voxel")
            if not reference.structurally_valid:
                raise ValueError(f"Invalid source SWC forest: {reference.validation}")
            selection = select_analysis_swc(reference, spacing_xyz_um=config.spacing_xyz_um,
                                            volume_shape_zyx=None, analysis_component_id=config.analysis_component_id)
            analysis = selection.analysis_swc
            item["original_selection"] = selection.summary
            item["component_records"] = selection.component_records
            sample_output = output_dir / "vascularmd" / record.source_stem
            physical = sample_output / "input" / (record.source_stem + "_analysis_um.swc")
            rows = np.column_stack((analysis.node_ids, analysis.type_codes, analysis.points_um,
                                    analysis.radius_raw_um, analysis.parent_ids))
            write_rows(physical, rows, f"Original mouse analysis component in physical micrometres\nSource: {record.swc_path}\nSource SHA256: {item['source_sha256']}\nXYZ multiplied by {list(config.spacing_xyz_um)}; radius unchanged; original IDs/TYPE/PARENT preserved")
            native = process_file(physical, sample_output / "native", options)
            item["native_qc"] = str(sample_output / "native" / (physical.stem + "_vmd_qc.json"))
            if native["status"] != "success":
                item["failed_components"] = native["failed_components"]
                raise ValueError("Native modeling/export did not succeed; no replacement SWC published")
            smooth = read_source(Path(native["outputs"]["swc"]))
            combined, mapping = replace_analysis_component(reference, analysis.node_ids, smooth.rows, config.spacing_xyz_um)
            target = output_dir / record.swc_path.relative_to(config.input_dir)
            # Validate the combined forest before placing it in the dataset's SWC directory.
            staged = sample_output / "staged_voxel_forest.swc"
            write_rows(staged, combined, f"VascularMD derived mouse SWC; only analysis component optimized\nModel mode: {model_mode}\nOriginal source: {record.swc_path}\nOriginal SHA256: {item['source_sha256']}\nXYZ in original voxel coordinates; radius in um\nOther {reference.component_count-1} components are unchanged REFERENCE_ONLY data\nNative QC: {item['native_qc']}")
            result = load_swc(staged, spacing_xyz_um=config.spacing_xyz_um, volume_shape_zyx=None)
            if not result.structurally_valid or result.component_count != reference.component_count:
                raise ValueError("Combined output changed reference component count or produced an invalid forest")
            selected = select_analysis_swc(result, spacing_xyz_um=config.spacing_xyz_um, volume_shape_zyx=None,
                                           analysis_component_id=config.analysis_component_id)
            if set(selected.analysis_swc.node_ids) != set(mapping.values()):
                raise ValueError("The existing mouse selection rule no longer selects the modeled component")
            np.testing.assert_allclose(selected.analysis_swc.points_um, smooth.rows[:, 2:5], rtol=1e-14, atol=1e-12)
            np.testing.assert_array_equal(selected.analysis_swc.radius_raw_um, smooth.rows[:, 5])
            reference_ids = set(selection.reference_only_node_ids)
            old_reference = rows_from_ids(reference, reference_ids)
            new_reference = rows_from_ids(result, reference_ids)
            np.testing.assert_array_equal(old_reference, new_reference)
            if sha256(record.swc_path) != item["source_sha256"]:
                raise ValueError("Original source changed during processing")
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists():
                raise FileExistsError(target)
            staged.rename(target)
            # Keep image/Mask alignment and data; do not resample or regenerate volumes.
            auxiliary = {}
            for name, source_path in (("image", record.image_path), ("mask", record.mask_path)):
                if source_path is not None:
                    link = output_dir / source_path.relative_to(config.input_dir)
                    link.parent.mkdir(parents=True, exist_ok=True)
                    link.symlink_to(source_path.resolve())
                    auxiliary[name] = {"source": str(source_path), "link": str(link), "sha256": sha256(source_path)}
            if config.cohort == "train" and record.split:
                split_file = output_dir / "train_data" / (record.split + ".txt")
                with split_file.open("a", encoding="utf-8") as handle:
                    handle.write(record.source_stem + "\n")
            item.update({"status": "success", "output_swc": str(target), "output_sha256": sha256(target),
                         "original_preserved": True, "reference_only_nodes_preserved": len(reference_ids),
                         "reference_component_count_preserved": result.component_count,
                         "optimized_node_count": len(smooth.rows), "optimized_node_id_mapping": mapping,
                         "topology": native["topology"], "comparison": native["comparison"], "auxiliary": auxiliary,
                         "output_analysis_selection": selected.summary})
        except Exception as exc:
            item["error"] = f"{type(exc).__name__}: {exc}"
    report["status"] = "success" if all(s["status"] == "success" for s in report["samples"]) else "failed"
    (output_dir / "vascularmd_mouse_manifest.json").write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    return report


def rows_from_ids(swc, ids):
    mask = np.isin(swc.node_ids, list(ids))
    return np.column_stack((swc.node_ids[mask], swc.type_codes[mask], swc.points_voxel_xyz[mask],
                            swc.radius_raw_um[mask], swc.parent_ids[mask]))
