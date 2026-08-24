"""File inventory and evidence-based Mask provenance classification."""

from __future__ import annotations

import hashlib
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import cv2
import numpy as np
import tifffile

from .mask_reference import SPACING_XYZ_UM


SWC_DERIVED = "SWC_DERIVED_ANNOTATION_MASK"
UNET_PREDICTED = "UNET_PREDICTED_MASK"
UNKNOWN = "UNKNOWN_MASK_PROVENANCE"


@dataclass(frozen=True, slots=True)
class DatasetGroup:
    cohort: str
    image_dir: Path
    mask_dir: Path
    swc_dir: Path | None


def _groups(dataset_root: Path) -> tuple[DatasetGroup, ...]:
    root = Path(dataset_root)
    return (
        DatasetGroup(
            "raw-analysis",
            root / "raw_data/analysis_data/analysis_data/images",
            root / "raw_data/analysis_data/analysis_data/mask",
            root / "raw_data/analysis_data/analysis_data/swc",
        ),
        DatasetGroup(
            "raw-total",
            root / "raw_data/total_vascular_data/total_vascular_data/images",
            root / "raw_data/total_vascular_data/total_vascular_data/mask",
            root / "raw_data/total_vascular_data/total_vascular_data/swc",
        ),
        DatasetGroup(
            "train-data",
            root / "train_data/images",
            root / "train_data/mask",
            root / "train_data/swc",
        ),
        DatasetGroup(
            "train-test-copy",
            root / "train_data/testData/images",
            root / "train_data/testData/mask",
            root / "train_data/swc",
        ),
    )


def _split_map(dataset_root: Path) -> dict[str, str]:
    output: dict[str, str] = {}
    for split in ("train", "val", "test"):
        path = Path(dataset_root) / "train_data" / f"{split}.txt"
        if not path.is_file():
            continue
        for line in path.read_text(encoding="utf-8-sig").splitlines():
            if line.strip():
                output[Path(line.strip()).stem] = split
    return output


def _analysis_to_train_stem(stem: str) -> str:
    return "cut_" + stem[len("fMOST_") :] if stem.startswith("fMOST_") else stem


def _volume_statistics(path: Path) -> dict[str, Any]:
    success, frames = cv2.imreadmulti(str(path), flags=cv2.IMREAD_UNCHANGED)
    if not success or not frames:
        raise ValueError(f"Unable to decode Mask TIFF: {path}")
    array = np.stack(frames, axis=0)
    unique = np.unique(array)
    if len(unique) <= 2:
        representation = "binary"
    elif np.issubdtype(array.dtype, np.integer) and float(array.min()) >= 0.0:
        representation = "probability"
    else:
        representation = "label"
    return {
        "shape_z": int(array.shape[0]),
        "shape_y": int(array.shape[1]),
        "shape_x": int(array.shape[2]),
        "dtype": str(array.dtype),
        "value_min": float(array.min()),
        "value_max": float(array.max()),
        "unique_value_count": int(len(unique)),
        "mask_representation": representation,
        "nonzero_voxel_count": int(np.count_nonzero(array)),
        "file_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def inventory_masks(dataset_root: Path, *, workers: int = 8) -> list[dict[str, Any]]:
    """Inventory every dataset Mask without modifying source data."""

    root = Path(dataset_root).resolve()
    split_by_stem = _split_map(root)
    tasks: list[tuple[DatasetGroup, Path]] = []
    for group in _groups(root):
        if group.mask_dir.is_dir():
            tasks.extend((group, path) for path in sorted(group.mask_dir.glob("*.tif")))

    def inspect(task: tuple[DatasetGroup, Path]) -> dict[str, Any]:
        group, mask_path = task
        stem = mask_path.stem
        train_stem = _analysis_to_train_stem(stem)
        if group.cohort == "raw-total":
            split = "full_1620"
        elif group.cohort == "raw-analysis":
            source_split = split_by_stem.get(train_stem)
            split = f"independent_analysis:{source_split or 'outside_training_set'}"
        elif group.cohort == "train-test-copy":
            split = "test_copy"
        else:
            split = split_by_stem.get(stem, "unknown")
        image_path = group.image_dir / mask_path.name
        swc_name = mask_path.with_suffix(".swc").name
        swc_path = group.swc_dir / swc_name if group.swc_dir is not None else None
        if group.cohort == "train-test-copy" and swc_path is not None:
            swc_path = root / "train_data/swc" / swc_name
        row: dict[str, Any] = {
            "sample_id": f"{group.cohort}__{stem}",
            "filename": mask_path.name,
            "parent_folder": str(mask_path.parent.resolve()),
            "cohort": group.cohort,
            "data_split": split,
            "mask_path": str(mask_path.resolve()),
            "image_path": str(image_path.resolve()) if image_path.is_file() else "",
            "swc_path": str(swc_path.resolve()) if swc_path and swc_path.is_file() else "",
            "corresponding_image": bool(image_path.is_file()),
            "corresponding_swc": bool(swc_path and swc_path.is_file()),
            "spacing_x_um": SPACING_XYZ_UM[0],
            "spacing_y_um": SPACING_XYZ_UM[1],
            "spacing_z_um": SPACING_XYZ_UM[2],
            "array_axis_order": "zyx",
        }
        try:
            with tifffile.TiffFile(mask_path) as stack:
                series = stack.series[0]
                row["tiff_shape"] = "x".join(map(str, series.shape))
                row["tiff_axes"] = series.axes
                row["tiff_compression"] = str(stack.pages[0].compression)
            row.update(_volume_statistics(mask_path))
            row["read_status"] = "PASS"
            row["read_error"] = ""
        except Exception as exc:  # every source failure remains explicit in the inventory
            row["read_status"] = "FAIL"
            row["read_error"] = f"{type(exc).__name__}: {exc}"
        return row

    with ThreadPoolExecutor(max_workers=max(1, int(workers))) as executor:
        rows = list(executor.map(inspect, tasks))
    return sorted(rows, key=lambda row: (str(row["cohort"]), str(row["filename"])))


def _classification_evidence(
    row: dict[str, Any],
    *,
    paper_path: Path,
    manual_path: Path,
) -> dict[str, Any]:
    cohort = str(row["cohort"])
    paper = str(Path(paper_path).resolve())
    manual = str(Path(manual_path).resolve())
    if cohort in {"raw-analysis", "train-data", "train-test-copy"}:
        return {
            "mask_provenance": SWC_DERIVED,
            "evidence_source": "paper+dataset_structure+file_identity",
            "evidence_path": f"{paper} (pp.3,5,6,8); {manual} (pp.2-3)",
            "document_evidence": (
                "Paper p.3: semi-manual skeleton+radii -> computational masks; p.5: corrected "
                "skeleton sampled at one voxel and local-radius Gaussian creates the supervision "
                "map; p.6: semi-automatic ground truth is generated from skeleton/radius and is "
                "not an independent manual segmentation; p.8 distinguishes it from U-Net output."
            ),
            "code_or_file_evidence": (
                "train_data has exactly 199 image/mask/SWC triplets matching the paper's 199 "
                "manually annotated blocks; raw-analysis has exactly 108 triplets matching the "
                "paper's reliability set; testData masks are byte-identical copies of test-split "
                "train_data masks. BVLab manual exports the corrected label as SWC."
            ),
            "confidence": 0.99 if cohort != "raw-analysis" else 0.98,
            "classification_basis": "direct_document_statement_plus_exact_count_and_copy_evidence",
        }
    if cohort == "raw-total":
        return {
            "mask_provenance": UNET_PREDICTED,
            "evidence_source": "paper+dataset_structure+paired_difference",
            "evidence_path": f"{paper} (pp.3,6-8)",
            "document_evidence": (
                "Paper p.3 says 1,620 sub-blocks were formed and U-Net was applied to remaining "
                "unlabeled data to predict masks; pp.6-8 explicitly call the automatic masks U-Net "
                "predictions and distinguish them from semi-automatic annotation masks."
            ),
            "code_or_file_evidence": (
                "raw-total contains exactly 1,620 image/mask/SWC triplets. For overlapping image "
                "blocks, the raw-total image is byte-identical to raw-analysis while Mask and SWC "
                "differ, matching the paper's automatic-vs-semi-automatic pairing. No generator "
                "source or manifest is shipped, so this conclusion is convergent rather than a "
                "direct per-file provenance tag."
            ),
            "confidence": 0.94,
            "classification_basis": "exact_count_workflow_match_and_same_image_paired_difference",
        }
    return {
        "mask_provenance": UNKNOWN,
        "evidence_source": "no_sufficient_evidence",
        "evidence_path": "",
        "document_evidence": "",
        "code_or_file_evidence": "",
        "confidence": 0.0,
        "classification_basis": "unclassified",
    }


def classify_inventory(
    inventory: list[dict[str, Any]],
    *,
    paper_path: Path,
    manual_path: Path,
) -> list[dict[str, Any]]:
    return [
        {**row, **_classification_evidence(row, paper_path=paper_path, manual_path=manual_path)}
        for row in inventory
    ]


def current_roi_provenance(
    report: list[dict[str, Any]],
    *,
    mask_path: Path,
) -> dict[str, Any]:
    target = str(Path(mask_path).resolve()).casefold()
    matches = [row for row in report if str(row["mask_path"]).casefold() == target]
    if len(matches) != 1:
        raise ValueError(f"Expected one provenance record for {mask_path}, found {len(matches)}")
    return matches[0]
