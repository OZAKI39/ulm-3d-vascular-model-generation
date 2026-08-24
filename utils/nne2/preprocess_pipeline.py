"""Reusable Step 1-2 artifact writing and loading for NNE2 stacks."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from ..io import write_csv, write_json
from .export import write_stack_nifti
from .manifest import resolve_and_verify_artifacts, write_stack_manifest
from .segmentation import SegmentationResult


@dataclass(frozen=True, slots=True)
class LoadedStep2Stack:
    stack_name: str
    normalized_zyx: np.ndarray
    candidate_mask_zyx: np.ndarray
    mask_zyx: np.ndarray
    removed_mask_zyx: np.ndarray
    skeleton_zyx: np.ndarray
    radius_zyx_um: np.ndarray
    spacing_xyz_um: tuple[float, float, float]
    shape_zyx: tuple[int, int, int]
    source_artifacts: dict[str, Path]


def write_step1_step2_artifacts(
    *,
    run_root: Path,
    output_dir: Path,
    stack_name: str,
    segmentation: SegmentationResult,
    skeleton_zyx: np.ndarray,
    radius_zyx_um: np.ndarray,
    spacing_xyz_um: tuple[float, float, float],
    stack_metadata_report: dict[str, Any],
    skeleton_report: dict[str, Any],
    cache_key: str,
    write_nifti: bool,
) -> tuple[list[Path], Path]:
    """Write the stable Step 1-2 contract consumed by later Step 3 runs."""
    volumes = output_dir / "volumes"
    reports = output_dir / "reports"
    tables = output_dir / "tables"
    volumes.mkdir(parents=True, exist_ok=True)
    reports.mkdir(parents=True, exist_ok=True)
    tables.mkdir(parents=True, exist_ok=True)

    arrays_path = volumes / "preprocess_arrays.npz"
    np.savez_compressed(
        arrays_path,
        normalized_zyx=np.asarray(segmentation.normalized_zyx, dtype=np.float16),
        candidate_mask_zyx=np.asarray(segmentation.candidate_mask_zyx, dtype=np.uint8),
        cleaned_mask_zyx=np.asarray(segmentation.mask_zyx, dtype=np.uint8),
        removed_mask_zyx=np.asarray(segmentation.removed_mask_zyx, dtype=np.uint8),
        skeleton_zyx=np.asarray(skeleton_zyx, dtype=np.uint8),
        radius_zyx_um=np.asarray(radius_zyx_um, dtype=np.float16),
    )
    metadata_path = reports / "stack_metadata.json"
    step1_path = reports / "step1_segmentation_report.json"
    step2_path = reports / "step2_centerline_report.json"
    decisions_path = tables / "component_decisions.csv"
    write_json(stack_metadata_report, metadata_path)
    write_json(segmentation.report(), step1_path)
    write_json(skeleton_report, step2_path)
    write_csv(segmentation.component_decisions, decisions_path)

    artifacts: dict[str, Path] = {
        "preprocess_arrays": arrays_path,
        "stack_metadata": metadata_path,
        "step1_report": step1_path,
        "step2_report": step2_path,
        "component_decisions": decisions_path,
    }
    if write_nifti:
        nifti_specs = (
            ("candidate_mask", segmentation.candidate_mask_zyx, np.uint8,
             "NNE2 Step 1 candidate vessel mask before component cleanup"),
            ("cleaned_mask", segmentation.mask_zyx, np.uint8,
             "NNE2 Step 1 cleaned vessel mask; stack XYZ coordinates"),
            ("removed_islands", segmentation.removed_mask_zyx, np.uint8,
             "NNE2 Step 1 removed small disconnected components"),
            ("coarse_centerline", skeleton_zyx, np.uint8,
             "NNE2 Step 2 coarse centerline; stack XYZ coordinates"),
            ("coarse_radius", radius_zyx_um, np.float32,
             "NNE2 Step 2 coarse radius estimate in micrometers"),
        )
        for name, array, dtype, description in nifti_specs:
            artifacts[name] = write_stack_nifti(
                array,
                volumes / f"{name}.nii.gz",
                spacing_xyz_um,
                description,
                dtype=dtype,
            )

    manifest_path = reports / "step1_step2_manifest.json"
    write_stack_manifest(
        manifest_path,
        run_root=run_root,
        stack_name=stack_name,
        shape_zyx=tuple(int(value) for value in segmentation.mask_zyx.shape),
        spacing_xyz_um=spacing_xyz_um,
        artifacts=artifacts,
        cache_key=cache_key,
    )
    return [*artifacts.values(), manifest_path], manifest_path


def load_step2_stack(
    source_run: Path, manifest: dict[str, Any]
) -> LoadedStep2Stack:
    artifacts = resolve_and_verify_artifacts(source_run, manifest)
    with np.load(artifacts["preprocess_arrays"], allow_pickle=False) as arrays:
        normalized = np.asarray(arrays["normalized_zyx"], dtype=np.float32)
        candidate = np.asarray(arrays["candidate_mask_zyx"], dtype=bool)
        cleaned = np.asarray(arrays["cleaned_mask_zyx"], dtype=bool)
        removed = np.asarray(arrays["removed_mask_zyx"], dtype=bool)
        skeleton = np.asarray(arrays["skeleton_zyx"], dtype=bool)
        radius = np.asarray(arrays["radius_zyx_um"], dtype=np.float32)
    expected_shape = tuple(int(value) for value in manifest["shape_zyx"])
    arrays_to_check = (normalized, candidate, cleaned, removed, skeleton, radius)
    if any(array.shape != expected_shape for array in arrays_to_check):
        raise ValueError(f"Step 1-2 arrays do not match manifest shape for {manifest['stack_name']}")
    if np.any(skeleton & ~cleaned):
        raise ValueError(f"Step 2 skeleton leaves cleaned mask for {manifest['stack_name']}")
    if np.any(cleaned & removed):
        raise ValueError(f"Kept and removed masks overlap for {manifest['stack_name']}")
    return LoadedStep2Stack(
        stack_name=str(manifest["stack_name"]),
        normalized_zyx=normalized,
        candidate_mask_zyx=candidate,
        mask_zyx=cleaned,
        removed_mask_zyx=removed,
        skeleton_zyx=skeleton,
        radius_zyx_um=radius,
        spacing_xyz_um=tuple(float(value) for value in manifest["spacing_xyz_um"]),
        shape_zyx=expected_shape,
        source_artifacts=artifacts,
    )
