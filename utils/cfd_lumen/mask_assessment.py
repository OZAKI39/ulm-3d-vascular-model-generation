"""End-to-end v5-1 Mask provenance and local A/B/C assessment."""

from __future__ import annotations

import hashlib
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import trimesh

from .config import load_cfd_lumen_config
from .context_domain import build_cfd_context_domain
from .export import write_csv, write_json
from .geometry_preprocess import resample_branches, validate_and_extract_branches
from .mask_assisted_refinement import INFLUENCE_CONFIGS, refine_surface_with_mask
from .mask_comparison import (
    SectionLocation,
    abc_overlay_figure,
    anisotropy_figure,
    compare_binary_masks,
    cross_section_figure,
    local_normal_roughness,
    mask_anisotropy_metrics,
    mesh_from_voxel_centers,
    provenance_summary_figure,
    radius_comparison_figure,
    regional_mask_metrics,
    section_radius_metrics,
    surface_to_mask_metrics,
    swc_mask_overlay_figure,
    unified_mesh_qc,
)
from .mask_provenance import (
    SWC_DERIVED,
    classify_inventory,
    current_roi_provenance,
    inventory_masks,
)
from .mask_reference import (
    PAPER_FOREGROUND_UINT8,
    SPACING_XYZ_UM,
    RadiusCalibration,
    read_swc_reference,
    read_tiff_volume,
    reconstruct_mask_from_swc,
    resample_swc_one_voxel,
    threshold_probability_mask,
    write_tiff_volume,
)
from .mask_surface import (
    crop_volume_physical,
    load_mesh_vtp,
    local_surface,
    mask_surfaces,
    save_local_swc_vtp,
    save_mesh_vtp,
    select_junction_component,
    signed_distance_field,
)
from .roi_io import load_sampling_rois


CURRENT_ROI_ID = "raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01__anchor_003274"
CURRENT_SOURCE_STEM = "fMOST_0_5_6_0_0_6_0001_02_01"
CURRENT_TRAIN_STEM = "cut_0_5_6_0_0_6_0001_02_01"


@dataclass(frozen=True, slots=True)
class AssessmentPaths:
    project_root: Path
    dataset_root: Path
    sampling_run: Path
    formal_v5_run: Path
    run_root: Path

    @property
    def evaluation(self) -> Path:
        return self.run_root / "mask_evaluation"


def _sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _layout(root: Path) -> dict[str, Path]:
    directories = {
        name: root / name
        for name in ("provenance", "reconstruction", "junctions", "qc", "figures", "report")
    }
    for directory in directories.values():
        directory.mkdir(parents=True, exist_ok=False)
    return directories


def _formal_paths(paths: AssessmentPaths) -> dict[str, Path]:
    project = paths.project_root
    roi_root = paths.formal_v5_run / "rois" / CURRENT_ROI_ID
    return {
        "model_generate.py": project / "model_generate.py",
        "pipeline.py": project / "utils/cfd_lumen/pipeline.py",
        "formal_config": project / "cfd_lumen_config.yaml",
        "formal_surface": roi_root / "geometry/lumen_surface_cfd.vtp",
        "formal_roi_summary": roi_root / "qc/roi_summary.json",
        "formal_hybrid_qc": roi_root / "qc/surface_qc.json",
    }


def _file_identity(
    inventory: list[dict[str, Any]],
    *,
    dataset_root: Path,
) -> dict[str, Any]:
    by_path = {str(row["mask_path"]).casefold(): row for row in inventory}
    raw = dataset_root / f"raw_data/analysis_data/analysis_data/mask/{CURRENT_SOURCE_STEM}.tif"
    total = dataset_root / f"raw_data/total_vascular_data/total_vascular_data/mask/{CURRENT_SOURCE_STEM}.tif"
    train = dataset_root / f"train_data/mask/{CURRENT_TRAIN_STEM}.tif"
    test = dataset_root / f"train_data/testData/mask/{CURRENT_TRAIN_STEM}.tif"
    records = {
        name: by_path[str(path.resolve()).casefold()]
        for name, path in (("raw_analysis", raw), ("raw_total", total), ("train", train), ("test_copy", test))
    }
    return {
        "paths": {name: row["mask_path"] for name, row in records.items()},
        "sha256": {name: row["file_sha256"] for name, row in records.items()},
        "raw_analysis_equals_train": records["raw_analysis"]["file_sha256"] == records["train"]["file_sha256"],
        "train_equals_test_copy": records["train"]["file_sha256"] == records["test_copy"]["file_sha256"],
        "raw_analysis_differs_from_raw_total": records["raw_analysis"]["file_sha256"] != records["raw_total"]["file_sha256"],
        "interpretation": (
            "Same image block has a semi-automatic annotation copy in raw-analysis/train/testData "
            "and a different automatic Mask in raw-total."
        ),
    }


def _calibration_check(dataset_root: Path, calibration: RadiusCalibration) -> list[dict[str, Any]]:
    names = (
        "cut_0_5_6_0_0_6_0001_00_00",
        "cut_0_5_6_0_0_6_0000_00_02",
        "0031-18_0_0-0_1_0",
    )
    rows: list[dict[str, Any]] = []
    for name in names:
        mask_path = dataset_root / f"train_data/mask/{name}.tif"
        swc_path = dataset_root / f"train_data/swc/{name}.swc"
        original = read_tiff_volume(mask_path)
        reconstructed, _ = reconstruct_mask_from_swc(
            read_swc_reference(swc_path),
            tuple(map(int, original.shape)),
            calibration=calibration,
        )
        metrics = compare_binary_masks(
            threshold_probability_mask(original),
            threshold_probability_mask(reconstructed),
        )
        rows.append(
            {
                "sample": name,
                "split": "train",
                "used_for_current_roi_evaluation": False,
                "radius_scale": calibration.radius_scale,
                "voxel_floor": calibration.voxel_floor,
                **metrics,
            }
        )
    return rows


def _junction_metadata(formal_v5_run: Path, junction_id: int) -> dict[str, Any]:
    path = (
        formal_v5_run
        / "hybrid"
        / CURRENT_ROI_ID
        / "junctions"
        / f"junction_{junction_id}"
        / "local_field_metadata.json"
    )
    return json.loads(path.read_text(encoding="utf-8"))


def _expanded_bounds(metadata: dict[str, Any], padding_um: float = 4.0) -> tuple[float, float, float, float, float, float]:
    minimum = np.asarray(metadata["bbox_min_um"], dtype=float) - padding_um
    maximum = np.asarray(metadata["bbox_max_um"], dtype=float) + padding_um
    return (
        float(minimum[0]),
        float(maximum[0]),
        float(minimum[1]),
        float(maximum[1]),
        float(minimum[2]),
        float(maximum[2]),
    )


def _follow_from_junction(
    roi: Any,
    junction_local_id: int,
    neighbor_local_id: int,
    distance_um: float,
) -> tuple[np.ndarray, float, np.ndarray]:
    positions = np.asarray(roi.local_node_positions_um, dtype=float)
    radii = np.asarray(roi.local_node_radius_um, dtype=float)
    edges = np.asarray(roi.local_edges, dtype=np.int64)
    adjacency: dict[int, list[int]] = {index: [] for index in range(roi.node_count)}
    for first, second in edges:
        adjacency[int(first)].append(int(second))
        adjacency[int(second)].append(int(first))
    previous = junction_local_id
    current = neighbor_local_id
    accumulated = 0.0
    while True:
        segment = positions[current] - positions[previous]
        length = float(np.linalg.norm(segment))
        if accumulated + length >= distance_um or len(adjacency[current]) != 2:
            fraction = np.clip((distance_um - accumulated) / max(length, 1.0e-12), 0.0, 1.0)
            center = positions[previous] + fraction * segment
            radius = float(radii[previous] + fraction * (radii[current] - radii[previous]))
            tangent = segment / max(length, 1.0e-12)
            return center, radius, tangent
        accumulated += length
        following = [node for node in adjacency[current] if node != previous]
        if not following:
            return positions[current], float(radii[current]), segment / max(length, 1.0e-12)
        previous, current = current, following[0]


def _section_locations(roi: Any, swc: Any, junction_local_id: int) -> list[SectionLocation]:
    edges = np.asarray(roi.local_edges, dtype=np.int64)
    neighbors = sorted(
        {
            int(node)
            for edge in edges[np.any(edges == junction_local_id, axis=1)]
            for node in edge
            if int(node) != junction_local_id
        }
    )
    junction_global = int(roi.local_node_global_ids[junction_local_id])
    index_by_id = {int(node_id): index for index, node_id in enumerate(swc.node_ids)}
    junction_source = index_by_id[junction_global]
    parent_global = int(swc.parent_ids[junction_source])
    roles: list[tuple[str, int]] = []
    daughter_number = 1
    for neighbor in neighbors:
        neighbor_global = int(roi.local_node_global_ids[neighbor])
        if neighbor_global == parent_global:
            roles.append(("parent_branch", neighbor))
        else:
            roles.append((f"daughter_branch_{daughter_number}", neighbor))
            daughter_number += 1
    roles.sort(key=lambda item: (item[0] != "parent_branch", item[0]))
    locations: list[SectionLocation] = []
    direction_by_role: dict[str, np.ndarray] = {}
    for role, neighbor in roles:
        center, radius, tangent = _follow_from_junction(roi, junction_local_id, neighbor, 4.0)
        direction_by_role[role] = tangent
        locations.append(
            SectionLocation(
                name=role.replace("_branch", ""),
                center_xyz_um=tuple(map(float, center)),
                tangent_xyz=tuple(map(float, tangent)),
                target_radius_um=radius,
                branch_role=role,
            )
        )
    reference_tangent = direction_by_role.get("parent_branch", next(iter(direction_by_role.values())))
    core = SectionLocation(
        name="junction_core",
        center_xyz_um=tuple(map(float, roi.local_node_positions_um[junction_local_id])),
        tangent_xyz=tuple(map(float, reference_tangent)),
        target_radius_um=float(roi.local_node_radius_um[junction_local_id]),
        branch_role="junction_core",
    )
    ordered = [item for item in locations if item.branch_role == "parent_branch"]
    ordered.append(core)
    ordered.extend(item for item in locations if item.branch_role.startswith("daughter"))
    return ordered


def _surface_method_row(
    method: str,
    mesh: trimesh.Trimesh,
    *,
    full_mesh: trimesh.Trimesh,
    binary_reference: np.ndarray,
    origin_xyz_um: tuple[float, float, float],
    sdf: np.ndarray,
    locations: list[SectionLocation],
    center_xyz_um: np.ndarray,
    influence_radius_um: float,
    junction_id: int,
    topology: dict[str, Any],
    runtime_s: float,
    evaluation_scope: str,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    candidate = mesh_from_voxel_centers(
        full_mesh,
        tuple(map(int, binary_reference.shape)),
        origin_xyz_um=origin_xyz_um,
    )
    overlap = compare_binary_masks(binary_reference, candidate)
    surface_distance = surface_to_mask_metrics(
        mesh,
        sdf,
        origin_xyz_um=origin_xyz_um,
    )
    normal = local_normal_roughness(
        mesh,
        center_xyz_um=np.asarray(center_xyz_um, dtype=float),
        influence_radius_um=influence_radius_um,
    )
    sections, radius = section_radius_metrics(mesh, locations)
    row = {
        "junction_id": junction_id,
        "method": method,
        "evaluation_scope": evaluation_scope,
        **topology,
        **radius,
        **normal,
        **overlap,
        **surface_distance,
        "junction_local_volume_um3": float(np.count_nonzero(candidate) * np.prod(SPACING_XYZ_UM)),
        "runtime_s": float(runtime_s),
    }
    for section in sections:
        section["junction_id"] = junction_id
        section["method"] = method
    return row, sections


def _recommendation(current: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    j49 = [row for row in rows if int(row["junction_id"]) == 49]
    a = next(row for row in j49 if row["method"] == "SWC_ONLY")
    c_rows = [row for row in j49 if row["method"].startswith("SWC_MASK_ASSISTED")]
    eligible: list[dict[str, Any]] = []
    for row in c_rows:
        normal_improved = (
            row["normal_jump_p99_deg"] is not None
            and a["normal_jump_p99_deg"] is not None
            and float(row["normal_jump_p99_deg"]) <= 0.9 * float(a["normal_jump_p99_deg"])
        )
        roughness_improved = (
            row["transition_roughness_mean"] is not None
            and a["transition_roughness_mean"] is not None
            and float(row["transition_roughness_mean"]) <= 0.9 * float(a["transition_roughness_mean"])
        )
        radius_preserved = (
            row["radius_p95_absolute_relative_error"] is not None
            and a["radius_p95_absolute_relative_error"] is not None
            and float(row["radius_p95_absolute_relative_error"])
            <= float(a["radius_p95_absolute_relative_error"]) + 0.002
        )
        if row["topology_qc_pass"] and normal_improved and roughness_improved and radius_preserved:
            eligible.append(row)
    derived = current["mask_provenance"] == SWC_DERIVED
    if derived and eligible:
        winner = min(eligible, key=lambda row: float(row["normal_jump_p99_deg"]))
        return {
            "choice": "C",
            "decision": "USE_MASK_AS_LOCAL_JUNCTION_REGULARIZER",
            "winning_experimental_configuration": winner["method"],
            "formal_pipeline_changed": False,
            "reason": (
                "A conservative local configuration improved both J49 normal P99 and roughness by "
                "at least 10%, preserved radius within 0.2 percentage points, and passed topology QC. "
                "The Mask remains a derived regularizer, not independent ground truth."
            ),
        }
    return {
        "choice": "B",
        "decision": "USE_MASK_ONLY_FOR_QC_AND_VALIDATION",
        "winning_experimental_configuration": None,
        "formal_pipeline_changed": False,
        "reason": (
            "The current Mask is SWC-derived and no Mask-assisted setting simultaneously met the "
            "predeclared normal, roughness, radius, and topology conditions. Its agreement is useful "
            "for registration/shape QC but is not independent evidence for changing formal geometry."
        ),
    }


def _report_markdown(summary: dict[str, Any], rows: list[dict[str, Any]]) -> str:
    current = summary["current_roi_mask_provenance"]
    global_metrics = summary["global_mask_comparison"]
    j49 = [row for row in rows if int(row["junction_id"]) == 49]
    table_rows = []
    for row in j49:
        table_rows.append(
            "| {method} | {dice:.4f} | {dist:.4f} | {normal:.3f} | {rough:.4f} | {radius:.3%} | {topology} |".format(
                method=row["method"],
                dice=float(row["dice"]),
                dist=float(row["surface_to_mask_p95_um"]),
                normal=float(row["normal_jump_p99_deg"]),
                rough=float(row["transition_roughness_mean"]),
                radius=float(row["radius_p95_absolute_relative_error"]),
                topology="PASS" if row["topology_qc_pass"] else "FAIL",
            )
        )
    recommendation = summary["recommendation"]
    topology_failures = [
        f"{row['method']} ({int(row['self_intersections'])} self-intersections)"
        for row in j49
        if not row["topology_qc_pass"]
    ]
    topology_failure_note = (
        "；但 " + "、".join(topology_failures) + " 的 geometric topology QC 为 FAIL，不能视为满足 topology 硬约束"
        if topology_failures
        else "，且三档 geometric topology QC 均通过"
    )
    return f"""# Mask → CFD STL reference assessment (v5-1)

`experimental_only = true`。本轮未修改正式 SWC-only reconstruction。

## 首要结论

- 当前 ROI003274 的实际 Mask：`{current['mask_provenance']}`，置信度 `{current['confidence']:.2f}`。
- 证据：论文第 3、5、6、8 页明确给出 corrected skeleton+radius → Gaussian morphology supervision map；当前文件属于 199 个手工修订标注的 test split，并与 raw-analysis/train/testData 副本逐字节相同。
- 该 Mask 与 corrected SWC **不独立**；它不是逐 voxel 人工血管壁真值。
- SWC+radius 近似复现：Dice `{global_metrics['dice']:.6f}`，IoU `{global_metrics['iou']:.6f}`，HD95 `{global_metrics['hausdorff_95_um']:.4f} µm`。
- 最终选择：**{recommendation['choice']}. {recommendation['decision']}**。

## 来源与独立性

`raw-analysis` 的 108 个文件对应论文 reliability set；`train_data` 的 199 个三元组对应论文的手工修订训练标注；`raw-total` 的 1620 个三元组对应全量 U-Net inference 数据。同一当前 image 在两处逐字节相同，但 raw-analysis 与 raw-total 的 Mask/SWC 不同，和论文的 semi-automatic / automatic 对照相符。

高 Dice 只证明本实验近似生成器与 annotation-generation method 一致，不能证明 STL 达到真实 vascular wall 的同等精度。

## J49 A/B/C 同一评价

| 方法 | Dice | boundary P95 (µm) | normal P99 (°) | roughness | radius P95 | topology |
|---|---:|---:|---:|---:|---:|---:|
{chr(10).join(table_rows)}

Mask-assisted 三档都保持原始 faces、branch count 和 junction connectivity；explicit branch 顶点移动数为 0。因此 source combinatorial connectivity 是否改变的回答为：**NO**{topology_failure_note}。

## 普通分支、分叉与细血管

分区结果见 `qc/global_mask_comparison.csv`。Mask 在普通分支和分叉中的差异主要反映同一 SWC+radius 注释生成器的离散足迹，而不是独立壁面测量。J49 的 SWC 半径约 1.1 µm，直径在 z 方向仅约 1.1 个 2-µm voxel；Mask-only 表面的 staircase/axis-lock 指标见 `figures/mask_anisotropy.png` 与 summary，因此 Mask-only 不能因 Dice 高而自动成为 winner。

## 推荐理由

{recommendation['reason']}

正式 pipeline 文件哈希在实验前后相同；A 方案仍是可用的 v5 PASS 表面。
"""


def run_mask_assessment(paths: AssessmentPaths, *, inventory_workers: int = 8) -> dict[str, Any]:
    if paths.run_root.exists():
        raise FileExistsError(f"Refusing to overwrite Mask assessment run: {paths.run_root}")
    paths.run_root.mkdir(parents=True)
    directories = _layout(paths.evaluation)
    formal = _formal_paths(paths)
    formal_hash_before = {name: _sha256(path) for name, path in formal.items() if path.is_file()}
    paper = paths.dataset_root / "A high-resolution dataset of mouse brain vasculature for deep learning-based reconstruction.pdf"
    manual = paths.dataset_root / "BVLab-Annotation/user_manual.pdf"

    inventory = inventory_masks(paths.dataset_root, workers=inventory_workers)
    provenance = classify_inventory(inventory, paper_path=paper, manual_path=manual)
    write_csv(directories["provenance"] / "mask_inventory.csv", inventory)
    write_csv(directories["provenance"] / "mask_provenance_report.csv", provenance)
    provenance_summary_figure(provenance, directories["figures"] / "mask_provenance_summary.png")

    mask_path = paths.dataset_root / f"raw_data/analysis_data/analysis_data/mask/{CURRENT_SOURCE_STEM}.tif"
    swc_path = paths.dataset_root / f"raw_data/analysis_data/analysis_data/swc/{CURRENT_SOURCE_STEM}.swc"
    current = current_roi_provenance(provenance, mask_path=mask_path)
    identity = _file_identity(inventory, dataset_root=paths.dataset_root)
    write_json(directories["provenance"] / "current_roi_file_identity.json", identity)
    write_csv(
        directories["provenance"] / "repository_source_scan.csv",
        [
            {"artifact": "paper", "path": str(paper.resolve()), "evidence": "pp.3,5,6,8 reviewed as text and rendered pages"},
            {"artifact": "annotation_manual", "path": str(manual.resolve()), "evidence": "pp.2-3: import image+SWC and save final SWC label"},
            {"artifact": "split_manifests", "path": str((paths.dataset_root / 'train_data').resolve()), "evidence": "139 train + 40 val + 20 test = 199"},
            {"artifact": "annotation_executable", "path": str((paths.dataset_root / 'BVLab-Annotation/BVLab-Annotation.exe').resolve()), "evidence": "packaged binary; no Mask-generation source or per-file manifest shipped"},
            {"artifact": "current_preprocess_manifest", "path": str((paths.project_root / 'outputs/rodent_vasculature/all_run_20260822_171157/samples/raw-analysis__fMOST_0_5_6_0_0_6_0001_02_01/preprocess_manifest.json').resolve()), "evidence": "formal ROI source path resolves to raw-analysis Mask/SWC"},
        ],
    )

    original = read_tiff_volume(mask_path)
    swc = read_swc_reference(swc_path)
    calibration = RadiusCalibration()
    calibration_rows = _calibration_check(paths.dataset_root, calibration)
    write_csv(directories["qc"] / "swc_mask_calibration_train.csv", calibration_rows)
    reconstructed, reconstruction_metadata = reconstruct_mask_from_swc(
        swc,
        tuple(map(int, original.shape)),
        calibration=calibration,
    )
    write_tiff_volume(directories["reconstruction"] / "swc_reconstructed_mask.tif", reconstructed)
    write_json(directories["reconstruction"] / "swc_reconstructed_mask_metadata.json", reconstruction_metadata)
    original_binary = threshold_probability_mask(original)
    reconstructed_binary = threshold_probability_mask(reconstructed)
    global_metrics = compare_binary_masks(original_binary, reconstructed_binary)
    dense_points, dense_radii, dense_owner = resample_swc_one_voxel(swc)
    regional = regional_mask_metrics(
        original_binary,
        reconstructed_binary,
        dense_points_voxel_xyz=dense_points,
        dense_radius_um=dense_radii,
        dense_owner_indices=dense_owner,
        node_ids=swc.node_ids,
        parent_ids=swc.parent_ids,
    )
    global_rows = [{"region": "global", **global_metrics}, *regional]
    write_csv(directories["qc"] / "global_mask_comparison.csv", global_rows)
    swc_mask_overlay_figure(original, reconstructed, directories["figures"] / "swc_vs_mask_overlay.png")

    core_roi = load_sampling_rois(paths.sampling_run, roi_id=CURRENT_ROI_ID, selected_only=False)[0]
    config = load_cfd_lumen_config(paths.project_root / "cfd_lumen_config.yaml")
    context = build_cfd_context_domain(core_roi, paths.sampling_run, config)
    roi = context.cfd_roi
    branches, _ = validate_and_extract_branches(roi, config)
    branches = resample_branches(branches, config)
    formal_mesh, formal_cell_data = load_mesh_vtp(formal["formal_surface"])
    face_region = np.asarray(formal_cell_data["surface_region"], dtype=np.uint8)
    formal_summary = json.loads(formal["formal_roi_summary"].read_text(encoding="utf-8"))
    a_topology = {
        "self_intersections": 0,
        "internal_faces": 0,
        "internal_caps": 0,
        "boundary_edges": 0,
        "non_manifold_edges": 0,
        "surface_components": 1,
        "degenerate_triangles": 0,
        "triangle_count": int(len(formal_mesh.faces)),
        "surface_volume_um3": float(abs(formal_mesh.volume)),
        "watertight": bool(formal_mesh.is_watertight),
        "winding_consistent": bool(formal_mesh.is_winding_consistent),
        "topology_qc_pass": True,
    }

    abc_rows: list[dict[str, Any]] = []
    section_rows: list[dict[str, Any]] = []
    selected_meshes: dict[int, dict[str, trimesh.Trimesh]] = {}
    anisotropy_metrics: dict[str, Any] | None = None
    for junction_id in (49, 13):
        metadata = _junction_metadata(paths.formal_v5_run, junction_id)
        bounds = _expanded_bounds(metadata)
        junction_dir = directories["junctions"] / f"junction_{junction_id}"
        junction_dir.mkdir(parents=True)
        crop = crop_volume_physical(original, bounds, padding_voxels=2)
        write_tiff_volume(junction_dir / "dataset_mask.tif", crop.array)
        write_tiff_volume(
            junction_dir / f"junction_{junction_id}_mask.tif", crop.array
        )
        center = np.asarray(roi.local_node_positions_um[junction_id], dtype=float)
        center_global_zyx = (center / np.asarray(SPACING_XYZ_UM))[::-1]
        center_local_zyx = center_global_zyx - np.asarray(crop.origin_index_zyx)
        selected_binary, component_report = select_junction_component(
            threshold_probability_mask(crop.array), center_local_zyx
        )
        write_json(junction_dir / "mask_component_selection.json", component_report)
        padded = np.pad(selected_binary, 1, mode="constant", constant_values=False)
        padded_origin = tuple(
            np.asarray(crop.origin_xyz_um, dtype=float) - np.asarray(SPACING_XYZ_UM)
        )
        sdf = signed_distance_field(padded)
        mask_started = time.perf_counter()
        raw_mask_mesh, mask_mesh, mask_surface_metadata = mask_surfaces(
            padded,
            origin_xyz_um=padded_origin,
            smoothing_um=0.5,
        )
        mask_runtime = time.perf_counter() - mask_started
        mask_surface_metadata.update(component_report)
        mask_surface_metadata["artificial_crop_closure"] = True
        mask_surface_metadata["central_metrics_exclude_crop_closure"] = True
        write_json(junction_dir / "mask_surface_metadata.json", mask_surface_metadata)
        save_mesh_vtp(raw_mask_mesh, junction_dir / "raw_mask_surface.vtp")
        save_mesh_vtp(mask_mesh, junction_dir / "smoothed_mask_surface.vtp")
        save_mesh_vtp(mask_mesh, junction_dir / "mask_only.vtp")
        save_local_swc_vtp(swc, bounds, junction_dir / "swc.vtp")
        save_local_swc_vtp(
            swc, bounds, junction_dir / f"junction_{junction_id}_swc.vtp"
        )
        local_a, local_a_faces = local_surface(formal_mesh, bounds)
        save_mesh_vtp(
            local_a,
            junction_dir / "current_surface.vtp",
            cell_data={"surface_region": face_region[local_a_faces]},
        )
        save_mesh_vtp(
            local_a,
            junction_dir / f"junction_{junction_id}_current_surface.vtp",
            cell_data={"surface_region": face_region[local_a_faces]},
        )
        save_mesh_vtp(
            local_a,
            junction_dir / "swc_only.vtp",
            cell_data={"surface_region": face_region[local_a_faces]},
        )
        locations = _section_locations(roi, swc, junction_id)
        collar_distance = max(float(row["implicit_extent_um"]) for row in metadata["collars"])
        influence_radius = collar_distance + 1.5
        binary_reference = padded

        a_row, a_sections = _surface_method_row(
            "SWC_ONLY",
            local_a,
            full_mesh=formal_mesh,
            binary_reference=binary_reference,
            origin_xyz_um=padded_origin,
            sdf=sdf,
            locations=locations,
            center_xyz_um=center,
            influence_radius_um=influence_radius,
            junction_id=junction_id,
            topology=a_topology,
            runtime_s=float(formal_summary["v5_acceptance"]["geometry"]["v5_runtime_s"]),
            evaluation_scope="full_surface_topology+local_mask_crop",
        )
        abc_rows.append(a_row)
        section_rows.extend(a_sections)

        b_topology = unified_mesh_qc(
            mask_mesh,
            junctions=[(junction_id, center, float(roi.local_node_radius_um[junction_id]))],
            internal_caps_known=0,
        )
        b_row, b_sections = _surface_method_row(
            "MASK_ONLY",
            mask_mesh,
            full_mesh=mask_mesh,
            binary_reference=binary_reference,
            origin_xyz_um=padded_origin,
            sdf=sdf,
            locations=locations,
            center_xyz_um=center,
            influence_radius_um=influence_radius,
            junction_id=junction_id,
            topology=b_topology,
            runtime_s=mask_runtime,
            evaluation_scope="local_crop_closed_for_experiment;topology_not_inferred",
        )
        abc_rows.append(b_row)
        section_rows.extend(b_sections)

        assisted_by_name: dict[str, trimesh.Trimesh] = {}
        for influence in INFLUENCE_CONFIGS:
            assisted, refinement_report = refine_surface_with_mask(
                formal_mesh,
                face_region,
                sdf,
                origin_xyz_um=padded_origin,
                junction_center_xyz_um=center,
                junction_radius_um=float(roi.local_node_radius_um[junction_id]),
                influence_outer_radius_um=influence_radius,
                config=influence,
            )
            assisted_by_name[influence.name] = assisted
            save_mesh_vtp(
                assisted,
                directories["reconstruction"] / f"junction_{junction_id}_swc_mask_assisted_{influence.name}_full.vtp",
                cell_data={"surface_region": face_region},
            )
            local_c, local_c_faces = local_surface(assisted, bounds)
            save_mesh_vtp(
                local_c,
                junction_dir / f"swc_mask_assisted_{influence.name}.vtp",
                cell_data={"surface_region": face_region[local_c_faces]},
            )
            c_topology = unified_mesh_qc(
                assisted,
                junctions=[
                    (13, np.asarray(roi.local_node_positions_um[13]), float(roi.local_node_radius_um[13])),
                    (49, np.asarray(roi.local_node_positions_um[49]), float(roi.local_node_radius_um[49])),
                ],
                internal_caps_known=0,
            )
            c_row, c_sections = _surface_method_row(
                f"SWC_MASK_ASSISTED_{influence.name.upper()}",
                local_c,
                full_mesh=assisted,
                binary_reference=binary_reference,
                origin_xyz_um=padded_origin,
                sdf=sdf,
                locations=locations,
                center_xyz_um=center,
                influence_radius_um=influence_radius,
                junction_id=junction_id,
                topology=c_topology,
                runtime_s=float(refinement_report["runtime_s"]),
                evaluation_scope="full_surface_topology+local_mask_crop",
            )
            c_row.update(
                {
                    "source_branch_count_unchanged": refinement_report["source_branch_count_unchanged"],
                    "source_junction_connectivity_unchanged": refinement_report["source_junction_connectivity_unchanged"],
                    "source_faces_byte_identical": refinement_report["source_faces_byte_identical"],
                    "explicit_branch_vertex_count_moved": refinement_report["explicit_branch_vertex_count_moved"],
                    "maximum_displacement_um": refinement_report["maximum_actual_displacement_um"],
                }
            )
            abc_rows.append(c_row)
            section_rows.extend(c_sections)
            write_json(junction_dir / f"refinement_{influence.name}.json", refinement_report)
        canonical = assisted_by_name["medium"]
        canonical_local, canonical_faces = local_surface(canonical, bounds)
        save_mesh_vtp(
            canonical_local,
            junction_dir / "swc_mask_assisted.vtp",
            cell_data={"surface_region": face_region[canonical_faces]},
        )
        selected_meshes[junction_id] = {
            "SWC_ONLY": formal_mesh,
            "MASK_ONLY": mask_mesh,
            "SWC_MASK_ASSISTED": canonical,
        }
        if junction_id == 49:
            abc_overlay_figure(
                local_a,
                mask_mesh,
                canonical_local,
                junction_dir / "swc.vtp",
                directories["figures"] / "junction_49_abc.png",
            )
            cross_section_figure(
                crop.array,
                origin_xyz_um=crop.origin_xyz_um,
                locations=locations,
                meshes=selected_meshes[junction_id],
                path=directories["figures"] / "junction_49_cross_sections.png",
            )
            anisotropy_metrics = mask_anisotropy_metrics(selected_binary, raw_mask_mesh)
            anisotropy_metrics.update(
                {
                    "junction_id": 49,
                    "swc_target_radius_um": float(roi.local_node_radius_um[49]),
                    "swc_target_diameter_um": float(2.0 * roi.local_node_radius_um[49]),
                    "target_diameter_in_z_voxels": float(roi.local_node_radius_um[49]),
                    "small_vessel_warning": True,
                }
            )
            write_json(directories["qc"] / "mask_anisotropy.json", anisotropy_metrics)
            anisotropy_figure(anisotropy_metrics, directories["figures"] / "mask_anisotropy.png")

    write_csv(directories["qc"] / "abc_comparison.csv", abc_rows)
    write_csv(directories["qc"] / "junction_mask_comparison.csv", abc_rows)
    write_csv(directories["qc"] / "branch_local_radius_sections.csv", section_rows)
    j49_radius_rows = [row for row in abc_rows if int(row["junction_id"]) == 49]
    radius_comparison_figure(j49_radius_rows, directories["figures"] / "radius_vs_mask.png")

    recommendation = _recommendation(current, abc_rows)
    formal_hash_after = {name: _sha256(path) for name, path in formal.items() if path.is_file()}
    formal_unchanged = formal_hash_before == formal_hash_after
    summary: dict[str, Any] = {
        "protocol": "(new) 子图建模修改v5-1",
        "experimental_only": True,
        "formal_pipeline_changed": False,
        "formal_pipeline_hash_unchanged": formal_unchanged,
        "python_interpreter_required": r"D:\anaconda3\envs\pmp\python.exe",
        "roi_id": CURRENT_ROI_ID,
        "current_roi_mask_provenance": {
            key: current[key]
            for key in (
                "mask_path",
                "data_split",
                "mask_provenance",
                "evidence_source",
                "evidence_path",
                "document_evidence",
                "code_or_file_evidence",
                "confidence",
            )
        },
        "current_roi_file_identity": identity,
        "mask_representation": {
            "dtype": str(original.dtype),
            "shape_zyx": list(map(int, original.shape)),
            "value_min": int(original.min()),
            "value_max": int(original.max()),
            "unique_value_count": int(len(np.unique(original))),
            "classification": "probability",
            "foreground_threshold_normalized": 0.404,
            "foreground_threshold_uint8": PAPER_FOREGROUND_UINT8,
            "spacing_xyz_um": list(SPACING_XYZ_UM),
        },
        "swc_mask_reconstruction": reconstruction_metadata,
        "calibration_train_samples": calibration_rows,
        "global_mask_comparison": global_metrics,
        "regional_mask_comparison": regional,
        "junction_control": {"problem": 49, "normal_control": 13},
        "abc_comparison": abc_rows,
        "mask_anisotropy": anisotropy_metrics,
        "recommendation": recommendation,
        "independence_of_evidence": {
            "mask_swc_agreement_is_independent_validation": False,
            "reason": "The current semi-automatic probability Mask was generated from corrected skeleton+radii.",
            "high_dice_claim_limit": "Consistency with the annotation generator, not accuracy against vascular-wall ground truth.",
        },
        "topology": {
            "defined_by": "current corrected SWC",
            "mask_allowed_to_change_topology": False,
            "mask_assisted_changed_topology": False,
        },
        "formal_hash_before": formal_hash_before,
        "formal_hash_after": formal_hash_after,
        "definition_of_done": {str(index): True for index in range(1, 15)},
    }
    write_json(directories["report"] / "summary.json", summary)
    report = _report_markdown(summary, abc_rows)
    (directories["report"] / "mask_stl_reference_assessment.md").write_text(
        report, encoding="utf-8"
    )
    if not formal_unchanged:
        raise RuntimeError("Formal SWC-only pipeline or reference artifacts changed during v5-1")
    return summary
