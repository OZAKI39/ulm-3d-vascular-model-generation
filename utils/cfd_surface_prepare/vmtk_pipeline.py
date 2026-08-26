"""Formal VMTK TPS boundary-normal cross-seam active-collar pipeline."""

from __future__ import annotations

import copy
import json
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv

from utils.cfd_lumen.model_yaml_config import load_swc_stl_yaml_config
from utils.cfd_lumen.roi_io import load_sampling_rois
from utils.cfd_lumen.ultraliser_qc import evaluate_radius_fidelity, validate_source_roi

from .config import SurfacePrepareConfig
from .export import write_csv
from .guarded_remesh import (
    assign_cross_seam_active_entities,
    assign_guarded_remesh_entities,
    cross_seam_entity_preservation_qc,
    cross_seam_intersection_qc,
    crossseam_active_mesh_quality,
    cut_seam_topology_qc,
    diagnose_previous_entityremesh,
    guarded_entity_preservation_qc,
    guarded_intersection_qc,
    guarded_region_mesh_quality,
    locked_entity_preservation_qc,
    previous_tail_fraction_inside_guard,
    seam_local_normal_diagnostic,
)
from .io import (
    SurfacePrepareError,
    load_original_surface,
    load_surface_inputs,
    read_json,
    sha256_file,
)
from .local_cut import local_plane_cut
from .mesh_quality import measure_local_original_mesh
from .types import TaggedSurface
from .vmtk_qc import (
    assign_remesh_entities,
    core_exact_preservation_qc,
    core_symmetric_distance_qc,
    entity_remesh_preservation_qc,
    extension_geometry_qc,
    extension_mesh_quality_from_surface,
    geometry_pressure_correction,
    interface_smoothness_from_raw,
    meter_scale_qc,
    normal_consistency_qc,
    open_profile_qc,
    polydata_mesh,
    previous_global_remesh_diagnostics,
    raw_core_exact_copy_qc,
    tag_and_export_final_surface,
    topology_qc,
    write_json,
)
from .vmtk_runner import (
    cap_official_vmtk,
    entity_remesh_official_vmtk,
    exchange_paths,
    official_source_provenance,
    parameter_mapping,
    run_official_vmtk,
)
from .vmtk_visualization import (
    save_crossseam_review_figures,
    save_entityremesh_review_figures,
    save_guarded_open_figures,
    save_guarded_three_way_comparison,
    save_previous_collision_closeup,
)


EXCLUSION_FAILURE = "VMTK_CROSS_SEAM_ENTITY_ASSIGNMENT_FAILED"
ASSIGNMENT_FAILURE = "VMTK_CROSS_SEAM_ENTITY_ASSIGNMENT_FAILED"
SUSPICIOUS_GUARD_FAILURE = "VMTK_CROSS_SEAM_ENTITY_ASSIGNMENT_FAILED"
CORE_FAILURE = "VMTK_CROSS_SEAM_FAR_CORE_MODIFIED"
GUARD_FAILURE = "VMTK_CROSS_SEAM_FAR_CORE_MODIFIED"
NO_EFFECT_FAILURE = "VMTK_CROSS_SEAM_GEOMETRY_FAILED"
GEOMETRY_FAILURE = "VMTK_CROSS_SEAM_GEOMETRY_FAILED"
TOPOLOGY_FAILURE = "VMTK_CROSS_SEAM_TOPOLOGY_FAILED"
BOUNDARY_FAILURE = "VMTK_CROSS_SEAM_GEOMETRY_FAILED"
RADIUS_FAILURE = "VMTK_CROSS_SEAM_GEOMETRY_FAILED"
COLLAR_FAILURE = "VMTK_CROSS_SEAM_COLLAR_GEOMETRY_FAILED"
RING_FAILURE = "VMTK_CROSS_SEAM_RING_TOPOLOGY_PRESERVED"
METER_FAILURE = "VMTK_METER_SCALE_SERIALIZATION_FAILED"
SUCCESS_STATUS = "VMTK_TPS_BOUNDARY_NORMAL_CROSS_SEAM_REMESH_PASS_PENDING_MANUAL_REVIEW"
TAIL_WARNING_STATUS = (
    "VMTK_TPS_BOUNDARY_NORMAL_CROSS_SEAM_REMESH_MESH_TAIL_WARNING_"
    "PENDING_MANUAL_REVIEW"
)
PASS_STATUSES = {SUCCESS_STATUS, TAIL_WARNING_STATUS}
SUCCESS_NEXT = "MANUALLY REVIEW CROSS-SEAM REMESHED CFD SURFACE"
FAIL_NEXT = "REVIEW CROSS-SEAM REMESH FAILURE"
PREVIOUS_GLOBAL_REMESH_RUN_ID = "vmtk_tps_boundarynormal_anchor003274_20260826_153709"
PREVIOUS_CAPONLY_RUN_ID = (
    "vmtk_tps_boundarynormal_caponly_anchor003274_20260826_165659"
)
PREVIOUS_ENTITY_REMESH_RUN_ID = (
    "vmtk_tps_boundarynormal_entityremesh_anchor003274_20260826_180737"
)
GUARDED_ENTITY_REMESH_RUN_ID = (
    "vmtk_tps_boundarynormal_guarded_entityremesh_anchor003274_20260826_185752"
)


@dataclass(frozen=True, slots=True)
class VmtkLayout:
    root: Path
    input: Path
    vmtk: Path
    geometry: Path
    boundaries: Path
    bc: Path
    qc: Path
    figures: Path
    report: Path


@dataclass(frozen=True, slots=True)
class VmtkSurfacePrepareResult:
    status: str
    next_stage: str
    run_root: Path
    runtime: dict[str, Any]
    assignment_qc: dict[str, Any]
    far_core_qc: dict[str, Any]
    boundary_lock_qc: dict[str, Any]
    collar_qc: dict[str, Any]
    extension_qc: dict[str, Any]
    seam_topology_qc: dict[str, Any]
    intersection_qc: dict[str, Any]
    boundaries: tuple[dict[str, Any], ...]
    mesh_quality: dict[str, Any]
    seam_quality: tuple[dict[str, Any], ...]
    open_topology: dict[str, Any]
    final_topology: dict[str, Any]
    radius_qc: dict[str, Any]
    meter_qc: dict[str, Any]
    collision_count: int
    pressure_rows: tuple[dict[str, Any], ...]
    output_paths: dict[str, Any]
    figures: tuple[Path, ...]


def _layout(output_root: Path, run_id: str) -> VmtkLayout:
    root = Path(output_root) / run_id
    if root.exists():
        raise FileExistsError(f"Refusing to overwrite VMTK run: {root}")
    directories = [
        root / name
        for name in (
            "input",
            "vmtk",
            "geometry",
            "boundaries",
            "bc",
            "qc",
            "figures",
            "report",
        )
    ]
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=False)
    return VmtkLayout(root.resolve(), *[path.resolve() for path in directories])


def _anchor_id(roi_id: str) -> str:
    marker = "__anchor_"
    return roi_id.rsplit(marker, 1)[1].split("__", 1)[0] if marker in roi_id else "unknown"


def _open_surface(surface: TaggedSurface, path: Path) -> None:
    surface.compact()
    faces = np.column_stack(
        (np.full(len(surface.faces), 3, dtype=np.int64), surface.faces)
    ).ravel()
    data = pv.PolyData(np.asarray(surface.vertices, dtype=float), faces)
    data.point_data["source_vertex_index"] = surface.source_vertex_index
    data.cell_data["surface_role"] = np.zeros(len(surface.faces), dtype=np.uint8)
    data.save(path, binary=True)


def _reference_paths(output_root: Path) -> dict[str, Path]:
    global_root = Path(output_root) / PREVIOUS_GLOBAL_REMESH_RUN_ID
    caponly_root = Path(output_root) / PREVIOUS_CAPONLY_RUN_ID
    previous_entity_root = Path(output_root) / PREVIOUS_ENTITY_REMESH_RUN_ID
    paths = {
        "global_root": global_root,
        "global_open_vtp": global_root / "input" / "open_surface_um.vtp",
        "global_raw_vtp": global_root / "geometry" / "vmtk_boundarynormal_raw_um.vtp",
        "global_final_vtp": (
            global_root / "geometry" / "cfd_surface_vmtk_tps_boundarynormal_um.vtp"
        ),
        "global_radius_qc": global_root / "qc" / "radius_fidelity.json",
        "global_summary": global_root / "qc" / "run_summary.json",
        "caponly_root": caponly_root,
        "caponly_raw_vtp": caponly_root / "geometry" / "vmtk_boundarynormal_raw_um.vtp",
        "caponly_final_vtp": (
            caponly_root
            / "geometry"
            / "cfd_surface_vmtk_tps_boundarynormal_caponly_um.vtp"
        ),
        "caponly_summary": caponly_root / "qc" / "run_summary.json",
        "previous_entity_root": previous_entity_root,
        "previous_entity_raw_vtp": (
            previous_entity_root / "geometry" / "vmtk_boundarynormal_raw_um.vtp"
        ),
        "previous_entity_remeshed_vtp": (
            previous_entity_root
            / "geometry"
            / "vmtk_boundarynormal_extension_remeshed_open_um.vtp"
        ),
        "previous_entity_summary": previous_entity_root / "qc" / "run_summary.json",
    }
    missing = [
        str(path)
        for key, path in paths.items()
        if not key.endswith("root") and not path.is_file()
    ]
    if missing:
        raise SurfacePrepareError("Missing read-only VMTK reference: " + ", ".join(missing))
    return {key: path.resolve() for key, path in paths.items()}


def _crossseam_reference_paths(output_root: Path) -> dict[str, Path]:
    """Resolve only the frozen GLOBAL and GUARDED comparison references."""

    global_root = Path(output_root) / PREVIOUS_GLOBAL_REMESH_RUN_ID
    guarded_root = Path(output_root) / GUARDED_ENTITY_REMESH_RUN_ID
    paths = {
        "global_root": global_root,
        "global_final_vtp": (
            global_root / "geometry" / "cfd_surface_vmtk_tps_boundarynormal_um.vtp"
        ),
        "global_summary": global_root / "qc" / "run_summary.json",
        "guarded_root": guarded_root,
        "guarded_final_vtp": (
            guarded_root
            / "geometry"
            / "cfd_surface_vmtk_tps_boundarynormal_guarded_entityremesh_um.vtp"
        ),
        "guarded_summary": guarded_root / "qc" / "run_summary.json",
    }
    missing = [
        str(path)
        for key, path in paths.items()
        if not key.endswith("root") and not path.is_file()
    ]
    if missing:
        raise SurfacePrepareError(
            "Missing read-only cross-seam comparison reference: "
            + ", ".join(missing)
        )
    return {key: path.resolve() for key, path in paths.items()}


def _manifest_records(boundaries: tuple[Any, ...]) -> tuple[dict[str, Any], ...]:
    return tuple(
        {
            "boundary_index": boundary.index,
            "port_id": boundary.port_id,
            "boundary_origin": boundary.boundary_origin,
            "role": boundary.role,
            "global_node_id": boundary.global_node_id,
            "global_edge_id": boundary.global_edge_id,
            "center_um": boundary.center_um.tolist(),
            "source_radius_um": boundary.source_radius_um,
            "planned_extension_length_um": boundary.extension_length_um,
            "centerline_usage": "NOT_GENERATED_BOUNDARY_NORMAL_MODE",
        }
        for boundary in boundaries
    )


def _write_manifest(
    path: Path,
    records: tuple[dict[str, Any], ...],
    profile: dict[str, Any],
) -> None:
    profile_by_id = {row["port_id"]: row for row in profile["boundaries"]}
    write_csv(
        path,
        [
            {
                **record,
                "open_profile_area_um2": profile_by_id[record["port_id"]]["area_um2"],
                "open_profile_point_count": profile_by_id[record["port_id"]]["point_count"],
            }
            for record in records
        ],
    )


def _radius_qc(*, final_mesh: Any, inputs: Any, project_root: Path) -> dict[str, Any]:
    model_root = Path(str(inputs.geometry_reference["run_root"]))
    source_config = model_root / "input" / "source_swc_stl_model_generate.yaml"
    model_config = load_swc_stl_yaml_config(source_config, project_root=project_root)
    if model_config.sampling_run is None:
        raise SurfacePrepareError("VMTK_SURFACE_QC_FAILED:missing_sampling_run")
    roi = load_sampling_rois(
        model_config.sampling_run,
        roi_id=str(inputs.geometry_reference["roi_id"]),
        selected_only=False,
    )[0]
    branches, source_report = validate_source_roi(roi)
    samples, report = evaluate_radius_fidelity(final_mesh, branches, roi, model_config.lumen)
    report["source_qc"] = source_report
    report["samples"] = [sample.report() for sample in samples]
    report["required_p95_maximum"] = 0.05
    report["status"] = (
        "PASS"
        if report["p95_absolute_relative_error"] is not None
        and float(report["p95_absolute_relative_error"]) <= 0.05
        else "FAIL"
    )
    return report


def _boundary_conditions(
    original: dict[str, Any], pressure_rows: list[dict[str, Any]]
) -> dict[str, Any]:
    output = copy.deepcopy(original)
    rows = {row["port_id"]: row for row in pressure_rows}
    inlet = output["inlet"]
    inlet_row = rows[inlet["port_id"]]
    inlet.update(
        {
            "Q_solver_m3_s": inlet_row["Q_solver_m3_s"],
            "profile": "PARABOLIC",
            "extension_pressure_correction_applied": False,
        }
    )
    for outlet in output["outlets"]:
        row = rows[outlet["port_id"]]
        outlet.update(
            {
                "P_original_1D_pa": row["P_original_1D_pa"],
                "Q_expected_1D_m3_s": row["Q_expected_1D_m3_s"],
                "extension_resistance_pa_s_m3": row["extension_resistance_pa_s_m3"],
                "predicted_extension_pressure_drop_pa": row[
                    "predicted_extension_pressure_drop_pa"
                ],
                "P_solver_boundary_pa": row["P_solver_boundary_pa"],
                "extension_pressure_correction_role": row["pressure_correction_role"],
            }
        )
    output["surface_geometry_method"] = (
        "OFFICIAL_VMTK_TPS_FLOWEXTENSIONS_ENTITY_AWARE_EXTENSION_REMESH"
    )
    output["pressure_correction_geometry"] = (
        "20 cross sections on final entity-aware VMTK geometry"
    )
    return output


def boundary_plane_alignment_pass(profile: dict[str, Any]) -> bool:
    return all(
        float(row["boundary_plane_normal_abs_dot_expected_outward"]) >= 0.999
        for row in profile["boundaries"]
    )


def raw_geometry_hard_gate_pass(*reports: dict[str, Any]) -> bool:
    return all(report.get("status") == "PASS" for report in reports)


def should_promote_raw_candidate(raw_hard_qc_pass: bool) -> bool:
    return raw_hard_qc_pass


def should_generate_manual_review_figures(
    *, raw_hard_qc_pass: bool, interface_status: str = "DIAGNOSTIC"
) -> bool:
    del interface_status
    return raw_hard_qc_pass


def final_candidate_status(
    *, final_hard_qc_pass: bool, interface_status: str = "DIAGNOSTIC"
) -> str:
    del interface_status
    return SUCCESS_STATUS if final_hard_qc_pass else TOPOLOGY_FAILURE


def resolve_entity_failure_status(
    *,
    topology: dict[str, Any],
    boundary_mapping: dict[str, Any],
    core_exact: dict[str, Any],
    core_distance: dict[str, Any],
    normal: dict[str, Any],
    radius: dict[str, Any],
    meter: dict[str, Any],
    pressure: dict[str, Any],
    integrity: dict[str, Any],
) -> str:
    if topology.get("status") != "PASS":
        return TOPOLOGY_FAILURE
    if boundary_mapping.get("status") != "PASS":
        return BOUNDARY_FAILURE
    if (
        core_exact.get("status") != "PASS"
        or core_distance.get("status") != "PASS"
        or normal.get("status") != "PASS"
    ):
        return CORE_FAILURE
    if radius.get("status") != "PASS":
        return RADIUS_FAILURE
    if meter.get("status") != "PASS":
        return METER_FAILURE
    if pressure.get("status") != "PASS" or integrity.get("status") != "PASS":
        return GEOMETRY_FAILURE
    return SUCCESS_STATUS


def _write_report(path: Path, summary: dict[str, Any]) -> None:
    runtime = summary.get("runtime", {})
    lines = [
        "# VMTK TPS boundary-normal cross-seam active-collar remesh experiment",
        "",
        f"- Status: `{summary['status']}`",
        "- Official remesher: `vmtkscripts.vmtkSurfaceRemeshing`",
        "- Cell entity array: `RemeshEntityId`",
        "- Remesh entities: `FAR_CORE=1`, `CROSS_SEAM_ACTIVE=2`",
        "- Excluded entity ids: `[1]`",
        "- Active entity ids: `[2]`",
        "- Global surface remeshing performed: `false`",
        f"- VMTK: `{runtime.get('vmtk_package_version', 'NOT_RUN')}`",
        f"- VTK: `{runtime.get('vtk_version', 'NOT_RUN')}`",
        "- Visual acceptance: `MANUAL_REVIEW_REQUIRED`",
        "- Volume mesh created: `false`",
        "- CFD run: `false`",
        "",
        f"NEXT: `{summary['next']}`",
        "",
    ]
    path.write_text("\n".join(lines), encoding="utf-8")


def _write_failure(
    layout: VmtkLayout,
    *,
    status: str,
    run_id: str,
    runtime: dict[str, Any],
    outputs: dict[str, Any],
    evidence: dict[str, Any],
) -> None:
    summary = {
        "status": status,
        "next": FAIL_NEXT,
        "run_id": run_id,
        "run_root": str(layout.root),
        "runtime": runtime,
        "outputs": outputs,
        "evidence": evidence,
        "single_primary_candidate_count": 1,
        "second_candidate_run": False,
        "automatic_fallback": False,
        "global_surface_remeshing_performed": False,
        "volume_mesh_created": False,
        "cfd_run": False,
        "microbubble_simulation_run": False,
    }
    write_json(layout.qc / "run_summary.json", summary)
    _write_report(
        layout.report / "vmtk_tps_boundarynormal_crossseam_report.md", summary
    )


def _synthetic_connected_entities(path: Path) -> None:
    """Create a curved CORE-GUARD-BODY tube for installed-VMTK safety proof."""

    count = 28
    levels = np.linspace(0.0, 7.0, 8)
    angle = np.linspace(0.0, 2.0 * np.pi, count, endpoint=False)
    rings: list[np.ndarray] = []
    for level in levels:
        center = np.asarray((0.08 * level**2, 0.0, level))
        tangent = np.asarray((0.16 * level, 0.0, 1.0))
        tangent /= np.linalg.norm(tangent)
        first = np.asarray((0.0, 1.0, 0.0))
        second = np.cross(tangent, first)
        second /= np.linalg.norm(second)
        rings.append(
            center
            + np.cos(angle)[:, None] * first
            + np.sin(angle)[:, None] * second
        )
    points = np.vstack(rings)
    faces: list[list[int]] = []
    entities: list[int] = []
    for layer in range(len(levels) - 1):
        for point in range(count):
            following = (point + 1) % count
            lower = layer * count
            upper = (layer + 1) * count
            faces.extend(
                (
                    [lower + point, lower + following, upper + following],
                    [lower + point, upper + following, upper + point],
                )
            )
            if layer == 0:
                entity = 1
            elif layer < 3:
                entity = 2
            else:
                entity = 3
            entities.extend((entity, entity))
    vtk_faces = np.column_stack(
        (np.full(len(faces), 3, dtype=np.int64), np.asarray(faces, dtype=np.int64))
    ).ravel()
    data = pv.PolyData(points, vtk_faces)
    entity_array = np.asarray(entities, dtype=np.int32)
    data.cell_data["RemeshEntityId"] = entity_array
    data.cell_data["SurfaceRegionId"] = np.where(entity_array == 1, 0, 1).astype(
        np.uint8
    )
    data.cell_data["SurfaceRegion"] = np.where(
        entity_array == 1, "CORE", "EXTENSION"
    )
    data.save(path, binary=True)


def run_synthetic_entity_exclusion_preflight(
    config: SurfacePrepareConfig,
    *,
    root: Path,
    tool_script: Path,
) -> dict[str, Any]:
    input_directory = root / "input"
    vmtk_directory = root / "vmtk"
    geometry_directory = root / "geometry"
    for directory in (input_directory, vmtk_directory, geometry_directory):
        directory.mkdir(parents=True, exist_ok=True)
    paths = exchange_paths(
        input_directory=input_directory,
        vmtk_directory=vmtk_directory,
        geometry_directory=geometry_directory,
        extension_mode="boundarynormal",
    )
    _synthetic_connected_entities(paths.raw_vtp)
    invocation = entity_remesh_official_vmtk(
        config=config.vmtk, paths=paths, tool_script=tool_script
    )
    core, guard, boundaries, body = guarded_entity_preservation_qc(
        paths.raw_vtp,
        paths.remeshed_open_vtp,
        entity_array_name=config.vmtk.entity_remesh.entity_array_name,
        core_entity_id=config.vmtk.entity_remesh.core_entity_id,
        guard_entity_id=config.vmtk.entity_remesh.guard_entity_id,
        body_entity_id=config.vmtk.entity_remesh.extension_body_entity_id,
    )
    _, mesh = polydata_mesh(paths.remeshed_open_vtp)
    topology, _ = topology_qc(
        mesh, expected_open_profile_count=2, allow_degenerate=False
    )
    intersections, _ = guarded_intersection_qc(paths.remeshed_open_vtp)
    topology["checks"]["zero_self_intersections"] = (
        intersections["true_self_intersection_count"] == 0
    )
    topology["self_intersection_count"] = intersections[
        "true_self_intersection_count"
    ]
    topology["status"] = (
        "PASS" if all(topology["checks"].values()) else "FAIL"
    )
    passed = all(
        report["status"] == "PASS"
        for report in (core, guard, boundaries, body, topology, intersections)
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "purpose": "curved CORE-GUARD-BODY ExcludeEntityIds safety proof",
        "runtime": invocation.runtime,
        "core": core,
        "guard": guard,
        "entity_boundaries": boundaries,
        "body": body,
        "topology": topology,
        "intersections": intersections,
    }


def _flatten_quality(report: dict[str, Any]) -> list[dict[str, Any]]:
    fields = {
        "triangle_count": "triangle_count",
        "edge_median_um": "edge_length_median_um",
        "edge_P95_um": "edge_length_P95_um",
        "minimum_angle_deg": "minimum_angle_deg",
        "angle_P05_deg": "angle_P05_deg",
        "aspect_P95": "aspect_ratio_P95",
        "aspect_max": "aspect_ratio_max",
        "neighbor_area_ratio_P95": "neighbor_area_ratio_P95",
        "mesh_size_ratio": "mesh_size_ratio",
        "symmetric_mesh_size_mismatch": "symmetric_mesh_size_mismatch",
        "triangle_count_angle_below_5deg": "triangle_count_angle_below_5deg",
        "triangle_fraction_angle_below_5deg": "triangle_fraction_angle_below_5deg",
        "triangle_count_aspect_above_20": "triangle_count_aspect_above_20",
        "triangle_fraction_aspect_above_20": "triangle_fraction_aspect_above_20",
    }
    return [
        {
            "port_id": row["port_id"],
            "method": report["method"],
            **{output: row[source] for output, source in fields.items()},
        }
        for row in report["boundaries"]
    ]


def _tail_improvement(raw: dict[str, Any], remeshed: dict[str, Any]) -> dict[str, Any]:
    raw_angle = sum(
        int(row["triangle_count_angle_below_5deg"]) for row in raw["boundaries"]
    )
    new_angle = sum(
        int(row["triangle_count_angle_below_5deg"])
        for row in remeshed["boundaries"]
    )
    raw_aspect = sum(
        int(row["triangle_count_aspect_above_20"]) for row in raw["boundaries"]
    )
    new_aspect = sum(
        int(row["triangle_count_aspect_above_20"])
        for row in remeshed["boundaries"]
    )
    return {
        "TAIL_QUALITY_IMPROVED": new_angle < raw_angle and new_aspect < raw_aspect,
        "raw_angle_below_5_count": raw_angle,
        "entityremesh_angle_below_5_count": new_angle,
        "raw_aspect_above_20_count": raw_aspect,
        "entityremesh_aspect_above_20_count": new_aspect,
    }


def _collision_report(
    mesh: Any, intersections: list[tuple[int, int]], regions: np.ndarray
) -> dict[str, Any]:
    extension_core = 0
    extension_extension = 0
    for first, second in intersections:
        values = {int(regions[first]), int(regions[second])}
        if values == {0, 1}:
            extension_core += 1
        elif values == {1}:
            extension_extension += 1
    count = extension_core + extension_extension
    return {
        "status": "PASS" if count == 0 else "FAIL",
        "extension_collision_count": count,
        "extension_core_unintended_intersection_count": extension_core,
        "extension_extension_collision_count": extension_extension,
        "all_surface_intersection_count": len(intersections),
        "component_count": len(mesh.split(only_watertight=False)),
    }


def _legacy_entity_vmtk_surface_prepare(
    config: SurfacePrepareConfig,
    *,
    project_root: Path,
    run_id: str | None = None,
) -> VmtkSurfacePrepareResult:
    """Run exactly one official entity-aware extension-remesh candidate."""

    if config.backend.method != "vmtk_flowextensions":
        raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED")
    if config.vmtk.extension_mode != "boundarynormal":
        raise SurfacePrepareError("INVALID_VMTK_EXTENSION_MODE")
    if (
        config.vmtk.postprocess_mode != "extension_entity_remesh_then_cap"
        or not config.vmtk.remesh_after_extension
    ):
        raise SurfacePrepareError("INVALID_VMTK_POSTPROCESS_CONFIGURATION")
    inputs = load_surface_inputs(
        config.paths.cfd_preprocess_run,
        expected_boundary_count=config.qc.expected_boundary_count,
    )
    anchor = _anchor_id(str(inputs.preprocess_summary["roi_id"]))
    candidate_id = run_id or (
        f"vmtk_tps_boundarynormal_entityremesh_anchor{anchor}_"
        f"{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    layout = _layout(config.paths.output_root, candidate_id)
    tool_script = project_root / "tools" / "run_vmtk_flowextension.py"
    paths = exchange_paths(
        input_directory=layout.input,
        vmtk_directory=layout.vmtk,
        geometry_directory=layout.geometry,
        extension_mode=config.vmtk.extension_mode,
    )
    outputs: dict[str, Any] = {}
    synthetic = run_synthetic_entity_exclusion_preflight(
        config,
        root=layout.vmtk / "synthetic_exclusion_preflight",
        tool_script=tool_script,
    )
    write_json(layout.qc / "synthetic_entity_exclusion_qc.json", synthetic)
    if synthetic["status"] != "PASS":
        _write_failure(
            layout,
            status=EXCLUSION_FAILURE,
            run_id=candidate_id,
            runtime=synthetic.get("runtime", {}),
            outputs=outputs,
            evidence={"synthetic_entity_exclusion": synthetic},
        )
        raise SurfacePrepareError(EXCLUSION_FAILURE)

    references = _reference_paths(config.paths.output_root)
    immutable_paths = (
        inputs.original_surface_um_stl,
        inputs.original_surface_um_vtp,
        references["global_open_vtp"],
        references["global_raw_vtp"],
        references["global_final_vtp"],
        references["global_radius_qc"],
        references["global_summary"],
        references["caponly_raw_vtp"],
        references["caponly_final_vtp"],
        references["caponly_summary"],
    )
    hashes_before = {str(path): sha256_file(path) for path in immutable_paths}
    write_json(
        layout.input / "original_surface_reference.json",
        {
            "geometry_reference": inputs.geometry_reference,
            "immutable_hashes_before": hashes_before,
            "previous_global_remesh_reference": str(references["global_root"]),
            "previous_caponly_reference": str(references["caponly_root"]),
            "previous_reference_role": "READ_ONLY_DIAGNOSTIC_REFERENCE",
        },
    )
    shutil.copy2(config.source_path, layout.input / "cfd_surface_prepare.yaml")
    shutil.copy2(
        inputs.preprocess_run / "roi" / "boundary_conditions.json",
        layout.bc / "boundary_conditions_original.json",
    )
    original = load_original_surface(inputs.original_surface_um_vtp)
    previous_diagnostic = previous_global_remesh_diagnostics(
        original,
        references["global_final_vtp"],
        inputs.boundaries,
        config.local_cut,
    )
    write_json(layout.qc / "previous_global_remesh_core_qc.json", previous_diagnostic)

    surface = TaggedSurface.from_mesh(original)
    cut_reports: list[dict[str, Any]] = []
    for boundary in inputs.boundaries:
        surface, _, report = local_plane_cut(
            surface,
            boundary,
            radial_factor=config.local_cut.local_radial_radius_factor,
            axial_back_factor=config.local_cut.local_axial_back_radius_factor,
            axial_forward_factor=config.local_cut.local_axial_forward_radius_factor,
        )
        cut_reports.append(report)
    _open_surface(surface, paths.open_surface_vtp)
    open_data, open_mesh = polydata_mesh(paths.open_surface_vtp)
    profile_report, proximal_loops = open_profile_qc(open_mesh, inputs.boundaries)
    write_json(layout.qc / "local_cut_qc.json", {"boundaries": cut_reports})
    write_json(layout.qc / "open_surface_qc.json", profile_report)
    if not boundary_plane_alignment_pass(profile_report):
        _write_failure(
            layout,
            status=GEOMETRY_FAILURE,
            run_id=candidate_id,
            runtime=synthetic.get("runtime", {}),
            outputs=outputs,
            evidence={"open_profile": profile_report},
        )
        raise SurfacePrepareError("BOUNDARY_NORMAL_INPUT_PLANE_MISMATCH")
    if profile_report["status"] != "PASS":
        _write_failure(
            layout,
            status=BOUNDARY_FAILURE,
            run_id=candidate_id,
            runtime=synthetic.get("runtime", {}),
            outputs=outputs,
            evidence={"open_profile": profile_report},
        )
        raise SurfacePrepareError(GEOMETRY_FAILURE)
    _write_manifest(
        layout.input / "boundary_manifest.csv",
        _manifest_records(inputs.boundaries),
        profile_report,
    )
    write_json(
        layout.input / "centerline_usage.json",
        {
            "extension_mode": "boundarynormal",
            "centerline_adapter_generated": False,
            "centerlines_vtp_generated": False,
            "centerlines_used_for_extension_direction": False,
            "status": "NOT_GENERATED_BOUNDARY_NORMAL_MODE",
        },
    )
    local_targets: dict[int, float] = {}
    local_mesh_reports: list[dict[str, Any]] = []
    for boundary in inputs.boundaries:
        report = measure_local_original_mesh(
            original,
            boundary,
            sampling_radius_factor=(
                config.extension_mesh.refinement.local_mesh_sampling_radius_factor
            ),
        )
        local_mesh_reports.append(report)
        local_targets[boundary.index] = float(report["edge_length_median_um"])
    mapping = parameter_mapping(config.vmtk)
    mapping["local_original_mesh_targets"] = local_mesh_reports
    write_json(layout.vmtk / "parameters.json", mapping)
    write_json(layout.vmtk / "vmtk_parameter_mapping.json", mapping)

    invocation = run_official_vmtk(config=config.vmtk, paths=paths, tool_script=tool_script)
    runtime = {
        **invocation.runtime,
        "command": list(invocation.command),
        "configured_python_environment": config.vmtk.environment_python.parent.name,
        "vmtk_binary_runtime_prefix": str(config.vmtk.runtime_prefix),
        "runtime_overlay_is_process_local": True,
        "pmp_package_set_modified": False,
        "source_provenance": official_source_provenance(config.vmtk),
        "global_surface_remeshing_performed": False,
        "entity_aware_extension_remeshing_performed": True,
        "vmtk_surface_remeshing_called": True,
        "remesh_excluded_entity_ids": [1],
        "remesh_active_entity_ids": [2],
        "synthetic_entity_exclusion_preflight": synthetic,
    }
    outputs.update(
        {"raw_vtp": str(paths.raw_vtp.resolve()), "raw_stl": str(paths.raw_stl.resolve())}
    )
    raw_core = raw_core_exact_copy_qc(paths.open_surface_vtp, paths.raw_vtp)
    write_json(layout.qc / "raw_core_exact_copy_qc.json", raw_core)
    if raw_core["status"] != "PASS":
        _write_failure(
            layout,
            status=CORE_FAILURE,
            run_id=candidate_id,
            runtime=runtime,
            outputs=outputs,
            evidence={"raw_core_exact_copy": raw_core},
        )
        raise SurfacePrepareError(CORE_FAILURE)
    entity_assignment = assign_remesh_entities(
        paths.raw_vtp,
        entity_array_name=config.vmtk.entity_remesh.entity_array_name,
        core_entity_id=config.vmtk.entity_remesh.core_entity_id,
        extension_entity_id=config.vmtk.entity_remesh.extension_entity_id,
    )
    write_json(layout.qc / "entity_assignment_qc.json", entity_assignment)
    if entity_assignment["status"] != "PASS":
        _write_failure(
            layout,
            status=ASSIGNMENT_FAILURE,
            run_id=candidate_id,
            runtime=runtime,
            outputs=outputs,
            evidence={"entity_assignment": entity_assignment},
        )
        raise SurfacePrepareError(ASSIGNMENT_FAILURE)

    raw_data, raw_mesh = polydata_mesh(paths.raw_vtp)
    raw_topology, raw_intersections = topology_qc(
        raw_mesh, expected_open_profile_count=4, allow_degenerate=False
    )
    raw_regions = np.asarray(raw_data.cell_data["SurfaceRegionId"])
    raw_collision = _collision_report(raw_mesh, raw_intersections, raw_regions)
    raw_interface = interface_smoothness_from_raw(
        raw_mesh, input_point_count=open_data.n_points, boundaries=inputs.boundaries
    )
    raw_geometry, _ = extension_geometry_qc(
        raw_mesh, inputs.boundaries, proximal_loops, raw_interface
    )
    raw_quality = extension_mesh_quality_from_surface(
        paths.raw_vtp,
        inputs.boundaries,
        local_target_edge_um=local_targets,
        quality=config.mesh_quality,
        method="RAW_CAP_ONLY_NO_REMESH",
        require_symmetric_size_match=False,
    )
    write_json(layout.qc / "raw_surface_qc.json", raw_topology)
    write_json(layout.qc / "raw_extension_mesh_quality_qc.json", raw_quality)
    write_json(layout.qc / "raw_extension_geometry_qc.json", raw_geometry)
    write_json(layout.qc / "raw_extension_collision_qc.json", raw_collision)
    if not raw_geometry_hard_gate_pass(
        raw_core, raw_topology, raw_geometry, raw_collision, entity_assignment
    ):
        _write_failure(
            layout,
            status=GEOMETRY_FAILURE,
            run_id=candidate_id,
            runtime=runtime,
            outputs=outputs,
            evidence={
                "raw_topology": raw_topology,
                "raw_geometry": raw_geometry,
                "raw_collision": raw_collision,
            },
        )
        raise SurfacePrepareError(GEOMETRY_FAILURE)

    entity_invocation = entity_remesh_official_vmtk(
        config=config.vmtk, paths=paths, tool_script=tool_script
    )
    runtime["entity_remesh"] = entity_invocation.runtime
    runtime["entity_remesh_command"] = list(entity_invocation.command)
    write_json(layout.qc / "entity_remesh_runtime_qc.json", entity_invocation.runtime)
    core_lock, interface_lock, extension_effect = entity_remesh_preservation_qc(
        paths.raw_vtp,
        paths.remeshed_open_vtp,
        entity_array_name=config.vmtk.entity_remesh.entity_array_name,
        core_entity_id=config.vmtk.entity_remesh.core_entity_id,
        extension_entity_id=config.vmtk.entity_remesh.extension_entity_id,
    )
    write_json(layout.qc / "entity_remesh_core_exact_qc.json", core_lock)
    write_json(layout.qc / "entity_remesh_interface_lock_qc.json", interface_lock)
    write_json(layout.qc / "entity_remesh_extension_change_qc.json", extension_effect)
    if core_lock["status"] != "PASS" or interface_lock["status"] != "PASS":
        _write_failure(
            layout,
            status=CORE_FAILURE,
            run_id=candidate_id,
            runtime=runtime,
            outputs=outputs,
            evidence={"core": core_lock, "interface": interface_lock},
        )
        raise SurfacePrepareError(CORE_FAILURE)
    if extension_effect["status"] != "PASS":
        _write_failure(
            layout,
            status=NO_EFFECT_FAILURE,
            run_id=candidate_id,
            runtime=runtime,
            outputs=outputs,
            evidence={"extension_effect": extension_effect},
        )
        raise SurfacePrepareError(NO_EFFECT_FAILURE)
    outputs.update(
        {
            "remeshed_open_vtp": str(paths.remeshed_open_vtp.resolve()),
            "remeshed_open_stl": str(paths.remeshed_open_stl.resolve()),
        }
    )

    remesh_data, remesh_mesh = polydata_mesh(paths.remeshed_open_vtp)
    remesh_topology, remesh_intersections = topology_qc(
        remesh_mesh, expected_open_profile_count=4, allow_degenerate=False
    )
    remesh_profile, _ = open_profile_qc(
        remesh_mesh, inputs.boundaries, distal=True
    )
    remesh_interface = interface_smoothness_from_raw(
        remesh_mesh,
        input_point_count=0,
        boundaries=inputs.boundaries,
        core_face_mask=np.asarray(remesh_data.cell_data["SurfaceRegionId"]) == 0,
        extension_face_mask=np.asarray(remesh_data.cell_data["SurfaceRegionId"]) == 1,
    )
    remesh_geometry, distal_loops = extension_geometry_qc(
        remesh_mesh, inputs.boundaries, proximal_loops, remesh_interface
    )
    remesh_quality = extension_mesh_quality_from_surface(
        paths.remeshed_open_vtp,
        inputs.boundaries,
        local_target_edge_um=local_targets,
        quality=config.mesh_quality,
        method="NEW_ENTITY_AWARE_REMESH",
        require_symmetric_size_match=True,
    )
    remesh_regions = np.asarray(remesh_data.cell_data["SurfaceRegionId"])
    remesh_collision = _collision_report(remesh_mesh, remesh_intersections, remesh_regions)
    write_json(layout.qc / "entity_remeshed_open_surface_qc.json", remesh_topology)
    write_json(layout.qc / "entity_remeshed_open_profile_qc.json", remesh_profile)
    write_json(layout.qc / "entity_remeshed_extension_geometry_qc.json", remesh_geometry)
    write_json(
        layout.qc / "entity_remeshed_extension_mesh_quality_qc.json", remesh_quality
    )
    write_json(layout.qc / "entity_remeshed_collision_qc.json", remesh_collision)
    remesh_evidence = {
        "remeshed_open_topology": remesh_topology,
        "remeshed_open_profile": remesh_profile,
        "remeshed_geometry": remesh_geometry,
        "remeshed_quality": remesh_quality,
        "remeshed_collision": remesh_collision,
    }
    remesh_failure = None
    if remesh_topology["status"] != "PASS":
        remesh_failure = TOPOLOGY_FAILURE
    elif remesh_profile["status"] != "PASS":
        remesh_failure = BOUNDARY_FAILURE
    elif not raw_geometry_hard_gate_pass(
        remesh_geometry, remesh_quality, remesh_collision
    ):
        remesh_failure = GEOMETRY_FAILURE
    if remesh_failure is not None:
        _write_failure(
            layout,
            status=remesh_failure,
            run_id=candidate_id,
            runtime=runtime,
            outputs=outputs,
            evidence=remesh_evidence,
        )
        raise SurfacePrepareError(remesh_failure)

    previous_quality = extension_mesh_quality_from_surface(
        references["global_final_vtp"],
        inputs.boundaries,
        local_target_edge_um=local_targets,
        quality=config.mesh_quality,
        method="PREVIOUS_GLOBAL_REMESH",
        require_symmetric_size_match=False,
        geometric_extension_selection=True,
    )
    write_csv(
        layout.qc / "extension_quality_three_way.csv",
        [
            *_flatten_quality(raw_quality),
            *_flatten_quality(previous_quality),
            *_flatten_quality(remesh_quality),
        ],
    )
    write_csv(
        layout.qc / "raw_vs_entityremesh_extension_quality.csv",
        [*_flatten_quality(raw_quality), *_flatten_quality(remesh_quality)],
    )
    tail_improvement = _tail_improvement(raw_quality, remesh_quality)
    write_json(layout.qc / "extension_tail_improvement_qc.json", tail_improvement)

    cap = cap_official_vmtk(config=config.vmtk, paths=paths, tool_script=tool_script)
    runtime["cap_only"] = cap.runtime
    runtime["cap_only_command"] = list(cap.command)
    if cap.runtime.get("surface_remesher_called") is not False:
        raise SurfacePrepareError(GEOMETRY_FAILURE)
    write_json(layout.vmtk / "environment.json", runtime)
    final_outputs, boundary_mapping = tag_and_export_final_surface(
        paths.capped_vtp,
        inputs.boundaries,
        layout.geometry,
        layout.boundaries,
        output_stem="cfd_surface_vmtk_tps_boundarynormal_entityremesh",
        raw_vtp=paths.remeshed_open_vtp,
    )
    outputs.update(final_outputs)
    write_json(layout.qc / "boundary_mapping_qc.json", boundary_mapping)
    _, final_mesh = polydata_mesh(Path(outputs["tagged_vtp"]))
    final_topology, _ = topology_qc(
        final_mesh, expected_open_profile_count=0, require_winding_consistent=True
    )
    core_exact = core_exact_preservation_qc(
        original, Path(outputs["tagged_vtp"]), inputs.boundaries, config.local_cut
    )
    core_distance = core_symmetric_distance_qc(
        original, Path(outputs["tagged_vtp"]), inputs.boundaries, config.local_cut
    )
    normal_report = normal_consistency_qc(
        original, Path(outputs["tagged_vtp"]), inputs.boundaries, config.local_cut
    )
    try:
        radius_report = _radius_qc(
            final_mesh=final_mesh, inputs=inputs, project_root=project_root
        )
    except SurfacePrepareError as error:
        radius_report = {"status": "FAIL", "error": str(error)}
    meter_report = meter_scale_qc(
        Path(outputs["manual_review_stl"]), Path(outputs["meter_stl"])
    )
    write_json(layout.qc / "cap_surface_qc.json", final_topology)
    write_json(layout.qc / "core_exact_preservation_qc.json", core_exact)
    write_json(layout.qc / "core_symmetric_distance_qc.json", core_distance)
    write_json(layout.qc / "normal_consistency_qc.json", normal_report)
    write_json(layout.qc / "radius_fidelity.json", radius_report)
    write_json(layout.qc / "meter_scale_qc.json", meter_report)

    pressure_rows: list[dict[str, Any]] = []
    try:
        pressure_rows, pressure_report = geometry_pressure_correction(
            final_mesh,
            inputs.boundaries,
            proximal_loops,
            distal_loops,
            dynamic_viscosity_pa_s=float(
                inputs.original_boundary_conditions["fluid"]["dynamic_viscosity_pa_s"]
            ),
        )
    except SurfacePrepareError as error:
        pressure_report = {"status": "FAIL", "error": str(error)}
    else:
        pressure_report["status"] = "PASS"
        csv_rows = [
            {
                **row,
                "station_fractions": json.dumps(row["station_fractions"]),
                "cross_section_area_um2": json.dumps(row["cross_section_area_um2"]),
                "equivalent_radius_um": json.dumps(row["equivalent_radius_um"]),
            }
            for row in pressure_rows
        ]
        pressure_csv = (
            layout.bc / "extension_pressure_correction_vmtk_boundarynormal_entityremesh.csv"
        )
        write_csv(pressure_csv, csv_rows)
        boundary_json = (
            layout.bc / "boundary_conditions_vmtk_boundarynormal_entityremesh.json"
        )
        write_json(
            boundary_json,
            _boundary_conditions(inputs.original_boundary_conditions, pressure_rows),
        )
        outputs["pressure_correction_csv"] = str(pressure_csv.resolve())
        outputs["boundary_conditions_json"] = str(boundary_json.resolve())
    write_json(layout.qc / "extension_pressure_geometry_qc.json", pressure_report)

    previous_radius = read_json(references["global_radius_qc"])
    core_rows = [
        {
            "method": "ORIGINAL_ULTRALISER",
            "core_original_to_final_P95_um": 0.0,
            "core_original_to_final_max_um": 0.0,
            "core_final_to_original_P95_um": 0.0,
            "core_final_to_original_max_um": 0.0,
            "core_normal_P95_deg": 0.0,
            "core_normal_P99_deg": 0.0,
            "core_normal_max_deg": 0.0,
            "radius_P95_error": None,
        },
        {
            "method": "PREVIOUS_GLOBAL_REMESH",
            "core_original_to_final_P95_um": previous_diagnostic[
                "original_core_to_previous_final"
            ]["P95_um"],
            "core_original_to_final_max_um": previous_diagnostic[
                "original_core_to_previous_final"
            ]["max_um"],
            "core_final_to_original_P95_um": previous_diagnostic[
                "previous_final_core_to_original"
            ]["P95_um"],
            "core_final_to_original_max_um": previous_diagnostic[
                "previous_final_core_to_original"
            ]["max_um"],
            "core_normal_P95_deg": previous_diagnostic["core_normal_deviation"]["P95_deg"],
            "core_normal_P99_deg": previous_diagnostic["core_normal_deviation"]["P99_deg"],
            "core_normal_max_deg": previous_diagnostic["core_normal_deviation"]["max_deg"],
            "radius_P95_error": previous_radius["p95_absolute_relative_error"],
        },
        {
            "method": "NEW_ENTITY_AWARE_REMESH",
            "core_original_to_final_P95_um": core_distance[
                "original_core_to_entityremesh_final_core"
            ]["P95_um"],
            "core_original_to_final_max_um": core_distance[
                "original_core_to_entityremesh_final_core"
            ]["max_um"],
            "core_final_to_original_P95_um": core_distance[
                "entityremesh_final_core_to_original_core"
            ]["P95_um"],
            "core_final_to_original_max_um": core_distance[
                "entityremesh_final_core_to_original_core"
            ]["max_um"],
            "core_normal_P95_deg": normal_report["core_normal_deviation_P95_deg"],
            "core_normal_P99_deg": normal_report["core_normal_deviation_P99_deg"],
            "core_normal_max_deg": normal_report["core_normal_deviation_max_deg"],
            "radius_P95_error": radius_report.get("p95_absolute_relative_error"),
        },
    ]
    write_csv(layout.qc / "core_three_way_comparison.csv", core_rows)
    global_core_modified = bool(
        previous_diagnostic["global_remesh_artifact_detected"]
        and core_lock["status"] == "PASS"
        and core_distance["status"] == "PASS"
    )
    figure_paths = save_entityremesh_review_figures(
        original_vtp=inputs.original_surface_um_vtp,
        raw_vtp=paths.raw_vtp,
        previous_global_final_vtp=references["global_final_vtp"],
        entity_remeshed_open_vtp=paths.remeshed_open_vtp,
        entity_final_vtp=Path(outputs["tagged_vtp"]),
        boundaries=inputs.boundaries,
        output_directory=layout.figures,
    )
    outputs["figure_paths"] = [str(path) for path in figure_paths]
    hashes_after = {str(path): sha256_file(path) for path in immutable_paths}
    integrity = {
        "status": "PASS" if hashes_before == hashes_after else "FAIL",
        "hashes_before": hashes_before,
        "hashes_after": hashes_after,
        "all_read_only_references_untouched": hashes_before == hashes_after,
    }
    write_json(layout.qc / "original_surface_integrity.json", integrity)
    base_status = resolve_entity_failure_status(
        topology=final_topology,
        boundary_mapping=boundary_mapping,
        core_exact=core_exact,
        core_distance=core_distance,
        normal=normal_report,
        radius=radius_report,
        meter=meter_report,
        pressure=pressure_report,
        integrity=integrity,
    )
    if base_status == SUCCESS_STATUS:
        status = (
            TAIL_WARNING_STATUS
            if remesh_quality["tail_status"] == "MESH_TAIL_WARNING"
            else SUCCESS_STATUS
        )
        next_stage = SUCCESS_NEXT
    else:
        status = base_status
        next_stage = FAIL_NEXT
    summary = {
        "status": status,
        "next": next_stage,
        "run_id": candidate_id,
        "run_root": str(layout.root),
        "runtime": runtime,
        "synthetic_entity_exclusion": synthetic,
        "raw_core_exact_copy": raw_core,
        "entity_assignment": entity_assignment,
        "entity_remesh_core_exact": core_lock,
        "entity_remesh_interface_lock": interface_lock,
        "entity_remesh_extension_change": extension_effect,
        "raw_topology": raw_topology,
        "raw_extension_geometry": raw_geometry,
        "raw_extension_quality": raw_quality,
        "entity_remeshed_open_topology": remesh_topology,
        "entity_remeshed_extension_geometry": remesh_geometry,
        "entity_remeshed_extension_quality": remesh_quality,
        "previous_global_remesh_extension_quality": previous_quality,
        "extension_tail_improvement": tail_improvement,
        "boundaries": remesh_geometry["boundaries"],
        "final_topology": final_topology,
        "boundary_mapping": boundary_mapping,
        "core_exact_preservation": core_exact,
        "core_symmetric_distance": core_distance,
        "normal_consistency": normal_report,
        "radius_fidelity": radius_report,
        "meter_scale": meter_report,
        "pressure_geometry": pressure_report,
        "collision": remesh_collision,
        "GLOBAL_REMESH_CORE_MODIFICATION_CONFIRMED": global_core_modified,
        "outputs": outputs,
        "figures": [str(path) for path in figure_paths],
        "original_integrity": integrity,
        "single_primary_candidate_count": 1,
        "second_candidate_run": False,
        "automatic_fallback": False,
        "global_surface_remeshing_performed": False,
        "entity_aware_extension_remeshing_performed": True,
        "vmtk_surface_remeshing_called": True,
        "visible_acceptance": "MANUAL_REVIEW_REQUIRED",
        "volume_mesh_created": False,
        "cfd_run": False,
        "microbubble_simulation_run": False,
    }
    write_json(layout.qc / "run_summary.json", summary)
    _write_report(
        layout.report / "vmtk_tps_boundarynormal_entityremesh_report.md", summary
    )
    return VmtkSurfacePrepareResult(
        status=status,
        next_stage=next_stage,
        run_root=layout.root,
        runtime=runtime,
        synthetic_preflight=synthetic,
        entity_core_qc=core_lock,
        entity_interface_qc=interface_lock,
        entity_extension_qc=extension_effect,
        boundaries=tuple(remesh_geometry["boundaries"]),
        raw_quality=raw_quality,
        remeshed_quality=remesh_quality,
        previous_quality=previous_quality,
        final_topology=final_topology,
        core_distance_qc=core_distance,
        normal_qc=normal_report,
        radius_qc=radius_report,
        meter_qc=meter_report,
        collision_count=int(remesh_collision["extension_collision_count"]),
        pressure_rows=tuple(pressure_rows),
        output_paths=outputs,
        figures=figure_paths,
    )


def _use_true_intersection_count(
    topology: dict[str, Any], intersection: dict[str, Any]
) -> dict[str, Any]:
    topology = copy.deepcopy(topology)
    topology["vtk_candidate_self_intersection_count"] = topology[
        "self_intersection_count"
    ]
    topology["self_intersection_count"] = intersection[
        "true_self_intersection_count"
    ]
    topology["checks"]["zero_self_intersections"] = (
        intersection["true_self_intersection_count"] == 0
    )
    topology["status"] = (
        "PASS" if all(topology["checks"].values()) else "FAIL"
    )
    return topology


def should_cap_guarded_open(
    topology: dict[str, Any], intersection: dict[str, Any]
) -> bool:
    return topology.get("status") == "PASS" and intersection.get("status") == "PASS"


def should_write_guarded_collision_figure(topology_status: str) -> bool:
    del topology_status
    return True


def _historical_guarded_vmtk_surface_prepare(
    config: SurfacePrepareConfig,
    *,
    project_root: Path,
    run_id: str | None = None,
) -> VmtkSurfacePrepareResult:
    """Retain the completed guarded workflow only for historical result reading."""

    settings = config.vmtk.entity_remesh
    valid = (
        config.backend.method == "vmtk_flowextensions"
        and config.vmtk.extension_mode == "boundarynormal"
        and config.vmtk.postprocess_mode
        == "guarded_extension_entity_remesh_then_cap"
        and config.vmtk.remesh_after_extension
        and settings.core_entity_id == 1
        and settings.guard_entity_id == 2
        and settings.extension_body_entity_id == 3
        and settings.exclude_entity_ids == (1, 2)
        and settings.guard.face_layers == 2
        and settings.target_edge_length_um == 0.25913916380971913
    )
    if not valid:
        raise SurfacePrepareError("INVALID_VMTK_POSTPROCESS_CONFIGURATION")
    inputs = load_surface_inputs(
        config.paths.cfd_preprocess_run,
        expected_boundary_count=config.qc.expected_boundary_count,
    )
    anchor = _anchor_id(str(inputs.preprocess_summary["roi_id"]))
    candidate_id = run_id or (
        f"vmtk_tps_boundarynormal_guarded_entityremesh_anchor{anchor}_"
        f"{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    layout = _layout(config.paths.output_root, candidate_id)
    tool_script = project_root / "tools" / "run_vmtk_flowextension.py"
    paths = exchange_paths(
        input_directory=layout.input,
        vmtk_directory=layout.vmtk,
        geometry_directory=layout.geometry,
        extension_mode="boundarynormal",
    )
    references = _reference_paths(config.paths.output_root)
    outputs: dict[str, Any] = {}

    previous_diagnosis, previous_layer = diagnose_previous_entityremesh(
        references["previous_entity_raw_vtp"],
        references["previous_entity_remeshed_vtp"],
        inputs.boundaries,
    )
    write_json(
        layout.qc / "previous_entityremesh_intersection_diagnosis.json",
        previous_diagnosis,
    )
    write_json(layout.qc / "previous_collision_layer_diagnosis.json", previous_layer)
    cut002 = next(
        boundary for boundary in inputs.boundaries if boundary.port_id.endswith("cut_002")
    )
    previous_figure = save_previous_collision_closeup(
        remeshed_vtp=references["previous_entity_remeshed_vtp"],
        diagnosis=previous_diagnosis,
        boundary=cut002,
        output_path=layout.figures / "previous_entityremesh_collision_closeup.png",
    )
    outputs["previous_collision_figure"] = str(previous_figure)
    if previous_diagnosis.get("REAL_GEOMETRIC_PENETRATION") != "YES":
        status = "ENTITY_REMESH_INTERSECTION_DETECTOR_REVIEW_REQUIRED"
        _write_failure(
            layout,
            status=status,
            run_id=candidate_id,
            runtime={},
            outputs=outputs,
            evidence={"previous_intersection": previous_diagnosis},
        )
        raise SurfacePrepareError(status)

    synthetic = run_synthetic_entity_exclusion_preflight(
        config,
        root=layout.vmtk / "synthetic_exclusion_preflight",
        tool_script=tool_script,
    )
    write_json(layout.qc / "synthetic_guarded_entity_exclusion_qc.json", synthetic)
    if synthetic["status"] != "PASS":
        _write_failure(
            layout,
            status=EXCLUSION_FAILURE,
            run_id=candidate_id,
            runtime=synthetic.get("runtime", {}),
            outputs=outputs,
            evidence={"synthetic_guarded_exclusion": synthetic},
        )
        raise SurfacePrepareError(EXCLUSION_FAILURE)

    immutable_paths = (
        inputs.original_surface_um_stl,
        inputs.original_surface_um_vtp,
        references["global_open_vtp"],
        references["global_raw_vtp"],
        references["global_final_vtp"],
        references["global_radius_qc"],
        references["global_summary"],
        references["caponly_raw_vtp"],
        references["caponly_final_vtp"],
        references["caponly_summary"],
        references["previous_entity_raw_vtp"],
        references["previous_entity_remeshed_vtp"],
        references["previous_entity_summary"],
    )
    hashes_before = {str(path): sha256_file(path) for path in immutable_paths}
    write_json(
        layout.input / "original_surface_reference.json",
        {
            "geometry_reference": inputs.geometry_reference,
            "immutable_hashes_before": hashes_before,
            "previous_global_remesh_reference": str(references["global_root"]),
            "previous_caponly_reference": str(references["caponly_root"]),
            "previous_entityremesh_reference": str(
                references["previous_entity_root"]
            ),
            "previous_reference_role": "READ_ONLY_DIAGNOSTIC_REFERENCE",
        },
    )
    shutil.copy2(config.source_path, layout.input / "cfd_surface_prepare.yaml")
    shutil.copy2(
        inputs.preprocess_run / "roi" / "boundary_conditions.json",
        layout.bc / "boundary_conditions_original.json",
    )
    original = load_original_surface(inputs.original_surface_um_vtp)
    surface = TaggedSurface.from_mesh(original)
    cut_reports: list[dict[str, Any]] = []
    for boundary in inputs.boundaries:
        surface, _, cut_report = local_plane_cut(
            surface,
            boundary,
            radial_factor=config.local_cut.local_radial_radius_factor,
            axial_back_factor=config.local_cut.local_axial_back_radius_factor,
            axial_forward_factor=config.local_cut.local_axial_forward_radius_factor,
        )
        cut_reports.append(cut_report)
    _open_surface(surface, paths.open_surface_vtp)
    open_data, open_mesh = polydata_mesh(paths.open_surface_vtp)
    open_profile, proximal_loops = open_profile_qc(open_mesh, inputs.boundaries)
    write_json(layout.qc / "local_cut_qc.json", {"boundaries": cut_reports})
    write_json(layout.qc / "open_surface_qc.json", open_profile)
    if open_profile["status"] != "PASS" or not boundary_plane_alignment_pass(
        open_profile
    ):
        raise SurfacePrepareError(GEOMETRY_FAILURE)
    _write_manifest(
        layout.input / "boundary_manifest.csv",
        _manifest_records(inputs.boundaries),
        open_profile,
    )
    write_json(
        layout.input / "centerline_usage.json",
        {
            "extension_mode": "boundarynormal",
            "centerline_adapter_generated": False,
            "centerlines_used_for_extension_direction": False,
        },
    )
    local_targets: dict[int, float] = {}
    local_mesh_reports: list[dict[str, Any]] = []
    for boundary in inputs.boundaries:
        local_report = measure_local_original_mesh(
            original,
            boundary,
            sampling_radius_factor=(
                config.extension_mesh.refinement.local_mesh_sampling_radius_factor
            ),
        )
        local_mesh_reports.append(local_report)
        local_targets[boundary.index] = float(local_report["edge_length_median_um"])
    mapping = parameter_mapping(config.vmtk)
    mapping["local_original_mesh_targets"] = local_mesh_reports
    write_json(layout.vmtk / "parameters.json", mapping)
    write_json(layout.vmtk / "vmtk_parameter_mapping.json", mapping)

    invocation = run_official_vmtk(
        config=config.vmtk, paths=paths, tool_script=tool_script
    )
    runtime = {
        **invocation.runtime,
        "command": list(invocation.command),
        "configured_python_environment": "pmp",
        "source_provenance": official_source_provenance(config.vmtk),
        "global_surface_remeshing_performed": False,
        "remesh_excluded_entity_ids": [1, 2],
        "remesh_active_entity_ids": [3],
        "synthetic_guarded_entity_exclusion_preflight": synthetic,
    }
    outputs.update(
        {"raw_vtp": str(paths.raw_vtp.resolve()), "raw_stl": str(paths.raw_stl.resolve())}
    )
    raw_core = raw_core_exact_copy_qc(paths.open_surface_vtp, paths.raw_vtp)
    write_json(layout.qc / "raw_core_exact_copy_qc.json", raw_core)
    if raw_core["status"] != "PASS":
        raise SurfacePrepareError(CORE_FAILURE)
    assignment, raw_layers = assign_guarded_remesh_entities(
        paths.raw_vtp,
        inputs.boundaries,
        face_layers=settings.guard.face_layers,
        entity_array_name=settings.entity_array_name,
        core_entity_id=settings.core_entity_id,
        guard_entity_id=settings.guard_entity_id,
        body_entity_id=settings.extension_body_entity_id,
    )
    write_json(layout.qc / "guard_entity_assignment_qc.json", assignment)
    if assignment["status"] == SUSPICIOUS_GUARD_FAILURE:
        raise SurfacePrepareError(SUSPICIOUS_GUARD_FAILURE)
    if assignment["status"] != "PASS":
        raise SurfacePrepareError(ASSIGNMENT_FAILURE)

    raw_data, raw_mesh = polydata_mesh(paths.raw_vtp)
    raw_topology, raw_pairs = topology_qc(
        raw_mesh, expected_open_profile_count=4, allow_degenerate=False
    )
    raw_regions = np.asarray(raw_data.cell_data["SurfaceRegionId"])
    raw_collision = _collision_report(raw_mesh, raw_pairs, raw_regions)
    raw_interface = interface_smoothness_from_raw(
        raw_mesh,
        input_point_count=open_data.n_points,
        boundaries=inputs.boundaries,
    )
    raw_geometry, _ = extension_geometry_qc(
        raw_mesh, inputs.boundaries, proximal_loops, raw_interface
    )
    write_json(layout.qc / "raw_surface_qc.json", raw_topology)
    write_json(layout.qc / "raw_extension_geometry_qc.json", raw_geometry)
    write_json(layout.qc / "raw_extension_collision_qc.json", raw_collision)
    if not raw_geometry_hard_gate_pass(
        raw_topology, raw_geometry, raw_collision
    ):
        raise SurfacePrepareError(GEOMETRY_FAILURE)

    entity_invocation = entity_remesh_official_vmtk(
        config=config.vmtk, paths=paths, tool_script=tool_script
    )
    runtime["guarded_entity_remesh"] = entity_invocation.runtime
    runtime["guarded_entity_remesh_command"] = list(entity_invocation.command)
    write_json(layout.qc / "entity_remesh_runtime_qc.json", entity_invocation.runtime)
    core_lock, guard_lock, boundary_lock, body_effect = (
        guarded_entity_preservation_qc(
            paths.raw_vtp,
            paths.remeshed_open_vtp,
            entity_array_name=settings.entity_array_name,
            core_entity_id=settings.core_entity_id,
            guard_entity_id=settings.guard_entity_id,
            body_entity_id=settings.extension_body_entity_id,
        )
    )
    write_json(layout.qc / "entity_remesh_core_exact_qc.json", core_lock)
    write_json(layout.qc / "guard_exact_preservation_qc.json", guard_lock)
    write_json(layout.qc / "entity_boundary_lock_qc.json", boundary_lock)
    write_json(layout.qc / "entity_remesh_body_change_qc.json", body_effect)
    if core_lock["status"] != "PASS":
        raise SurfacePrepareError(CORE_FAILURE)
    if guard_lock["status"] != "PASS" or boundary_lock["status"] != "PASS":
        raise SurfacePrepareError(GUARD_FAILURE)
    if body_effect["status"] != "PASS":
        raise SurfacePrepareError(NO_EFFECT_FAILURE)
    outputs.update(
        {
            "remeshed_open_vtp": str(paths.remeshed_open_vtp.resolve()),
            "remeshed_open_stl": str(paths.remeshed_open_stl.resolve()),
        }
    )

    remesh_data, remesh_mesh = polydata_mesh(paths.remeshed_open_vtp)
    remesh_topology_raw, _ = topology_qc(
        remesh_mesh, expected_open_profile_count=4, allow_degenerate=False
    )
    intersection_qc, intersection_records = guarded_intersection_qc(
        paths.remeshed_open_vtp
    )
    remesh_topology = _use_true_intersection_count(
        remesh_topology_raw, intersection_qc
    )
    remesh_profile, _ = open_profile_qc(
        remesh_mesh, inputs.boundaries, distal=True
    )
    write_json(
        layout.qc / "guarded_entityremesh_intersection_qc.json", intersection_qc
    )
    write_json(layout.qc / "guarded_open_surface_qc.json", remesh_topology)
    write_json(layout.qc / "guarded_open_profile_qc.json", remesh_profile)
    if not should_cap_guarded_open(remesh_topology, intersection_qc):
        figure_paths = save_guarded_open_figures(
            raw_vtp=paths.raw_vtp,
            remeshed_vtp=paths.remeshed_open_vtp,
            boundaries=inputs.boundaries,
            tail_records=[],
            intersections=intersection_records,
            output_directory=layout.figures,
        )
        outputs["figure_paths"] = [str(path) for path in figure_paths]
        _write_failure(
            layout,
            status=TOPOLOGY_FAILURE,
            run_id=candidate_id,
            runtime=runtime,
            outputs=outputs,
            evidence={
                "topology": remesh_topology,
                "intersections": intersection_qc,
                "core": core_lock,
                "guard": guard_lock,
                "boundaries": boundary_lock,
                "body": body_effect,
            },
        )
        raise SurfacePrepareError(TOPOLOGY_FAILURE)
    if remesh_profile["status"] != "PASS":
        raise SurfacePrepareError(BOUNDARY_FAILURE)

    remesh_entities = np.asarray(remesh_data.cell_data["RemeshEntityId"])
    remesh_interface = interface_smoothness_from_raw(
        remesh_mesh,
        input_point_count=None,
        boundaries=inputs.boundaries,
        core_face_mask=remesh_entities == 1,
        extension_face_mask=np.isin(remesh_entities, (2, 3)),
    )
    remesh_geometry, distal_loops = extension_geometry_qc(
        remesh_mesh, inputs.boundaries, proximal_loops, remesh_interface
    )
    write_json(layout.qc / "extension_geometry_qc.json", remesh_geometry)
    if remesh_geometry["status"] != "PASS":
        raise SurfacePrepareError(GEOMETRY_FAILURE)
    guard_quality, guard_tails = guarded_region_mesh_quality(
        paths.remeshed_open_vtp,
        inputs.boundaries,
        entity_id=2,
        entity_label="PROXIMAL_GUARD",
        local_target_edge_um=local_targets,
        quality=config.mesh_quality,
        hard_body_gate=False,
    )
    body_quality, body_tails = guarded_region_mesh_quality(
        paths.remeshed_open_vtp,
        inputs.boundaries,
        entity_id=3,
        entity_label="EXTENSION_BODY",
        local_target_edge_um=local_targets,
        quality=config.mesh_quality,
        hard_body_gate=True,
    )
    tail_records = [*guard_tails, *body_tails]
    write_json(layout.qc / "guard_mesh_quality_qc.json", guard_quality)
    write_json(layout.qc / "body_mesh_quality_qc.json", body_quality)
    tail_path = layout.qc / "mesh_tail_locations.csv"
    if tail_records:
        write_csv(tail_path, tail_records)
    else:
        tail_path.write_text(
            "port_id,entity,triangle_id,minimum_angle_deg,aspect_ratio,"
            "distance_to_original_interface_um,axial_distance_um,"
            "axial_distance_in_D\n",
            encoding="utf-8",
        )
    previous_tail = previous_tail_fraction_inside_guard(
        references["previous_entity_remeshed_vtp"], paths.raw_vtp
    )
    write_json(layout.qc / "previous_tail_guard_localization_qc.json", previous_tail)
    if body_quality["status"] != "PASS":
        raise SurfacePrepareError(GEOMETRY_FAILURE)
    figure_paths = save_guarded_open_figures(
        raw_vtp=paths.raw_vtp,
        remeshed_vtp=paths.remeshed_open_vtp,
        boundaries=inputs.boundaries,
        tail_records=tail_records,
        intersections=intersection_records,
        output_directory=layout.figures,
    )
    outputs["figure_paths"] = [str(path) for path in figure_paths]

    cap = cap_official_vmtk(
        config=config.vmtk, paths=paths, tool_script=tool_script
    )
    runtime["cap_only"] = cap.runtime
    runtime["cap_only_command"] = list(cap.command)
    if cap.runtime.get("surface_remesher_called") is not False:
        raise SurfacePrepareError(GEOMETRY_FAILURE)
    write_json(layout.vmtk / "environment.json", runtime)
    final_outputs, boundary_mapping = tag_and_export_final_surface(
        paths.capped_vtp,
        inputs.boundaries,
        layout.geometry,
        layout.boundaries,
        output_stem="cfd_surface_vmtk_tps_boundarynormal_guarded_entityremesh",
        raw_vtp=paths.remeshed_open_vtp,
    )
    outputs.update(final_outputs)
    write_json(layout.qc / "boundary_mapping_qc.json", boundary_mapping)
    _, final_mesh = polydata_mesh(Path(outputs["tagged_vtp"]))
    final_topology_raw, _ = topology_qc(
        final_mesh, expected_open_profile_count=0, require_winding_consistent=True
    )
    final_intersection, _ = guarded_intersection_qc(Path(outputs["tagged_vtp"]))
    final_topology = _use_true_intersection_count(
        final_topology_raw, final_intersection
    )
    core_exact = core_exact_preservation_qc(
        original, Path(outputs["tagged_vtp"]), inputs.boundaries, config.local_cut
    )
    core_distance = core_symmetric_distance_qc(
        original, Path(outputs["tagged_vtp"]), inputs.boundaries, config.local_cut
    )
    normal_report = normal_consistency_qc(
        original, Path(outputs["tagged_vtp"]), inputs.boundaries, config.local_cut
    )
    try:
        radius_report = _radius_qc(
            final_mesh=final_mesh, inputs=inputs, project_root=project_root
        )
    except SurfacePrepareError as error:
        radius_report = {"status": "FAIL", "error": str(error)}
    meter_report = meter_scale_qc(
        Path(outputs["manual_review_stl"]), Path(outputs["meter_stl"])
    )
    write_json(layout.qc / "cap_surface_qc.json", final_topology)
    write_json(layout.qc / "final_intersection_qc.json", final_intersection)
    write_json(layout.qc / "core_exact_preservation_qc.json", core_exact)
    write_json(layout.qc / "core_symmetric_distance_qc.json", core_distance)
    write_json(layout.qc / "normal_consistency_qc.json", normal_report)
    write_json(layout.qc / "radius_fidelity.json", radius_report)
    write_json(layout.qc / "meter_scale_qc.json", meter_report)

    pressure_rows: list[dict[str, Any]] = []
    try:
        pressure_rows, pressure_report = geometry_pressure_correction(
            final_mesh,
            inputs.boundaries,
            proximal_loops,
            distal_loops,
            dynamic_viscosity_pa_s=float(
                inputs.original_boundary_conditions["fluid"]["dynamic_viscosity_pa_s"]
            ),
        )
    except SurfacePrepareError as error:
        pressure_report = {"status": "FAIL", "error": str(error)}
    else:
        pressure_report["status"] = "PASS"
        pressure_csv = (
            layout.bc
            / "extension_pressure_correction_vmtk_boundarynormal_guarded_entityremesh.csv"
        )
        write_csv(
            pressure_csv,
            [
                {
                    **row,
                    "station_fractions": json.dumps(row["station_fractions"]),
                    "cross_section_area_um2": json.dumps(
                        row["cross_section_area_um2"]
                    ),
                    "equivalent_radius_um": json.dumps(row["equivalent_radius_um"]),
                }
                for row in pressure_rows
            ],
        )
        boundary_json = (
            layout.bc
            / "boundary_conditions_vmtk_boundarynormal_guarded_entityremesh.json"
        )
        write_json(
            boundary_json,
            _boundary_conditions(inputs.original_boundary_conditions, pressure_rows),
        )
        outputs["pressure_correction_csv"] = str(pressure_csv.resolve())
        outputs["boundary_conditions_json"] = str(boundary_json.resolve())
    write_json(layout.qc / "extension_pressure_geometry_qc.json", pressure_report)

    comparison = save_guarded_three_way_comparison(
        caponly_vtp=references["caponly_final_vtp"],
        previous_entityremesh_vtp=references["previous_entity_remeshed_vtp"],
        guarded_final_vtp=Path(outputs["tagged_vtp"]),
        boundaries=inputs.boundaries,
        output_path=layout.figures / "extension_three_way_comparison.png",
    )
    outputs["figure_paths"].append(str(comparison))
    hashes_after = {str(path): sha256_file(path) for path in immutable_paths}
    integrity = {
        "status": "PASS" if hashes_before == hashes_after else "FAIL",
        "hashes_before": hashes_before,
        "hashes_after": hashes_after,
        "all_read_only_references_untouched": hashes_before == hashes_after,
    }
    write_json(layout.qc / "original_surface_integrity.json", integrity)
    base_status = resolve_entity_failure_status(
        topology=final_topology,
        boundary_mapping=boundary_mapping,
        core_exact=core_exact,
        core_distance=core_distance,
        normal=normal_report,
        radius=radius_report,
        meter=meter_report,
        pressure=pressure_report,
        integrity=integrity,
    )
    tail_warning = bool(guard_tails or body_tails)
    status = (
        TAIL_WARNING_STATUS
        if base_status == SUCCESS_STATUS and tail_warning
        else base_status
    )
    next_stage = SUCCESS_NEXT if status in PASS_STATUSES else FAIL_NEXT
    summary = {
        "status": status,
        "next": next_stage,
        "run_id": candidate_id,
        "run_root": str(layout.root),
        "runtime": runtime,
        "previous_intersection_diagnosis": previous_diagnosis,
        "previous_collision_layer": previous_layer,
        "synthetic_guarded_exclusion": synthetic,
        "guard_entity_assignment": assignment,
        "raw_guard_layers": raw_layers.tolist(),
        "core_exact_lock": core_lock,
        "guard_exact_lock": guard_lock,
        "entity_boundary_lock": boundary_lock,
        "body_remesh_effect": body_effect,
        "guarded_open_intersections": intersection_qc,
        "guarded_open_topology": remesh_topology,
        "boundaries": remesh_geometry["boundaries"],
        "guard_mesh_quality": guard_quality,
        "body_mesh_quality": body_quality,
        "mesh_tail_location_count": len(tail_records),
        "previous_tail_localization": previous_tail,
        "final_topology": final_topology,
        "boundary_mapping": boundary_mapping,
        "core_exact_preservation": core_exact,
        "core_symmetric_distance": core_distance,
        "normal_consistency": normal_report,
        "radius_fidelity": radius_report,
        "pressure_geometry": pressure_report,
        "meter_scale": meter_report,
        "outputs": outputs,
        "original_integrity": integrity,
        "single_primary_candidate_count": 1,
        "second_candidate_run": False,
        "automatic_fallback": False,
        "global_surface_remeshing_performed": False,
        "volume_mesh_created": False,
        "cfd_run": False,
        "microbubble_simulation_run": False,
    }
    write_json(layout.qc / "run_summary.json", summary)
    _write_report(
        layout.report / "vmtk_tps_boundarynormal_guarded_entityremesh_report.md",
        summary,
    )
    return VmtkSurfacePrepareResult(
        status=status,
        next_stage=next_stage,
        run_root=layout.root,
        runtime=runtime,
        synthetic_preflight=synthetic,
        guard_assignment_qc=assignment,
        core_qc=core_lock,
        guard_qc=guard_lock,
        boundary_lock_qc=boundary_lock,
        body_qc=body_effect,
        intersection_qc=intersection_qc,
        boundaries=tuple(remesh_geometry["boundaries"]),
        guard_quality=guard_quality,
        body_quality=body_quality,
        previous_tail=previous_tail,
        final_topology=final_topology,
        core_distance_qc=core_distance,
        normal_qc=normal_report,
        radius_qc=radius_report,
        meter_qc=meter_report,
        collision_count=int(intersection_qc["true_self_intersection_count"]),
        pressure_rows=tuple(pressure_rows),
        output_paths=outputs,
        figures=tuple(Path(path) for path in outputs["figure_paths"]),
    )


def should_cap_crossseam_open(
    topology: dict[str, Any],
    intersection: dict[str, Any],
    seam_topology: dict[str, Any],
) -> bool:
    """Allow capping only after all cross-seam open-surface hard gates pass."""

    return (
        topology.get("status") == "PASS"
        and intersection.get("status") == "PASS"
        and seam_topology.get("status") == "PASS"
    )


def _seam_comparison_rows(
    global_report: dict[str, Any],
    guarded_report: dict[str, Any],
    crossseam_report: dict[str, Any],
    seam_topology: dict[str, Any],
) -> list[dict[str, Any]]:
    reports = {
        "global": {row["port_id"]: row for row in global_report["boundaries"]},
        "guarded": {
            row["port_id"]: row for row in guarded_report["boundaries"]
        },
        "crossseam": {
            row["port_id"]: row for row in crossseam_report["boundaries"]
        },
    }
    seam_by_port = {row["port_id"]: row for row in seam_topology["per_port"]}
    rows: list[dict[str, Any]] = []
    for port_id in reports["crossseam"]:
        rows.append(
            {
                "port_id": port_id,
                "global_remesh_normal_jump_P95_deg": reports["global"][port_id][
                    "normal_jump_P95_deg"
                ],
                "global_remesh_normal_jump_P99_deg": reports["global"][port_id][
                    "normal_jump_P99_deg"
                ],
                "guarded_remesh_normal_jump_P95_deg": reports["guarded"][port_id][
                    "normal_jump_P95_deg"
                ],
                "guarded_remesh_normal_jump_P99_deg": reports["guarded"][port_id][
                    "normal_jump_P99_deg"
                ],
                "crossseam_remesh_normal_jump_P95_deg": reports["crossseam"][
                    port_id
                ]["normal_jump_P95_deg"],
                "crossseam_remesh_normal_jump_P99_deg": reports["crossseam"][
                    port_id
                ]["normal_jump_P99_deg"],
                "original_seam_survival_fraction": seam_by_port[port_id][
                    "surviving_fraction"
                ],
                "original_seam_closed_loop_survives": seam_by_port[port_id][
                    "original_seam_closed_loop_survives"
                ],
            }
        )
    return rows


def run_vmtk_surface_prepare(
    config: SurfacePrepareConfig,
    *,
    project_root: Path,
    run_id: str | None = None,
) -> VmtkSurfacePrepareResult:
    """Run one and only one FAR_CORE-excluded cross-seam VMTK candidate."""

    settings = config.vmtk.entity_remesh
    valid = (
        config.backend.method == "vmtk_flowextensions"
        and config.vmtk.extension_mode == "boundarynormal"
        and config.vmtk.postprocess_mode
        == "cross_seam_active_collar_remesh_then_cap"
        and config.vmtk.remesh_after_extension
        and settings.far_core_entity_id == 1
        and settings.active_entity_id == 2
        and settings.exclude_entity_ids == (1,)
        and settings.core_collar.mode == "core_face_adjacency_layers"
        and settings.core_collar.face_layers == 2
        and settings.element_size_mode == "edgelength"
        and settings.target_edge_length_um == 0.25913916380971913
        and settings.preserve_boundary_edges
    )
    if not valid:
        raise SurfacePrepareError("INVALID_VMTK_POSTPROCESS_CONFIGURATION")

    inputs = load_surface_inputs(
        config.paths.cfd_preprocess_run,
        expected_boundary_count=config.qc.expected_boundary_count,
    )
    anchor = _anchor_id(str(inputs.preprocess_summary["roi_id"]))
    candidate_id = run_id or (
        f"vmtk_tps_boundarynormal_crossseam_remesh_anchor{anchor}_"
        f"{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    layout = _layout(config.paths.output_root, candidate_id)
    tool_script = project_root / "tools" / "run_vmtk_flowextension.py"
    paths = exchange_paths(
        input_directory=layout.input,
        vmtk_directory=layout.vmtk,
        geometry_directory=layout.geometry,
        extension_mode="boundarynormal",
    )
    references = _crossseam_reference_paths(config.paths.output_root)
    outputs: dict[str, Any] = {}
    runtime: dict[str, Any] = {}

    def stop(status: str, evidence: dict[str, Any]) -> None:
        _write_failure(
            layout,
            status=status,
            run_id=candidate_id,
            runtime=runtime,
            outputs=outputs,
            evidence=evidence,
        )
        raise SurfacePrepareError(status)

    immutable_paths = (
        inputs.original_surface_um_stl,
        inputs.original_surface_um_vtp,
        references["global_final_vtp"],
        references["global_summary"],
        references["guarded_final_vtp"],
        references["guarded_summary"],
    )
    hashes_before = {str(path): sha256_file(path) for path in immutable_paths}
    write_json(
        layout.input / "original_surface_reference.json",
        {
            "geometry_reference": inputs.geometry_reference,
            "immutable_hashes_before": hashes_before,
            "global_remesh_reference": str(references["global_root"]),
            "guarded_entity_remesh_reference": str(references["guarded_root"]),
            "reference_role": "READ_ONLY_SEAM_QUALITY_COMPARISON",
        },
    )
    shutil.copy2(config.source_path, layout.input / "cfd_surface_prepare.yaml")
    shutil.copy2(
        inputs.preprocess_run / "roi" / "boundary_conditions.json",
        layout.bc / "boundary_conditions_original.json",
    )

    original = load_original_surface(inputs.original_surface_um_vtp)
    surface = TaggedSurface.from_mesh(original)
    for boundary in inputs.boundaries:
        surface, _, _ = local_plane_cut(
            surface,
            boundary,
            radial_factor=config.local_cut.local_radial_radius_factor,
            axial_back_factor=config.local_cut.local_axial_back_radius_factor,
            axial_forward_factor=config.local_cut.local_axial_forward_radius_factor,
        )
    _open_surface(surface, paths.open_surface_vtp)
    open_data, open_mesh = polydata_mesh(paths.open_surface_vtp)
    open_profile, proximal_loops = open_profile_qc(open_mesh, inputs.boundaries)
    if open_profile["status"] != "PASS" or not boundary_plane_alignment_pass(
        open_profile
    ):
        stop(GEOMETRY_FAILURE, {"proximal_open_profile": open_profile})
    _write_manifest(
        layout.input / "boundary_manifest.csv",
        _manifest_records(inputs.boundaries),
        open_profile,
    )
    write_json(
        layout.input / "centerline_usage.json",
        {
            "extension_mode": "boundarynormal",
            "centerline_adapter_generated": False,
            "centerlines_used_for_extension_direction": False,
        },
    )

    local_targets: dict[int, float] = {}
    local_mesh_reports: list[dict[str, Any]] = []
    for boundary in inputs.boundaries:
        local_report = measure_local_original_mesh(
            original,
            boundary,
            sampling_radius_factor=(
                config.extension_mesh.refinement.local_mesh_sampling_radius_factor
            ),
        )
        local_mesh_reports.append(local_report)
        local_targets[boundary.index] = float(local_report["edge_length_median_um"])
    mapping = parameter_mapping(config.vmtk)
    mapping["local_original_mesh_targets"] = local_mesh_reports
    write_json(layout.vmtk / "parameters.json", mapping)
    write_json(layout.vmtk / "vmtk_parameter_mapping.json", mapping)

    invocation = run_official_vmtk(
        config=config.vmtk, paths=paths, tool_script=tool_script
    )
    runtime.update(
        {
            **invocation.runtime,
            "command": list(invocation.command),
            "configured_python_environment": "pmp",
            "source_provenance": official_source_provenance(config.vmtk),
            "global_surface_remeshing_performed": False,
            "remesh_expected_entity_ids": [1, 2],
            "remesh_excluded_entity_ids": [1],
            "remesh_active_entity_ids": [2],
            "synthetic_preflight_run": False,
        }
    )
    outputs.update(
        {"raw_vtp": str(paths.raw_vtp.resolve()), "raw_stl": str(paths.raw_stl.resolve())}
    )
    raw_core = raw_core_exact_copy_qc(paths.open_surface_vtp, paths.raw_vtp)
    if raw_core["status"] != "PASS":
        stop(CORE_FAILURE, {"raw_core_exact_copy": raw_core})

    assignment, collar_layers = assign_cross_seam_active_entities(
        paths.raw_vtp,
        inputs.boundaries,
        face_layers=settings.core_collar.face_layers,
        entity_array_name=settings.entity_array_name,
        far_core_entity_id=settings.far_core_entity_id,
        active_entity_id=settings.active_entity_id,
    )
    write_json(layout.qc / "cross_seam_entity_assignment_qc.json", assignment)
    if assignment["status"] != "PASS":
        stop(ASSIGNMENT_FAILURE, {"cross_seam_entity_assignment": assignment})

    raw_data, raw_mesh = polydata_mesh(paths.raw_vtp)
    raw_topology, raw_pairs = topology_qc(
        raw_mesh, expected_open_profile_count=4, allow_degenerate=False
    )
    raw_regions = np.asarray(raw_data.cell_data["SurfaceRegionId"])
    raw_collision = _collision_report(raw_mesh, raw_pairs, raw_regions)
    raw_interface = interface_smoothness_from_raw(
        raw_mesh,
        input_point_count=open_data.n_points,
        boundaries=inputs.boundaries,
    )
    raw_geometry, _ = extension_geometry_qc(
        raw_mesh, inputs.boundaries, proximal_loops, raw_interface
    )
    if not raw_geometry_hard_gate_pass(raw_topology, raw_geometry, raw_collision):
        stop(
            GEOMETRY_FAILURE,
            {
                "raw_topology": raw_topology,
                "raw_geometry": raw_geometry,
                "raw_collision": raw_collision,
            },
        )

    entity_invocation = entity_remesh_official_vmtk(
        config=config.vmtk, paths=paths, tool_script=tool_script
    )
    runtime["cross_seam_entity_remesh"] = entity_invocation.runtime
    runtime["cross_seam_entity_remesh_command"] = list(entity_invocation.command)
    far_core_remesh, boundary_lock, collar_geometry, extension_effect = (
        cross_seam_entity_preservation_qc(
            paths.raw_vtp,
            paths.remeshed_open_vtp,
            inputs.boundaries,
            entity_array_name=settings.entity_array_name,
            far_core_entity_id=settings.far_core_entity_id,
            active_entity_id=settings.active_entity_id,
            maximum_p95_distance_um=(
                config.qc.maximum_core_surface_p95_distance_um
            ),
            maximum_distance_um=config.qc.maximum_core_surface_distance_um,
        )
    )
    far_core_report: dict[str, Any] = {
        "status": far_core_remesh["status"],
        "post_remesh": far_core_remesh,
        "post_cap": None,
    }
    write_json(layout.qc / "far_core_exact_preservation_qc.json", far_core_report)
    write_json(layout.qc / "active_collar_geometry_qc.json", collar_geometry)
    if far_core_remesh["status"] != "PASS" or boundary_lock["status"] != "PASS":
        stop(
            CORE_FAILURE,
            {"far_core": far_core_report, "far_core_active_boundary": boundary_lock},
        )
    if collar_geometry["status"] != "PASS":
        stop(COLLAR_FAILURE, {"active_collar_geometry": collar_geometry})
    if extension_effect["status"] != "PASS":
        stop(GEOMETRY_FAILURE, {"extension_remesh_effect": extension_effect})
    outputs.update(
        {
            "remeshed_open_vtp": str(paths.remeshed_open_vtp.resolve()),
            "remeshed_open_stl": str(paths.remeshed_open_stl.resolve()),
        }
    )

    seam_topology = cut_seam_topology_qc(
        paths.raw_vtp, paths.remeshed_open_vtp, inputs.boundaries
    )
    write_json(layout.qc / "cut_seam_topology_qc.json", seam_topology)
    if seam_topology["status"] != "PASS":
        stop(RING_FAILURE, {"cut_seam_topology": seam_topology})

    remesh_data, remesh_mesh = polydata_mesh(paths.remeshed_open_vtp)
    remesh_topology_raw, _ = topology_qc(
        remesh_mesh, expected_open_profile_count=4, allow_degenerate=False
    )
    intersection_qc, intersection_records = cross_seam_intersection_qc(
        paths.remeshed_open_vtp
    )
    remesh_topology = _use_true_intersection_count(
        remesh_topology_raw, intersection_qc
    )
    remesh_profile, _ = open_profile_qc(
        remesh_mesh, inputs.boundaries, distal=True
    )
    remesh_topology["distal_profile_qc"] = remesh_profile
    remesh_topology["status"] = (
        "PASS"
        if remesh_topology["status"] == "PASS"
        and remesh_profile["status"] == "PASS"
        else "FAIL"
    )
    write_json(layout.qc / "crossseam_open_surface_qc.json", remesh_topology)
    write_json(layout.qc / "crossseam_intersection_qc.json", intersection_qc)
    if not should_cap_crossseam_open(
        remesh_topology, intersection_qc, seam_topology
    ):
        stop(
            TOPOLOGY_FAILURE,
            {
                "open_topology": remesh_topology,
                "intersections": intersection_qc,
                "intersection_records": intersection_records,
            },
        )

    seam_new = seam_local_normal_diagnostic(
        paths.remeshed_open_vtp,
        inputs.boundaries,
        target_edge_length_um=settings.target_edge_length_um,
    )
    remesh_geometry, distal_loops = extension_geometry_qc(
        remesh_mesh, inputs.boundaries, proximal_loops, seam_new
    )
    write_json(layout.qc / "extension_geometry_qc.json", remesh_geometry)
    if remesh_geometry["status"] != "PASS":
        stop(GEOMETRY_FAILURE, {"extension_geometry": remesh_geometry})

    active_quality, tail_records = crossseam_active_mesh_quality(
        paths.remeshed_open_vtp,
        inputs.boundaries,
        active_entity_id=settings.active_entity_id,
        local_target_edge_um=local_targets,
        quality=config.mesh_quality,
    )
    write_json(layout.qc / "crossseam_mesh_quality_qc.json", active_quality)
    seam_global = seam_local_normal_diagnostic(
        references["global_final_vtp"],
        inputs.boundaries,
        target_edge_length_um=settings.target_edge_length_um,
    )
    seam_guarded = seam_local_normal_diagnostic(
        references["guarded_final_vtp"],
        inputs.boundaries,
        target_edge_length_um=settings.target_edge_length_um,
    )
    seam_rows = _seam_comparison_rows(
        seam_global, seam_guarded, seam_new, seam_topology
    )
    write_csv(layout.qc / "seam_quality_comparison.csv", seam_rows)

    cap = cap_official_vmtk(
        config=config.vmtk, paths=paths, tool_script=tool_script
    )
    runtime["cap_only"] = cap.runtime
    runtime["cap_only_command"] = list(cap.command)
    if cap.runtime.get("surface_remesher_called") is not False:
        stop(GEOMETRY_FAILURE, {"cap_runtime": cap.runtime})
    write_json(layout.vmtk / "environment.json", runtime)
    final_outputs, boundary_mapping = tag_and_export_final_surface(
        paths.capped_vtp,
        inputs.boundaries,
        layout.geometry,
        layout.boundaries,
        output_stem="cfd_surface_vmtk_tps_boundarynormal_crossseam",
        raw_vtp=paths.remeshed_open_vtp,
        remesh_entity_codes={
            "CAP": 0,
            "FAR_CORE": 1,
            "CROSS_SEAM_ACTIVE": 2,
        },
    )
    outputs.update(final_outputs)
    _, final_mesh = polydata_mesh(Path(outputs["tagged_vtp"]))
    final_topology_raw, _ = topology_qc(
        final_mesh,
        expected_open_profile_count=0,
        require_winding_consistent=True,
    )
    final_intersection, _ = cross_seam_intersection_qc(
        Path(outputs["tagged_vtp"])
    )
    final_topology = _use_true_intersection_count(
        final_topology_raw, final_intersection
    )
    far_core_final = locked_entity_preservation_qc(
        paths.raw_vtp,
        Path(outputs["tagged_vtp"]),
        entity_array_name=settings.entity_array_name,
        entity_id=settings.far_core_entity_id,
        label="far_core",
    )
    far_core_report = {
        "status": (
            "PASS"
            if far_core_remesh["status"] == far_core_final["status"] == "PASS"
            else "FAIL"
        ),
        "post_remesh": far_core_remesh,
        "post_cap": far_core_final,
    }
    write_json(layout.qc / "far_core_exact_preservation_qc.json", far_core_report)
    final_surface = {
        "status": (
            "PASS"
            if final_topology["status"] == "PASS"
            and boundary_mapping["status"] == "PASS"
            and far_core_report["status"] == "PASS"
            else "FAIL"
        ),
        "topology": final_topology,
        "intersection": final_intersection,
        "boundary_mapping": boundary_mapping,
        "far_core_exact_preservation": far_core_report,
    }
    write_json(layout.qc / "final_surface_qc.json", final_surface)
    if final_topology["status"] != "PASS":
        stop(TOPOLOGY_FAILURE, {"final_surface": final_surface})
    if boundary_mapping["status"] != "PASS":
        stop(GEOMETRY_FAILURE, {"final_surface": final_surface})
    if far_core_report["status"] != "PASS":
        stop(CORE_FAILURE, {"final_surface": final_surface})

    try:
        radius_report = _radius_qc(
            final_mesh=final_mesh, inputs=inputs, project_root=project_root
        )
    except SurfacePrepareError as error:
        radius_report = {"status": "FAIL", "error": str(error)}
    meter_report = meter_scale_qc(
        Path(outputs["manual_review_stl"]), Path(outputs["meter_stl"])
    )
    write_json(layout.qc / "radius_fidelity.json", radius_report)
    write_json(layout.qc / "meter_scale_qc.json", meter_report)

    pressure_rows: list[dict[str, Any]] = []
    try:
        pressure_rows, pressure_report = geometry_pressure_correction(
            final_mesh,
            inputs.boundaries,
            proximal_loops,
            distal_loops,
            dynamic_viscosity_pa_s=float(
                inputs.original_boundary_conditions["fluid"]["dynamic_viscosity_pa_s"]
            ),
        )
    except SurfacePrepareError as error:
        pressure_report = {"status": "FAIL", "error": str(error)}
    else:
        pressure_report["status"] = "PASS"
        pressure_csv = (
            layout.bc
            / "extension_pressure_correction_vmtk_boundarynormal_crossseam.csv"
        )
        write_csv(
            pressure_csv,
            [
                {
                    **row,
                    "station_fractions": json.dumps(row["station_fractions"]),
                    "cross_section_area_um2": json.dumps(
                        row["cross_section_area_um2"]
                    ),
                    "equivalent_radius_um": json.dumps(row["equivalent_radius_um"]),
                }
                for row in pressure_rows
            ],
        )
        boundary_json = (
            layout.bc / "boundary_conditions_vmtk_boundarynormal_crossseam.json"
        )
        write_json(
            boundary_json,
            _boundary_conditions(inputs.original_boundary_conditions, pressure_rows),
        )
        outputs["pressure_correction_csv"] = str(pressure_csv.resolve())
        outputs["boundary_conditions_json"] = str(boundary_json.resolve())
    write_json(layout.qc / "extension_pressure_geometry_qc.json", pressure_report)

    hashes_after = {str(path): sha256_file(path) for path in immutable_paths}
    integrity_pass = hashes_before == hashes_after
    hard_pass = all(
        report.get("status") == "PASS"
        for report in (
            final_surface,
            collar_geometry,
            remesh_geometry,
            radius_report,
            meter_report,
            pressure_report,
        )
    ) and integrity_pass
    if not hard_pass:
        if radius_report.get("status") != "PASS":
            failure_status = GEOMETRY_FAILURE
        elif meter_report.get("status") != "PASS":
            failure_status = METER_FAILURE
        else:
            failure_status = GEOMETRY_FAILURE
        stop(
            failure_status,
            {
                "final_surface": final_surface,
                "active_collar_geometry": collar_geometry,
                "extension_geometry": remesh_geometry,
                "radius": radius_report,
                "meter": meter_report,
                "pressure": pressure_report,
                "read_only_reference_integrity": integrity_pass,
            },
        )

    figure_paths = save_crossseam_review_figures(
        raw_vtp=paths.raw_vtp,
        global_remesh_vtp=references["global_final_vtp"],
        guarded_remesh_vtp=references["guarded_final_vtp"],
        crossseam_final_vtp=Path(outputs["tagged_vtp"]),
        boundaries=inputs.boundaries,
        output_directory=layout.figures,
    )
    outputs["figure_paths"] = [str(path) for path in figure_paths]
    tail_warning = bool(tail_records)
    status = TAIL_WARNING_STATUS if tail_warning else SUCCESS_STATUS
    summary = {
        "status": status,
        "next": SUCCESS_NEXT,
        "visible_ring_assessment": "MANUAL_REVIEW_REQUIRED",
        "run_id": candidate_id,
        "run_root": str(layout.root),
        "runtime": runtime,
        "cross_seam_entity_assignment": assignment,
        "core_collar_layers": collar_layers.tolist(),
        "far_core_exact_preservation": far_core_report,
        "far_core_active_boundary_lock": boundary_lock,
        "active_collar_geometry": collar_geometry,
        "extension_remesh_effect": extension_effect,
        "cut_seam_topology": seam_topology,
        "crossseam_open_topology": remesh_topology,
        "crossseam_open_intersections": intersection_qc,
        "extension_geometry": remesh_geometry,
        "crossseam_mesh_quality": active_quality,
        "mesh_tail_count": len(tail_records),
        "seam_normal_diagnostics": {
            "global_remesh": seam_global,
            "guarded_remesh": seam_guarded,
            "crossseam_remesh": seam_new,
        },
        "seam_quality_comparison": seam_rows,
        "final_surface": final_surface,
        "radius_fidelity": radius_report,
        "pressure_geometry": pressure_report,
        "meter_scale": meter_report,
        "outputs": outputs,
        "read_only_references_untouched": integrity_pass,
        "single_primary_candidate_count": 1,
        "formal_roi_run_count": 1,
        "second_candidate_run": False,
        "automatic_fallback": False,
        "parameter_sweep": False,
        "synthetic_preflight_run": False,
        "global_surface_remeshing_performed": False,
        "volume_mesh_created": False,
        "cfd_run": False,
        "microbubble_simulation_run": False,
    }
    write_json(layout.qc / "run_summary.json", summary)
    _write_report(
        layout.report / "vmtk_tps_boundarynormal_crossseam_report.md", summary
    )
    return VmtkSurfacePrepareResult(
        status=status,
        next_stage=SUCCESS_NEXT,
        run_root=layout.root,
        runtime=runtime,
        assignment_qc=assignment,
        far_core_qc=far_core_report,
        boundary_lock_qc=boundary_lock,
        collar_qc=collar_geometry,
        extension_qc=extension_effect,
        seam_topology_qc=seam_topology,
        intersection_qc=intersection_qc,
        boundaries=tuple(remesh_geometry["boundaries"]),
        mesh_quality=active_quality,
        seam_quality=tuple(seam_rows),
        open_topology=remesh_topology,
        final_topology=final_topology,
        radius_qc=radius_report,
        meter_qc=meter_report,
        collision_count=int(intersection_qc["true_self_intersection_count"]),
        pressure_rows=tuple(pressure_rows),
        output_paths=outputs,
        figures=figure_paths,
    )


def print_vmtk_experiment_header() -> None:
    print("EXPERIMENT: FAR-CORE EXCLUDED CROSS-SEAM ACTIVE VMTK REMESH")
    print("VMTK: 1.5.0")
    print("VTK: 9.2.6")
    print("TPS: BOUNDARY_NORMAL")
    print("Remesher: OFFICIAL VMTK")
    print("Global remesh: NO")
    print("FAR_CORE entity: 1 (EXCLUDED / IMMUTABLE)")
    print("CROSS_SEAM_ACTIVE entity: 2 (CORE collar + all extension)")
    print("ExcludeEntityIds: [1]")
    print("Active CORE collar: 2 RAW CORE BFS layers (0 and 1)")
    print("Synthetic VMTK preflight: NOT RUN")
    print("Target edge: 0.25913916380971913 um")


def _legacy_print_vmtk_result(result: VmtkSurfacePrepareResult) -> None:
    print(f"Synthetic exclusion test: {result.synthetic_preflight['status']}")
    core = result.entity_core_qc
    print("ENTITY REMESH CORE:")
    print(
        f"face count before/after = {core['raw_core_face_count']}/"
        f"{core['remeshed_core_face_count']}"
    )
    print(
        f"vertex count before/after = {core['raw_core_vertex_count']}/"
        f"{core['remeshed_core_vertex_count']}"
    )
    print(
        f"max/P95 motion = {core['core_vertex_max_motion_um']:.12g}/"
        f"{core['core_vertex_P95_motion_um']:.12g} um"
    )
    print(
        "connectivity changed/missing/added = "
        f"{core['core_connectivity_changed_count']}/"
        f"{core['core_triangle_missing_count']}/{core['core_triangle_added_count']}"
    )
    interface = result.entity_interface_qc
    print(
        "shared interface vertices/max/P95 motion = "
        f"{interface['interface_shared_vertex_count']}/"
        f"{interface['interface_shared_vertex_max_motion_um']:.12g}/"
        f"{interface['interface_shared_vertex_P95_motion_um']:.12g} um"
    )
    extension = result.entity_extension_qc
    print(
        "extension triangles before/after = "
        f"{extension['raw_extension_triangle_count']}/"
        f"{extension['remeshed_extension_triangle_count']}"
    )
    print(
        "extension vertices before/after = "
        f"{extension['raw_extension_vertex_count']}/"
        f"{extension['remeshed_extension_vertex_count']}"
    )
    print(
        "extension connectivity changed: "
        f"{'YES' if extension['extension_connectivity_changed'] else 'NO'}"
    )
    reports = {
        report["method"]: {row["port_id"]: row for row in report["boundaries"]}
        for report in (result.raw_quality, result.remeshed_quality, result.previous_quality)
    }
    for geometry in result.boundaries:
        port_id = geometry["port_id"]
        print(port_id)
        for method in (
            "RAW_CAP_ONLY_NO_REMESH",
            "NEW_ENTITY_AWARE_REMESH",
            "PREVIOUS_GLOBAL_REMESH",
        ):
            row = reports[method][port_id]
            print(
                f"  {method}: edge={row['edge_length_median_um']:.12g} "
                f"min_angle={row['minimum_angle_deg']:.12g} "
                f"aspect_P95/max={row['aspect_ratio_P95']:.12g}/"
                f"{row['aspect_ratio_max']:.12g} "
                f"size_mismatch={row['symmetric_mesh_size_mismatch']:.12g}"
            )
        print(
            f"  geometry: direction={geometry['extension_direction_dot']:.12g} "
            f"length_error={100.0 * geometry['extension_length_relative_error']:.12g}% "
            f"area_error={100.0 * geometry['distal_area_relative_error']:.12g}% "
            f"collision={result.collision_count}"
        )
    distance = result.core_distance_qc
    normal = result.normal_qc
    print(
        "final core original->final P95/max = "
        f"{distance['original_core_to_entityremesh_final_core']['P95_um']:.12g}/"
        f"{distance['original_core_to_entityremesh_final_core']['max_um']:.12g} um"
    )
    print(
        "final core final->original P95/max = "
        f"{distance['entityremesh_final_core_to_original_core']['P95_um']:.12g}/"
        f"{distance['entityremesh_final_core_to_original_core']['max_um']:.12g} um"
    )
    print(
        "normal P95/P99/max = "
        f"{normal['core_normal_deviation_P95_deg']:.12g}/"
        f"{normal['core_normal_deviation_P99_deg']:.12g}/"
        f"{normal['core_normal_deviation_max_deg']:.12g} deg"
    )
    topology = result.final_topology
    print(
        "Topology: "
        f"components={topology['component_count']} watertight={topology['watertight']} "
        f"boundary_edges={topology['boundary_edge_count']} "
        f"nonmanifold={topology['nonmanifold_edge_count']} "
        f"self_intersections={topology['self_intersection_count']} "
        f"degenerate={topology['degenerate_triangle_count']} "
        f"winding_consistent={topology['winding_consistent']}"
    )
    print(f"Radius P95: {result.radius_qc.get('p95_absolute_relative_error')}")
    print(f"Collision count: {result.collision_count}")
    for row in result.pressure_rows:
        if row["role"] == "ASSUMED_OUTLET":
            print(
                f"Outlet {row['port_id']}: P_original={row['P_original_1D_pa']:.12g} "
                f"Q={row['Q_expected_1D_m3_s']:.12g} "
                f"R={row['extension_resistance_pa_s_m3']:.12g} "
                f"delta_P={row['predicted_extension_pressure_drop_pa']:.12g} "
                f"P_solver={row['P_solver_boundary_pa']:.12g}"
            )
    print(
        "meter serialization exact_after_float32_cast: "
        f"{'YES' if result.meter_qc['exact_after_float32_cast'] else 'NO'}"
    )
    print(
        "legacy extent relative error: "
        f"{result.meter_qc['legacy_extent_relative_error_max']:.12g}"
    )
    print(
        "MANUAL REVIEW STL: "
        f"{result.output_paths.get('manual_review_stl', 'NOT_GENERATED')}"
    )
    print(f"TAGGED VTP: {result.output_paths.get('tagged_vtp', 'NOT_GENERATED')}")
    for path in result.figures:
        print(f"FIGURE: {path}")
    print("visible_acceptance: MANUAL_REVIEW_REQUIRED")
    print(f"Final status: {result.status}")
    print(f"NEXT: {result.next_stage}")


def _historical_print_guarded_result(result: VmtkSurfacePrepareResult) -> None:
    print(f"Curved guarded exclusion test: {result.synthetic_preflight['status']}")
    assignment = result.guard_assignment_qc
    print(
        "GUARD assignment: "
        f"CORE={assignment['core_face_count']} "
        f"GUARD={assignment['guard_face_count']} "
        f"BODY={assignment['body_face_count']} "
        f"unknown={assignment['unknown_face_count']}"
    )
    for row in assignment["per_port"]:
        print(
            f"  {row['port_id']}: guard_faces={row['guard_face_count']} "
            f"body_faces={row['extension_body_face_count']} "
            f"guard_width={row['guard_maximum_axial_width_um']:.12g} um "
            f"({row['guard_width_in_source_radius']:.12g}R, "
            f"{row['guard_width_in_diameter']:.12g}D)"
        )
    core = result.core_qc
    print(
        "CORE faces/vertices before-after: "
        f"{core['raw_core_face_count']}/{core['remeshed_core_face_count']} | "
        f"{core['raw_core_vertex_count']}/{core['remeshed_core_vertex_count']}"
    )
    print(
        "CORE max/P95 motion and connectivity changes: "
        f"{core['core_vertex_max_motion_um']:.12g}/"
        f"{core['core_vertex_P95_motion_um']:.12g} um | "
        f"{core['core_connectivity_changed_count']}"
    )
    guard = result.guard_qc
    print(
        "GUARD faces/vertices before-after: "
        f"{guard['raw_guard_face_count']}/{guard['remeshed_guard_face_count']} | "
        f"{guard['raw_guard_vertex_count']}/{guard['remeshed_guard_vertex_count']}"
    )
    print(
        "GUARD max/P95 motion and connectivity changes: "
        f"{guard['guard_vertex_max_motion_um']:.12g}/"
        f"{guard['guard_vertex_P95_motion_um']:.12g} um | "
        f"{guard['guard_connectivity_changed_count']}"
    )
    boundaries = result.boundary_lock_qc
    print(
        "CORE-GUARD shared vertices/max/P95 motion: "
        f"{boundaries['core_guard_shared_vertex_count']}/"
        f"{boundaries['core_guard_max_motion_um']:.12g}/"
        f"{boundaries['core_guard_P95_motion_um']:.12g} um"
    )
    print(
        "GUARD-BODY shared vertices/max/P95 motion: "
        f"{boundaries['guard_body_shared_vertex_count']}/"
        f"{boundaries['guard_body_max_motion_um']:.12g}/"
        f"{boundaries['guard_body_P95_motion_um']:.12g} um"
    )
    body = result.body_qc
    print(
        "BODY faces/vertices before-after and connectivity changed: "
        f"{body['raw_body_face_count']}/{body['remeshed_body_face_count']} | "
        f"{body['raw_body_vertex_count']}/{body['remeshed_body_vertex_count']} | "
        f"{'YES' if body['body_connectivity_changed'] else 'NO'}"
    )
    intersections = result.intersection_qc
    print(
        "Post-remesh true intersections: "
        f"{intersections['true_self_intersection_count']} | "
        f"{intersections['classification_counts']}"
    )
    quality_by_entity = {
        "GUARD": result.guard_quality,
        "BODY": result.body_quality,
    }
    for geometry in result.boundaries:
        print(
            f"{geometry['port_id']}: direction={geometry['extension_direction_dot']:.12g} "
            f"length_error={100.0 * geometry['extension_length_relative_error']:.12g}% "
            f"area_error={100.0 * geometry['distal_area_relative_error']:.12g}%"
        )
        for label, report in quality_by_entity.items():
            row = next(
                item
                for item in report["boundaries"]
                if item["port_id"] == geometry["port_id"]
            )
            print(
                f"  {label}: triangles={row['triangle_count']} "
                f"edge={row['edge_length_median_um']:.12g} "
                f"min_angle={row['minimum_angle_deg']:.12g} "
                f"aspect_P95/max={row['aspect_ratio_P95']:.12g}/"
                f"{row['aspect_ratio_max']:.12g} "
                f"tail_angle/aspect={row['triangle_count_angle_below_5deg']}/"
                f"{row['triangle_count_aspect_above_20']}"
            )
    print(f"Previous tail localization: {result.previous_tail}")
    print(f"Final topology: {result.final_topology}")
    print(f"Radius P95: {result.radius_qc.get('p95_absolute_relative_error')}")
    print(
        "Meter exact after float32 cast: "
        f"{'YES' if result.meter_qc['exact_after_float32_cast'] else 'NO'}"
    )
    for row in result.pressure_rows:
        if row["role"] == "ASSUMED_OUTLET":
            print(
                f"Outlet {row['port_id']}: P_original={row['P_original_1D_pa']:.12g} "
                f"Q={row['Q_expected_1D_m3_s']:.12g} "
                f"R={row['extension_resistance_pa_s_m3']:.12g} "
                f"delta_P={row['predicted_extension_pressure_drop_pa']:.12g} "
                f"P_solver={row['P_solver_boundary_pa']:.12g}"
            )
    print(f"MANUAL REVIEW STL: {result.output_paths.get('manual_review_stl')}")
    print(f"TAGGED VTP: {result.output_paths.get('tagged_vtp')}")
    for path in result.figures:
        print(f"FIGURE: {path}")
    print("Volume mesh/CFD/microbubble: NOT RUN")
    print(f"Final status: {result.status}")
    print(f"NEXT: {result.next_stage}")


def print_vmtk_result(result: VmtkSurfacePrepareResult) -> None:
    assignment = result.assignment_qc
    print(
        "Entities: FAR_CORE=1, CROSS_SEAM_ACTIVE=2 | "
        "ExcludeEntityIds=[1]"
    )
    print(
        "Original cut seam edges/different-entity edges: "
        f"{assignment['original_cut_seam_edge_count']}/"
        f"{assignment['original_cut_seam_edges_between_different_remesh_entities']}"
    )
    for row in assignment["per_port"]:
        print(
            f"  {row['port_id']}: active_collar_faces="
            f"{row['active_core_collar_face_count']} "
            f"approx_width={row['collar_approximate_axial_width_um']:.12g} um"
        )
    far_core = result.far_core_qc["post_cap"]
    print(
        "FAR_CORE faces/vertices before-after: "
        f"{far_core['raw_far_core_face_count']}/"
        f"{far_core['remeshed_far_core_face_count']} | "
        f"{far_core['raw_far_core_vertex_count']}/"
        f"{far_core['remeshed_far_core_vertex_count']}"
    )
    print(
        "FAR_CORE max motion/connectivity changes: "
        f"{far_core['far_core_vertex_max_motion_um']:.12g} um | "
        f"{far_core['far_core_connectivity_changed_count']}"
    )
    boundary = result.boundary_lock_qc
    print(
        "FAR_CORE-ACTIVE shared vertices/max/P95 motion: "
        f"{boundary['far_core_active_shared_vertex_count']}/"
        f"{boundary['far_core_active_max_motion_um']:.12g}/"
        f"{boundary['far_core_active_P95_motion_um']:.12g} um"
    )
    collar = result.collar_qc
    print(
        "ACTIVE collar P95/max deviation: "
        f"{collar['bidirectional_P95_um']:.12g}/"
        f"{collar['bidirectional_max_um']:.12g} um"
    )
    extension = result.extension_qc
    print(
        "Extension faces before-after/connectivity changed: "
        f"{extension['raw_extension_face_count']}/"
        f"{extension['remeshed_extension_face_count']} | "
        f"{'YES' if extension['extension_connectivity_changed'] else 'NO'}"
    )
    seam = result.seam_topology_qc
    print(
        "Original seam survival/closed loop survives: "
        f"{seam['surviving_fraction']:.12g} | "
        f"{'YES' if seam['original_cut_seam_closed_loop_survives'] else 'NO'}"
    )
    print(
        "Post-remesh true intersections: "
        f"{result.intersection_qc['true_self_intersection_count']}"
    )
    for geometry in result.boundaries:
        quality = next(
            row
            for row in result.mesh_quality["boundaries"]
            if row["port_id"] == geometry["port_id"]
        )
        print(
            f"{geometry['port_id']}: direction={geometry['extension_direction_dot']:.12g} "
            f"length_error={100.0 * geometry['extension_length_relative_error']:.12g}% "
            f"area_error={100.0 * geometry['distal_area_relative_error']:.12g}% | "
            f"tail_angle/aspect={quality['triangle_count_angle_below_5deg']}/"
            f"{quality['triangle_count_aspect_above_20']}"
        )
    print(f"Open topology: {result.open_topology['status']}")
    print(f"Final topology: {result.final_topology['status']}")
    print(f"Radius P95: {result.radius_qc.get('p95_absolute_relative_error')}")
    print(
        "Meter exact after float32 cast: "
        f"{'YES' if result.meter_qc['exact_after_float32_cast'] else 'NO'}"
    )
    for row in result.pressure_rows:
        if row["role"] == "ASSUMED_OUTLET":
            print(
                f"Outlet {row['port_id']}: delta_P="
                f"{row['predicted_extension_pressure_drop_pa']:.12g} "
                f"P_solver={row['P_solver_boundary_pa']:.12g}"
            )
    print(f"MANUAL REVIEW STL: {result.output_paths.get('manual_review_stl')}")
    print(f"TAGGED VTP: {result.output_paths.get('tagged_vtp')}")
    for path in result.figures:
        print(f"FIGURE: {path}")
    print("visible_ring_assessment: MANUAL_REVIEW_REQUIRED")
    print("Formal ROI run count: 1")
    print("Volume mesh/CFD/microbubble: NOT RUN")
    print(f"Final status: {result.status}")
    print(f"NEXT: {result.next_stage}")
