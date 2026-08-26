"""Formal single-candidate VMTK TPS CFD-surface preparation pipeline."""

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
from utils.cfd_lumen.ultraliser_qc import (
    evaluate_radius_fidelity,
    validate_source_roi,
)

from .config import SurfacePrepareConfig
from .export import write_csv
from .io import (
    SurfacePrepareError,
    load_original_surface,
    load_surface_inputs,
    read_json,
    sha256_file,
)
from .local_cut import local_plane_cut
from .qc import core_surface_preservation_qc
from .types import TaggedSurface
from .vmtk_qc import (
    extension_geometry_qc,
    extension_mesh_metrics,
    geometry_pressure_correction,
    interface_smoothness_from_old_custom,
    interface_smoothness_from_raw,
    meter_scale_qc,
    open_profile_qc,
    polydata_mesh,
    tag_and_export_final_surface,
    topology_qc,
    write_json,
)
from .vmtk_runner import (
    exchange_paths,
    official_source_provenance,
    parameter_mapping,
    promote_official_vmtk,
    run_official_vmtk,
)
from .vmtk_visualization import save_vmtk_review_figures


RAW_FAILURE_STATUS = "VMTK_BOUNDARY_NORMAL_RAW_GEOMETRY_FAILED"
FINAL_FAILURE_STATUS = "VMTK_BOUNDARY_NORMAL_FINAL_SURFACE_FAILED"
WARNING_STATUS = (
    "VMTK_TPS_BOUNDARY_NORMAL_INTERFACE_WARNING_PENDING_MANUAL_REVIEW"
)
SUCCESS_STATUS = "VMTK_TPS_BOUNDARY_NORMAL_PASS_PENDING_MANUAL_REVIEW"
PASS_STATUSES = {WARNING_STATUS, SUCCESS_STATUS}
SUCCESS_NEXT = "MANUALLY REVIEW VMTK BOUNDARY-NORMAL TPS SURFACE"
FAIL_NEXT = "REVIEW VMTK BOUNDARY-NORMAL TPS FAILURE"
CENTERLINE_REFERENCE_RUN_ID = "vmtk_tps_anchor003274_20260826_144705"


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
    boundaries: tuple[dict[str, Any], ...]
    raw_topology: dict[str, Any]
    final_topology: dict[str, Any]
    radius_qc: dict[str, Any]
    core_qc: dict[str, Any]
    collision_count: int
    interface_comparison: dict[str, Any]
    cut002_diagnosis: dict[str, Any]
    output_paths: dict[str, Any]
    figures: tuple[Path, ...]


def _layout(output_root: Path, run_id: str) -> VmtkLayout:
    root = Path(output_root) / run_id
    if root.exists():
        raise FileExistsError(f"Refusing to overwrite VMTK run: {root}")
    directories = [root / name for name in ("input", "vmtk", "geometry", "boundaries", "bc", "qc", "figures", "report")]
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=False)
    return VmtkLayout(root.resolve(), *[directory.resolve() for directory in directories])


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


def _global_median_edge_length(mesh: Any) -> float:
    unique = np.unique(np.sort(np.asarray(mesh.edges, dtype=np.int64), axis=1), axis=0)
    lengths = np.linalg.norm(mesh.vertices[unique[:, 1]] - mesh.vertices[unique[:, 0]], axis=1)
    value = float(np.median(lengths))
    if not np.isfinite(value) or value <= 0.0:
        raise SurfacePrepareError("VMTK_EXTENSION_GEOMETRY_FAILED:invalid_remesh_target")
    return value


def _old_reference(root: Path) -> tuple[Path, Path]:
    vtp = root / "geometry" / "cfd_surface_refined_um.vtp"
    stl = root / "geometry" / "cfd_surface_refined_um.stl"
    if not vtp.is_file() or not stl.is_file():
        raise SurfacePrepareError(f"Missing OLD_CUSTOM_EXTENSION_REFERENCE: {root}")
    return vtp.resolve(), stl.resolve()


def _centerline_reference(output_root: Path) -> dict[str, Path]:
    root = Path(output_root) / CENTERLINE_REFERENCE_RUN_ID
    paths = {
        "root": root,
        "raw_vtp": root / "geometry" / "vmtk_flowextension_raw_um.vtp",
        "open_vtp": root / "input" / "open_surface_um.vtp",
        "geometry_qc": root / "qc" / "extension_geometry_qc.json",
        "interface_qc": root / "qc" / "interface_smoothness_qc.json",
    }
    missing = [str(path) for key, path in paths.items() if key != "root" and not path.is_file()]
    if missing:
        raise SurfacePrepareError(
            "Missing CENTERLINE_VMTK_REFERENCE: " + ", ".join(missing)
        )
    return {key: path.resolve() for key, path in paths.items()}


def _boundary_manifest_records(boundaries: tuple[Any, ...]) -> tuple[dict[str, Any], ...]:
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


def _write_boundary_manifest(path: Path, records: tuple[dict[str, Any], ...], profile: dict[str, Any]) -> None:
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


def _radius_qc(
    *,
    final_mesh: Any,
    inputs: Any,
    project_root: Path,
) -> dict[str, Any]:
    model_root = Path(str(inputs.geometry_reference["run_root"]))
    source_config = model_root / "input" / "source_swc_stl_model_generate.yaml"
    model_config = load_swc_stl_yaml_config(source_config, project_root=project_root)
    if model_config.sampling_run is None:
        raise SurfacePrepareError("VMTK_SURFACE_QC_FAILED:missing_sampling_run")
    rois = load_sampling_rois(
        model_config.sampling_run,
        roi_id=str(inputs.geometry_reference["roi_id"]),
        selected_only=False,
    )
    roi = rois[0]
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
                "predicted_extension_pressure_drop_pa": row["predicted_extension_pressure_drop_pa"],
                "P_solver_boundary_pa": row["P_solver_boundary_pa"],
                "extension_pressure_correction_role": row["pressure_correction_role"],
            }
        )
    output["surface_geometry_method"] = "OFFICIAL_VMTK_TPS_FLOWEXTENSIONS"
    output["pressure_correction_geometry"] = "20 cross sections on final VMTK geometry"
    return output


def boundary_plane_alignment_pass(profile: dict[str, Any]) -> bool:
    return all(
        float(row["boundary_plane_normal_abs_dot_expected_outward"]) >= 0.999
        for row in profile["boundaries"]
    )


def raw_geometry_hard_gate_pass(
    topology: dict[str, Any],
    geometry: dict[str, Any],
    collision: dict[str, Any],
) -> bool:
    return all(
        report["status"] == "PASS" for report in (topology, geometry, collision)
    )


def should_promote_raw_candidate(raw_hard_qc_pass: bool) -> bool:
    return raw_hard_qc_pass


def should_generate_manual_review_figures(
    *, raw_hard_qc_pass: bool, interface_status: str
) -> bool:
    del interface_status
    return raw_hard_qc_pass


def final_candidate_status(
    *, final_hard_qc_pass: bool, interface_status: str
) -> str:
    if not final_hard_qc_pass:
        return FINAL_FAILURE_STATUS
    return SUCCESS_STATUS if interface_status == "PASS" else WARNING_STATUS


def _extension_direction_comparison(
    custom_geometry: dict[str, Any],
    centerline_geometry: dict[str, Any],
    boundarynormal_geometry: dict[str, Any],
) -> list[dict[str, Any]]:
    custom = {row["port_id"]: row for row in custom_geometry["boundaries"]}
    centerline = {
        row["port_id"]: row for row in centerline_geometry["boundaries"]
    }
    boundarynormal = {
        row["port_id"]: row for row in boundarynormal_geometry["boundaries"]
    }
    rows: list[dict[str, Any]] = []
    for port_id in custom:
        previous = centerline[port_id]
        current = boundarynormal[port_id]
        rows.append(
            {
                "port_id": port_id,
                "custom_direction_dot": custom[port_id]["extension_axis_dot"],
                "centerline_vmtk_direction_dot": previous["extension_direction_dot"],
                "boundarynormal_vmtk_direction_dot": current["extension_direction_dot"],
                "centerline_actual_length_um": previous["actual_extension_length_um"],
                "boundarynormal_actual_length_um": current["actual_axial_length_um"],
                "planned_length_um": current["planned_extension_length_um"],
                "centerline_area_error": previous["distal_area_relative_error"],
                "boundarynormal_area_error": current["distal_area_relative_error"],
            }
        )
    return rows


def _interface_three_way_comparison(
    old_interface: dict[str, Any],
    centerline_interface: dict[str, Any],
    boundarynormal_interface: dict[str, Any],
) -> dict[str, Any]:
    old = {row["port_id"]: row for row in old_interface["boundaries"]}
    centerline = {
        row["port_id"]: row for row in centerline_interface["boundaries"]
    }
    boundarynormal = {
        row["port_id"]: row for row in boundarynormal_interface["boundaries"]
    }
    rows: list[dict[str, Any]] = []
    for port_id in old:
        row: dict[str, Any] = {"port_id": port_id}
        for label, values in (
            ("old_custom", old[port_id]),
            ("centerline_vmtk", centerline[port_id]),
            ("boundarynormal_vmtk_raw", boundarynormal[port_id]),
        ):
            for percentile in ("P50", "P95", "P99", "max"):
                row[f"{label}_{percentile}_deg"] = values[
                    f"normal_jump_{percentile}_deg"
                ]
        row.update(
            {
                "boundarynormal_vmtk_final_P50_deg": None,
                "boundarynormal_vmtk_final_P95_deg": None,
                "boundarynormal_vmtk_final_P99_deg": None,
                "boundarynormal_vmtk_final_max_deg": None,
                "boundarynormal_final_metric_status": (
                    "NOT_COMPARABLE_AFTER_GLOBAL_VMTK_REMESH_FACE_CLASSIFICATION_LOST"
                ),
                "boundarynormal_raw_P95_not_worse_than_old": (
                    boundarynormal[port_id]["normal_jump_P95_deg"]
                    <= old[port_id]["normal_jump_P95_deg"]
                ),
                "boundarynormal_raw_P99_not_worse_than_old": (
                    boundarynormal[port_id]["normal_jump_P99_deg"]
                    <= old[port_id]["normal_jump_P99_deg"]
                ),
            }
        )
        rows.append(row)
    not_worse = all(
        row["boundarynormal_raw_P95_not_worse_than_old"]
        and row["boundarynormal_raw_P99_not_worse_than_old"]
        for row in rows
    )
    return {
        "status": "PASS" if not_worse else "WARNING",
        "hard_gate": False,
        "visible_fold_assessment": "MANUAL_REVIEW_REQUIRED",
        "final_metric_limitation": (
            "RAW values are the comparable interface metrics because global VMTK "
            "remeshing does not preserve the core/extension face classification."
        ),
        "boundaries": rows,
    }


def _cut002_diagnosis(
    centerline_geometry: dict[str, Any], boundarynormal_geometry: dict[str, Any]
) -> dict[str, Any]:
    previous = next(
        row
        for row in centerline_geometry["boundaries"]
        if row["port_id"].endswith("__cut_002")
    )
    current = next(
        row
        for row in boundarynormal_geometry["boundaries"]
        if row["port_id"].endswith("__cut_002")
    )
    recovered = (
        current["extension_direction_dot"] >= 0.999
        and current["extension_length_relative_error"] <= 0.02
        and current["distal_area_relative_error"] <= 0.05
    )
    return {
        "port_id": current["port_id"],
        "previous_centerline_direction_dot": previous["extension_direction_dot"],
        "new_boundarynormal_direction_dot": current["extension_direction_dot"],
        "previous_axial_length_um": previous["actual_extension_length_um"],
        "new_axial_length_um": current["actual_axial_length_um"],
        "planned_length_um": current["planned_extension_length_um"],
        "previous_area_error": previous["distal_area_relative_error"],
        "new_area_error": current["distal_area_relative_error"],
        "previous_area_error_percent": 100.0
        * previous["distal_area_relative_error"],
        "new_area_error_percent": 100.0 * current["distal_area_relative_error"],
        "recovery_thresholds": {
            "minimum_signed_direction_dot": 0.999,
            "maximum_length_relative_error": 0.02,
            "maximum_area_relative_error": 0.05,
        },
        "recovered": recovered,
    }


def _write_report(path: Path, summary: dict[str, Any]) -> None:
    runtime = summary.get("runtime", {})
    geometry = summary.get("boundaries", [])
    lines = [
        "# VMTK TPS boundary-normal direction experiment",
        "",
        f"- Status: `{summary['status']}`",
        f"- VMTK runtime: `{runtime.get('vmtk_package_version', 'NOT_RUN')}`",
        f"- VTK runtime: `{runtime.get('vtk_version', 'NOT_RUN')}`",
        "- Official filter: `vtkvmtkPolyDataFlowExtensionsFilter`",
        "- Interpolation: `thinplatespline`",
        "- Custom TPS implementation: `false`",
        "- Previous extension mode: `centerlinedirection`",
        "- Current extension mode: `boundarynormal`",
        "- Official direction API: `SetExtensionModeToUseNormalToBoundary`",
        "- Centerlines used by VMTK: `false`",
        "- Transition ratio: `0.5`",
        "- Preserve cross-section shape: `false`",
        "- Parameter sweeps/fallbacks: `none`",
        "- Visible fold assessment: `MANUAL_REVIEW_REQUIRED`",
        "",
        "## Boundary measurements",
        "",
    ]
    for row in geometry:
        lines.append(
            f"- `{row['port_id']}`: axial length={row['actual_axial_length_um']:.6g} um; "
            f"area error={100.0 * row['distal_area_relative_error']:.4g}%; "
            f"direction dot={row['extension_direction_dot']:.9f}; "
            f"interface P95/P99={row['normal_jump_P95_deg']:.5g}/"
            f"{row['normal_jump_P99_deg']:.5g} deg"
        )
    lines.extend(
        [
            "",
            "## Final checks",
            "",
            f"- RAW topology: `{summary['raw_topology']['status']}`",
            f"- Remeshed/capped topology: `{summary['final_topology']['status']}`",
            f"- Radius P95: `{summary['radius_fidelity'].get('p95_absolute_relative_error', 'NOT_RUN')}`",
            f"- Core P95/max (um): `{summary['core_fidelity'].get('P95_core_surface_distance_um', 'NOT_RUN')}` / `{summary['core_fidelity'].get('max_core_surface_distance_um', 'NOT_RUN')}`",
            f"- Manual-review STL: `{summary['outputs'].get('manual_review_stl', 'NOT_GENERATED')}`",
            f"- Tagged VTP: `{summary['outputs'].get('tagged_vtp', 'NOT_GENERATED')}`",
            "- Volume mesh created: `false`",
            "- CFD run: `false`",
            "",
            f"NEXT: `{summary['next']}`",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def run_vmtk_surface_prepare(
    config: SurfacePrepareConfig,
    *,
    project_root: Path,
    run_id: str | None = None,
) -> VmtkSurfacePrepareResult:
    """Run exactly one official boundary-normal thin-plate-spline candidate."""

    if config.backend.method != "vmtk_flowextensions":
        raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED")
    if config.vmtk.extension_mode != "boundarynormal":
        raise SurfacePrepareError("INVALID_VMTK_EXTENSION_MODE")
    inputs = load_surface_inputs(
        config.paths.cfd_preprocess_run,
        expected_boundary_count=config.qc.expected_boundary_count,
    )
    anchor = _anchor_id(str(inputs.preprocess_summary["roi_id"]))
    candidate_id = run_id or (
        f"vmtk_tps_boundarynormal_anchor{anchor}_"
        f"{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    layout = _layout(config.paths.output_root, candidate_id)
    paths = exchange_paths(
        input_directory=layout.input,
        vmtk_directory=layout.vmtk,
        geometry_directory=layout.geometry,
        extension_mode=config.vmtk.extension_mode,
    )
    old_vtp, old_stl = _old_reference(config.manual_review.previous_surface_run)
    centerline_reference = _centerline_reference(config.paths.output_root)
    old_geometry_qc_path = (
        config.manual_review.previous_surface_run / "qc" / "boundary_geometry_qc.json"
    )
    if not old_geometry_qc_path.is_file():
        raise SurfacePrepareError(f"Missing OLD_CUSTOM_QC: {old_geometry_qc_path}")
    immutable_paths = (
        inputs.original_surface_um_stl,
        inputs.original_surface_um_vtp,
        old_stl,
        old_vtp,
        centerline_reference["raw_vtp"],
        centerline_reference["open_vtp"],
        centerline_reference["geometry_qc"],
        centerline_reference["interface_qc"],
    )
    hashes_before = {str(path): sha256_file(path) for path in immutable_paths}
    write_json(
        layout.input / "original_surface_reference.json",
        {
            "geometry_reference": inputs.geometry_reference,
            "immutable_hashes_before": hashes_before,
            "old_custom_extension_reference": str(
                config.manual_review.previous_surface_run
            ),
            "centerline_vmtk_failed_reference": str(centerline_reference["root"]),
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
        raise SurfacePrepareError("BOUNDARY_NORMAL_INPUT_PLANE_MISMATCH")
    if profile_report["status"] != "PASS":
        raise SurfacePrepareError("VMTK_BOUNDARY_NORMAL_RAW_GEOMETRY_FAILED")

    manifest_records = _boundary_manifest_records(inputs.boundaries)
    _write_boundary_manifest(
        layout.input / "boundary_manifest.csv", manifest_records, profile_report
    )
    write_json(
        layout.input / "centerline_usage.json",
        {
            "extension_mode": "boundarynormal",
            "centerline_adapter_generated": False,
            "centerlines_vtp_generated": False,
            "centerlines_used_for_extension_direction": False,
            "status": "NOT_GENERATED_BOUNDARY_NORMAL_MODE",
            "retained_adapter_module": "utils/cfd_surface_prepare/vmtk_adapter.py",
        },
    )
    target_edge = _global_median_edge_length(original)
    mapping = parameter_mapping(config.vmtk)
    mapping["global_original_median_edge_length_um"] = target_edge
    write_json(layout.vmtk / "parameters.json", mapping)
    write_json(layout.vmtk / "vmtk_parameter_mapping.json", mapping)

    invocation = run_official_vmtk(
        config=config.vmtk,
        paths=paths,
        tool_script=project_root / "tools" / "run_vmtk_flowextension.py",
    )
    runtime = {
        **invocation.runtime,
        "command": list(invocation.command),
        "configured_python_environment": config.vmtk.environment_python.parent.name,
        "vmtk_binary_runtime_prefix": str(config.vmtk.runtime_prefix),
        "runtime_overlay_is_process_local": True,
        "pmp_package_set_modified": False,
        "source_provenance": official_source_provenance(config.vmtk),
    }
    write_json(layout.vmtk / "environment.json", runtime)

    _, raw_mesh = polydata_mesh(paths.raw_vtp)
    raw_topology, raw_intersections = topology_qc(
        raw_mesh, expected_open_profile_count=4, allow_degenerate=True
    )
    raw_faces = np.asarray(raw_mesh.faces, dtype=np.int64)
    added_mask = np.any(raw_faces >= open_data.n_points, axis=1)
    collision_pairs = [
        [int(first), int(second)]
        for first, second in raw_intersections
        if added_mask[first] or added_mask[second]
    ]
    collision = {
        "status": "PASS" if not collision_pairs else "FAIL",
        "extension_collision_count": len(collision_pairs),
        "intersection_face_pairs": collision_pairs,
    }
    boundarynormal_interface = interface_smoothness_from_raw(
        raw_mesh,
        input_point_count=open_data.n_points,
        boundaries=inputs.boundaries,
    )
    geometry_report, distal_loops = extension_geometry_qc(
        raw_mesh, inputs.boundaries, proximal_loops, boundarynormal_interface
    )
    profile_by_id = {
        row["port_id"]: row for row in profile_report["boundaries"]
    }
    for row in geometry_report["boundaries"]:
        profile = profile_by_id[row["port_id"]]
        row["pre_cut_plane_abs_dot"] = profile[
            "boundary_plane_normal_abs_dot_expected_outward"
        ]
        row["pre_cut_center_distance_um"] = profile[
            "center_distance_to_expected_boundary_um"
        ]

    old_interface = interface_smoothness_from_old_custom(old_vtp, inputs.boundaries)
    centerline_open_data = pv.read(centerline_reference["open_vtp"]).triangulate()
    _, centerline_raw_mesh = polydata_mesh(centerline_reference["raw_vtp"])
    centerline_interface = interface_smoothness_from_raw(
        centerline_raw_mesh,
        input_point_count=centerline_open_data.n_points,
        boundaries=inputs.boundaries,
    )
    interface_comparison = _interface_three_way_comparison(
        old_interface, centerline_interface, boundarynormal_interface
    )
    centerline_geometry = read_json(centerline_reference["geometry_qc"])
    old_geometry = read_json(old_geometry_qc_path)
    direction_rows = _extension_direction_comparison(
        old_geometry, centerline_geometry, geometry_report
    )
    cut002 = _cut002_diagnosis(centerline_geometry, geometry_report)
    write_csv(layout.qc / "extension_direction_comparison.csv", direction_rows)
    write_csv(
        layout.qc / "interface_three_way_comparison.csv",
        interface_comparison["boundaries"],
    )
    write_json(layout.qc / "cut002_direction_diagnosis.json", cut002)
    write_json(layout.qc / "raw_surface_qc.json", raw_topology)
    write_json(layout.qc / "extension_collision_qc.json", collision)
    write_json(layout.qc / "extension_geometry_qc.json", geometry_report)
    write_json(
        layout.qc / "interface_smoothness_qc.json",
        {
            **interface_comparison,
            "old_custom": old_interface,
            "vmtk_centerline_raw": centerline_interface,
            "vmtk_boundarynormal_raw": boundarynormal_interface,
        },
    )
    _, old_mesh = polydata_mesh(old_vtp)
    old_metrics = extension_mesh_metrics(old_mesh, inputs.boundaries)
    centerline_metrics = extension_mesh_metrics(
        centerline_raw_mesh, inputs.boundaries
    )
    boundarynormal_metrics = extension_mesh_metrics(
        raw_mesh, inputs.boundaries, added_face_mask=added_mask
    )
    write_json(
        layout.qc / "three_way_raw_mesh_quality.json",
        {
            "old_custom": old_metrics,
            "vmtk_centerline_raw": centerline_metrics,
            "vmtk_boundarynormal_raw": boundarynormal_metrics,
        },
    )

    outputs: dict[str, Any] = {
        "raw_vtp": str(paths.raw_vtp.resolve()),
        "raw_stl": str(paths.raw_stl.resolve()),
        "extension_direction_comparison_csv": str(
            (layout.qc / "extension_direction_comparison.csv").resolve()
        ),
        "interface_three_way_comparison_csv": str(
            (layout.qc / "interface_three_way_comparison.csv").resolve()
        ),
        "cut002_direction_diagnosis_json": str(
            (layout.qc / "cut002_direction_diagnosis.json").resolve()
        ),
    }
    hashes_after_raw = {str(path): sha256_file(path) for path in immutable_paths}
    integrity = {
        "status": "PASS" if hashes_before == hashes_after_raw else "FAIL",
        "hashes_before": hashes_before,
        "hashes_after": hashes_after_raw,
        "original_ultraliser_stl_unchanged": hashes_before[
            str(inputs.original_surface_um_stl)
        ]
        == hashes_after_raw[str(inputs.original_surface_um_stl)],
        "original_ultraliser_vtp_unchanged": hashes_before[
            str(inputs.original_surface_um_vtp)
        ]
        == hashes_after_raw[str(inputs.original_surface_um_vtp)],
        "old_custom_reference_unchanged": hashes_before[str(old_stl)]
        == hashes_after_raw[str(old_stl)]
        and hashes_before[str(old_vtp)] == hashes_after_raw[str(old_vtp)],
        "centerline_vmtk_reference_unchanged": all(
            hashes_before[str(path)] == hashes_after_raw[str(path)]
            for path in (
                centerline_reference["raw_vtp"],
                centerline_reference["open_vtp"],
                centerline_reference["geometry_qc"],
                centerline_reference["interface_qc"],
            )
        ),
    }
    write_json(layout.qc / "original_surface_integrity.json", integrity)
    if integrity["status"] != "PASS":
        raise SurfacePrepareError("ORIGINAL_ULTRALISER_GEOMETRY_MODIFIED")

    boundary_rows = tuple(geometry_report["boundaries"])
    not_run = {"status": "NOT_RUN_RAW_HARD_GATE_FAILED"}
    raw_hard_pass = raw_geometry_hard_gate_pass(
        raw_topology, geometry_report, collision
    )
    if not should_promote_raw_candidate(raw_hard_pass):
        summary = {
            "status": RAW_FAILURE_STATUS,
            "next": FAIL_NEXT,
            "run_id": candidate_id,
            "run_root": str(layout.root),
            "backend": config.backend.method,
            "single_primary_candidate_count": 1,
            "second_candidate_run": False,
            "automatic_fallback": False,
            "custom_tps_implementation": False,
            "runtime": runtime,
            "boundaries": boundary_rows,
            "raw_topology": raw_topology,
            "remeshed_open_topology": not_run,
            "final_topology": not_run,
            "radius_fidelity": not_run,
            "core_fidelity": not_run,
            "collision": collision,
            "interface_comparison": interface_comparison,
            "cut002_recovery": cut002,
            "outputs": outputs,
            "figures": [],
            "original_integrity": integrity,
            "remesh_cap_promoted": False,
            "visible_fold_assessment": "MANUAL_REVIEW_REQUIRED",
            "volume_mesh_created": False,
            "cfd_run": False,
            "microbubble_simulation_run": False,
        }
        write_json(layout.qc / "run_summary.json", summary)
        _write_report(
            layout.report / "vmtk_tps_boundarynormal_report.md", summary
        )
        return VmtkSurfacePrepareResult(
            RAW_FAILURE_STATUS,
            FAIL_NEXT,
            layout.root,
            runtime,
            boundary_rows,
            raw_topology,
            not_run,
            not_run,
            not_run,
            int(collision["extension_collision_count"]),
            interface_comparison,
            cut002,
            outputs,
            (),
        )

    promotion = promote_official_vmtk(
        config=config.vmtk,
        paths=paths,
        tool_script=project_root / "tools" / "run_vmtk_flowextension.py",
        target_edge_length_um=target_edge,
    )
    runtime["promotion"] = promotion.runtime
    runtime["promotion_command"] = list(promotion.command)
    write_json(layout.vmtk / "environment.json", runtime)

    _, remeshed_open_mesh = polydata_mesh(paths.remeshed_open_vtp)
    remeshed_open_topology, _ = topology_qc(
        remeshed_open_mesh, expected_open_profile_count=4
    )
    write_json(layout.qc / "remeshed_open_surface_qc.json", remeshed_open_topology)
    final_outputs, boundary_mapping = tag_and_export_final_surface(
        paths.capped_vtp,
        inputs.boundaries,
        layout.geometry,
        layout.boundaries,
        output_stem="cfd_surface_vmtk_tps_boundarynormal",
    )
    outputs.update(final_outputs)
    write_json(layout.qc / "boundary_mapping_qc.json", boundary_mapping)

    if not should_generate_manual_review_figures(
        raw_hard_qc_pass=raw_hard_pass,
        interface_status=interface_comparison["status"],
    ):
        raise AssertionError("RAW PASS must always produce manual-review figures")
    figure_paths = save_vmtk_review_figures(
        old_custom_vtp=old_vtp,
        centerline_raw_vtp=centerline_reference["raw_vtp"],
        boundarynormal_raw_vtp=paths.raw_vtp,
        boundarynormal_remeshed_open_vtp=paths.remeshed_open_vtp,
        boundarynormal_final_vtp=Path(outputs["tagged_vtp"]),
        boundaries=inputs.boundaries,
        output_directory=layout.figures,
    )
    outputs["figure_paths"] = [str(path) for path in figure_paths]

    _, final_mesh = polydata_mesh(Path(outputs["tagged_vtp"]))
    final_topology, _ = topology_qc(final_mesh, expected_open_profile_count=0)
    core_report = core_surface_preservation_qc(
        original,
        final_mesh,
        inputs.boundaries,
        config.local_cut,
        config.qc,
    )
    try:
        radius_report = _radius_qc(
            final_mesh=final_mesh, inputs=inputs, project_root=project_root
        )
    except SurfacePrepareError as error:
        radius_report = {"status": "FAIL", "error": str(error)}
    scale_report = meter_scale_qc(
        Path(outputs["manual_review_stl"]), Path(outputs["meter_stl"])
    )
    write_json(layout.qc / "surface_qc.json", final_topology)
    write_json(layout.qc / "core_fidelity_qc.json", core_report)
    write_json(layout.qc / "radius_fidelity.json", radius_report)
    write_json(layout.qc / "meter_scale_qc.json", scale_report)

    hashes_after = {str(path): sha256_file(path) for path in immutable_paths}
    integrity["hashes_after"] = hashes_after
    integrity["status"] = "PASS" if hashes_before == hashes_after else "FAIL"
    write_json(layout.qc / "original_surface_integrity.json", integrity)
    final_hard_pass = all(
        report["status"] == "PASS"
        for report in (
            remeshed_open_topology,
            final_topology,
            core_report,
            radius_report,
            scale_report,
            integrity,
        )
    )

    if final_hard_pass:
        try:
            pressure_rows, pressure_report = geometry_pressure_correction(
                final_mesh,
                inputs.boundaries,
                proximal_loops,
                distal_loops,
                dynamic_viscosity_pa_s=float(
                    inputs.original_boundary_conditions["fluid"][
                        "dynamic_viscosity_pa_s"
                    ]
                ),
            )
            write_json(
                layout.qc / "extension_pressure_geometry_qc.json", pressure_report
            )
            csv_rows = [
                {
                    **row,
                    "station_fractions": json.dumps(row["station_fractions"]),
                    "cross_section_area_um2": json.dumps(
                        row["cross_section_area_um2"]
                    ),
                    "equivalent_radius_um": json.dumps(row["equivalent_radius_um"]),
                }
                for row in pressure_rows
            ]
            write_csv(
                layout.bc
                / "extension_pressure_correction_vmtk_boundarynormal.csv",
                csv_rows,
            )
            write_json(
                layout.bc / "boundary_conditions_vmtk_boundarynormal.json",
                _boundary_conditions(
                    inputs.original_boundary_conditions, pressure_rows
                ),
            )
        except SurfacePrepareError as error:
            final_hard_pass = False
            write_json(
                layout.qc / "extension_pressure_geometry_qc.json",
                {"status": "FAIL", "error": str(error)},
            )

    status = final_candidate_status(
        final_hard_qc_pass=final_hard_pass,
        interface_status=interface_comparison["status"],
    )
    next_stage = SUCCESS_NEXT if status in PASS_STATUSES else FAIL_NEXT
    summary = {
        "status": status,
        "next": next_stage,
        "run_id": candidate_id,
        "run_root": str(layout.root),
        "backend": config.backend.method,
        "single_primary_candidate_count": 1,
        "second_candidate_run": False,
        "automatic_fallback": False,
        "custom_tps_implementation": False,
        "runtime": runtime,
        "boundaries": boundary_rows,
        "raw_topology": raw_topology,
        "remeshed_open_topology": remeshed_open_topology,
        "final_topology": final_topology,
        "radius_fidelity": radius_report,
        "core_fidelity": core_report,
        "collision": collision,
        "interface_comparison": interface_comparison,
        "cut002_recovery": cut002,
        "outputs": outputs,
        "figures": [str(path) for path in figure_paths],
        "original_integrity": integrity,
        "remesh_cap_promoted": True,
        "visible_fold_assessment": "MANUAL_REVIEW_REQUIRED",
        "volume_mesh_created": False,
        "cfd_run": False,
        "microbubble_simulation_run": False,
    }
    write_json(layout.qc / "run_summary.json", summary)
    _write_report(layout.report / "vmtk_tps_boundarynormal_report.md", summary)
    return VmtkSurfacePrepareResult(
        status,
        next_stage,
        layout.root,
        runtime,
        boundary_rows,
        raw_topology,
        final_topology,
        radius_report,
        core_report,
        int(collision["extension_collision_count"]),
        interface_comparison,
        cut002,
        outputs,
        figure_paths,
    )


def print_vmtk_result(result: VmtkSurfacePrepareResult) -> None:
    print("EXPERIMENT: VMTK TPS DIRECTION SOURCE TEST")
    print("Previous mode: CENTERLINE_DIRECTION")
    print("Current mode: BOUNDARY_NORMAL")
    print(f"VMTK: {result.runtime['vmtk_package_version']}")
    print(f"VTK: {result.runtime['vtk_version']}")
    print("TPS: YES")
    print("Sigma: 1.0")
    print("Transition ratio: 0.5")
    print("Extension ratio: 10.0")
    print("Centerlines used by VMTK: NO")
    for row in result.boundaries:
        print(
            f"{row['port_id']}: "
            f"pre_cut_plane_abs_dot={row['pre_cut_plane_abs_dot']:.12g} "
            f"signed_direction_dot={row['extension_direction_dot']:.12g} "
            f"planned_length_um={row['planned_extension_length_um']:.12g} "
            f"actual_axial_length_um={row['actual_axial_length_um']:.12g} "
            f"actual_vector_norm_um={row['actual_center_to_center_norm_um']:.12g} "
            f"length_error_percent={100.0 * row['extension_length_relative_error']:.12g} "
            f"proximal_area_um2={row['proximal_area_um2']:.12g} "
            f"distal_area_um2={row['distal_area_um2']:.12g} "
            f"area_error_percent={100.0 * row['distal_area_relative_error']:.12g} "
            f"RAW_interface_P95={row['normal_jump_P95_deg']:.12g} "
            f"RAW_interface_P99={row['normal_jump_P99_deg']:.12g} "
            f"collision={result.collision_count}"
        )
    cut002 = result.cut002_diagnosis
    print("CUT_002 RECOVERY CHECK")
    print(
        "CENTERLINE_DIRECTION: "
        f"dot={cut002['previous_centerline_direction_dot']:.12g} "
        f"axial_length_um={cut002['previous_axial_length_um']:.12g} "
        f"area_error_percent={cut002['previous_area_error_percent']:.12g}"
    )
    print(
        "BOUNDARY_NORMAL: "
        f"dot={cut002['new_boundarynormal_direction_dot']:.12g} "
        f"axial_length_um={cut002['new_axial_length_um']:.12g} "
        f"area_error_percent={cut002['new_area_error_percent']:.12g}"
    )
    print(f"RECOVERED: {'YES' if cut002['recovered'] else 'NO'}")
    print(f"RAW topology: {result.raw_topology['status']}")
    print(f"Final topology: {result.final_topology['status']}")
    core_p95 = result.core_qc.get("P95_core_surface_distance_um", "NOT_RUN")
    core_max = result.core_qc.get("max_core_surface_distance_um", "NOT_RUN")
    radius_p95 = result.radius_qc.get(
        "p95_absolute_relative_error", "NOT_RUN"
    )
    print(f"Core P95/max: {core_p95} / {core_max} um")
    print(f"Radius P95: {radius_p95}")
    print(f"Collision count: {result.collision_count}")
    print(
        "MANUAL REVIEW STL: "
        f"{result.output_paths.get('manual_review_stl', 'NOT_GENERATED')}"
    )
    print(f"TAGGED VTP: {result.output_paths.get('tagged_vtp', 'NOT_GENERATED')}")
    figures = {path.name: path for path in result.figures}
    print(
        "THREE-WAY SURFACE FIGURE: "
        f"{figures.get('old_custom_vs_centerline_vs_boundarynormal_surface.png', 'NOT_GENERATED')}"
    )
    print(
        "THREE-WAY INTERFACE FIGURE: "
        f"{figures.get('three_way_interface_closeups.png', 'NOT_GENERATED')}"
    )
    print(
        "THREE-WAY WIREFRAME: "
        f"{figures.get('three_way_interface_wireframe.png', 'NOT_GENERATED')}"
    )
    print(
        "CUT002 COMPARISON: "
        f"{figures.get('cut002_centerline_vs_boundarynormal.png', 'NOT_GENERATED')}"
    )
    print("visible_fold_assessment: MANUAL_REVIEW_REQUIRED")
    print(f"Final status: {result.status}")
    print(f"NEXT: {result.next_stage}")
