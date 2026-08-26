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
    sha256_file,
)
from .local_cut import local_plane_cut
from .qc import core_surface_preservation_qc
from .types import TaggedSurface
from .vmtk_adapter import build_centerline_adapter
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
    run_official_vmtk,
)
from .vmtk_visualization import save_vmtk_review_figures


SUCCESS_STATUS = "VMTK_TPS_EXTENSION_PASS_PENDING_MANUAL_REVIEW"
SUCCESS_NEXT = "MANUALLY REVIEW VMTK TPS CFD SURFACE"
FAIL_NEXT = "REVIEW VMTK TPS FAILURE"


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


def _comparison(
    old_interface: dict[str, Any],
    new_interface: dict[str, Any],
    old_mesh: dict[str, Any],
    new_mesh: dict[str, Any],
) -> dict[str, Any]:
    old_i = {row["port_id"]: row for row in old_interface["boundaries"]}
    new_i = {row["port_id"]: row for row in new_interface["boundaries"]}
    old_m = {row["port_id"]: row for row in old_mesh["boundaries"]}
    new_m = {row["port_id"]: row for row in new_mesh["boundaries"]}
    rows: list[dict[str, Any]] = []
    for port_id in old_i:
        p95_improved = new_i[port_id]["normal_jump_P95_deg"] < old_i[port_id]["normal_jump_P95_deg"]
        p99_improved = new_i[port_id]["normal_jump_P99_deg"] < old_i[port_id]["normal_jump_P99_deg"]
        rows.append(
            {
                "port_id": port_id,
                "old_interface_P95_deg": old_i[port_id]["normal_jump_P95_deg"],
                "old_interface_P99_deg": old_i[port_id]["normal_jump_P99_deg"],
                "vmtk_interface_P95_deg": new_i[port_id]["normal_jump_P95_deg"],
                "vmtk_interface_P99_deg": new_i[port_id]["normal_jump_P99_deg"],
                "interface_P95_improved": p95_improved,
                "interface_P99_improved": p99_improved,
                "old_aspect_ratio_P95": old_m[port_id]["aspect_ratio_P95"],
                "vmtk_raw_aspect_ratio_P95": new_m[port_id]["aspect_ratio_P95"],
                "old_minimum_angle_deg": old_m[port_id]["minimum_angle_deg"],
                "vmtk_raw_minimum_angle_deg": new_m[port_id]["minimum_angle_deg"],
            }
        )
    return {
        "status": "PASS" if all(row["interface_P95_improved"] and row["interface_P99_improved"] for row in rows) else "FAIL",
        "visible_fold_assessment": "PENDING_MANUAL_REVIEW",
        "boundaries": rows,
    }


def _write_report(
    path: Path,
    *,
    runtime: dict[str, Any],
    geometry: dict[str, Any],
    raw_topology: dict[str, Any],
    final_topology: dict[str, Any],
    core: dict[str, Any],
    radius: dict[str, Any],
    outputs: dict[str, Any],
    status: str,
) -> None:
    lines = [
        "# VMTK TPS flow-extension validation",
        "",
        f"- Status: `{status}`",
        f"- VMTK runtime: `{runtime['vmtk_package_version']}`",
        f"- VTK runtime: `{runtime['vtk_version']}`",
        "- Official filter: `vtkvmtkPolyDataFlowExtensionsFilter`",
        "- Interpolation: `thinplatespline`",
        "- Custom TPS implementation: `false`",
        "- Extension mode: `centerlinedirection`",
        "- Transition ratio: `0.5`",
        "- Preserve cross-section shape: `false`",
        "- Parameter sweeps/fallbacks: `none`",
        "",
        "## Boundary measurements",
        "",
    ]
    for row in geometry["boundaries"]:
        lines.append(
            f"- `{row['port_id']}`: length={row['actual_extension_length_um']:.6g} um; "
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
            f"- RAW topology: `{raw_topology['status']}`",
            f"- Remeshed/capped topology: `{final_topology['status']}`",
            f"- Radius P95: `{radius['p95_absolute_relative_error']}`",
            f"- Core P95/max (um): `{core['P95_core_surface_distance_um']}` / `{core['max_core_surface_distance_um']}`",
            f"- Manual-review STL: `{outputs['manual_review_stl']}`",
            f"- Tagged VTP: `{outputs['tagged_vtp']}`",
            "- Volume mesh created: `false`",
            "- CFD run: `false`",
            "",
            f"NEXT: `{SUCCESS_NEXT if status == SUCCESS_STATUS else FAIL_NEXT}`",
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
    """Run one, and only one, official VMTK thin-plate-spline candidate."""

    if config.backend.method != "vmtk_flowextensions":
        raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED")
    inputs = load_surface_inputs(
        config.paths.cfd_preprocess_run,
        expected_boundary_count=config.qc.expected_boundary_count,
    )
    anchor = _anchor_id(str(inputs.preprocess_summary["roi_id"]))
    candidate_id = run_id or f"vmtk_tps_anchor{anchor}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    layout = _layout(config.paths.output_root, candidate_id)
    paths = exchange_paths(
        input_directory=layout.input,
        vmtk_directory=layout.vmtk,
        geometry_directory=layout.geometry,
    )
    old_vtp, old_stl = _old_reference(config.manual_review.previous_surface_run)
    immutable_paths = (
        inputs.original_surface_um_stl,
        inputs.original_surface_um_vtp,
        old_stl,
        old_vtp,
    )
    hashes_before = {str(path): sha256_file(path) for path in immutable_paths}
    write_json(
        layout.input / "original_surface_reference.json",
        {
            "geometry_reference": inputs.geometry_reference,
            "immutable_hashes_before": hashes_before,
            "old_custom_extension_reference": str(config.manual_review.previous_surface_run),
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
    if profile_report["status"] != "PASS":
        raise SurfacePrepareError("VMTK_EXTENSION_GEOMETRY_FAILED:open_profiles")

    adapter = build_centerline_adapter(
        inputs.preprocess_run, inputs.boundaries, paths.centerlines_vtp
    )
    _write_boundary_manifest(layout.input / "boundary_manifest.csv", adapter.records, profile_report)
    target_edge = _global_median_edge_length(original)
    mapping = parameter_mapping(config.vmtk)
    mapping["global_original_median_edge_length_um"] = target_edge
    write_json(layout.vmtk / "parameters.json", mapping)
    write_json(layout.vmtk / "vmtk_parameter_mapping.json", mapping)

    invocation = run_official_vmtk(
        config=config.vmtk,
        paths=paths,
        tool_script=project_root / "tools" / "run_vmtk_flowextension.py",
        target_edge_length_um=target_edge,
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
    new_interface = interface_smoothness_from_raw(
        raw_mesh,
        input_point_count=open_data.n_points,
        boundaries=inputs.boundaries,
    )
    geometry_report, distal_loops = extension_geometry_qc(
        raw_mesh, inputs.boundaries, proximal_loops, new_interface
    )
    old_interface = interface_smoothness_from_old_custom(old_vtp, inputs.boundaries)
    old_data, old_mesh = polydata_mesh(old_vtp)
    old_metrics = extension_mesh_metrics(old_mesh, inputs.boundaries)
    new_metrics = extension_mesh_metrics(
        raw_mesh, inputs.boundaries, added_face_mask=added_mask
    )
    comparison = _comparison(old_interface, new_interface, old_metrics, new_metrics)
    interface_report = {
        "status": comparison["status"],
        "old_custom": old_interface,
        "vmtk_tps": new_interface,
        "comparison": comparison,
    }
    write_json(layout.qc / "raw_surface_qc.json", raw_topology)
    write_json(layout.qc / "extension_collision_qc.json", collision)
    write_json(layout.qc / "extension_geometry_qc.json", geometry_report)
    write_json(layout.qc / "interface_smoothness_qc.json", interface_report)
    write_json(
        layout.qc / "old_custom_vs_vmtk_mesh_quality.json",
        {"old_custom": old_metrics, "vmtk_tps_raw": new_metrics},
    )
    if raw_topology["status"] != "PASS" or geometry_report["status"] != "PASS" or collision["status"] != "PASS":
        raise SurfacePrepareError("VMTK_EXTENSION_GEOMETRY_FAILED")
    if comparison["status"] != "PASS":
        raise SurfacePrepareError("VMTK_TPS_EXTENSION_FAILED")

    _, remeshed_open_mesh = polydata_mesh(paths.remeshed_open_vtp)
    remeshed_open_topology, _ = topology_qc(
        remeshed_open_mesh, expected_open_profile_count=4
    )
    write_json(layout.qc / "remeshed_open_surface_qc.json", remeshed_open_topology)
    outputs, boundary_mapping = tag_and_export_final_surface(
        paths.capped_vtp, inputs.boundaries, layout.geometry, layout.boundaries
    )
    write_json(layout.qc / "boundary_mapping_qc.json", boundary_mapping)
    _, final_mesh = polydata_mesh(Path(outputs["tagged_vtp"]))
    final_topology, _ = topology_qc(final_mesh, expected_open_profile_count=0)
    core_report = core_surface_preservation_qc(
        original,
        final_mesh,
        inputs.boundaries,
        config.local_cut,
        config.qc,
    )
    radius_report = _radius_qc(
        final_mesh=final_mesh, inputs=inputs, project_root=project_root
    )
    scale_report = meter_scale_qc(
        Path(outputs["manual_review_stl"]), Path(outputs["meter_stl"])
    )
    write_json(layout.qc / "surface_qc.json", final_topology)
    write_json(layout.qc / "core_fidelity_qc.json", core_report)
    write_json(layout.qc / "radius_fidelity.json", radius_report)
    write_json(layout.qc / "meter_scale_qc.json", scale_report)
    if core_report["status"] != "PASS":
        raise SurfacePrepareError("VMTK_REMESH_CORE_FIDELITY_FAILED")
    if final_topology["status"] != "PASS" or radius_report["status"] != "PASS" or scale_report["status"] != "PASS":
        raise SurfacePrepareError("VMTK_SURFACE_QC_FAILED")

    pressure_rows, pressure_report = geometry_pressure_correction(
        final_mesh,
        inputs.boundaries,
        proximal_loops,
        distal_loops,
        dynamic_viscosity_pa_s=float(
            inputs.original_boundary_conditions["fluid"]["dynamic_viscosity_pa_s"]
        ),
    )
    write_json(layout.qc / "extension_pressure_geometry_qc.json", pressure_report)
    csv_rows = [
        {
            **row,
            "station_fractions": json.dumps(row["station_fractions"]),
            "cross_section_area_um2": json.dumps(row["cross_section_area_um2"]),
            "equivalent_radius_um": json.dumps(row["equivalent_radius_um"]),
        }
        for row in pressure_rows
    ]
    write_csv(layout.bc / "extension_pressure_correction_vmtk.csv", csv_rows)
    write_json(
        layout.bc / "boundary_conditions_vmtk.json",
        _boundary_conditions(inputs.original_boundary_conditions, pressure_rows),
    )

    figure_paths = save_vmtk_review_figures(
        old_custom_vtp=old_vtp,
        raw_vtp=paths.raw_vtp,
        final_vtp=Path(outputs["tagged_vtp"]),
        boundaries=inputs.boundaries,
        output_directory=layout.figures,
    )
    hashes_after = {str(path): sha256_file(path) for path in immutable_paths}
    integrity = {
        "status": "PASS" if hashes_before == hashes_after else "FAIL",
        "hashes_before": hashes_before,
        "hashes_after": hashes_after,
        "original_ultraliser_stl_unchanged": hashes_before[str(inputs.original_surface_um_stl)] == hashes_after[str(inputs.original_surface_um_stl)],
        "original_ultraliser_vtp_unchanged": hashes_before[str(inputs.original_surface_um_vtp)] == hashes_after[str(inputs.original_surface_um_vtp)],
        "old_custom_reference_unchanged": hashes_before[str(old_stl)] == hashes_after[str(old_stl)] and hashes_before[str(old_vtp)] == hashes_after[str(old_vtp)],
    }
    write_json(layout.qc / "original_surface_integrity.json", integrity)
    if integrity["status"] != "PASS":
        raise SurfacePrepareError("ORIGINAL_ULTRALISER_GEOMETRY_MODIFIED")

    boundary_rows = tuple(geometry_report["boundaries"])
    summary = {
        "status": SUCCESS_STATUS,
        "next": SUCCESS_NEXT,
        "run_id": candidate_id,
        "run_root": str(layout.root),
        "backend": config.backend.method,
        "single_primary_candidate_count": 1,
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
        "interface_comparison": comparison,
        "outputs": outputs,
        "figures": [str(path) for path in figure_paths],
        "original_integrity": integrity,
        "volume_mesh_created": False,
        "cfd_run": False,
        "microbubble_simulation_run": False,
    }
    write_json(layout.qc / "run_summary.json", summary)
    _write_report(
        layout.report / "vmtk_tps_extension_report.md",
        runtime=runtime,
        geometry=geometry_report,
        raw_topology=raw_topology,
        final_topology=final_topology,
        core=core_report,
        radius=radius_report,
        outputs=outputs,
        status=SUCCESS_STATUS,
    )
    return VmtkSurfacePrepareResult(
        SUCCESS_STATUS,
        SUCCESS_NEXT,
        layout.root,
        runtime,
        boundary_rows,
        raw_topology,
        final_topology,
        radius_report,
        core_report,
        int(collision["extension_collision_count"]),
        outputs,
        figure_paths,
    )


def print_vmtk_result(result: VmtkSurfacePrepareResult) -> None:
    print(f"VMTK version: {result.runtime['vmtk_package_version']}")
    print(f"VMTK environment: {result.runtime['python_executable']}")
    print("Interpolation: THIN_PLATE_SPLINE")
    print("Transition ratio: 0.5")
    print("Extension mode: CENTERLINE_DIRECTION")
    print("Preserve shape: FALSE")
    for row in result.boundaries:
        print(
            f"{row['port_id']}: length={row['actual_extension_length_um']:.12g} "
            f"area_error={row['distal_area_relative_error']:.12g} "
            f"direction_dot={row['extension_direction_dot']:.12g} "
            f"interface_P95={row['normal_jump_P95_deg']:.12g} "
            f"interface_P99={row['normal_jump_P99_deg']:.12g}"
        )
    print(f"RAW topology: {result.raw_topology['status']}")
    print(f"Remeshed topology: {result.final_topology['status']}")
    print(f"Radius P95: {result.radius_qc['p95_absolute_relative_error']:.12g}")
    print(
        f"Core P95/max: {result.core_qc['P95_core_surface_distance_um']:.12g} / "
        f"{result.core_qc['max_core_surface_distance_um']:.12g} um"
    )
    print(f"Collision count: {result.collision_count}")
    print(f"MANUAL REVIEW STL PATH: {result.output_paths['manual_review_stl']}")
    print(f"TAGGED VTP PATH: {result.output_paths['tagged_vtp']}")
    print(f"Final status: {result.status}")
    print(f"NEXT: {result.next_stage}")
