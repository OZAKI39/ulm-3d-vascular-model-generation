"""Formal single-candidate VMTK TPS boundary-normal cap-only pipeline."""

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
    core_exact_preservation_qc,
    core_symmetric_distance_qc,
    extension_geometry_qc,
    extension_mesh_quality_from_raw,
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
    exchange_paths,
    official_source_provenance,
    parameter_mapping,
    run_official_vmtk,
)
from .vmtk_visualization import save_caponly_review_figures


RAW_GEOMETRY_FAILURE = "VMTK_BOUNDARY_NORMAL_RAW_GEOMETRY_FAILED"
RAW_CORE_FAILURE = "VMTK_RAW_CORE_NOT_EXACT_COPY"
RAW_MESH_FAILURE = "VMTK_RAW_EXTENSION_MESH_QUALITY_FAILED"
CAP_FAILURE = "VMTK_RAW_DIRECT_CAP_FAILED"
TOPOLOGY_FAILURE = "VMTK_CAPONLY_TOPOLOGY_FAILED"
CORE_FAILURE = "VMTK_CAPONLY_CORE_PRESERVATION_FAILED"
RADIUS_FAILURE = "VMTK_CAPONLY_RADIUS_FIDELITY_FAILED"
BOUNDARY_FAILURE = "VMTK_CAPONLY_BOUNDARY_MAPPING_FAILED"
SUCCESS_STATUS = "VMTK_TPS_BOUNDARY_NORMAL_CAPONLY_PASS_PENDING_MANUAL_REVIEW"
PASS_STATUSES = {SUCCESS_STATUS}
SUCCESS_NEXT = "MANUALLY REVIEW CAP-ONLY CFD SURFACE"
FAIL_NEXT = "REVIEW CAP-ONLY SURFACE FAILURE"
PREVIOUS_GLOBAL_REMESH_RUN_ID = "vmtk_tps_boundarynormal_anchor003274_20260826_153709"


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
    raw_core_qc: dict[str, Any]
    raw_topology: dict[str, Any]
    extension_mesh_quality: dict[str, Any]
    final_topology: dict[str, Any]
    core_exact_qc: dict[str, Any]
    core_distance_qc: dict[str, Any]
    normal_qc: dict[str, Any]
    previous_global_qc: dict[str, Any]
    artifact_confirmed: bool
    radius_qc: dict[str, Any]
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


def _previous_reference(output_root: Path) -> dict[str, Path]:
    root = Path(output_root) / PREVIOUS_GLOBAL_REMESH_RUN_ID
    paths = {
        "root": root,
        "open_vtp": root / "input" / "open_surface_um.vtp",
        "raw_vtp": root / "geometry" / "vmtk_boundarynormal_raw_um.vtp",
        "final_vtp": root / "geometry" / "cfd_surface_vmtk_tps_boundarynormal_um.vtp",
        "radius_qc": root / "qc" / "radius_fidelity.json",
        "summary": root / "qc" / "run_summary.json",
    }
    missing = [str(path) for key, path in paths.items() if key != "root" and not path.is_file()]
    if missing:
        raise SurfacePrepareError(
            "Missing PREVIOUS_GLOBAL_REMESH_REFERENCE: " + ", ".join(missing)
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
                "predicted_extension_pressure_drop_pa": row[
                    "predicted_extension_pressure_drop_pa"
                ],
                "P_solver_boundary_pa": row["P_solver_boundary_pa"],
                "extension_pressure_correction_role": row["pressure_correction_role"],
            }
        )
    output["surface_geometry_method"] = "OFFICIAL_VMTK_TPS_FLOWEXTENSIONS_DIRECT_CAP_ONLY"
    output["pressure_correction_geometry"] = (
        "20 cross sections on final VMTK direct-cap geometry"
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


def _write_report(path: Path, summary: dict[str, Any]) -> None:
    runtime = summary.get("runtime", {})
    lines = [
        "# VMTK TPS boundary-normal direct-cap experiment",
        "",
        f"- Status: `{summary['status']}`",
        "- Extension: `VMTK TPS BOUNDARY_NORMAL`",
        "- Postprocess: `DIRECT VMTK CAP ONLY`",
        "- Global surface remeshing performed: `false`",
        "- VMTK surface remesher called: `false`",
        "- Official capper: `vmtkscripts.vmtkSurfaceCapper`",
        f"- VMTK: `{runtime.get('vmtk_package_version', 'NOT_RUN')}`",
        f"- VTK: `{runtime.get('vtk_version', 'NOT_RUN')}`",
        f"- Global-remesh artifact confirmed: `{summary.get('GLOBAL_REMESH_ARTIFACT_CONFIRMED', False)}`",
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
        "vmtk_surface_remeshing_called": False,
        "volume_mesh_created": False,
        "cfd_run": False,
        "microbubble_simulation_run": False,
    }
    write_json(layout.qc / "run_summary.json", summary)
    _write_report(layout.report / "vmtk_tps_boundarynormal_caponly_report.md", summary)


def run_vmtk_surface_prepare(
    config: SurfacePrepareConfig,
    *,
    project_root: Path,
    run_id: str | None = None,
) -> VmtkSurfacePrepareResult:
    """Run one official boundary-normal TPS candidate and cap it directly."""

    if config.backend.method != "vmtk_flowextensions":
        raise SurfacePrepareError("VMTK_ENVIRONMENT_BLOCKED")
    if config.vmtk.extension_mode != "boundarynormal":
        raise SurfacePrepareError("INVALID_VMTK_EXTENSION_MODE")
    if config.vmtk.postprocess_mode != "cap_only" or config.vmtk.remesh_after_extension:
        raise SurfacePrepareError("INVALID_VMTK_POSTPROCESS_CONFIGURATION")

    inputs = load_surface_inputs(
        config.paths.cfd_preprocess_run,
        expected_boundary_count=config.qc.expected_boundary_count,
    )
    anchor = _anchor_id(str(inputs.preprocess_summary["roi_id"]))
    candidate_id = run_id or (
        f"vmtk_tps_boundarynormal_caponly_anchor{anchor}_"
        f"{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    layout = _layout(config.paths.output_root, candidate_id)
    paths = exchange_paths(
        input_directory=layout.input,
        vmtk_directory=layout.vmtk,
        geometry_directory=layout.geometry,
        extension_mode=config.vmtk.extension_mode,
    )
    previous = _previous_reference(config.paths.output_root)
    immutable_paths = (
        inputs.original_surface_um_stl,
        inputs.original_surface_um_vtp,
        previous["open_vtp"],
        previous["raw_vtp"],
        previous["final_vtp"],
        previous["radius_qc"],
        previous["summary"],
    )
    hashes_before = {str(path): sha256_file(path) for path in immutable_paths}
    write_json(
        layout.input / "original_surface_reference.json",
        {
            "geometry_reference": inputs.geometry_reference,
            "immutable_hashes_before": hashes_before,
            "previous_global_remesh_reference": str(previous["root"]),
            "previous_reference_role": "DIAGNOSTIC_REFERENCE_ONLY",
        },
    )
    shutil.copy2(config.source_path, layout.input / "cfd_surface_prepare.yaml")
    shutil.copy2(
        inputs.preprocess_run / "roi" / "boundary_conditions.json",
        layout.bc / "boundary_conditions_original.json",
    )

    original = load_original_surface(inputs.original_surface_um_vtp)
    previous_raw_exact = raw_core_exact_copy_qc(
        previous["open_vtp"], previous["raw_vtp"], tag_regions=False
    )
    previous_diagnostic = previous_global_remesh_diagnostics(
        original,
        previous["final_vtp"],
        inputs.boundaries,
        config.local_cut,
    )
    write_json(layout.qc / "previous_global_remesh_core_qc.json", previous_diagnostic)
    write_csv(
        layout.qc / "previous_global_remesh_hotspots.csv",
        previous_diagnostic["hotspots"],
    )

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
        raise SurfacePrepareError(RAW_GEOMETRY_FAILURE)

    manifest = _manifest_records(inputs.boundaries)
    _write_manifest(layout.input / "boundary_manifest.csv", manifest, profile_report)
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
        "global_surface_remeshing_performed": False,
        "vmtk_surface_remeshing_called": False,
    }
    write_json(layout.vmtk / "environment.json", runtime)
    outputs: dict[str, Any] = {
        "raw_vtp": str(paths.raw_vtp.resolve()),
        "raw_stl": str(paths.raw_stl.resolve()),
    }

    raw_core = raw_core_exact_copy_qc(paths.open_surface_vtp, paths.raw_vtp)
    write_json(layout.qc / "raw_core_exact_copy_qc.json", raw_core)
    if raw_core["status"] != "PASS":
        _write_failure(
            layout,
            status=RAW_CORE_FAILURE,
            run_id=candidate_id,
            runtime=runtime,
            outputs=outputs,
            evidence={"raw_core_exact_copy": raw_core},
        )
        raise SurfacePrepareError(RAW_CORE_FAILURE)

    _, raw_mesh = polydata_mesh(paths.raw_vtp)
    raw_topology, raw_intersections = topology_qc(
        raw_mesh, expected_open_profile_count=4, allow_degenerate=False
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
    interface = interface_smoothness_from_raw(
        raw_mesh,
        input_point_count=open_data.n_points,
        boundaries=inputs.boundaries,
    )
    geometry_report, distal_loops = extension_geometry_qc(
        raw_mesh, inputs.boundaries, proximal_loops, interface
    )
    profile_by_id = {row["port_id"]: row for row in profile_report["boundaries"]}
    for row in geometry_report["boundaries"]:
        profile = profile_by_id[row["port_id"]]
        row["pre_cut_plane_abs_dot"] = profile[
            "boundary_plane_normal_abs_dot_expected_outward"
        ]
        row["pre_cut_center_distance_um"] = profile[
            "center_distance_to_expected_boundary_um"
        ]
    extension_quality = extension_mesh_quality_from_raw(
        paths.raw_vtp,
        inputs.boundaries,
        local_target_edge_um=local_targets,
        quality=config.mesh_quality,
        input_point_count=open_data.n_points,
    )
    write_json(layout.qc / "raw_surface_qc.json", raw_topology)
    write_json(layout.qc / "extension_geometry_qc.json", geometry_report)
    write_json(layout.qc / "extension_mesh_quality_qc.json", extension_quality)
    write_json(layout.qc / "extension_collision_qc.json", collision)
    write_json(
        layout.qc / "interface_smoothness_qc.json",
        {
            "status": "DIAGNOSTIC_ONLY",
            "hard_gate": False,
            "vmtk_boundarynormal_raw": interface,
        },
    )
    raw_hard_pass = raw_geometry_hard_gate_pass(
        raw_topology, geometry_report, extension_quality, collision, raw_core
    )
    if not raw_hard_pass:
        status = (
            RAW_MESH_FAILURE
            if extension_quality["status"] != "PASS"
            else RAW_GEOMETRY_FAILURE
        )
        _write_failure(
            layout,
            status=status,
            run_id=candidate_id,
            runtime=runtime,
            outputs=outputs,
            evidence={
                "raw_topology": raw_topology,
                "extension_geometry": geometry_report,
                "extension_mesh_quality": extension_quality,
                "collision": collision,
            },
        )
        raise SurfacePrepareError(status)

    try:
        cap = cap_official_vmtk(
            config=config.vmtk,
            paths=paths,
            tool_script=project_root / "tools" / "run_vmtk_flowextension.py",
        )
    except SurfacePrepareError:
        _write_failure(
            layout,
            status=CAP_FAILURE,
            run_id=candidate_id,
            runtime=runtime,
            outputs=outputs,
            evidence={"raw_hard_gate": "PASS"},
        )
        raise
    runtime["cap_only"] = cap.runtime
    runtime["cap_only_command"] = list(cap.command)
    runtime["global_surface_remeshing_performed"] = False
    runtime["vmtk_surface_remeshing_called"] = False
    write_json(layout.vmtk / "environment.json", runtime)

    try:
        final_outputs, boundary_mapping = tag_and_export_final_surface(
            paths.capped_vtp,
            inputs.boundaries,
            layout.geometry,
            layout.boundaries,
            output_stem="cfd_surface_vmtk_tps_boundarynormal_caponly",
            raw_vtp=paths.raw_vtp,
        )
    except SurfacePrepareError:
        _write_failure(
            layout,
            status=CORE_FAILURE,
            run_id=candidate_id,
            runtime=runtime,
            outputs=outputs,
            evidence={"cap_runtime": cap.runtime},
        )
        raise
    outputs.update(final_outputs)
    write_json(layout.qc / "boundary_mapping_qc.json", boundary_mapping)

    _, final_mesh = polydata_mesh(Path(outputs["tagged_vtp"]))
    final_topology, _ = topology_qc(
        final_mesh,
        expected_open_profile_count=0,
        require_winding_consistent=True,
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
    scale_report = meter_scale_qc(
        Path(outputs["manual_review_stl"]), Path(outputs["meter_stl"])
    )
    write_json(layout.qc / "caponly_surface_qc.json", final_topology)
    write_json(layout.qc / "core_exact_preservation_qc.json", core_exact)
    write_json(layout.qc / "core_symmetric_distance_qc.json", core_distance)
    write_json(layout.qc / "normal_consistency_qc.json", normal_report)
    write_json(layout.qc / "radius_fidelity.json", radius_report)
    write_json(layout.qc / "meter_scale_qc.json", scale_report)

    pressure_rows: list[dict[str, Any]] = []
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
            layout.bc
            / "extension_pressure_correction_vmtk_boundarynormal_caponly.csv"
        )
        write_csv(pressure_csv, csv_rows)
        boundary_json = (
            layout.bc / "boundary_conditions_vmtk_boundarynormal_caponly.json"
        )
        write_json(
            boundary_json,
            _boundary_conditions(inputs.original_boundary_conditions, pressure_rows),
        )
        outputs["pressure_correction_csv"] = str(pressure_csv.resolve())
        outputs["boundary_conditions_json"] = str(boundary_json.resolve())
    write_json(layout.qc / "extension_pressure_geometry_qc.json", pressure_report)

    previous_radius = read_json(previous["radius_qc"])
    comparison_rows = [
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
            "core_normal_P95_deg": previous_diagnostic["core_normal_deviation"][
                "P95_deg"
            ],
            "core_normal_P99_deg": previous_diagnostic["core_normal_deviation"][
                "P99_deg"
            ],
            "core_normal_max_deg": previous_diagnostic["core_normal_deviation"][
                "max_deg"
            ],
            "winding_consistent": previous_diagnostic["winding_consistent"],
            "radius_P95_error": previous_radius["p95_absolute_relative_error"],
        },
        {
            "method": "NEW_CAP_ONLY",
            "core_original_to_final_P95_um": core_distance[
                "original_core_to_caponly_final_core"
            ]["P95_um"],
            "core_original_to_final_max_um": core_distance[
                "original_core_to_caponly_final_core"
            ]["max_um"],
            "core_final_to_original_P95_um": core_distance[
                "caponly_final_core_to_original_core"
            ]["P95_um"],
            "core_final_to_original_max_um": core_distance[
                "caponly_final_core_to_original_core"
            ]["max_um"],
            "core_normal_P95_deg": normal_report[
                "core_normal_deviation_P95_deg"
            ],
            "core_normal_P99_deg": normal_report[
                "core_normal_deviation_P99_deg"
            ],
            "core_normal_max_deg": normal_report[
                "core_normal_deviation_max_deg"
            ],
            "winding_consistent": normal_report["winding_consistent"],
            "radius_P95_error": radius_report.get("p95_absolute_relative_error"),
        },
    ]
    comparison_csv = layout.qc / "global_remesh_vs_caponly_core_comparison.csv"
    write_csv(comparison_csv, comparison_rows)
    outputs["global_remesh_vs_caponly_core_comparison_csv"] = str(
        comparison_csv.resolve()
    )
    artifact_confirmed = bool(
        previous_raw_exact["status"] == "PASS"
        and previous_diagnostic["global_remesh_artifact_detected"]
        and core_exact["status"] == "PASS"
        and core_distance["status"] == "PASS"
        and normal_report["status"] == "PASS"
    )

    figure_paths = save_caponly_review_figures(
        original_vtp=inputs.original_surface_um_vtp,
        raw_vtp=paths.raw_vtp,
        previous_global_final_vtp=previous["final_vtp"],
        caponly_final_vtp=Path(outputs["tagged_vtp"]),
        boundaries=inputs.boundaries,
        hotspots=previous_diagnostic["hotspots"],
        output_directory=layout.figures,
    )
    outputs["figure_paths"] = [str(path) for path in figure_paths]
    outputs["previous_global_remesh_artifact_hotspots_figure"] = str(
        (layout.figures / "previous_global_remesh_artifact_hotspots.png").resolve()
    )

    hashes_after = {str(path): sha256_file(path) for path in immutable_paths}
    integrity = {
        "status": "PASS" if hashes_before == hashes_after else "FAIL",
        "hashes_before": hashes_before,
        "hashes_after": hashes_after,
        "previous_global_remesh_reference_untouched": all(
            hashes_before[str(path)] == hashes_after[str(path)]
            for path in (
                previous["open_vtp"],
                previous["raw_vtp"],
                previous["final_vtp"],
                previous["radius_qc"],
                previous["summary"],
            )
        ),
    }
    write_json(layout.qc / "original_surface_integrity.json", integrity)
    reports = (
        final_topology,
        boundary_mapping,
        core_exact,
        core_distance,
        normal_report,
        radius_report,
        scale_report,
        pressure_report,
        integrity,
    )
    final_hard_pass = all(report.get("status") == "PASS" for report in reports)
    if final_hard_pass:
        status = SUCCESS_STATUS
    elif boundary_mapping.get("status") != "PASS":
        status = BOUNDARY_FAILURE
    elif radius_report.get("status") != "PASS":
        status = RADIUS_FAILURE
    elif core_exact.get("status") != "PASS" or core_distance.get("status") != "PASS":
        status = CORE_FAILURE
    else:
        status = TOPOLOGY_FAILURE
    next_stage = SUCCESS_NEXT if final_hard_pass else FAIL_NEXT
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
        "boundaries": geometry_report["boundaries"],
        "raw_core_exact_copy": raw_core,
        "raw_topology": raw_topology,
        "extension_mesh_quality": extension_quality,
        "final_topology": final_topology,
        "core_exact_preservation": core_exact,
        "core_symmetric_distance": core_distance,
        "normal_consistency": normal_report,
        "previous_global_remesh": previous_diagnostic,
        "GLOBAL_REMESH_ARTIFACT_CONFIRMED": artifact_confirmed,
        "radius_fidelity": radius_report,
        "collision": collision,
        "boundary_mapping": boundary_mapping,
        "pressure_geometry": pressure_report,
        "outputs": outputs,
        "figures": [str(path) for path in figure_paths],
        "original_integrity": integrity,
        "global_surface_remeshing_performed": False,
        "vmtk_surface_remeshing_called": False,
        "visible_acceptance": "MANUAL_REVIEW_REQUIRED",
        "volume_mesh_created": False,
        "cfd_run": False,
        "microbubble_simulation_run": False,
    }
    write_json(layout.qc / "run_summary.json", summary)
    _write_report(layout.report / "vmtk_tps_boundarynormal_caponly_report.md", summary)
    return VmtkSurfacePrepareResult(
        status=status,
        next_stage=next_stage,
        run_root=layout.root,
        runtime=runtime,
        boundaries=tuple(geometry_report["boundaries"]),
        raw_core_qc=raw_core,
        raw_topology=raw_topology,
        extension_mesh_quality=extension_quality,
        final_topology=final_topology,
        core_exact_qc=core_exact,
        core_distance_qc=core_distance,
        normal_qc=normal_report,
        previous_global_qc=previous_diagnostic,
        artifact_confirmed=artifact_confirmed,
        radius_qc=radius_report,
        collision_count=int(collision["extension_collision_count"]),
        pressure_rows=tuple(pressure_rows),
        output_paths=outputs,
        figures=figure_paths,
    )


def print_vmtk_result(result: VmtkSurfacePrepareResult) -> None:
    print("EXPERIMENT: REMOVE GLOBAL VMTK SURFACE REMESHING")
    print("Extension: VMTK TPS BOUNDARY_NORMAL")
    print("Postprocess: DIRECT VMTK CAP ONLY")
    print("Global remesh: NO")
    print(f"VMTK: {result.runtime['vmtk_package_version']}")
    print(f"VTK: {result.runtime['vtk_version']}")
    raw_core = result.raw_core_qc
    print("RAW CORE EXACT COPY:")
    print(
        "retained point max motion = "
        f"{raw_core['retained_input_point_max_motion_um']:.12g} um"
    )
    print(
        "connectivity changed count = "
        f"{raw_core['original_input_cell_connectivity_changed_count']}"
    )
    print(raw_core["status"])
    quality = {
        row["port_id"]: row for row in result.extension_mesh_quality["boundaries"]
    }
    for row in result.boundaries:
        mesh = quality[row["port_id"]]
        print(
            f"{row['port_id']}: direction_dot={row['extension_direction_dot']:.12g} "
            f"planned/actual_um={row['planned_extension_length_um']:.12g}/"
            f"{row['actual_axial_length_um']:.12g} "
            f"length_error_percent={100.0 * row['extension_length_relative_error']:.12g} "
            f"area_error_percent={100.0 * row['distal_area_relative_error']:.12g} "
            f"RAW_extension_aspect_P95={mesh['aspect_ratio_P95']:.12g} "
            f"RAW_minimum_angle={mesh['minimum_angle_deg']:.12g} "
            f"collision={result.collision_count}"
        )
    core = result.core_exact_qc
    distance = result.core_distance_qc
    normal = result.normal_qc
    print("CAP-ONLY CORE:")
    print(
        "retained original vertex max motion = "
        f"{core['retained_original_vertex_max_motion_um']:.12g} um"
    )
    print(
        "connectivity changed count = "
        f"{core['core_triangle_connectivity_changed_count']}"
    )
    print(
        "original->final P95/max = "
        f"{distance['original_core_to_caponly_final_core']['P95_um']:.12g}/"
        f"{distance['original_core_to_caponly_final_core']['max_um']:.12g} um"
    )
    print(
        "final->original P95/max = "
        f"{distance['caponly_final_core_to_original_core']['P95_um']:.12g}/"
        f"{distance['caponly_final_core_to_original_core']['max_um']:.12g} um"
    )
    print(
        "normal deviation P95/P99/max = "
        f"{normal['core_normal_deviation_P95_deg']:.12g}/"
        f"{normal['core_normal_deviation_P99_deg']:.12g}/"
        f"{normal['core_normal_deviation_max_deg']:.12g} deg"
    )
    print(f"winding consistent: {'YES' if normal['winding_consistent'] else 'NO'}")
    previous = result.previous_global_qc
    print("PREVIOUS GLOBAL REMESH:")
    print(
        "core distance P95/max = "
        f"{previous['original_core_to_previous_final']['P95_um']:.12g}/"
        f"{previous['original_core_to_previous_final']['max_um']:.12g} um"
    )
    print(
        "reverse core distance P95/max = "
        f"{previous['previous_final_core_to_original']['P95_um']:.12g}/"
        f"{previous['previous_final_core_to_original']['max_um']:.12g} um"
    )
    print(
        "normal deviation P95/P99/max = "
        f"{previous['core_normal_deviation']['P95_deg']:.12g}/"
        f"{previous['core_normal_deviation']['P99_deg']:.12g}/"
        f"{previous['core_normal_deviation']['max_deg']:.12g} deg"
    )
    print(
        "GLOBAL_REMESH_ARTIFACT_CONFIRMED: "
        f"{'YES' if result.artifact_confirmed else 'NO'}"
    )
    topology = result.final_topology
    print(
        "Final surface: "
        f"vertices={topology['vertex_count']} triangles={topology['triangle_count']} "
        f"components={topology['component_count']} watertight={topology['watertight']} "
        f"boundary_edges={topology['boundary_edge_count']} "
        f"nonmanifold={topology['nonmanifold_edge_count']} "
        f"self_intersections={topology['self_intersection_count']} "
        f"degenerate={topology['degenerate_triangle_count']} "
        f"winding_consistent={topology['winding_consistent']}"
    )
    print(
        "Radius P95: "
        f"{result.radius_qc.get('p95_absolute_relative_error', 'NOT_RUN')}"
    )
    print(f"Collision count: {result.collision_count}")
    print(
        "MANUAL REVIEW CAP-ONLY STL: "
        f"{result.output_paths.get('manual_review_stl', 'NOT_GENERATED')}"
    )
    print(f"TAGGED VTP: {result.output_paths.get('tagged_vtp', 'NOT_GENERATED')}")
    figures = {path.name: path for path in result.figures}
    for label, name in (
        ("BIFURCATION COMPARISON", "bifurcation_artifact_comparison.png"),
        ("GLOBAL REMESH VS CAPONLY", "global_remesh_vs_caponly_whole_surface.png"),
        ("EXTENSION CLOSEUPS", "extension_caponly_closeups.png"),
        ("WIREFRAME", "extension_caponly_wireframe.png"),
        ("NORMAL HOTSPOTS", "core_normal_deviation_hotspots.png"),
    ):
        print(f"{label}: {figures.get(name, 'NOT_GENERATED')}")
    print("visible_acceptance: MANUAL_REVIEW_REQUIRED")
    print(f"Final status: {result.status}")
    print(f"NEXT: {result.next_stage}")
