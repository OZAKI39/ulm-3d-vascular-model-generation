"""Formal VMTK TPS boundary-normal cross-seam surface pipeline."""

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
    active_collar_original_side_distance_qc,
    assign_cross_seam_active_entities,
    cross_seam_entity_preservation_qc,
    cross_seam_intersection_qc,
    crossseam_active_mesh_quality,
    cut_seam_topology_qc,
    locked_entity_preservation_qc,
    seam_local_normal_diagnostic,
)
from .io import (
    SurfacePrepareError,
    load_original_surface,
    load_surface_inputs,
    sha256_file,
)
from .local_cut import local_plane_cut
from .mesh_quality import measure_local_original_mesh
from .types import TaggedSurface
from .vmtk_qc import (
    active_collar_cross_section_fidelity_qc,
    extension_geometry_qc,
    geometry_pressure_correction,
    interface_smoothness_from_raw,
    meter_scale_qc,
    open_profile_qc,
    polydata_mesh,
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
from .vmtk_visualization import save_production_review_figures


ASSIGNMENT_FAILURE = "VMTK_CROSS_SEAM_ENTITY_ASSIGNMENT_FAILED"
CORE_FAILURE = "VMTK_CROSS_SEAM_FAR_CORE_MODIFIED"
COLLAR_FAILURE = "VMTK_CROSS_SEAM_COLLAR_GEOMETRY_FAILED"
RING_FAILURE = "VMTK_CROSS_SEAM_RING_TOPOLOGY_PRESERVED"
TOPOLOGY_FAILURE = "VMTK_CROSS_SEAM_TOPOLOGY_FAILED"
GEOMETRY_FAILURE = "VMTK_CROSS_SEAM_GEOMETRY_FAILED"
METER_FAILURE = "VMTK_METER_SCALE_SERIALIZATION_FAILED"
SUCCESS_STATUS = "CFD_SURFACE_CROSS_SEAM_FINAL_PASS_PENDING_MANUAL_REVIEW"
SUCCESS_NEXT = "MANUALLY REVIEW FINAL CAPPED CFD SURFACE"
FAIL_NEXT = "REVIEW CFD SURFACE PREPARATION FAILURE"


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


@dataclass(frozen=True, slots=True)
class VmtkSurfacePrepareResult:
    status: str
    next_stage: str
    run_root: Path
    summary: dict[str, Any]


def _layout(output_root: Path, run_id: str) -> VmtkLayout:
    root = Path(output_root) / run_id
    if root.exists():
        raise FileExistsError(f"Refusing to overwrite VMTK run: {root}")
    directories = [
        root / name
        for name in (
            "input", "vmtk", "geometry", "boundaries", "bc", "qc", "figures"
        )
    ]
    for directory in directories:
        directory.mkdir(parents=True, exist_ok=False)
    return VmtkLayout(root.resolve(), *[path.resolve() for path in directories])


def _anchor_id(roi_id: str) -> str:
    marker = "__anchor_"
    if marker not in roi_id:
        return "unknown"
    return roi_id.rsplit(marker, 1)[1].split("__", 1)[0]


def _open_surface(surface: TaggedSurface, path: Path) -> None:
    surface.compact()
    faces = np.column_stack(
        (np.full(len(surface.faces), 3, dtype=np.int64), surface.faces)
    ).ravel()
    data = pv.PolyData(np.asarray(surface.vertices, dtype=float), faces)
    data.point_data["source_vertex_index"] = surface.source_vertex_index
    data.cell_data["surface_role"] = np.zeros(len(surface.faces), dtype=np.uint8)
    data.save(path, binary=True)


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
                "open_profile_area_um2": profile_by_id[record["port_id"]][
                    "area_um2"
                ],
                "open_profile_point_count": profile_by_id[record["port_id"]][
                    "point_count"
                ],
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
    inlet.update(
        {
            "Q_solver_m3_s": rows[inlet["port_id"]]["Q_solver_m3_s"],
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
                "extension_resistance_pa_s_m3": row[
                    "extension_resistance_pa_s_m3"
                ],
                "predicted_extension_pressure_drop_pa": row[
                    "predicted_extension_pressure_drop_pa"
                ],
                "P_solver_boundary_pa": row["P_solver_boundary_pa"],
                "extension_pressure_correction_role": row[
                    "pressure_correction_role"
                ],
            }
        )
    output["surface_geometry_method"] = (
        "OFFICIAL_VMTK_TPS_BOUNDARY_NORMAL_CROSS_SEAM_REMESH"
    )
    output["pressure_correction_geometry"] = (
        "20 cross sections on final entity-aware VMTK geometry"
    )
    return output


def boundary_plane_alignment_pass(profile: dict[str, Any]) -> bool:
    """Keep the verified 0.999 boundary-normal alignment gate."""

    return all(
        float(row["boundary_plane_normal_abs_dot_expected_outward"]) >= 0.999
        for row in profile["boundaries"]
    )


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


def _use_true_intersection_count(
    topology: dict[str, Any], intersection: dict[str, Any]
) -> dict[str, Any]:
    result = copy.deepcopy(topology)
    result["vtk_candidate_self_intersection_count"] = result[
        "self_intersection_count"
    ]
    result["self_intersection_count"] = intersection[
        "true_self_intersection_count"
    ]
    result["checks"]["zero_self_intersections"] = (
        intersection["true_self_intersection_count"] == 0
    )
    result["status"] = "PASS" if all(result["checks"].values()) else "FAIL"
    return result


def _write_failure(
    layout: VmtkLayout,
    *,
    status: str,
    run_id: str,
    roi_id: str,
    runtime: dict[str, Any],
    outputs: dict[str, Any],
    evidence: dict[str, Any],
) -> None:
    write_json(
        layout.qc / "run_summary.json",
        {
            "status": status,
            "next": FAIL_NEXT,
            "run_id": run_id,
            "run_root": str(layout.root),
            "roi_id": roi_id,
            "runtime": runtime,
            "outputs": outputs,
            "evidence": evidence,
        },
    )


def run_vmtk_surface_prepare(
    config: SurfacePrepareConfig,
    *,
    project_root: Path,
    run_id: str | None = None,
) -> VmtkSurfacePrepareResult:
    """Run the single formal VMTK extension/remesh/cap production path."""

    settings = config.vmtk.entity_remesh
    inputs = load_surface_inputs(
        config.paths.cfd_preprocess_run,
        expected_boundary_count=config.qc.expected_boundary_count,
    )
    roi_id = str(inputs.preprocess_summary["roi_id"])
    anchor = _anchor_id(roi_id)
    candidate_id = run_id or (
        f"vmtk_tps_boundarynormal_crossseam_anchor{anchor}_"
        f"{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    layout = _layout(config.paths.output_root, candidate_id)
    paths = exchange_paths(
        input_directory=layout.input,
        vmtk_directory=layout.vmtk,
        geometry_directory=layout.geometry,
    )
    tool_script = project_root / "tools" / "run_vmtk_flowextension.py"
    outputs: dict[str, Any] = {}
    runtime: dict[str, Any] = {
        "source_provenance": official_source_provenance(config.vmtk)
    }

    def stop(status: str, evidence: dict[str, Any]) -> None:
        _write_failure(
            layout,
            status=status,
            run_id=candidate_id,
            roi_id=roi_id,
            runtime=runtime,
            outputs=outputs,
            evidence=evidence,
        )
        raise SurfacePrepareError(status)

    immutable_paths = (
        inputs.original_surface_um_stl,
        inputs.original_surface_um_vtp,
    )
    hashes_before = {str(path): sha256_file(path) for path in immutable_paths}
    write_json(
        layout.input / "original_surface_reference.json",
        {
            "geometry_reference": inputs.geometry_reference,
            "immutable_hashes_before": hashes_before,
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
    proximal_profile, proximal_loops = open_profile_qc(open_mesh, inputs.boundaries)
    write_json(layout.qc / "proximal_open_profile_qc.json", proximal_profile)
    if proximal_profile["status"] != "PASS" or not boundary_plane_alignment_pass(
        proximal_profile
    ):
        stop(GEOMETRY_FAILURE, {"proximal_open_profile": proximal_profile})
    _write_manifest(
        layout.input / "boundary_manifest.csv",
        _manifest_records(inputs.boundaries),
        proximal_profile,
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

    extension = run_official_vmtk(
        config=config.vmtk, paths=paths, tool_script=tool_script
    )
    runtime["extension"] = extension.runtime
    runtime["extension_command"] = list(extension.command)
    outputs.update(
        raw_vtp=str(paths.raw_vtp.resolve()),
        raw_stl=str(paths.raw_stl.resolve()),
    )

    raw_core = raw_core_exact_copy_qc(paths.open_surface_vtp, paths.raw_vtp)
    write_json(layout.qc / "raw_core_exact_copy_qc.json", raw_core)
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
    raw_collision = _collision_report(
        raw_mesh,
        raw_pairs,
        np.asarray(raw_data.cell_data["SurfaceRegionId"]),
    )
    raw_interface = interface_smoothness_from_raw(
        raw_mesh,
        input_point_count=open_data.n_points,
        boundaries=inputs.boundaries,
    )
    raw_geometry, _ = extension_geometry_qc(
        raw_mesh, inputs.boundaries, proximal_loops, raw_interface
    )
    write_json(layout.qc / "raw_surface_qc.json", raw_topology)
    write_json(layout.qc / "raw_extension_collision_qc.json", raw_collision)
    write_json(layout.qc / "raw_extension_geometry_qc.json", raw_geometry)
    if any(
        report["status"] != "PASS"
        for report in (raw_topology, raw_collision, raw_geometry)
    ):
        stop(
            GEOMETRY_FAILURE,
            {
                "raw_topology": raw_topology,
                "raw_collision": raw_collision,
                "raw_geometry": raw_geometry,
            },
        )

    remesh = entity_remesh_official_vmtk(
        config=config.vmtk, paths=paths, tool_script=tool_script
    )
    runtime["cross_seam_remesh"] = remesh.runtime
    runtime["cross_seam_remesh_command"] = list(remesh.command)
    outputs.update(
        remeshed_open_vtp=str(paths.remeshed_open_vtp.resolve()),
        remeshed_open_stl=str(paths.remeshed_open_stl.resolve()),
    )

    far_core_remesh, boundary_lock, extension_effect = (
        cross_seam_entity_preservation_qc(
            paths.raw_vtp,
            paths.remeshed_open_vtp,
            inputs.boundaries,
            entity_array_name=settings.entity_array_name,
            far_core_entity_id=settings.far_core_entity_id,
            active_entity_id=settings.active_entity_id,
        )
    )
    write_json(layout.qc / "far_core_post_remesh_qc.json", far_core_remesh)
    write_json(layout.qc / "far_core_active_boundary_lock_qc.json", boundary_lock)
    write_json(layout.qc / "extension_remesh_effect_qc.json", extension_effect)
    if far_core_remesh["status"] != "PASS" or boundary_lock["status"] != "PASS":
        stop(
            CORE_FAILURE,
            {"far_core": far_core_remesh, "boundary_lock": boundary_lock},
        )
    if extension_effect["status"] != "PASS":
        stop(GEOMETRY_FAILURE, {"extension_remesh_effect": extension_effect})

    collar_distance, _ = active_collar_original_side_distance_qc(
        paths.raw_vtp,
        paths.remeshed_open_vtp,
        inputs.boundaries,
        maximum_p95_distance_um=config.qc.maximum_core_surface_p95_distance_um,
        cross_seam_neighborhood_um=settings.target_edge_length_um,
    )
    cross_sections = active_collar_cross_section_fidelity_qc(
        paths.raw_vtp,
        paths.remeshed_open_vtp,
        inputs.boundaries,
        maximum_equivalent_radius_relative_error=(
            config.qc.maximum_equivalent_radius_relative_error
        ),
    )
    write_json(
        layout.qc / "active_collar_original_side_distance_qc.json",
        collar_distance,
    )
    write_json(
        layout.qc / "active_collar_cross_section_fidelity.json", cross_sections
    )
    if collar_distance["status"] != "PASS" or cross_sections["status"] != "PASS":
        stop(
            COLLAR_FAILURE,
            {"collar_distance": collar_distance, "cross_sections": cross_sections},
        )

    seam_topology = cut_seam_topology_qc(
        paths.raw_vtp, paths.remeshed_open_vtp, inputs.boundaries
    )
    write_json(layout.qc / "cut_seam_topology_qc.json", seam_topology)
    if seam_topology["status"] != "PASS":
        stop(RING_FAILURE, {"cut_seam_topology": seam_topology})

    _, remesh_mesh = polydata_mesh(paths.remeshed_open_vtp)
    open_topology_raw, _ = topology_qc(
        remesh_mesh,
        expected_open_profile_count=config.qc.expected_boundary_count,
        allow_degenerate=False,
    )
    intersections, intersection_records = cross_seam_intersection_qc(
        paths.remeshed_open_vtp
    )
    open_topology = _use_true_intersection_count(open_topology_raw, intersections)
    distal_profile, _ = open_profile_qc(
        remesh_mesh, inputs.boundaries, distal=True
    )
    open_topology["distal_profile_qc"] = distal_profile
    if distal_profile["status"] != "PASS":
        open_topology["status"] = "FAIL"
    write_json(layout.qc / "crossseam_open_surface_qc.json", open_topology)
    write_json(layout.qc / "crossseam_intersection_qc.json", intersections)

    seam_diagnostic = seam_local_normal_diagnostic(
        paths.remeshed_open_vtp,
        inputs.boundaries,
        target_edge_length_um=settings.target_edge_length_um,
    )
    write_json(layout.qc / "seam_normal_diagnostic.json", seam_diagnostic)
    extension_geometry, distal_loops = extension_geometry_qc(
        remesh_mesh, inputs.boundaries, proximal_loops, seam_diagnostic
    )
    write_json(layout.qc / "extension_geometry_qc.json", extension_geometry)
    active_quality, tail_records = crossseam_active_mesh_quality(
        paths.remeshed_open_vtp,
        inputs.boundaries,
        active_entity_id=settings.active_entity_id,
        local_target_edge_um=local_targets,
        quality=config.mesh_quality,
    )
    write_json(layout.qc / "crossseam_mesh_quality_qc.json", active_quality)
    if any(
        report["status"] != "PASS"
        for report in (
            open_topology,
            intersections,
            seam_topology,
            extension_geometry,
        )
    ):
        stop(
            TOPOLOGY_FAILURE,
            {
                "open_topology": open_topology,
                "intersections": intersections,
                "intersection_records": intersection_records,
                "seam_topology": seam_topology,
                "extension_geometry": extension_geometry,
            },
        )

    cap = cap_official_vmtk(
        config=config.vmtk, paths=paths, tool_script=tool_script
    )
    runtime["cap"] = cap.runtime
    runtime["cap_command"] = list(cap.command)
    write_json(layout.vmtk / "environment.json", runtime)

    final_outputs, boundary_mapping = tag_and_export_final_surface(
        paths.capped_vtp,
        inputs.boundaries,
        layout.geometry,
        layout.boundaries,
        output_stem="cfd_surface_vmtk_tps_boundarynormal_crossseam",
        raw_vtp=paths.remeshed_open_vtp,
        remesh_entity_codes={"CAP": 0, "FAR_CORE": 1, "CROSS_SEAM_ACTIVE": 2},
    )
    outputs.update(final_outputs)
    final_vtp = Path(outputs["tagged_vtp"])
    _, final_mesh = polydata_mesh(final_vtp)
    final_topology_raw, _ = topology_qc(
        final_mesh, expected_open_profile_count=0, require_winding_consistent=True
    )
    final_intersection, _ = cross_seam_intersection_qc(final_vtp)
    final_topology = _use_true_intersection_count(
        final_topology_raw, final_intersection
    )
    far_core_final = locked_entity_preservation_qc(
        paths.raw_vtp,
        final_vtp,
        entity_array_name=settings.entity_array_name,
        entity_id=settings.far_core_entity_id,
        label="far_core",
    )
    far_core = {
        "status": (
            "PASS"
            if far_core_remesh["status"] == far_core_final["status"] == "PASS"
            else "FAIL"
        ),
        "post_remesh": far_core_remesh,
        "post_cap": far_core_final,
    }
    final_surface = {
        "status": (
            "PASS"
            if final_topology["status"] == "PASS"
            and boundary_mapping["status"] == "PASS"
            and far_core["status"] == "PASS"
            else "FAIL"
        ),
        "topology": final_topology,
        "intersection": final_intersection,
        "boundary_mapping": boundary_mapping,
        "far_core_exact_preservation": far_core,
    }
    write_json(layout.qc / "far_core_exact_preservation_qc.json", far_core)
    write_json(layout.qc / "final_surface_qc.json", final_surface)
    if final_surface["status"] != "PASS":
        stop(TOPOLOGY_FAILURE, {"final_surface": final_surface})

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
                inputs.original_boundary_conditions["fluid"][
                    "dynamic_viscosity_pa_s"
                ]
            ),
        )
    except SurfacePrepareError as error:
        pressure_report = {"status": "FAIL", "error": str(error)}
    else:
        pressure_report["status"] = "PASS"
        pressure_csv = layout.bc / "extension_pressure_correction.csv"
        write_csv(
            pressure_csv,
            [
                {
                    **row,
                    "station_fractions": json.dumps(row["station_fractions"]),
                    "cross_section_area_um2": json.dumps(
                        row["cross_section_area_um2"]
                    ),
                    "equivalent_radius_um": json.dumps(
                        row["equivalent_radius_um"]
                    ),
                }
                for row in pressure_rows
            ],
        )
        boundary_json = layout.bc / "boundary_conditions.json"
        write_json(
            boundary_json,
            _boundary_conditions(inputs.original_boundary_conditions, pressure_rows),
        )
        outputs["pressure_correction_csv"] = str(pressure_csv.resolve())
        outputs["boundary_conditions_json"] = str(boundary_json.resolve())
    write_json(layout.qc / "extension_pressure_geometry_qc.json", pressure_report)

    hashes_after = {str(path): sha256_file(path) for path in immutable_paths}
    integrity_pass = hashes_before == hashes_after
    if not all(
        report["status"] == "PASS"
        for report in (radius_report, meter_report, pressure_report)
    ) or not integrity_pass:
        failure = (
            METER_FAILURE if meter_report["status"] != "PASS" else GEOMETRY_FAILURE
        )
        stop(
            failure,
            {
                "radius": radius_report,
                "meter": meter_report,
                "pressure": pressure_report,
                "source_integrity": integrity_pass,
            },
        )

    figure_paths = save_production_review_figures(
        raw_vtp=paths.raw_vtp,
        remeshed_open_vtp=paths.remeshed_open_vtp,
        final_vtp=final_vtp,
        boundaries=inputs.boundaries,
        output_directory=layout.figures,
    )
    outputs["figure_paths"] = [str(path) for path in figure_paths]
    summary = {
        "status": SUCCESS_STATUS,
        "next": SUCCESS_NEXT,
        "run_id": candidate_id,
        "run_root": str(layout.root),
        "roi_id": roi_id,
        "method": "VMTK TPS BOUNDARY_NORMAL + CROSS_SEAM REMESH",
        "runtime": runtime,
        "input_geometry_reference": inputs.geometry_reference,
        "cross_seam_entity_assignment": assignment,
        "core_collar_layers": collar_layers.tolist(),
        "far_core_exact_preservation": far_core,
        "active_collar_original_side_distance": collar_distance,
        "active_collar_cross_section_fidelity": cross_sections,
        "cut_seam_topology": seam_topology,
        "open_topology": open_topology,
        "open_intersections": intersections,
        "seam_normal_diagnostic": seam_diagnostic,
        "extension_geometry": extension_geometry,
        "mesh_quality": active_quality,
        "mesh_tail_count": len(tail_records),
        "final_surface": final_surface,
        "radius_fidelity": radius_report,
        "pressure_geometry": pressure_report,
        "meter_scale": meter_report,
        "source_geometry_unchanged": integrity_pass,
        "outputs": outputs,
    }
    write_json(layout.qc / "run_summary.json", summary)
    return VmtkSurfacePrepareResult(
        status=SUCCESS_STATUS,
        next_stage=SUCCESS_NEXT,
        run_root=layout.root,
        summary=summary,
    )


def print_vmtk_result(result: VmtkSurfacePrepareResult) -> None:
    summary = result.summary
    final = summary["final_surface"]
    mapping = final["boundary_mapping"]
    outputs = summary["outputs"]
    inlet_count = sum(
        row["role"] == "ASSUMED_INLET" for row in mapping["boundaries"]
    )
    outlet_count = sum(
        row["role"] == "ASSUMED_OUTLET" for row in mapping["boundaries"]
    )
    print("CFD SURFACE PREPARE")
    print(f"ROI: {summary['roi_id']}")
    print(f"Method: {summary['method']}")
    print(f"OPEN topology: {summary['open_topology']['status']}")
    print(f"FINAL topology: {final['topology']['status']}")
    print(
        "Radius P95: "
        f"{summary['radius_fidelity']['p95_absolute_relative_error']:.9g}"
    )
    print(
        "Boundary mapping: "
        f"{inlet_count} inlet / {outlet_count} outlets"
    )
    print(f"Final STL: {outputs['manual_review_stl']}")
    print(f"Final VTP: {outputs['tagged_vtp']}")
    print(f"STATUS: {result.status}")
    print(f"NEXT: {result.next_stage}")
