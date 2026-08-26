"""Orchestration for immutable-core, refined CFD extension surfaces."""

from __future__ import annotations

import re
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from .config import SurfacePrepareConfig
from .export import (
    OutputLayout,
    create_layout,
    export_geometry,
    meter_scale_qc,
    write_csv,
    write_json,
)
from .extension import extrude_and_cap
from .io import (
    SurfacePrepareError,
    load_original_surface,
    load_surface_inputs,
    read_json,
    sha256_file,
)
from .local_cut import local_plane_cut
from .mesh_quality import (
    ASPECT_RATIO_DEFINITION,
    compare_previous_extensions,
    extension_mesh_quality_qc,
    measure_local_original_mesh,
)
from .pressure_correction import (
    build_extended_boundary_conditions,
    calculate_pressure_corrections,
)
from .qc import (
    boundary_geometry_qc,
    core_surface_preservation_qc,
    extension_collision_qc,
    original_locked_vertex_motion_qc,
    surface_topology_qc,
)
from .types import TaggedSurface
from .visualization import (
    load_previous_tagged_surface,
    save_refinement_review_figures,
)


SUCCESS_STATUS = "CFD_EXTENSION_MESH_REFINED_PENDING_MANUAL_REVIEW"
NEXT_STAGE = "MANUAL_REFINED_CFD_SURFACE_REVIEW"


@dataclass(frozen=True, slots=True)
class SurfacePrepareResult:
    status: str
    run_root: Path
    input_run: Path
    original_stl: Path
    original_sha_before: str
    original_sha_after: str
    previous_reference_unchanged: bool
    boundaries: tuple[dict[str, Any], ...]
    surface_qc: dict[str, Any]
    core_qc: dict[str, Any]
    locked_vertex_qc: dict[str, Any]
    collision_qc: dict[str, Any]
    mesh_quality_qc: dict[str, Any]
    output_paths: dict[str, Any]
    figure_paths: tuple[Path, Path, Path]


@dataclass(frozen=True, slots=True)
class PreviousReference:
    root: Path
    stl: Path
    vtp: Path
    sha_before: str


def _anchor_id(roi_id: str) -> int:
    match = re.search(r"anchor_(\d+)", roi_id)
    if not match:
        raise SurfacePrepareError(f"Cannot identify anchor from ROI ID: {roi_id}")
    return int(match.group(1))


def _new_run_id(output_root: Path, anchor: int) -> str:
    base = f"cfd_surface_refined_anchor{anchor:06d}_{datetime.now():%Y%m%d_%H%M%S}"
    candidate = base
    suffix = 1
    while (output_root / candidate).exists():
        candidate = f"{base}_{suffix:02d}"
        suffix += 1
    return candidate


def _load_previous_reference(root: Path) -> PreviousReference:
    previous_root = Path(root).resolve()
    summary = read_json(previous_root / "qc" / "run_summary.json")
    if summary.get("status") != "CFD_SURFACE_PREPARE_PASS_PENDING_MANUAL_REVIEW":
        raise SurfacePrepareError("PREVIOUS_DIRECT_EXTRUSION_REFERENCE_NOT_PASS")
    stl = previous_root / "geometry" / "cfd_surface_extended_um.stl"
    vtp = previous_root / "geometry" / "cfd_surface_extended_um.vtp"
    if not stl.is_file() or not vtp.is_file():
        raise SurfacePrepareError("PREVIOUS_DIRECT_EXTRUSION_REFERENCE_MISSING")
    return PreviousReference(previous_root, stl, vtp, sha256_file(stl))


def _write_input_records(
    layout: OutputLayout,
    config: SurfacePrepareConfig,
    inputs: Any,
    previous: PreviousReference,
    sha_before: str,
) -> None:
    shutil.copy2(config.source_path, layout.input / config.source_path.name)
    artifacts = {
        "boundary_conditions.json": inputs.preprocess_run
        / "roi"
        / "boundary_conditions.json",
        "port_classification.csv": inputs.preprocess_run
        / "roi"
        / "port_classification.csv",
        "port_extension_plan.csv": inputs.preprocess_run
        / "roi"
        / "port_extension_plan.csv",
        "port_planes.vtp": inputs.preprocess_run / "roi" / "port_planes.vtp",
        "geometry_reference.json": inputs.preprocess_run
        / "input"
        / "geometry_reference.json",
    }
    write_json(
        layout.input / "source_manifest.json",
        {
            "configuration": config.source_path,
            "upstream_regeneration_performed": False,
            "whole_surface_remeshing_performed": False,
            "volume_mesh_created": False,
            "three_dimensional_cfd_run": False,
            "artifacts": {
                name: {"path": path.resolve(), "sha256": sha256_file(path)}
                for name, path in artifacts.items()
            },
        },
    )
    write_json(
        layout.input / "original_surface_reference.json",
        {
            "read_only": True,
            "lumen_surface_um_stl": inputs.original_surface_um_stl,
            "lumen_surface_um_vtp": inputs.original_surface_um_vtp,
            "lumen_surface_m_stl": inputs.original_surface_m_stl,
            "original_surface_sha256_before": sha_before,
            "source_geometry_reference": inputs.geometry_reference,
        },
    )
    write_json(
        layout.input / "cfd_preprocess_reference.json",
        {
            "run_root": inputs.preprocess_run,
            "status": inputs.preprocess_summary["status"],
            "run_summary": inputs.preprocess_summary,
        },
    )
    write_json(
        layout.input / "previous_direct_extrusion_reference.json",
        {
            "role": "PREVIOUS_DIRECT_EXTRUSION_REFERENCE",
            "read_only": True,
            "run_root": previous.root,
            "surface_stl": previous.stl,
            "surface_vtp": previous.vtp,
            "surface_stl_sha256_before": previous.sha_before,
        },
    )


def _write_report(
    path: Path,
    *,
    inputs: Any,
    integrity: dict[str, Any],
    boundary_rows: list[dict[str, Any]],
    surface_report: dict[str, Any],
    core_report: dict[str, Any],
    locked_report: dict[str, Any],
    collision_report: dict[str, Any],
    mesh_report: dict[str, Any],
    paths: dict[str, Any],
) -> None:
    quality_by_id = {item["port_id"]: item for item in mesh_report["boundaries"]}
    lines = [
        "# Refined CFD extension surface",
        "",
        "This run keeps the validated Ultraliser core and every proximal cut ring locked. Only "
        "the artificial extension surface is refined with multiple axial rings, a one-diameter "
        "shape transition, quality-driven quad diagonals, and constrained Taubin smoothing.",
        "",
        f"- Input preprocess run: `{inputs.preprocess_run}`",
        f"- Original STL unchanged: `{integrity['original_unchanged']}`",
        f"- Previous direct-extrusion reference unchanged: `{integrity['previous_reference_unchanged']}`",
        "- Whole-surface smoothing/remeshing/reconstruction: `False`",
        "- Volume mesh created: `False`",
        "- 3D CFD run: `False`",
        "",
        "## Mesh method and quality",
        "",
        f"Aspect ratio definition: `{ASPECT_RATIO_DEFINITION}`. The refined extension is intended "
        "to improve surface triangle quality before volume meshing.",
        "",
    ]
    for boundary in boundary_rows:
        quality = quality_by_id[boundary["port_id"]]
        lines.append(
            f"- `{boundary['port_id']}`: rings={boundary['ring_count']}, "
            f"target={boundary['target_edge_length_um']:.6g} um, triangles="
            f"{quality['extension_triangle_count']}, P95 aspect="
            f"{quality['aspect_ratio_p95']:.6g}, bad fraction="
            f"{quality['bad_triangle_fraction']:.6g}, P_solver="
            f"{boundary['P_solver_boundary_pa']}"
        )
    lines.extend(
        (
            "",
            "## Quality controls",
            "",
            f"- Extension mesh quality: `{mesh_report['status']}`",
            f"- Surface topology: `{surface_report['status']}`",
            f"- Core closest-point preservation: `{core_report['status']}`; max "
            f"{core_report['max_core_surface_distance_um']:.6g} um, P95 "
            f"{core_report['P95_core_surface_distance_um']:.6g} um",
            f"- Direct locked-original-vertex motion: `{locked_report['status']}`; max "
            f"{locked_report['original_locked_vertex_motion_max_um']:.6g} um",
            f"- Extension collision QC: `{collision_report['status']}`; count "
            f"{collision_report['extension_collision_count']}",
            "",
            "## Manual review",
            "",
            f"- Refined STL: `{paths['cfd_surface_refined_um_stl']}`",
            f"- Tagged refined VTP: `{paths['cfd_surface_refined_um_vtp']}`",
            f"- Status: `{SUCCESS_STATUS}`",
            f"- Next stage: `{NEXT_STAGE}`",
            "",
        )
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def _short_port_id(port_id: str) -> str:
    return port_id.rsplit("__", maxsplit=1)[-1]


def _print_centerline_qc(records: Any) -> None:
    print("INTERMEDIATE RING CENTERLINE QC")
    for boundary in records:
        print(
            f"{_short_port_id(str(boundary['port_id']))}: "
            f"axis dot={boundary['extension_axis_dot']:.12g} "
            f"ring count={boundary['ring_count']} "
            f"max drift um="
            f"{boundary['maximum_intermediate_ring_center_drift_um']:.12g} "
            f"P95 drift um="
            f"{boundary['P95_intermediate_ring_center_drift_um']:.12g} "
            f"mean drift um="
            f"{boundary['mean_intermediate_ring_center_drift_um']:.12g} "
            f"worst ring={boundary['worst_ring_index']} "
            f"worst station um={boundary['worst_ring_axial_station_um']:.12g}"
        )


def run_surface_prepare(
    config: SurfacePrepareConfig, *, project_root: Path
) -> SurfacePrepareResult:
    del project_root
    inputs = load_surface_inputs(
        config.paths.cfd_preprocess_run,
        expected_boundary_count=config.qc.expected_boundary_count,
    )
    previous = _load_previous_reference(config.manual_review.previous_surface_run)
    run_id = _new_run_id(
        config.paths.output_root,
        _anchor_id(str(inputs.preprocess_summary["roi_id"])),
    )
    layout = create_layout(config.paths.output_root, run_id)
    sha_before = sha256_file(inputs.original_surface_um_stl)
    _write_input_records(layout, config, inputs, previous, sha_before)

    original = load_original_surface(inputs.original_surface_um_vtp)
    previous_surface = load_previous_tagged_surface(previous.vtp, inputs.boundaries)
    local_statistics = [
        measure_local_original_mesh(
            original,
            boundary,
            sampling_radius_factor=config.extension_mesh.refinement.local_mesh_sampling_radius_factor,
        )
        for boundary in inputs.boundaries
    ]
    local_by_id = {item["port_id"]: item for item in local_statistics}
    write_csv(layout.qc / "local_original_mesh_statistics.csv", local_statistics)

    surface = TaggedSurface.from_mesh(original)
    cut_reports: list[dict[str, Any]] = []
    boundary_results = []
    for boundary in inputs.boundaries:
        surface, loop, cut_report = local_plane_cut(
            surface,
            boundary,
            radial_factor=config.local_cut.local_radial_radius_factor,
            axial_back_factor=config.local_cut.local_axial_back_radius_factor,
            axial_forward_factor=config.local_cut.local_axial_forward_radius_factor,
        )
        surface, boundary_result = extrude_and_cap(
            surface,
            loop,
            boundary,
            local_original_median_edge_length_um=local_by_id[boundary.port_id][
                "edge_length_median_um"
            ],
            mesh_config=config.extension_mesh,
        )
        cut_reports.append(cut_report)
        boundary_results.append(boundary_result)
    surface.compact()

    boundary_report = boundary_geometry_qc(
        inputs.boundaries, cut_reports, boundary_results, config.qc
    )
    mesh_report = extension_mesh_quality_qc(
        surface, inputs.boundaries, boundary_results, config.mesh_quality
    )
    for record in mesh_report["boundaries"]:
        record.update(
            {
                "edge_split_count": 0,
                "edge_split_safety_decision": (
                    "not applied; no topology-safe split was needed outside the separately "
                    "checked locked interface"
                ),
            }
        )
    comparison_rows = compare_previous_extensions(
        previous.vtp, surface, inputs.boundaries, boundary_results
    )
    quality_by_id = {item["port_id"]: item for item in mesh_report["boundaries"]}
    for row in comparison_rows:
        quality = quality_by_id[row["port_id"]]
        row.update(
            {
                "refined_interface_aspect_ratio_p95": quality[
                    "interface_triangle_aspect_ratio_p95"
                ],
                "refined_interface_minimum_angle_deg": quality[
                    "interface_minimum_angle_deg"
                ],
                "interface_p95_aspect_ratio_improved": quality[
                    "interface_triangle_aspect_ratio_p95"
                ]
                < row["previous_aspect_ratio_p95"],
                "interface_minimum_angle_improved": quality[
                    "interface_minimum_angle_deg"
                ]
                > row["previous_minimum_angle_deg"],
            }
        )
        row["mesh_quality_improved"] = bool(
            row["triangle_count_increased"]
            and row["p95_aspect_ratio_improved"]
            and row["minimum_angle_improved"]
            and row["interface_p95_aspect_ratio_improved"]
            and row["interface_minimum_angle_improved"]
            and quality["status"] == "PASS"
        )
    comparison_pass = all(
        row["mesh_quality_improved"] for row in comparison_rows
    )
    mesh_report["previous_vs_refined_checks"] = {
        "status": "PASS" if comparison_pass else "FAIL",
        "all_triangle_counts_increased": all(
            row["triangle_count_increased"] for row in comparison_rows
        ),
        "all_p95_aspect_ratios_improved": all(
            row["p95_aspect_ratio_improved"] for row in comparison_rows
        ),
        "all_minimum_angles_improved": all(
            row["minimum_angle_improved"] for row in comparison_rows
        ),
        "all_interface_p95_aspect_ratios_improved": all(
            row["interface_p95_aspect_ratio_improved"]
            for row in comparison_rows
        ),
        "all_interface_minimum_angles_improved": all(
            row["interface_minimum_angle_improved"] for row in comparison_rows
        ),
    }
    if not comparison_pass:
        mesh_report["status"] = "FAIL"
    surface_report, intersections = surface_topology_qc(surface, config.qc)
    core_report = core_surface_preservation_qc(
        original,
        surface.mesh(),
        inputs.boundaries,
        config.local_cut,
        config.qc,
    )
    locked_report = original_locked_vertex_motion_qc(
        original, surface, inputs.boundaries, config.local_cut
    )
    collision_report = extension_collision_qc(surface, intersections)
    sha_after = sha256_file(inputs.original_surface_um_stl)
    previous_sha_after = sha256_file(previous.stl)
    integrity = {
        "status": "PASS"
        if sha_before == sha_after and previous.sha_before == previous_sha_after
        else "FAIL",
        "original_surface_path": inputs.original_surface_um_stl,
        "original_surface_sha256_before": sha_before,
        "original_surface_sha256_after": sha_after,
        "original_unchanged": sha_before == sha_after,
        "previous_reference_path": previous.stl,
        "previous_reference_sha256_before": previous.sha_before,
        "previous_reference_sha256_after": previous_sha_after,
        "previous_reference_unchanged": previous.sha_before == previous_sha_after,
    }
    write_json(layout.qc / "original_surface_integrity.json", integrity)
    write_json(layout.qc / "local_cut_qc.json", {"boundaries": cut_reports})
    write_json(layout.qc / "boundary_geometry_qc.json", boundary_report)
    write_json(layout.qc / "extension_mesh_quality_qc.json", mesh_report)
    write_csv(layout.qc / "extension_mesh_before_after.csv", comparison_rows)
    write_json(layout.qc / "core_surface_preservation_qc.json", core_report)
    write_json(layout.qc / "original_locked_vertex_motion_qc.json", locked_report)
    write_json(layout.qc / "surface_qc.json", surface_report)
    write_json(layout.qc / "extension_collision_qc.json", collision_report)
    failures = [
        name
        for name, report in (
            ("surface_integrity", integrity),
            ("boundary_geometry", boundary_report),
            ("extension_mesh_quality", mesh_report),
            ("surface_topology", surface_report),
            ("core_surface_preservation", core_report),
            ("original_locked_vertex_motion", locked_report),
            ("extension_collision", collision_report),
        )
        if report["status"] != "PASS"
    ]
    if failures:
        write_json(
            layout.qc / "run_summary.json",
            {
                "status": "CFD_EXTENSION_MESH_REFINEMENT_FAILED",
                "run_id": run_id,
                "run_root": layout.root.resolve(),
                "failed_checks": failures,
                "formal_refined_surface_emitted": False,
                "volume_mesh_created": False,
                "three_dimensional_cfd_run": False,
            },
        )
        _print_centerline_qc(boundary_report["boundaries"])
        print("FAILED MANDATORY QC:")
        for name in failures:
            print(name)
        for boundary in boundary_report["boundaries"]:
            failed_boundary_checks = [
                name for name, passed in boundary["checks"].items() if not passed
            ]
            if failed_boundary_checks:
                print(
                    f"FAILED BOUNDARY: {boundary['port_id']} "
                    f"checks={','.join(failed_boundary_checks)} "
                    f"axis_dot={boundary['extension_axis_dot']:.12g} "
                    f"axis_threshold="
                    f"{boundary['minimum_extension_axis_dot_allowed']:.12g} "
                    f"max_center_drift_um="
                    f"{boundary['maximum_intermediate_ring_center_drift_um']:.12g} "
                    f"centerline_tolerance_um="
                    f"{boundary['intermediate_ring_centerline_tolerance_um']:.12g}"
                )
        raise SurfacePrepareError("; ".join(failures))

    pressure_rows = calculate_pressure_corrections(
        inputs.boundaries,
        boundary_results,
        dynamic_viscosity_pa_s=float(
            inputs.original_boundary_conditions["fluid"]["dynamic_viscosity_pa_s"]
        ),
        allow_negative_gauge_pressure=config.pressure_correction.allow_negative_gauge_pressure,
    )
    refined_bc = build_extended_boundary_conditions(
        inputs.original_boundary_conditions, inputs.boundaries, pressure_rows
    )
    shutil.copy2(
        inputs.preprocess_run / "roi" / "boundary_conditions.json",
        layout.bc / "boundary_conditions_original.json",
    )
    write_json(layout.bc / "boundary_conditions_refined.json", refined_bc)
    write_csv(
        layout.bc / "extension_pressure_correction_refined.csv", pressure_rows
    )

    output_paths = export_geometry(
        surface,
        inputs.boundaries,
        layout,
        create_meter_copy=config.geometry.create_meter_copy,
    )
    scale_report = meter_scale_qc(
        output_paths["cfd_surface_refined_um_stl"],
        output_paths["cfd_surface_refined_m_stl"],
    )
    write_json(layout.qc / "meter_scale_qc.json", scale_report)
    if scale_report["status"] != "PASS":
        raise SurfacePrepareError("METER_STL_SCALE_INVALID")
    figure_paths = save_refinement_review_figures(
        previous_surface, surface, inputs.boundaries, layout.figures
    )
    if len(figure_paths) != 3:
        raise SurfacePrepareError("REQUIRED_REFINEMENT_REVIEW_FIGURES_NOT_CREATED")

    comparison_by_id = {item["port_id"]: item for item in comparison_rows}
    boundary_rows = [
        {
            **result.report(),
            **quality_by_id[result.port_id],
            **comparison_by_id[result.port_id],
            **next(row for row in pressure_rows if row["port_id"] == result.port_id),
        }
        for result in boundary_results
    ]
    _write_report(
        layout.report / "cfd_surface_prepare_report.md",
        inputs=inputs,
        integrity=integrity,
        boundary_rows=boundary_rows,
        surface_report=surface_report,
        core_report=core_report,
        locked_report=locked_report,
        collision_report=collision_report,
        mesh_report=mesh_report,
        paths=output_paths,
    )
    summary = {
        "status": SUCCESS_STATUS,
        "next_stage": NEXT_STAGE,
        "run_id": run_id,
        "run_root": layout.root.resolve(),
        "configuration": config.source_path,
        "input_cfd_preprocess_run": inputs.preprocess_run,
        "previous_direct_extrusion_run": previous.root,
        "original_surface_path": inputs.original_surface_um_stl,
        "original_surface_sha256_before": sha_before,
        "original_surface_sha256_after": sha_after,
        "original_unchanged": True,
        "previous_reference_unchanged": True,
        "boundary_count": len(inputs.boundaries),
        "boundaries": boundary_rows,
        "extension_mesh_quality_qc": mesh_report,
        "surface_qc": surface_report,
        "core_surface_preservation_qc": core_report,
        "original_locked_vertex_motion_qc": locked_report,
        "extension_collision_qc": collision_report,
        "meter_scale_qc": scale_report,
        "outputs": output_paths,
        "figures": figure_paths,
        "upstream_regeneration_performed": False,
        "whole_surface_smoothing_performed": False,
        "whole_surface_remeshing_performed": False,
        "volume_mesh_created": False,
        "three_dimensional_cfd_run": False,
    }
    write_json(layout.qc / "run_summary.json", summary)
    return SurfacePrepareResult(
        status=SUCCESS_STATUS,
        run_root=layout.root.resolve(),
        input_run=inputs.preprocess_run,
        original_stl=inputs.original_surface_um_stl,
        original_sha_before=sha_before,
        original_sha_after=sha_after,
        previous_reference_unchanged=bool(integrity["previous_reference_unchanged"]),
        boundaries=tuple(boundary_rows),
        surface_qc=surface_report,
        core_qc=core_report,
        locked_vertex_qc=locked_report,
        collision_qc=collision_report,
        mesh_quality_qc=mesh_report,
        output_paths=output_paths,
        figure_paths=figure_paths,
    )


def print_result(result: SurfacePrepareResult) -> None:
    _print_centerline_qc(result.boundaries)
    print("MESH QUALITY COMPARISON")
    for boundary in result.boundaries:
        print(
            f"{_short_port_id(str(boundary['port_id']))}: "
            f"old P95 aspect={boundary['previous_aspect_ratio_p95']:.12g} "
            f"new P95 aspect={boundary['refined_aspect_ratio_p95']:.12g} "
            f"old min angle={boundary['previous_minimum_angle_deg']:.12g} "
            f"new min angle={boundary['refined_minimum_angle_deg']:.12g} "
            f"old triangles={boundary['previous_triangle_count']} "
            f"new triangles={boundary['refined_triangle_count']} "
            f"mesh quality improved="
            f"{'YES' if boundary['mesh_quality_improved'] else 'NO'}"
        )
    surface = result.surface_qc
    print("GEOMETRY QC")
    print(
        f"Original STL unchanged: "
        f"{'YES' if result.original_sha_before == result.original_sha_after else 'NO'}"
    )
    print(
        f"Old reference STL unchanged: "
        f"{'YES' if result.previous_reference_unchanged else 'NO'}"
    )
    print(
        "Locked original vertex max motion: "
        f"{result.locked_vertex_qc['original_locked_vertex_motion_max_um']:.12g}"
    )
    print(
        f"Core preservation max: "
        f"{result.core_qc['max_core_surface_distance_um']:.12g}"
    )
    print(
        f"Core preservation P95: "
        f"{result.core_qc['P95_core_surface_distance_um']:.12g}"
    )
    print(f"Components: {surface['component_count']}")
    print(f"Watertight: {'YES' if surface['watertight'] else 'NO'}")
    print(f"Boundary edges: {surface['boundary_edge_count']}")
    print(f"Nonmanifold: {surface['nonmanifold_edge_count']}")
    print(f"Self intersections: {surface['self_intersection_count']}")
    print(f"Degenerate: {surface['degenerate_triangle_count']}")
    print(
        f"Extension collisions: {result.collision_qc['extension_collision_count']}"
    )
    print("PRESSURE RESULTS")
    for boundary in result.boundaries:
        if boundary["role"] != "ASSUMED_OUTLET":
            continue
        print(
            f"{boundary['port_id']}: "
            f"P_original_1D={boundary['P_original_1D_pa']:.12g} "
            f"Q_expected_1D={boundary['Q_expected_1D_m3_s']:.12g} "
            f"final cap equivalent radius="
            f"{boundary['equivalent_radius_um']:.12g} "
            f"predicted pressure drop="
            f"{boundary['predicted_extension_pressure_drop_pa']:.12g} "
            f"P_solver_refined={boundary['P_solver_boundary_pa']:.12g}"
        )
    print("MANUAL REVIEW REFINED STL:")
    print(result.output_paths["cfd_surface_refined_um_stl"])
    print("TAGGED REFINED VTP:")
    print(result.output_paths["cfd_surface_refined_um_vtp"])
    print("BOUNDARY STL DIRECTORY:")
    print(result.output_paths["boundary_stl_directory"])
    print("SURFACE BEFORE/AFTER:")
    print(result.figure_paths[0])
    print("INTERFACE BEFORE/AFTER:")
    print(result.figure_paths[1])
    print("WIREFRAME BEFORE/AFTER:")
    print(result.figure_paths[2])
    print(f"Final status: {result.status}")
    print("NEXT:")
    print("MANUALLY REVIEW REFINED CFD STL SURFACE")
