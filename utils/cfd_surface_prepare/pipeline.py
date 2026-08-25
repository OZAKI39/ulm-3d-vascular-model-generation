"""Orchestration for immutable-input, local CFD surface preparation."""

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
    sha256_file,
)
from .local_cut import local_plane_cut
from .pressure_correction import (
    build_extended_boundary_conditions,
    calculate_pressure_corrections,
)
from .qc import (
    boundary_geometry_qc,
    core_surface_preservation_qc,
    extension_collision_qc,
    surface_topology_qc,
)
from .types import TaggedSurface
from .visualization import save_review_figures


SUCCESS_STATUS = "CFD_SURFACE_PREPARE_PASS_PENDING_MANUAL_REVIEW"
NEXT_STAGE = "MANUAL_CFD_SURFACE_REVIEW"


@dataclass(frozen=True, slots=True)
class SurfacePrepareResult:
    status: str
    run_root: Path
    input_run: Path
    original_stl: Path
    original_sha_before: str
    original_sha_after: str
    boundaries: tuple[dict[str, Any], ...]
    surface_qc: dict[str, Any]
    core_qc: dict[str, Any]
    collision_qc: dict[str, Any]
    output_paths: dict[str, Any]
    figure_paths: tuple[Path, Path]


def _anchor_id(roi_id: str) -> int:
    match = re.search(r"anchor_(\d+)", roi_id)
    if not match:
        raise SurfacePrepareError(f"Cannot identify anchor from ROI ID: {roi_id}")
    return int(match.group(1))


def _new_run_id(output_root: Path, anchor: int) -> str:
    base = f"cfd_surface_anchor{anchor:06d}_{datetime.now():%Y%m%d_%H%M%S}"
    candidate = base
    suffix = 1
    while (output_root / candidate).exists():
        candidate = f"{base}_{suffix:02d}"
        suffix += 1
    return candidate


def _write_input_records(
    layout: OutputLayout,
    config: SurfacePrepareConfig,
    inputs: Any,
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
            "configuration": str(config.source_path),
            "upstream_regeneration_performed": False,
            "volume_mesh_created": False,
            "three_dimensional_cfd_run": False,
            "artifacts": {
                name: {"path": str(path.resolve()), "sha256": sha256_file(path)}
                for name, path in artifacts.items()
            },
        },
    )
    write_json(
        layout.input / "original_surface_reference.json",
        {
            "read_only": True,
            "lumen_surface_um_stl": str(inputs.original_surface_um_stl),
            "lumen_surface_um_vtp": str(inputs.original_surface_um_vtp),
            "lumen_surface_m_stl": str(inputs.original_surface_m_stl),
            "original_surface_sha256_before": sha_before,
            "source_geometry_reference": inputs.geometry_reference,
        },
    )
    write_json(
        layout.input / "cfd_preprocess_reference.json",
        {
            "run_root": str(inputs.preprocess_run),
            "status": inputs.preprocess_summary["status"],
            "run_summary": inputs.preprocess_summary,
        },
    )


def _write_report(
    path: Path,
    *,
    inputs: Any,
    integrity: dict[str, Any],
    boundary_report: dict[str, Any],
    pressure_rows: list[dict[str, Any]],
    surface_report: dict[str, Any],
    core_report: dict[str, Any],
    collision_report: dict[str, Any],
    paths: dict[str, Any],
) -> None:
    pressure_by_id = {row["port_id"]: row for row in pressure_rows}
    lines = [
        "# CFD surface preparation report",
        "",
        "This run derives a CFD-only surface from the validated Ultraliser lumen. The source "
        "STL was treated as immutable; only the four saved boundary neighborhoods were cut, "
        "extruded, and capped.",
        "",
        f"- Input preprocess run: `{inputs.preprocess_run}`",
        f"- Original STL unchanged: `{integrity['original_unchanged']}`",
        "- Whole-surface reconstruction performed: `False`",
        "- Volume mesh created: `False`",
        "- 3D CFD run: `False`",
        "",
        "## Boundary construction and numerical pressure correction",
        "",
        "The outlet correction accounts only for the predicted loss in each artificial straight "
        "extension. It is not a physiological outlet model, resistance boundary condition, or "
        "Windkessel model. Negative gauge pressure is valid where produced.",
        "",
    ]
    for boundary in boundary_report["boundaries"]:
        pressure = pressure_by_id[boundary["port_id"]]
        lines.extend(
            (
                f"- `{boundary['port_id']}`: {boundary['role']}, "
                f"{boundary['boundary_origin']}, L={boundary['extension_length_um']:.6g} um, "
                f"r_eq={boundary['equivalent_radius_um']:.6g} um, "
                f"P_original={pressure['P_original_1D_pa']:.9g} Pa, "
                f"Q_expected={pressure['Q_expected_1D_m3_s']:.9g} m3/s, "
                f"P_solver={pressure['P_solver_boundary_pa']}",
            )
        )
    lines.extend(
        (
            "",
            "## Quality control",
            "",
            f"- Surface topology: `{surface_report['status']}`",
            f"- Core preservation: `{core_report['status']}`; max "
            f"{core_report['max_core_surface_distance_um']:.6g} um, P95 "
            f"{core_report['P95_core_surface_distance_um']:.6g} um",
            f"- Extension collision QC: `{collision_report['status']}`; count "
            f"{collision_report['extension_collision_count']}",
            "",
            "## Manual review",
            "",
            f"- STL: `{paths['cfd_surface_extended_um_stl']}`",
            f"- Tagged VTP: `{paths['cfd_surface_extended_um_vtp']}`",
            f"- Status: `{SUCCESS_STATUS}`",
            f"- Next stage: `{NEXT_STAGE}`",
            "",
        )
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def run_surface_prepare(
    config: SurfacePrepareConfig, *, project_root: Path
) -> SurfacePrepareResult:
    """Prepare one derived surface from the frozen saved PASS artifacts."""

    del project_root  # Paths are already made absolute by the strict loader.
    inputs = load_surface_inputs(
        config.paths.cfd_preprocess_run,
        expected_boundary_count=config.qc.expected_boundary_count,
    )
    roi_id = str(inputs.preprocess_summary["roi_id"])
    run_id = _new_run_id(config.paths.output_root, _anchor_id(roi_id))
    layout = create_layout(config.paths.output_root, run_id)
    sha_before = sha256_file(inputs.original_surface_um_stl)
    _write_input_records(layout, config, inputs, sha_before)
    original = load_original_surface(inputs.original_surface_um_vtp)
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
        surface, boundary_result = extrude_and_cap(surface, loop, boundary)
        cut_reports.append(cut_report)
        boundary_results.append(boundary_result)
    surface.compact()

    boundary_report = boundary_geometry_qc(
        inputs.boundaries, cut_reports, boundary_results, config.qc
    )
    surface_report, intersections = surface_topology_qc(surface, config.qc)
    core_report = core_surface_preservation_qc(
        original,
        surface.mesh(),
        inputs.boundaries,
        config.local_cut,
        config.qc,
    )
    collision_report = extension_collision_qc(surface, intersections)
    sha_after = sha256_file(inputs.original_surface_um_stl)
    integrity = {
        "status": "PASS" if sha_before == sha_after else "FAIL",
        "original_surface_path": str(inputs.original_surface_um_stl),
        "original_surface_sha256_before": sha_before,
        "original_surface_sha256_after": sha_after,
        "original_unchanged": sha_before == sha_after,
    }
    write_json(layout.qc / "original_surface_integrity.json", integrity)
    write_json(layout.qc / "local_cut_qc.json", {"boundaries": cut_reports})
    write_json(layout.qc / "boundary_geometry_qc.json", boundary_report)
    write_json(layout.qc / "core_surface_preservation_qc.json", core_report)
    write_json(layout.qc / "surface_qc.json", surface_report)
    write_json(layout.qc / "extension_collision_qc.json", collision_report)
    failures = [
        name
        for name, report in (
            ("original_surface_integrity", integrity),
            ("boundary_geometry", boundary_report),
            ("surface_topology", surface_report),
            ("core_surface_preservation", core_report),
            ("extension_collision", collision_report),
        )
        if report["status"] != "PASS"
    ]
    if failures:
        write_json(
            layout.qc / "run_summary.json",
            {
                "status": "CFD_SURFACE_PREPARE_FAILED",
                "run_id": run_id,
                "run_root": layout.root.resolve(),
                "failed_checks": failures,
                "original_surface_sha256_before": sha_before,
                "original_surface_sha256_after": sha_after,
                "original_unchanged": sha_before == sha_after,
                "formal_cfd_surface_emitted": False,
                "volume_mesh_created": False,
                "three_dimensional_cfd_run": False,
            },
        )
        raise SurfacePrepareError("; ".join(failures))

    dynamic_viscosity = float(
        inputs.original_boundary_conditions["fluid"]["dynamic_viscosity_pa_s"]
    )
    pressure_rows = calculate_pressure_corrections(
        inputs.boundaries,
        boundary_results,
        dynamic_viscosity_pa_s=dynamic_viscosity,
        allow_negative_gauge_pressure=config.pressure_correction.allow_negative_gauge_pressure,
    )
    extended_bc = build_extended_boundary_conditions(
        inputs.original_boundary_conditions, inputs.boundaries, pressure_rows
    )
    shutil.copy2(
        inputs.preprocess_run / "roi" / "boundary_conditions.json",
        layout.bc / "boundary_conditions_original.json",
    )
    write_json(layout.bc / "boundary_conditions_extended.json", extended_bc)
    write_csv(layout.bc / "extension_pressure_correction.csv", pressure_rows)

    output_paths = export_geometry(
        surface,
        inputs.boundaries,
        layout,
        create_meter_copy=config.geometry.create_meter_copy,
    )
    scale_report = meter_scale_qc(
        output_paths["cfd_surface_extended_um_stl"],
        output_paths["cfd_surface_extended_m_stl"],
    )
    write_json(layout.qc / "meter_scale_qc.json", scale_report)
    if scale_report["status"] != "PASS":
        raise SurfacePrepareError("METER_STL_SCALE_INVALID")
    figure_paths = (
        save_review_figures(
            original, surface, inputs.boundaries, layout.figures
        )
        if config.outputs.save_figures
        else tuple()
    )
    if len(figure_paths) != 2:
        raise SurfacePrepareError("REQUIRED_MANUAL_REVIEW_FIGURES_NOT_CREATED")
    _write_report(
        layout.report / "cfd_surface_prepare_report.md",
        inputs=inputs,
        integrity=integrity,
        boundary_report=boundary_report,
        pressure_rows=pressure_rows,
        surface_report=surface_report,
        core_report=core_report,
        collision_report=collision_report,
        paths=output_paths,
    )
    summary = {
        "status": SUCCESS_STATUS,
        "next_stage": NEXT_STAGE,
        "run_id": run_id,
        "run_root": layout.root.resolve(),
        "configuration": config.source_path,
        "input_cfd_preprocess_run": inputs.preprocess_run,
        "input_cfd_preprocess_status": inputs.preprocess_summary["status"],
        "original_surface_path": inputs.original_surface_um_stl,
        "original_surface_sha256_before": sha_before,
        "original_surface_sha256_after": sha_after,
        "original_unchanged": True,
        "boundary_count": len(inputs.boundaries),
        "boundaries": [
            {
                **result.report(),
                **next(
                    row for row in pressure_rows if row["port_id"] == result.port_id
                ),
            }
            for result in boundary_results
        ],
        "surface_qc": surface_report,
        "core_surface_preservation_qc": core_report,
        "extension_collision_qc": collision_report,
        "meter_scale_qc": scale_report,
        "outputs": output_paths,
        "figures": figure_paths,
        "upstream_regeneration_performed": False,
        "whole_surface_reconstruction_performed": False,
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
        boundaries=tuple(summary["boundaries"]),
        surface_qc=surface_report,
        core_qc=core_report,
        collision_qc=collision_report,
        output_paths=output_paths,
        figure_paths=figure_paths,
    )


def print_result(result: SurfacePrepareResult) -> None:
    print(f"Input CFD preprocess run: {result.input_run}")
    print(f"Original STL: {result.original_stl}")
    print(f"Original STL SHA before: {result.original_sha_before}")
    print(f"Original STL SHA after:  {result.original_sha_after}")
    print(
        "Original unchanged: "
        + ("YES" if result.original_sha_before == result.original_sha_after else "NO")
    )
    print(f"Boundary count: {len(result.boundaries)}")
    for boundary in result.boundaries:
        print(
            "BOUNDARY: "
            f"port_id={boundary['port_id']} origin={boundary['boundary_origin']} "
            f"role={boundary['role']} extension_length_um={boundary['extension_length_um']:.9g} "
            f"actual_cap_radius_um={boundary['equivalent_radius_um']:.9g} "
            f"P_original_pa={boundary['P_original_1D_pa']:.9g} "
            f"Q_expected_m3_s={boundary['Q_expected_1D_m3_s']:.9g} "
            f"P_solver_pa={boundary['P_solver_boundary_pa']}"
        )
    surface = result.surface_qc
    print(f"single_component = {surface['component_count'] == 1}")
    print(f"watertight = {surface['watertight']}")
    print(f"boundary_edges = {surface['boundary_edge_count']}")
    print(f"nonmanifold = {surface['nonmanifold_edge_count']}")
    print(f"self_intersections = {surface['self_intersection_count']}")
    print(f"degenerate_triangles = {surface['degenerate_triangle_count']}")
    print(
        "core_surface_max_change_um = "
        f"{result.core_qc['max_core_surface_distance_um']:.9g}"
    )
    print(
        "core_surface_p95_change_um = "
        f"{result.core_qc['P95_core_surface_distance_um']:.9g}"
    )
    print(
        "extension_collision_count = "
        f"{result.collision_qc['extension_collision_count']}"
    )
    print("MANUAL REVIEW STL:")
    print(result.output_paths["cfd_surface_extended_um_stl"])
    print("TAGGED VTP:")
    print(result.output_paths["cfd_surface_extended_um_vtp"])
    print("BOUNDARY STL FILES:")
    for path in result.output_paths["boundary_stl_paths"]:
        print(path)
    print(f"Final status: {result.status}")
    print("NEXT:")
    print("MANUALLY REVIEW CFD STL SURFACE")
