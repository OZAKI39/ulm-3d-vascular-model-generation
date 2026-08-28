"""Orchestration for formal global-to-ROI CFD boundary preprocessing."""

from __future__ import annotations

import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any

from .config import CFDPreprocessConfig
from .export import (
    OutputLayout,
    create_layout,
    write_boundary_package,
    write_extension_plan,
    write_global_tables,
    write_global_vtp,
    write_json,
    write_port_classification,
    write_port_vtp,
)
from .io import (
    load_matching_global_model,
    resolve_model_run,
    resolve_rodent_run,
    resolve_sampling_run,
    select_roi,
    verify_global_edge_manifest,
)
from .one_d_flow import GlobalFlowResult, solve_global_flow
from .port_transfer import PortTransfer, transfer_all_boundaries
from .qc import global_graph_qc, roi_readiness
from .visualization import global_pressure_figure, roi_boundary_figure


@dataclass(frozen=True, slots=True)
class CFDPreprocessResult:
    layout: OutputLayout
    status: str
    summary: dict[str, Any]
    transfers: list[PortTransfer]
    boundary_conditions_path: Path | None
    port_planes_path: Path | None


def _unique_run_id(output_root: Path, anchor: int) -> str:
    base = (
        f"global_to_roi_anchor{anchor:06d}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    candidate = base
    revision = 1
    while (output_root / candidate).exists():
        revision += 1
        candidate = f"{base}_r{revision}"
    return candidate


def next_stage_for_status(status: str) -> str:
    """Return the only permitted workflow transition for a final status."""

    return {
        "CFD_PREPROCESS_BASELINE_PASS": "CFD_SURFACE_VOLUME_MESH_PREPARATION",
        "CFD_ROI_NOT_READY": "REVIEW_CFD_PREPROCESS_FAILURE",
        "GLOBAL_1D_FLOW_FAILED": "REVIEW_GLOBAL_1D_FLOW",
        "GLOBAL_TO_ROI_TRANSFER_FAILED": "REVIEW_GLOBAL_TO_ROI_BOUNDARY_TRANSFER",
        "CFD_GEOMETRY_REFERENCE_INVALID": "REVIEW_CFD_PREPROCESS_FAILURE",
        "CFD_PREPROCESS_INTERNAL_ERROR": "REVIEW_CFD_PREPROCESS_FAILURE",
    }.get(status, "REVIEW_CFD_PREPROCESS_FAILURE")


def next_stage_display(next_stage: str) -> str:
    return {
        "CFD_SURFACE_VOLUME_MESH_PREPARATION": (
            "CFD SURFACE / VOLUME MESH PREPARATION"
        ),
        "REVIEW_CFD_PREPROCESS_FAILURE": "REVIEW CFD PREPROCESS FAILURE REASONS",
        "REVIEW_GLOBAL_1D_FLOW": "REVIEW GLOBAL 1D FLOW",
        "REVIEW_GLOBAL_TO_ROI_BOUNDARY_TRANSFER": (
            "REVIEW GLOBAL-TO-ROI BOUNDARY TRANSFER"
        ),
    }[next_stage]


def _solver_summary(
    flow: GlobalFlowResult, *, leaf_pressure_pa: float
) -> dict[str, Any]:
    return {
        "status": "PASS",
        "backend": "scipy_sparse_spsolve",
        "simulation_direction_basis": "SWC_PARENT_TO_CURRENT",
        "simulation_direction_is_measured": False,
        "simulation_direction_is_physiological_ground_truth": False,
        "structural_root_node_id": flow.root_node_id,
        "structural_leaf_count": len(flow.leaf_node_ids),
        "root_mean_velocity_mm_s": flow.root_mean_velocity_mm_s,
        "global_root_mean_velocity_is_measured": False,
        "global_root_mean_velocity_role": "LITERATURE_DERIVED_BASELINE_ASSUMPTION",
        "root_flow_rate_m3_s": flow.root_flow_m3_s,
        "leaf_flow_sum_m3_s": flow.leaf_flow_sum_m3_s,
        "relative_mass_error": flow.relative_mass_error,
        "maximum_internal_relative_residual": flow.maximum_internal_relative_residual,
        "reverse_flow_count": flow.reverse_flow_count,
        "leaf_pressure_pa": leaf_pressure_pa,
        "leaf_pressure_role": "GAUGE_PRESSURE_REFERENCE_NOT_PHYSIOLOGICAL_ZERO",
        "radius_source": "ORIGINAL_ANALYSIS_SWC_SOURCE_RADIUS_NOT_RADIUS_SCALE_0.91",
    }


def _boundary_qc(
    transfers: list[PortTransfer],
    boundary_mass_error: float,
    *,
    mass_tolerance: float,
    maximum_position_error_um: float,
    maximum_radius_relative_error: float,
) -> dict[str, Any]:
    position_pass = all(
        item.position_error_um <= maximum_position_error_um for item in transfers
    )
    radius_pass = all(
        item.radius_relative_error <= maximum_radius_relative_error
        for item in transfers
    )
    mass_pass = boundary_mass_error <= mass_tolerance
    return {
        "status": "PASS" if position_pass and radius_pass and mass_pass else "FAIL",
        "all_global_edges_found": True,
        "all_boundary_positions_match": position_pass,
        "all_boundary_radii_match": radius_pass,
        "pressure_interpolation": "EXACT_PARTIAL_LINEAR_RADIUS_HYDRAULIC_RESISTANCE",
        "relative_boundary_mass_error": boundary_mass_error,
        "relative_boundary_mass_tolerance": mass_tolerance,
        "maximum_allowed_position_error_um": maximum_position_error_um,
        "maximum_allowed_radius_relative_error": maximum_radius_relative_error,
        "boundaries": [
            {
                "port_id": item.port_id,
                "role": item.role,
                "boundary_origin": item.boundary_origin,
                "global_node_id": item.global_node_id,
                "global_edge_id": item.global_edge_id,
                "alpha_on_global_edge": item.alpha_on_global_edge,
                "position_error_um": item.position_error_um,
                "radius_relative_error": item.radius_relative_error,
                "pressure_pa": item.pressure_pa,
                "flow_rate_m3_s": item.role_flow_m3_s,
                "signed_parent_to_child_flow_m3_s": item.signed_parent_to_child_flow_m3_s,
            }
            for item in transfers
        ],
    }


def run_cfd_preprocess(
    config: CFDPreprocessConfig,
    *,
    project_root: Path,
) -> CFDPreprocessResult:
    """Run the complete saved-input preprocessing workflow once."""

    project_root = Path(project_root).resolve()
    sampling_run = resolve_sampling_run(
        config.paths.sampling_run, project_root=project_root
    )
    roi = select_roi(
        sampling_run,
        anchor=config.selection.roi_anchor,
        roi_id=config.selection.roi_id,
    )
    run_id = _unique_run_id(config.paths.output_root, roi.anchor_id)
    layout = create_layout(config.paths.output_root, run_id)
    shutil.copy2(config.source_path, layout.config / "source_cfd_preprocess.yaml")

    rodent_root = project_root / "outputs" / "rodent_vasculature"
    rodent_run, rodent_auto = resolve_rodent_run(
        config.paths.rodent_run,
        output_root=rodent_root,
        source_model_id=roi.source_model_id,
    )
    model = load_matching_global_model(rodent_run, roi.source_model_id)
    verify_global_edge_manifest(model, sampling_run / "manifests" / "global_edges.csv")

    mu = config.fluid.dynamic_viscosity_pa_s
    flow = solve_global_flow(
        model,
        mu_pa_s=mu,
        leaf_pressure_pa=config.global_1d.leaf_pressure_pa,
        boundary_type=config.global_1d.inlet.boundary_type,
        mean_velocity_mm_s=config.global_1d.inlet.mean_velocity_mm_s,
        prescribed_flow_m3_s=config.global_1d.inlet.flow_rate_m3_s,
        relative_mass_tolerance=config.global_1d.solver.relative_mass_tolerance,
        relative_node_residual_tolerance=(
            config.global_1d.solver.relative_node_residual_tolerance
        ),
        reverse_flow_tolerance_m3_s=(
            config.global_1d.solver.reverse_flow_tolerance_m3_s
        ),
    )
    write_global_tables(layout, model, flow)
    solver_summary = _solver_summary(
        flow, leaf_pressure_pa=config.global_1d.leaf_pressure_pa
    )
    solver_summary["dynamic_viscosity_pa_s"] = mu
    solver_summary["node_count"] = model.node_count
    solver_summary["edge_count"] = model.edge_count
    write_json(layout.global_1d / "solver_summary.json", solver_summary)
    global_qc = global_graph_qc(model, flow)
    global_qc["relative_mass_tolerance"] = (
        config.global_1d.solver.relative_mass_tolerance
    )
    global_qc["relative_node_residual_tolerance"] = (
        config.global_1d.solver.relative_node_residual_tolerance
    )
    write_json(layout.qc / "global_1d_qc.json", global_qc)
    if config.outputs.save_global_vtp:
        write_global_vtp(layout.global_1d / "global_1d_flow.vtp", model, flow)

    transfers, boundary_mass_error = transfer_all_boundaries(
        roi,
        model,
        flow,
        mu_pa_s=mu,
        maximum_position_error_um=config.transfer.maximum_cut_position_error_um,
        maximum_radius_relative_error=config.transfer.maximum_cut_radius_relative_error,
        inlet_length_diameters=config.extension.inlet_length_diameters,
        outlet_length_diameters=config.extension.outlet_length_diameters,
    )
    write_port_classification(layout.roi / "port_classification.csv", transfers)
    boundary_qc = _boundary_qc(
        transfers,
        boundary_mass_error,
        mass_tolerance=config.transfer.relative_port_mass_tolerance,
        maximum_position_error_um=config.transfer.maximum_cut_position_error_um,
        maximum_radius_relative_error=(
            config.transfer.maximum_cut_radius_relative_error
        ),
    )
    write_json(layout.qc / "port_transfer_qc.json", boundary_qc)

    geometry, model_auto = resolve_model_run(
        config.paths.model_run,
        model_output_root=config.paths.model_output_root,
        roi_id=roi.roi_id,
    )
    geometry_report = geometry.report()
    write_json(layout.input / "geometry_reference.json", geometry_report)
    write_json(
        layout.input / "input_manifest.json",
        {
            "configuration": str(config.source_path),
            "sampling_run": str(sampling_run),
            "rodent_run": str(rodent_run),
            "rodent_run_auto_resolved": rodent_auto,
            "model_run": str(geometry.run_root),
            "model_run_auto_resolved": model_auto,
            "roi_id": roi.roi_id,
            "source_model_id": roi.source_model_id,
            "analysis_swc_radius_role": "SOURCE_RADIUS_FOR_1D_FLOW_NO_0.91_SCALING",
        },
    )

    readiness = roi_readiness(
        roi,
        transfers,
        boundary_mass_error,
        config.readiness,
        boundary_mass_tolerance=config.transfer.relative_port_mass_tolerance,
        maximum_position_error_um=config.transfer.maximum_cut_position_error_um,
        maximum_radius_relative_error=(
            config.transfer.maximum_cut_radius_relative_error
        ),
    )
    readiness["geometry_reference_pass"] = True
    readiness["global_qc_pass"] = global_qc["status"] == "PASS"
    final_status = (
        "CFD_PREPROCESS_BASELINE_PASS"
        if readiness["status"] == "PASS"
        else "CFD_ROI_NOT_READY"
    )
    next_stage = next_stage_for_status(final_status)
    readiness["final_status"] = final_status
    readiness["next_stage"] = next_stage
    write_json(layout.qc / "cfd_readiness.json", readiness)

    boundary_path: Path | None = None
    port_vtp_path: Path | None = None
    if final_status == "CFD_PREPROCESS_BASELINE_PASS":
        _, boundary_path, _ = write_boundary_package(layout, transfers, config)
        if config.extension.enabled:
            write_extension_plan(layout.roi / "port_extension_plan.csv", transfers)
        if config.outputs.save_port_vtp:
            port_vtp_path = write_port_vtp(layout.roi / "port_planes.vtp", transfers)

    if config.outputs.save_figures:
        global_pressure_figure(layout.figures / "global_1d_pressure.png", model, flow)
        roi_boundary_figure(
            layout.figures / "roi_boundary_conditions.png",
            geometry.surface_um_vtp,
            transfers,
            final_status=final_status,
        )
    summary = {
        "status": final_status,
        "next_stage": next_stage,
        "run_id": run_id,
        "run_root": str(layout.run_root.resolve()),
        "configuration": str(config.source_path),
        "roi_id": roi.roi_id,
        "source_model_id": roi.source_model_id,
        "global_node_count": model.node_count,
        "global_edge_count": model.edge_count,
        "structural_root_node_id": flow.root_node_id,
        "structural_leaf_count": len(flow.leaf_node_ids),
        "dynamic_viscosity_pa_s": mu,
        "root_mean_velocity_mm_s": flow.root_mean_velocity_mm_s,
        "root_flow_m3_s": flow.root_flow_m3_s,
        "global_relative_mass_error": flow.relative_mass_error,
        "cut_port_count": roi.cut_port_count,
        "true_terminal_count": roi.true_terminal_count,
        "total_boundary_count": len(transfers),
        "assumed_inlet_count": readiness["assumed_inlet_count"],
        "assumed_outlet_count": readiness["assumed_outlet_count"],
        "cut_port_outlet_count": readiness["cut_port_outlet_count"],
        "true_terminal_outlet_count": readiness["true_terminal_outlet_count"],
        "relative_boundary_mass_error": boundary_mass_error,
        "boundaries": boundary_qc["boundaries"],
        "geometry": geometry_report,
        "boundary_conditions_json": str(boundary_path) if boundary_path else None,
        "port_planes_vtp": str(port_vtp_path) if port_vtp_path else None,
    }
    write_json(layout.qc / "run_summary.json", summary)
    return CFDPreprocessResult(
        layout=layout,
        status=final_status,
        summary=summary,
        transfers=transfers,
        boundary_conditions_path=boundary_path,
        port_planes_path=port_vtp_path,
    )


def print_result(result: CFDPreprocessResult) -> None:
    summary = result.summary
    print("CFD PREPROCESS")
    print(f"ROI: {summary['roi_id']}")
    print(
        "Global: "
        f"nodes={summary['global_node_count']} | edges={summary['global_edge_count']} | "
        f"mass={summary['global_relative_mass_error']:.6e}"
    )
    print(
        f"Boundaries: inlet={summary['assumed_inlet_count']} | "
        f"outlet={summary['assumed_outlet_count']} | "
        f"mass={summary['relative_boundary_mass_error']:.6e}"
    )
    print(f"Boundary conditions: {summary['boundary_conditions_json']}")
    print(f"STATUS: {result.status}")
    print(f"NEXT: {next_stage_display(summary['next_stage'])}")
