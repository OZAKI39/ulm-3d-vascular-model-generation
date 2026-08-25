"""Orchestration for formal global-to-ROI CFD boundary preprocessing."""

from __future__ import annotations

import logging
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
from .port_transfer import PortTransfer, transfer_all_ports
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


def _logger(path: Path, verbose: bool) -> logging.Logger:
    logger = logging.getLogger(f"cfd_preprocess.{path.parent.parent.name}")
    logger.handlers.clear()
    logger.setLevel(logging.DEBUG if verbose else logging.INFO)
    handler = logging.FileHandler(path, encoding="utf-8")
    handler.setFormatter(logging.Formatter("%(asctime)s | %(levelname)s | %(message)s"))
    logger.addHandler(handler)
    logger.propagate = False
    return logger


def _solver_summary(flow: GlobalFlowResult) -> dict[str, Any]:
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
        "leaf_pressure_pa": 0.0,
        "leaf_pressure_role": "GAUGE_PRESSURE_REFERENCE_NOT_PHYSIOLOGICAL_ZERO",
        "radius_source": "ORIGINAL_ANALYSIS_SWC_SOURCE_RADIUS_NOT_RADIUS_SCALE_0.91",
    }


def _port_qc(
    transfers: list[PortTransfer], port_mass_error: float, tolerance: float
) -> dict[str, Any]:
    return {
        "status": "PASS" if port_mass_error <= tolerance else "FAIL",
        "all_global_edges_found": True,
        "all_positions_match": True,
        "all_radii_match": True,
        "pressure_interpolation": "EXACT_PARTIAL_LINEAR_RADIUS_HYDRAULIC_RESISTANCE",
        "relative_port_mass_error": port_mass_error,
        "relative_port_mass_tolerance": tolerance,
        "ports": [
            {
                "port_id": item.port_id,
                "role": item.role,
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


def _write_report(
    path: Path,
    *,
    status: str,
    roi_id: str,
    source_model_id: str,
    flow: GlobalFlowResult,
    transfers: list[PortTransfer],
    readiness: dict[str, Any],
    geometry: dict[str, Any],
) -> Path:
    port_lines = "\n".join(
        f"- `{item.port_id}` — {item.role}; P={item.pressure_pa:.9g} Pa; "
        f"Q={item.role_flow_m3_s:.9g} m³/s"
        for item in transfers
    )
    reasons = readiness.get("failure_reasons") or []
    reason_text = ", ".join(reasons) if reasons else "none"
    text = f"""# CFD boundary preprocessing report

Final status: **{status}**

ROI: `{roi_id}`

Global source model: `{source_model_id}`

## Formal baseline and its limits

SWC parent→current direction is used only as the **simulation direction**. It is not an
experimentally measured blood-flow direction and is not physiological ground truth. The single
structural root is temporarily treated as `ASSUMED_GLOBAL_INLET`. Its configured
{flow.root_mean_velocity_mm_s:.9g} mm/s velocity is a literature-derived baseline assumption,
not a measurement from the current fMOST mouse.

All structural leaves are assigned 0 Pa gauge pressure as the global baseline reference. This
does not mean that physiological blood pressure is zero. The 1D calculation uses original
analysis-SWC source radii; the Ultraliser radius scale of 0.91 is only a surface-reconstruction
compensation and is not used in hydraulic resistance.

ROI inlet flow and outlet pressure are transferred from the global 1D model. The intended outlet
condition is **DIRECT 1D PRESSURE**, not a resistance or Windkessel condition. Version-2
physiological refinements are disabled.

## Global solution

- Structural root: {flow.root_node_id}
- Structural leaves: {len(flow.leaf_node_ids)}
- Root flow: {flow.root_flow_m3_s:.12g} m³/s
- Global relative mass error: {flow.relative_mass_error:.6e}
- Maximum internal relative residual: {flow.maximum_internal_relative_residual:.6e}

## CUT_PORT transfer

{port_lines}

ROI relative port mass error: {readiness["relative_port_mass_error"]:.6e}

Readiness failure reasons: {reason_text}

## Geometry reference

- Model run: `{geometry["run_id"]}`
- Geometry status: {geometry["status"]}
- Radius scale: {geometry["radius_scale"]}
- Radius P95 absolute relative error: {geometry["radius_p95_absolute_relative_error"]:.9g}

No surface was modified, no volume mesh was created, and no 3D CFD or microbubble simulation was
run in this stage.
"""
    path.write_text(text, encoding="utf-8")
    return path


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
    logger = _logger(layout.logs / "cfd_preprocess.log", config.verbose)
    logger.info("Configuration: %s", config.source_path)
    logger.info("Sampling run: %s", sampling_run)
    logger.info("ROI: %s", roi.roi_id)

    rodent_root = project_root / "outputs" / "rodent_vasculature"
    rodent_run, rodent_auto = resolve_rodent_run(
        config.paths.rodent_run,
        output_root=rodent_root,
        source_model_id=roi.source_model_id,
    )
    if rodent_auto:
        logger.info("AUTO_RESOLVED_RODENT_RUN: %s", rodent_run)
    model = load_matching_global_model(rodent_run, roi.source_model_id)
    verify_global_edge_manifest(model, sampling_run / "manifests" / "global_edges.csv")
    logger.info("GLOBAL_EDGE_MAPPING_MATCH: %d edges", model.edge_count)

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
    logger.info(
        "Global solve PASS: nodes=%d edges=%d mass_error=%.3e",
        model.node_count,
        model.edge_count,
        flow.relative_mass_error,
    )
    write_global_tables(layout, model, flow)
    solver_summary = _solver_summary(flow)
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

    transfers, port_mass_error = transfer_all_ports(
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
    port_qc = _port_qc(
        transfers,
        port_mass_error,
        config.transfer.relative_port_mass_tolerance,
    )
    write_json(layout.qc / "port_transfer_qc.json", port_qc)

    geometry, model_auto = resolve_model_run(
        config.paths.model_run,
        model_output_root=config.paths.model_output_root,
        roi_id=roi.roi_id,
    )
    if model_auto:
        logger.info("AUTO_RESOLVED_MODEL_RUN: %s", geometry.run_root)
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
        port_mass_error,
        config.readiness,
        port_mass_tolerance=config.transfer.relative_port_mass_tolerance,
    )
    readiness["geometry_reference_pass"] = True
    readiness["global_qc_pass"] = global_qc["status"] == "PASS"
    final_status = (
        "CFD_PREPROCESS_BASELINE_PASS"
        if readiness["status"] == "PASS"
        else "CFD_ROI_NOT_READY"
    )
    readiness["final_status"] = final_status
    write_json(layout.qc / "cfd_readiness.json", readiness)

    boundary_path: Path | None = None
    port_vtp_path: Path | None = None
    if final_status == "CFD_PREPROCESS_BASELINE_PASS":
        _, boundary_path, _ = write_boundary_package(layout, transfers, config)
        if config.extension.enabled:
            write_extension_plan(layout.roi / "port_extension_plan.csv", transfers)
        if config.outputs.save_port_vtp:
            port_vtp_path = write_port_vtp(layout.roi / "port_planes.vtp", transfers)
    else:
        logger.warning(
            "CFD_ROI_NOT_READY: boundary package not generated; reasons=%s",
            ", ".join(readiness["failure_reasons"]),
        )

    if config.outputs.save_figures:
        global_pressure_figure(layout.figures / "global_1d_pressure.png", model, flow)
        roi_boundary_figure(
            layout.figures / "roi_boundary_conditions.png",
            geometry.surface_um_vtp,
            transfers,
            final_status=final_status,
        )
    _write_report(
        layout.report / "cfd_preprocess_report.md",
        status=final_status,
        roi_id=roi.roi_id,
        source_model_id=roi.source_model_id,
        flow=flow,
        transfers=transfers,
        readiness=readiness,
        geometry=geometry_report,
    )
    summary = {
        "status": final_status,
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
        "assumed_inlet_count": readiness["assumed_inlet_count"],
        "assumed_outlet_count": readiness["assumed_outlet_count"],
        "relative_port_mass_error": port_mass_error,
        "ports": port_qc["ports"],
        "geometry": geometry_report,
        "boundary_conditions_json": str(boundary_path) if boundary_path else None,
        "port_planes_vtp": str(port_vtp_path) if port_vtp_path else None,
        "upstream_regeneration_performed": False,
        "surface_modified": False,
        "volume_mesh_created": False,
        "three_dimensional_cfd_run": False,
    }
    write_json(layout.qc / "run_summary.json", summary)
    logger.info("Final status: %s", final_status)
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
    print(f"Configuration: {summary['configuration']}")
    print(f"ROI: {summary['roi_id']}")
    print(f"Global source model: {summary['source_model_id']}")
    print(f"Global nodes: {summary['global_node_count']}")
    print(f"Global edges: {summary['global_edge_count']}")
    print(f"Structural root: {summary['structural_root_node_id']}")
    print(f"Global leaves: {summary['structural_leaf_count']}")
    print(
        f"Global inlet baseline: {summary['root_mean_velocity_mm_s']:.9g} mm/s (NOT MEASURED)"
    )
    print(f"Q_root: {summary['root_flow_m3_s']:.12g} m3/s")
    print(f"Global mass error: {summary['global_relative_mass_error']:.6e}")
    print(f"CUT_PORT count: {summary['cut_port_count']}")
    print(f"ASSUMED_INLET count: {summary['assumed_inlet_count']}")
    print(f"ASSUMED_OUTLET count: {summary['assumed_outlet_count']}")
    print(f"TRUE_TERMINAL count: {summary['true_terminal_count']}")
    for port in summary["ports"]:
        print(
            f"Port: {port['port_id']} | {port['role']} | "
            f"P_1D={port['pressure_pa']:.9g} Pa | Q_1D={port['flow_rate_m3_s']:.9g} m3/s"
        )
    print(f"ROI port mass error: {summary['relative_port_mass_error']:.6e}")
    geometry = summary["geometry"]
    print(
        "Ultraliser geometry: PASS | "
        f"radius scale={geometry['radius_scale']} | "
        f"radius P95={geometry['radius_p95_absolute_relative_error']:.9g}"
    )
    print(f"Final status: {result.status}")
