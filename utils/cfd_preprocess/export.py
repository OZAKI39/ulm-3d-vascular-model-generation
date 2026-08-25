"""CSV, JSON, and VTP exports for the solver-independent baseline package."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pyvista as pv

from utils.sampling.sampling_types import GlobalVascularModel

from .config import CFDPreprocessConfig
from .one_d_flow import GlobalFlowResult
from .port_transfer import PortTransfer


@dataclass(frozen=True, slots=True)
class OutputLayout:
    run_root: Path
    config: Path
    input: Path
    logs: Path
    global_1d: Path
    roi: Path
    qc: Path
    figures: Path
    report: Path


def create_layout(output_root: Path, run_id: str) -> OutputLayout:
    root = Path(output_root) / run_id
    folders = {
        name: root / name
        for name in (
            "config",
            "input",
            "logs",
            "global_1d",
            "roi",
            "qc",
            "figures",
            "report",
        )
    }
    for folder in folders.values():
        folder.mkdir(parents=True, exist_ok=False)
    return OutputLayout(root, **folders)


def write_json(path: Path, payload: Any) -> Path:
    path.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return path


def write_csv(
    path: Path,
    rows: Iterable[dict[str, Any]],
    *,
    fieldnames: Iterable[str] | None = None,
) -> Path:
    materialized = list(rows)
    names = list(fieldnames or (materialized[0].keys() if materialized else []))
    with path.open("w", encoding="utf-8-sig", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=names, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(materialized)
    return path


def write_global_tables(
    layout: OutputLayout,
    model: GlobalVascularModel,
    flow: GlobalFlowResult,
) -> tuple[Path, Path]:
    root_set = {flow.root_node_id}
    leaf_set = set(map(int, flow.leaf_node_ids))
    nodes = []
    for index, node_id in enumerate(model.node_ids):
        x, y, z = model.node_positions_um[index]
        pressure = flow.pressures_pa[index]
        nodes.append(
            {
                "node_id": int(node_id),
                "parent_id": int(model.parent_ids[index]),
                "x_um": x,
                "y_um": y,
                "z_um": z,
                "radius_um": model.node_radius_um[index],
                "pressure_pa": pressure,
                "pressure_mmhg": pressure / 133.32236842105263,
                "is_root": int(int(node_id) in root_set),
                "is_leaf": int(int(node_id) in leaf_set),
            }
        )
    edges = []
    for edge in model.edges:
        edge_id = edge.edge_id
        flow_rate = flow.flows_m3_s[edge_id]
        parent_area = np.pi * (edge.upstream_radius_um * 1.0e-6) ** 2
        child_area = np.pi * (edge.downstream_radius_um * 1.0e-6) ** 2
        edges.append(
            {
                "global_edge_id": edge_id,
                "parent_node_id": edge.upstream_node_id,
                "child_node_id": edge.downstream_node_id,
                "length_um": flow.edge_lengths_um[edge_id],
                "radius_parent_um": edge.upstream_radius_um,
                "radius_child_um": edge.downstream_radius_um,
                "resistance_pa_s_m3": flow.resistances_pa_s_m3[edge_id],
                "flow_rate_m3_s": flow_rate,
                "flow_rate_pL_s": flow_rate * 1.0e15,
                "mean_velocity_parent_mm_s": flow_rate / parent_area * 1.0e3,
                "mean_velocity_child_mm_s": flow_rate / child_area * 1.0e3,
                "flow_sign_parent_to_child": int(np.sign(flow_rate)),
            }
        )
    return (
        write_csv(layout.global_1d / "nodes.csv", nodes),
        write_csv(layout.global_1d / "edges.csv", edges),
    )


def write_global_vtp(
    path: Path,
    model: GlobalVascularModel,
    flow: GlobalFlowResult,
) -> Path:
    lines = np.asarray(
        [
            (
                2,
                model.node_index_by_id[edge.upstream_node_id],
                model.node_index_by_id[edge.downstream_node_id],
            )
            for edge in model.edges
        ],
        dtype=np.int64,
    ).ravel()
    mesh = pv.PolyData(np.asarray(model.node_positions_um, dtype=float), lines=lines)
    leaf_set = set(map(int, flow.leaf_node_ids))
    mesh.point_data["node_id"] = np.asarray(model.node_ids, dtype=np.int64)
    mesh.point_data["pressure_pa"] = flow.pressures_pa
    mesh.point_data["radius_um"] = np.asarray(model.node_radius_um, dtype=float)
    mesh.point_data["is_root"] = (model.node_ids == flow.root_node_id).astype(np.uint8)
    mesh.point_data["is_leaf"] = np.asarray(
        [int(int(node_id) in leaf_set) for node_id in model.node_ids], dtype=np.uint8
    )
    mesh.cell_data["global_edge_id"] = flow.edge_ids
    mesh.cell_data["resistance_pa_s_m3"] = flow.resistances_pa_s_m3
    mesh.cell_data["flow_rate_m3_s"] = flow.flows_m3_s
    mesh.cell_data["flow_rate_pL_s"] = flow.flows_m3_s * 1.0e15
    mesh.save(path)
    return path


def port_rows(transfers: list[PortTransfer]) -> list[dict[str, Any]]:
    rows = []
    for item in transfers:
        x, y, z = item.center_um
        tangent = item.geometry.simulation_tangent
        normal = item.geometry.outward_normal
        end = item.geometry.extension_end_um
        nominal_area = np.pi * (item.source_radius_um * 1.0e-6) ** 2
        mean_velocity = item.role_flow_m3_s / nominal_area * 1.0e3
        bc_type = (
            "VOLUMETRIC_FLOW_RATE"
            if item.role == "ASSUMED_INLET"
            else "PRESSURE_DIRICHLET"
        )
        bc_value = (
            item.role_flow_m3_s if item.role == "ASSUMED_INLET" else item.pressure_pa
        )
        rows.append(
            {
                "port_id": item.port_id,
                "role": item.role,
                "boundary_origin": item.boundary_origin,
                "local_node_id": item.local_node_id,
                "global_node_id": item.global_node_id,
                "global_edge_id": item.global_edge_id,
                "x_um": x,
                "y_um": y,
                "z_um": z,
                "radius_um": item.source_radius_um,
                "diameter_um": 2.0 * item.source_radius_um,
                "alpha_on_global_edge": item.alpha_on_global_edge,
                "position_error_um": item.position_error_um,
                "radius_relative_error": item.radius_relative_error,
                "P_1D_pa": item.pressure_pa,
                "P_1D_mmHg": item.pressure_mmhg,
                "Q_1D_m3_s": item.role_flow_m3_s,
                "Q_1D_pL_s": item.flow_pl_s,
                "signed_parent_to_child_flow_m3_s": item.signed_parent_to_child_flow_m3_s,
                "bc_type": bc_type,
                "bc_value_si": bc_value,
                "nominal_cross_section_area_m2": nominal_area,
                "nominal_mean_velocity_mm_s": mean_velocity,
                "simulation_tangent_x": tangent[0],
                "simulation_tangent_y": tangent[1],
                "simulation_tangent_z": tangent[2],
                "outward_normal_x": normal[0],
                "outward_normal_y": normal[1],
                "outward_normal_z": normal[2],
                "extension_length_um": item.geometry.extension_length_um,
                "extension_end_x_um": end[0],
                "extension_end_y_um": end[1],
                "extension_end_z_um": end[2],
            }
        )
    return rows


def write_port_classification(path: Path, transfers: list[PortTransfer]) -> Path:
    return write_csv(path, port_rows(transfers))


def write_boundary_package(
    layout: OutputLayout,
    transfers: list[PortTransfer],
    config: CFDPreprocessConfig,
) -> tuple[Path, Path, Path]:
    rows = port_rows(transfers)
    boundary_csv = write_csv(layout.roi / "boundary_conditions.csv", rows)
    inlet = next(item for item in transfers if item.role == "ASSUMED_INLET")
    outlets = [item for item in transfers if item.role == "ASSUMED_OUTLET"]
    inlet_area = np.pi * (inlet.source_radius_um * 1.0e-6) ** 2
    inlet_mean = inlet.role_flow_m3_s / inlet_area
    boundary_json = write_json(
        layout.roi / "boundary_conditions.json",
        {
            "method": "GLOBAL_1D_TO_ROI_DIRECT_PRESSURE_BASELINE",
            "simulation_direction": {
                "basis": "SWC_PARENT_TO_CURRENT",
                "is_measured": False,
                "is_physiological_ground_truth": False,
            },
            "fluid": {
                "model": "NEWTONIAN",
                "density_kg_m3": config.fluid.density_kg_m3,
                "kinematic_viscosity_m2_s": config.fluid.kinematic_viscosity_m2_s,
                "dynamic_viscosity_pa_s": config.fluid.dynamic_viscosity_pa_s,
            },
            "flow_model": {
                "steady": True,
                "incompressible": True,
                "rigid_wall": True,
                "no_slip_wall": True,
                "turbulence": False,
                "pulsatility": False,
            },
            "inlet": {
                "port_id": inlet.port_id,
                "role": inlet.role,
                "boundary_origin": inlet.boundary_origin,
                "global_node_id": inlet.global_node_id,
                "global_edge_id": inlet.global_edge_id,
                "type": "VOLUMETRIC_FLOW_RATE",
                "flow_rate_m3_s": inlet.role_flow_m3_s,
                "profile": "PARABOLIC",
                "nominal_area_m2": inlet_area,
                "nominal_mean_velocity_m_s": inlet_mean,
                "nominal_centerline_max_velocity_m_s": 2.0 * inlet_mean,
                "actual_mesh_area_required_before_solver_setup": True,
            },
            "outlets": [
                {
                    "port_id": item.port_id,
                    "role": item.role,
                    "boundary_origin": item.boundary_origin,
                    "global_node_id": item.global_node_id,
                    "global_edge_id": item.global_edge_id,
                    "type": "PRESSURE_DIRICHLET",
                    "pressure_pa": item.pressure_pa,
                    "expected_1d_flow_m3_s": item.role_flow_m3_s,
                }
                for item in outlets
            ],
            "wall": {"type": "NO_SLIP"},
            "pressure_reference": "GLOBAL_STRUCTURAL_LEAVES_ZERO_GAUGE",
        },
    )
    physics = write_json(
        layout.roi / "physics_baseline.json",
        {
            "steady": True,
            "incompressible": True,
            "blood_model": "newtonian",
            "density_kg_m3": config.fluid.density_kg_m3,
            "kinematic_viscosity_m2_s": config.fluid.kinematic_viscosity_m2_s,
            "dynamic_viscosity_pa_s": config.fluid.dynamic_viscosity_pa_s,
            "wall": "rigid_no_slip",
            "inlet": "Q_1D_plus_parabolic_profile",
            "outlet": "direct_P_1D",
            "turbulence": "none",
            "pulsatility": False,
            "version2_physiology": False,
        },
    )
    return boundary_csv, boundary_json, physics


def write_extension_plan(path: Path, transfers: list[PortTransfer]) -> Path:
    fields = [
        "port_id",
        "role",
        "boundary_origin",
        "global_node_id",
        "global_edge_id",
        "extension_length_um",
        "extension_end_x_um",
        "extension_end_y_um",
        "extension_end_z_um",
        "surface_modified",
    ]
    rows = []
    for item in transfers:
        end = item.geometry.extension_end_um
        rows.append(
            {
                "port_id": item.port_id,
                "role": item.role,
                "boundary_origin": item.boundary_origin,
                "global_node_id": item.global_node_id,
                "global_edge_id": item.global_edge_id,
                "extension_length_um": item.geometry.extension_length_um,
                "extension_end_x_um": end[0],
                "extension_end_y_um": end[1],
                "extension_end_z_um": end[2],
                "surface_modified": False,
            }
        )
    return write_csv(path, rows, fieldnames=fields)


def write_port_vtp(
    path: Path, transfers: list[PortTransfer], *, resolution: int = 48
) -> Path:
    all_points: list[np.ndarray] = []
    faces: list[int] = []
    cell_port_id: list[str] = []
    cell_role: list[int] = []
    cell_origin: list[int] = []
    cell_global_node: list[int] = []
    cell_global_edge: list[int] = []
    cell_pressure: list[float] = []
    cell_flow: list[float] = []
    cell_radius: list[float] = []
    offset = 0
    for item in transfers:
        normal = item.geometry.outward_normal
        helper = np.array([1.0, 0.0, 0.0])
        if abs(float(np.dot(helper, normal))) > 0.9:
            helper = np.array([0.0, 1.0, 0.0])
        basis1 = np.cross(normal, helper)
        basis1 /= np.linalg.norm(basis1)
        basis2 = np.cross(normal, basis1)
        angles = np.linspace(0.0, 2.0 * np.pi, resolution, endpoint=False)
        ring = item.center_um + item.source_radius_um * (
            np.cos(angles)[:, None] * basis1 + np.sin(angles)[:, None] * basis2
        )
        points = np.vstack((item.center_um, ring))
        all_points.append(points)
        for index in range(resolution):
            faces.extend(
                (3, offset, offset + 1 + index, offset + 1 + (index + 1) % resolution)
            )
            cell_port_id.append(item.port_id)
            cell_role.append(1 if item.role == "ASSUMED_INLET" else 2)
            cell_origin.append(1 if item.boundary_origin == "CUT_PORT" else 2)
            cell_global_node.append(
                item.global_node_id if item.global_node_id is not None else -1
            )
            cell_global_edge.append(item.global_edge_id)
            cell_pressure.append(item.pressure_pa)
            cell_flow.append(item.role_flow_m3_s)
            cell_radius.append(item.source_radius_um)
        offset += len(points)
    mesh = pv.PolyData(np.vstack(all_points), faces=np.asarray(faces, dtype=np.int64))
    mesh.cell_data["port_id"] = np.asarray(cell_port_id)
    mesh.cell_data["role_code"] = np.asarray(cell_role, dtype=np.uint8)
    mesh.cell_data["boundary_origin_code"] = np.asarray(cell_origin, dtype=np.uint8)
    mesh.cell_data["global_node_id"] = np.asarray(cell_global_node, dtype=np.int64)
    mesh.cell_data["global_edge_id"] = np.asarray(cell_global_edge, dtype=np.int64)
    mesh.cell_data["pressure_pa"] = np.asarray(cell_pressure, dtype=float)
    mesh.cell_data["flow_rate_m3_s"] = np.asarray(cell_flow, dtype=float)
    mesh.cell_data["radius_um"] = np.asarray(cell_radius, dtype=float)
    mesh.save(path)
    return path
