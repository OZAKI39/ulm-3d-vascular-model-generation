"""Pure QC calculations and PROTEUS-facing VTU validation."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
import pyvista as pv
import trimesh

from .apes import BoundaryConditions
from .geometry import SurfacePartition
from .io import FlowError, write_json


def percentile_summary(values: np.ndarray, percentiles: tuple[int, ...]) -> dict[str, float]:
    array = np.asarray(values, dtype=float).reshape(-1)
    result = {"min": float(np.min(array)), "max": float(np.max(array))}
    for percentile in percentiles:
        result[f"p{percentile:02d}"] = float(np.percentile(array, percentile))
    return result


def evaluate_mass_conservation(q_in: float, q_out: tuple[float, ...]) -> dict[str, Any]:
    if q_in == 0.0:
        error = math.inf
    else:
        error = abs(q_in - sum(q_out)) / abs(q_in)
    return {
        "q_in_m3_s": float(q_in),
        "q_out_m3_s": [float(value) for value in q_out],
        "q_out_sum_m3_s": float(sum(q_out)),
        "relative_error": float(error),
        "flow_signs_pass": bool(q_in > 0.0 and all(value > 0.0 for value in q_out)),
    }


def numerical_port_fluxes(
    grid: pv.DataSet,
    partition: SurfacePartition,
    dx_m: float,
) -> tuple[dict[str, float], dict[str, float]]:
    """Integrate fields on each real cap translated inward by exactly 2dx.

    Three symmetric barycentric quadrature points are used per source cap
    triangle.  Cell data are deterministically converted to point data before
    VTK interpolation, avoiding the area bias of selecting a one-cell slab.
    """

    interpolator = grid.cell_data_to_point_data(pass_cell_data=True)
    fluxes: dict[str, float] = {}
    pressures: dict[str, float] = {}
    for patch in partition.patches:
        if patch.label == "wall":
            continue
        triangles = (
            partition.points_um[partition.faces[patch.face_indices]] * 1.0e-6
            - patch.outward_normal[None, None, :] * (2.0 * dx_m)
        )
        cross = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
        areas = 0.5 * np.linalg.norm(cross, axis=1)
        barycentric = np.asarray(
            (
                (2.0 / 3.0, 1.0 / 6.0, 1.0 / 6.0),
                (1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0),
                (1.0 / 6.0, 1.0 / 6.0, 2.0 / 3.0),
            )
        )
        sample_points = np.einsum("qv,tvc->tqc", barycentric, triangles).reshape(-1, 3)
        sampled = pv.PolyData(sample_points).sample(interpolator)
        valid = np.asarray(sampled.point_data.get("vtkValidPointMask", np.ones(len(sample_points))), dtype=bool)
        if not np.all(valid):
            raise FlowError(
                "CFD_FLOW_OUTPUT_INVALID",
                f"{patch.label} internal plane has {int(np.count_nonzero(~valid))} invalid samples",
            )
        velocity = np.asarray(sampled.point_data["velocity_phy"], dtype=float).reshape(-1, 3, 3)
        pressure = np.asarray(sampled.point_data["pressure_gauge_pa"], dtype=float).reshape(-1, 3)
        triangle_axial_velocity = np.mean(velocity @ patch.outward_normal, axis=1)
        outward_flux = float(np.sum(triangle_axial_velocity * areas))
        fluxes[patch.label] = -outward_flux if patch.label == "inlet" else outward_flux
        pressures[patch.label] = float(np.average(np.mean(pressure, axis=1), weights=areas))
    return fluxes, pressures


def fluid_domain_geometry_qc(
    grid: pv.UnstructuredGrid,
    partition: SurfacePartition,
    dx_m: float,
    tolerance_cells: float,
) -> dict[str, Any]:
    """Verify one connected fluid domain contained by the frozen lumen."""

    connected = grid.connectivity(extraction_mode="all")
    region_ids = np.asarray(connected.cell_data["RegionId"], dtype=np.int64)
    region_count = int(len(np.unique(region_ids)))

    centers_um = np.asarray(grid.cell_centers().points, dtype=float) * 1.0e6
    proximity = trimesh.proximity.ProximityQuery(partition.mesh_um)
    signed_chunks = []
    for start in range(0, len(centers_um), 10_000):
        signed_chunks.append(proximity.signed_distance(centers_um[start : start + 10_000]))
    signed_um = np.concatenate(signed_chunks) if signed_chunks else np.empty(0)
    tolerance_um = tolerance_cells * dx_m * 1.0e6
    beyond = signed_um < -tolerance_um
    strictly_inside = signed_um >= 0.0
    return {
        "status": "PASS" if region_count == 1 and not np.any(beyond) else "FAIL",
        "single_fluid_domain": region_count == 1,
        "connected_region_count": region_count,
        "fluid_cell_center_count": int(len(centers_um)),
        "strictly_inside_or_on_count": int(np.count_nonzero(strictly_inside)),
        "boundary_tolerance_count": int(np.count_nonzero((signed_um < 0.0) & ~beyond)),
        "outside_beyond_tolerance_count": int(np.count_nonzero(beyond)),
        "tolerance_cells": float(tolerance_cells),
        "tolerance_um": float(tolerance_um),
        "minimum_signed_distance_um": float(np.min(signed_um)) if len(signed_um) else None,
        "signed_distance_convention": "positive inside; negative outside",
    }


def validate_and_convert_flow_vtu(
    source: Path,
    destination: Path,
    *,
    pressure_reference_pa: float,
) -> tuple[pv.UnstructuredGrid, dict[str, Any]]:
    """Create a single meter-coordinate VTU with required cellData fields."""

    data = pv.read(source)
    conversion = "NONE_ALREADY_CELL_DATA"
    required = {"velocity_phy", "pressure_phy"}
    if not required.issubset(data.cell_data.keys()) and required.issubset(data.point_data.keys()):
        data = data.point_data_to_cell_data(pass_point_data=False)
        conversion = "DETERMINISTIC_POINT_TO_CELL_AVERAGE"
    if not required.issubset(data.cell_data.keys()):
        raise FlowError("CFD_FLOW_OUTPUT_INVALID", "Musubi output lacks pressure_phy/velocity_phy")
    velocity = np.asarray(data.cell_data["velocity_phy"], dtype=float)
    pressure = np.asarray(data.cell_data["pressure_phy"], dtype=float).reshape(-1)
    if velocity.ndim != 2 or velocity.shape != (data.n_cells, 3):
        raise FlowError("CFD_FLOW_OUTPUT_INVALID", "velocity_phy is not 3-component cellData")
    if pressure.shape != (data.n_cells,):
        raise FlowError("CFD_FLOW_OUTPUT_INVALID", "pressure_phy is not scalar cellData")
    finite = bool(np.all(np.isfinite(velocity)) and np.all(np.isfinite(pressure)))
    if not finite:
        raise FlowError("CFD_FLOW_OUTPUT_INVALID", "Flow field contains NaN or Inf")
    if float(np.max(np.abs(np.asarray(data.bounds)))) >= 1.0e-2:
        raise FlowError("CFD_FLOW_OUTPUT_INVALID", "VTU coordinates are not in meters")
    data.cell_data["pressure_gauge_pa"] = pressure - pressure_reference_pa
    grid = data.cast_to_unstructured_grid()
    celltypes = np.unique(np.asarray(grid.celltypes, dtype=np.uint8))
    allowed = {11, 12}
    topology_pass = set(int(value) for value in celltypes).issubset(allowed)
    if not topology_pass:
        raise FlowError("CFD_FLOW_OUTPUT_INVALID", f"Non-Cartesian cell types: {celltypes.tolist()}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    grid.save(destination, binary=True)
    speed = np.linalg.norm(velocity, axis=1)
    gauge = pressure - pressure_reference_pa
    qc = {
        "status": "PASS",
        "source": str(source),
        "output": str(destination),
        "coordinate_unit": "m",
        "point_to_cell_conversion": conversion,
        "cell_count": int(grid.n_cells),
        "point_count": int(grid.n_points),
        "cell_types": [int(value) for value in celltypes],
        "cell_topology": "UNIFORM_CARTESIAN_HEXAHEDRAL",
        "velocity_phy_cell_data": True,
        "velocity_phy_components": 3,
        "pressure_phy_cell_data": True,
        "pressure_gauge_pa_cell_data": True,
        "finite": finite,
        "nan_velocity": int(np.isnan(velocity).sum()),
        "inf_velocity": int(np.isinf(velocity).sum()),
        "nan_pressure": int(np.isnan(pressure).sum()),
        "inf_pressure": int(np.isinf(pressure).sum()),
        "velocity_m_s": percentile_summary(speed, (50, 95, 99)),
        "pressure_gauge_pa": percentile_summary(gauge, (1, 50, 99)),
    }
    return grid, qc


def reynolds_diagnostics(
    fluxes: dict[str, float], partition: SurfacePartition, nu_m2_s: float
) -> dict[str, float]:
    result: dict[str, float] = {}
    for patch in partition.patches:
        if patch.label == "wall":
            continue
        area = patch.area_um2 * 1.0e-12
        diameter = 2.0 * patch.equivalent_radius_um * 1.0e-6
        result[patch.label] = abs(fluxes[patch.label]) / area * diameter / nu_m2_s
    return result


def write_proteus_metadata(
    path: Path,
    *,
    inlet_equivalent_diameter_m: float,
    source_flow_vtu: Path,
) -> dict[str, Any]:
    metadata = {
        "status": "PROTEUS_IMPORT_READY",
        "lengthUnit": 1.0,
        "velocityUnit": 1.0,
        "velocityField": "velocity_phy",
        "pressureField": "pressure_gauge_pa",
        "flowFieldVTU": str(source_flow_vtu),
        "inletEquivalentDiameter": float(inlet_equivalent_diameter_m),
        "inletNormal": None,
        "inletNormalSource": "AUTO_DETECT_BY_BACKPROPAGATION",
        "surfaceGeometryModified": False,
        "microbubbleSimulationRun": False,
        "backpropagationRun": False,
    }
    write_json(path, metadata)
    return metadata


def boundary_condition_qc(bc: BoundaryConditions) -> dict[str, Any]:
    return {
        "status": "PASS",
        "source_of_truth": "boundary_conditions_vmtk_boundarynormal_crossseam.json",
        "fluid": {
            "density_kg_m3": bc.density_kg_m3,
            "kinematic_viscosity_m2_s": bc.kinematic_viscosity_m2_s,
            "dynamic_viscosity_pa_s": bc.dynamic_viscosity_pa_s,
        },
        "inlet": {
            "port_id": bc.inlet_port_id,
            "kind": "mfr_eq",
            "requested_flow_m3_s": bc.inlet_flow_m3_s,
            "mass_flow_kg_s": bc.density_kg_m3 * bc.inlet_flow_m3_s,
            "requested_profile": bc.inlet_profile_requested,
            "effective_profile": "MFR_EQ_NATIVE",
            "profile_exactly_preserved": False,
            "flow_rate_exactly_requested": True,
            "priority": "EXACT_VOLUMETRIC_FLOW",
        },
        "outlets": [
            {"port_id": port_id, "label": f"outlet_{index:02d}", "kind": "pressure_eq", "P_solver_boundary_pa_gauge": pressure}
            for index, (port_id, pressure) in enumerate(zip(bc.outlet_port_ids, bc.outlet_gauge_pressures_pa, strict=True), start=1)
        ],
        "wall": {"kind": "wall_libb", "rigid": True, "no_slip": True},
    }
