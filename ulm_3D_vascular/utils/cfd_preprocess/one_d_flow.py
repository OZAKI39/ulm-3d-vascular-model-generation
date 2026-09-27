"""Sparse Newtonian resistor-network solver over source SWC edges."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve

from utils.sampling.sampling_types import GlobalVascularModel


class GlobalFlowError(RuntimeError):
    """A global baseline network or sparse-solve failure."""


@dataclass(frozen=True, slots=True)
class GlobalFlowResult:
    node_ids: np.ndarray
    parent_ids: np.ndarray
    pressures_pa: np.ndarray
    edge_ids: np.ndarray
    edge_lengths_um: np.ndarray
    resistances_pa_s_m3: np.ndarray
    flows_m3_s: np.ndarray
    root_node_id: int
    leaf_node_ids: np.ndarray
    root_flow_m3_s: float
    root_mean_velocity_mm_s: float
    leaf_flow_sum_m3_s: float
    relative_mass_error: float
    maximum_internal_relative_residual: float
    reverse_flow_count: int

    @property
    def pressure_by_node_id(self) -> dict[int, float]:
        return dict(zip(self.node_ids.tolist(), self.pressures_pa.tolist()))

    @property
    def flow_by_edge_id(self) -> dict[int, float]:
        return dict(zip(self.edge_ids.tolist(), self.flows_m3_s.tolist()))


def edge_resistance(
    length_m: float, radius0_m: float, radius1_m: float, mu_pa_s: float
) -> float:
    """Exact resistance for an edge whose radius varies linearly in arc length."""

    values = np.asarray((length_m, radius0_m, radius1_m, mu_pa_s), dtype=float)
    if (
        not np.all(np.isfinite(values))
        or length_m < 0
        or min(radius0_m, radius1_m, mu_pa_s) <= 0
    ):
        raise ValueError(
            "Length must be non-negative and radii/viscosity finite positive"
        )
    if length_m == 0:
        return 0.0
    if np.isclose(radius0_m, radius1_m, rtol=1.0e-10, atol=0.0):
        integral = length_m / radius0_m**4
    else:
        integral = (
            length_m
            / (3.0 * (radius1_m - radius0_m))
            * (1.0 / radius0_m**3 - 1.0 / radius1_m**3)
        )
    result = 8.0 * mu_pa_s / np.pi * integral
    if not np.isfinite(result) or result <= 0:
        raise ValueError("Edge resistance is not finite positive")
    return float(result)


def _root_and_leaves(model: GlobalVascularModel) -> tuple[int, np.ndarray]:
    roots = model.node_ids[np.asarray(model.parent_ids, dtype=np.int64) == -1]
    if len(roots) != 1:
        raise GlobalFlowError("GLOBAL_BASELINE_REQUIRES_SINGLE_ROOT")
    out_degree = {int(node_id): 0 for node_id in model.node_ids}
    for edge in model.edges:
        out_degree[edge.upstream_node_id] += 1
    root = int(roots[0])
    leaves = np.asarray(
        [
            node_id
            for node_id, degree in out_degree.items()
            if degree == 0 and node_id != root
        ],
        dtype=np.int64,
    )
    if not len(leaves):
        raise GlobalFlowError("GLOBAL_1D_FLOW_FAILED: no structural leaf outlets")
    return root, leaves


def solve_global_flow(
    model: GlobalVascularModel,
    *,
    mu_pa_s: float,
    leaf_pressure_pa: float,
    boundary_type: str,
    mean_velocity_mm_s: float | None,
    prescribed_flow_m3_s: float | None,
    relative_mass_tolerance: float,
    relative_node_residual_tolerance: float,
    reverse_flow_tolerance_m3_s: float,
) -> GlobalFlowResult:
    """Solve source-edge pressures and parent-to-current flows using SciPy sparse."""

    root, leaves = _root_and_leaves(model)
    if model.edge_count != model.node_count - 1:
        raise GlobalFlowError("GLOBAL_1D_FLOW_FAILED: global graph is not cycle-free")
    root_index = model.node_index_by_id[root]
    root_radius_m = float(model.node_radius_um[root_index]) * 1.0e-6
    if boundary_type == "mean_velocity":
        assert mean_velocity_mm_s is not None and prescribed_flow_m3_s is None
        root_flow = mean_velocity_mm_s * 1.0e-3 * np.pi * root_radius_m**2
        root_velocity = mean_velocity_mm_s
    elif boundary_type == "flow_rate":
        assert prescribed_flow_m3_s is not None and mean_velocity_mm_s is None
        root_flow = prescribed_flow_m3_s
        root_velocity = root_flow / (np.pi * root_radius_m**2) * 1.0e3
    else:
        raise GlobalFlowError("GLOBAL_1D_FLOW_FAILED: unsupported root boundary type")
    if not np.isfinite(root_flow) or root_flow <= 0:
        raise GlobalFlowError("GLOBAL_1D_FLOW_FAILED: root flow is not finite positive")

    edge_lengths_um = np.empty(model.edge_count, dtype=float)
    resistances = np.empty(model.edge_count, dtype=float)
    for edge in model.edges:
        length_um = float(
            np.linalg.norm(edge.downstream_position_um - edge.upstream_position_um)
        )
        edge_lengths_um[edge.edge_id] = length_um
        resistances[edge.edge_id] = edge_resistance(
            length_um * 1.0e-6,
            edge.upstream_radius_um * 1.0e-6,
            edge.downstream_radius_um * 1.0e-6,
            mu_pa_s,
        )

    leaf_set = set(map(int, leaves))
    unknown_ids = [
        int(node_id) for node_id in model.node_ids if int(node_id) not in leaf_set
    ]
    equation_index = {node_id: index for index, node_id in enumerate(unknown_ids)}
    rows: list[int] = []
    columns: list[int] = []
    data: list[float] = []
    rhs = np.zeros(len(unknown_ids), dtype=float)
    rhs[equation_index[root]] = root_flow
    for edge in model.edges:
        a = edge.upstream_node_id
        b = edge.downstream_node_id
        conductance = 1.0 / resistances[edge.edge_id]
        for node, neighbor in ((a, b), (b, a)):
            if node in leaf_set:
                continue
            row = equation_index[node]
            rows.append(row)
            columns.append(row)
            data.append(conductance)
            if neighbor in leaf_set:
                rhs[row] += conductance * leaf_pressure_pa
            else:
                rows.append(row)
                columns.append(equation_index[neighbor])
                data.append(-conductance)
    matrix = coo_matrix((data, (rows, columns)), shape=(len(unknown_ids),) * 2).tocsr()
    pressures_unknown = np.asarray(spsolve(matrix, rhs), dtype=float)
    pressures = np.full(model.node_count, leaf_pressure_pa, dtype=float)
    for node_id, pressure in zip(unknown_ids, pressures_unknown):
        pressures[model.node_index_by_id[node_id]] = pressure
    if not np.all(np.isfinite(pressures)):
        raise GlobalFlowError(
            "GLOBAL_1D_FLOW_FAILED: sparse solver returned non-finite pressure"
        )

    flows = np.empty(model.edge_count, dtype=float)
    node_flux = np.zeros(model.node_count, dtype=float)
    leaf_inflows = 0.0
    for edge in model.edges:
        parent_index = model.node_index_by_id[edge.upstream_node_id]
        child_index = model.node_index_by_id[edge.downstream_node_id]
        flow = (pressures[parent_index] - pressures[child_index]) / resistances[
            edge.edge_id
        ]
        flows[edge.edge_id] = flow
        node_flux[parent_index] += flow
        node_flux[child_index] -= flow
        if edge.downstream_node_id in leaf_set:
            leaf_inflows += flow
    mass_error = abs(root_flow - leaf_inflows) / abs(root_flow)
    internal_mask = np.ones(model.node_count, dtype=bool)
    internal_mask[root_index] = False
    internal_mask[[model.node_index_by_id[int(node)] for node in leaves]] = False
    maximum_residual = (
        float(np.max(np.abs(node_flux[internal_mask])) / abs(root_flow))
        if np.any(internal_mask)
        else 0.0
    )
    reverse_count = int(np.count_nonzero(flows < -reverse_flow_tolerance_m3_s))
    if mass_error > relative_mass_tolerance:
        raise GlobalFlowError(f"GLOBAL_1D_FLOW_FAILED: mass error {mass_error:.3e}")
    if maximum_residual > relative_node_residual_tolerance:
        raise GlobalFlowError(
            f"GLOBAL_1D_FLOW_FAILED: internal residual {maximum_residual:.3e}"
        )
    if reverse_count:
        raise GlobalFlowError(
            f"GLOBAL_1D_FLOW_FAILED: reverse edge count {reverse_count}"
        )
    return GlobalFlowResult(
        node_ids=np.asarray(model.node_ids, dtype=np.int64),
        parent_ids=np.asarray(model.parent_ids, dtype=np.int64),
        pressures_pa=pressures,
        edge_ids=np.asarray([edge.edge_id for edge in model.edges], dtype=np.int64),
        edge_lengths_um=edge_lengths_um,
        resistances_pa_s_m3=resistances,
        flows_m3_s=flows,
        root_node_id=root,
        leaf_node_ids=leaves,
        root_flow_m3_s=float(root_flow),
        root_mean_velocity_mm_s=float(root_velocity),
        leaf_flow_sum_m3_s=float(leaf_inflows),
        relative_mass_error=float(mass_error),
        maximum_internal_relative_residual=maximum_residual,
        reverse_flow_count=reverse_count,
    )
