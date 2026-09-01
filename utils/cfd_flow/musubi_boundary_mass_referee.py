"""Source-proven Musubi one-step boundary mass referee.

This module is deliberately separate from the production CFD pipeline.  It
replays the pinned D3Q19/PULL boundary writes on existing restart PDFs. The
complete identity is owned by :mod:`full_timestep_mass_referee`; boundary-only
accounting remains diagnostic support and is never a final acceptance referee.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from .exact_link_flux import (
    CS2,
    INVERSE_DIRECTIONS,
    BoundaryReconstruction,
    build_coordinate_lookup,
    equilibrium_pdf,
    pull_fetch_pdfs,
    reconstruct_boundary,
    velocity_from_pdf,
)
from .port_flux_audit import (
    extract_boundary_property_indices,
    parse_bnd_header,
    parse_boundary_property_header,
    read_boundary_ids,
)
from .restart_decode import (
    D3Q19_DIRECTIONS,
    read_treelm_elemlist,
    tree_ids_to_ijk,
)
from .validated_contract import (
    BASE_DX_M as DX_M,
    BOUNDARY_WINDOW_CLOSURE_GATE,
    FULL_TIMESTEP_IDENTITY_GATE,
    INLET_GATE,
    MASS_GATE,
    MAXIMUM_LATTICE_SPEED,
    PRESSURE_GATE,
    RHO0_KG_M3 as RHO0,
    TARGET_MASS_FLOW_KG_S as TARGET_MASS_FLOW,
    VELOCITY_GATE,
    ValidatedTau1Contract,
)


REFEREE_REVISION_OLD = "ENDPOINT_INLET_MINUS_OUTLETS_V1_DIAGNOSTIC_ONLY"
REFEREE_REVISION_NEW = "MUSUBI_ONE_STEP_DISCRETE_MASS_IDENTITY_V2"
RUNTIME_ORDER = ("wall", "outlet_02", "outlet_03", "inlet", "outlet_01")
PORTS = ("inlet", "outlet_01", "outlet_02", "outlet_03")
PRESSURE_NORMALS = {
    "inlet": np.asarray((0, 0, -1), dtype=np.int64),
    "outlet_01": np.asarray((0, -1, -1), dtype=np.int64),
    "outlet_02": np.asarray((1, 0, 1), dtype=np.int64),
    "outlet_03": np.asarray((-1, 0, 0), dtype=np.int64),
}

DT_S = ValidatedTau1Contract().dt_s
OUTLET_PRESSURES_PA = ValidatedTau1Contract().outlet_absolute_pressures_pa


@dataclass(frozen=True, slots=True)
class MeshContract:
    mesh: Path
    tree_ids: np.ndarray
    cell_ijk: np.ndarray
    lookup: dict[tuple[int, int, int], int]
    boundaries: Mapping[str, BoundaryReconstruction]
    boundary_labels: tuple[str, ...]
    qvalues_by_cell: np.ndarray


def kg_s_per_lattice_population(
    density_kg_m3: float = RHO0, dx_m: float = DX_M, dt_s: float = DT_S
) -> float:
    return float(density_kg_m3) * float(dx_m) ** 3 / float(dt_s)


def lattice_delta_to_kg_s(value: float) -> float:
    return float(value) * kg_s_per_lattice_population()


def conservation_identity_residual(
    predicted_delta: float, actual_delta: float, target_lattice_flux: float
) -> float:
    return abs(float(predicted_delta) - float(actual_delta)) / max(
        abs(float(target_lattice_flux)), np.finfo(np.float64).tiny
    )


def boundary_window_closure(
    predicted_mass_change_kg: float,
    observed_mass_change_kg: float,
    inlet_mass_over_window_kg: float,
) -> float:
    return abs(float(predicted_mass_change_kg) - float(observed_mass_change_kg)) / max(
        abs(float(inlet_mass_over_window_kg)), np.finfo(np.float64).tiny
    )


def trapezoidal_integral(iterations: Sequence[int], values: Sequence[float]) -> float:
    x = np.asarray(iterations, dtype=np.float64)
    y = np.asarray(values, dtype=np.float64)
    if x.ndim != 1 or y.ndim != 1 or len(x) != len(y) or len(x) < 2:
        raise ValueError("trapezoidal integration needs matching 1-D samples")
    if np.any(np.diff(x) <= 0):
        raise ValueError("iterations must be strictly increasing")
    return float(np.trapezoid(y, x=x))


def significant_time_averaged_backflow(
    mean_inlet: float, mean_outlets: Iterable[float]
) -> bool:
    limit = 0.05 * abs(float(mean_inlet))
    return any(float(value) < 0.0 and abs(float(value)) > limit for value in mean_outlets)


def fetch_storage_slot(
    *,
    target_cell: int,
    direction: int,
    cell_ijk: np.ndarray,
    coordinate_lookup: Mapping[tuple[int, int, int], int],
) -> tuple[int, int, str]:
    """Return the AOS slot selected by pinned PULL ``FETCH``.

    Musubi connectivity maps a missing source link to the target element's
    inverse PDF.  Otherwise the same PDF direction is read at the source
    element in ``target-c_i``.
    """

    coord = np.asarray(cell_ijk[int(target_cell)], dtype=np.int64)
    source_coord = tuple(
        int(value) for value in coord - D3Q19_DIRECTIONS[int(direction)]
    )
    source = coordinate_lookup.get(source_coord)
    if source is None:
        return int(target_cell), int(INVERSE_DIRECTIONS[int(direction)]), "MISSING_NEIGHBOR_LOCAL_INVERSE"
    return int(source), int(direction), "PULL_SOURCE_SAME_DIRECTION"


def boundary_write_slots(
    boundary: BoundaryReconstruction, selected_rows: np.ndarray | None = None
) -> set[tuple[int, int]]:
    rows = (
        np.arange(len(boundary.cell_indices), dtype=np.int64)
        if selected_rows is None
        else np.asarray(selected_rows, dtype=np.int64)
    )
    result: set[tuple[int, int]] = set()
    for row in rows:
        cell = int(boundary.cell_indices[int(row)])
        # assignBCList stores inverse(raw Seeder side) in globBC bitmask.
        for active_dir in np.flatnonzero(boundary.incoming_masks[int(row)]):
            result.add((cell, int(INVERSE_DIRECTIONS[int(active_dir)])))
    return result


def detect_slot_overlaps(
    slots: Mapping[str, set[tuple[int, int]]]
) -> dict[tuple[str, str], set[tuple[int, int]]]:
    names = list(slots)
    return {
        (names[i], names[j]): slots[names[i]] & slots[names[j]]
        for i in range(len(names))
        for j in range(i + 1, len(names))
    }


def replay_sequential_writes(
    state: np.ndarray,
    operations: Sequence[tuple[str, int, int, float]],
) -> tuple[np.ndarray, dict[str, float]]:
    replay = np.asarray(state, dtype=np.float64).copy()
    deltas: dict[str, float] = {}
    for label, cell, direction, replacement in operations:
        old = float(replay[int(cell), int(direction)])
        replay[int(cell), int(direction)] = float(replacement)
        deltas[label] = deltas.get(label, 0.0) + float(replacement) - old
    return replay, deltas


def _property_definition(header_text: str, label: str) -> tuple[int, int]:
    match = re.search(r"property\s*=\s*\{(.*)\}\s*$", header_text, re.DOTALL)
    if not match:
        raise ValueError("mesh property table is missing")
    for block in re.findall(r"\{(.*?)\}", match.group(1), re.DOTALL):
        found = re.search(r"label\s*=\s*['\"]([^'\"]+)['\"]", block)
        if found and found.group(1).strip() == label:
            bit = re.search(r"bitpos\s*=\s*(\d+)", block)
            count = re.search(r"nElems\s*=\s*(\d+)", block)
            if bit and count:
                return int(bit.group(1)), int(count.group(1))
    raise ValueError(f"mesh property is missing: {label}")


def read_qvalue_rows(
    path: Path,
    *,
    element_count: int,
    side_count: int,
) -> np.ndarray:
    """Read the source-proven TreElm qVal property record layout.

    ``qval.lsb`` is a little-endian sequence of double-precision, element-major
    records.  Within each record the values follow TreElm ``qOffset(1:qQQQ)``.
    """

    values = np.fromfile(Path(path), dtype="<f8")
    expected = int(element_count) * int(side_count)
    if values.size != expected:
        raise ValueError(
            f"qVal property has {values.size} values; expected {expected} "
            f"({element_count} elements x {side_count} sides)"
        )
    return values.reshape(int(element_count), int(side_count))


def scatter_qvalue_rows(
    rows: np.ndarray,
    *,
    cell_count: int,
    property_cells: np.ndarray,
) -> np.ndarray:
    """Scatter qVal rows onto ascending tree elements carrying the property."""

    values = np.asarray(rows, dtype=np.float64)
    cells = np.asarray(property_cells, dtype=np.int64)
    if values.ndim != 2 or values.shape[0] != len(cells):
        raise ValueError("qVal rows and property-cell ordering do not match")
    result = np.full((int(cell_count), values.shape[1]), np.nan, dtype=np.float64)
    result[cells] = values
    return result


def load_mesh_contract(
    mesh: Path,
    *,
    expected_cells: int | None = None,
    allow_zero_normals: bool = False,
    require_runtime_order: bool = True,
) -> MeshContract:
    mesh = Path(mesh).resolve()
    header_text = (mesh / "header.lua").read_text(encoding="utf-8")
    count_match = re.search(r"(?m)^\s*nElems\s*=\s*(\d+)", header_text)
    if not count_match:
        raise ValueError("mesh header is missing nElems")
    cell_count = int(count_match.group(1))
    if expected_cells is not None and cell_count != int(expected_cells):
        raise ValueError(
            f"mesh cell count {cell_count} != expected {int(expected_cells)}"
        )
    boundary_property = parse_boundary_property_header(header_text)
    boundary_header = parse_bnd_header((mesh / "bnd.lua").read_text(encoding="utf-8"))
    if require_runtime_order and tuple(boundary_header.labels) != RUNTIME_ORDER:
        raise ValueError(f"runtime boundary order changed: {boundary_header.labels}")
    tree_ids, property_bits, _ = read_treelm_elemlist(
        mesh / "elemlist.lsb", n_elems=cell_count
    )
    cell_ijk = tree_ids_to_ijk(tree_ids)
    lookup = build_coordinate_lookup(cell_ijk)
    boundary_property_cells = extract_boundary_property_indices(
        property_bits, boundary_property.bit_position
    )
    ids = read_boundary_ids(
        mesh / "bnd.lsb",
        element_count=boundary_property.element_count,
        side_count=boundary_header.side_count,
    )
    label_to_id = {label: index for index, label in enumerate(boundary_header.labels, 1)}
    boundaries = {
        label: reconstruct_boundary(
            ids,
            boundary_property_cells,
            label=label,
            boundary_id=label_to_id[label],
            allow_zero_normals=allow_zero_normals,
        )
        for label in boundary_header.labels
    }

    q_bit, q_count = _property_definition(header_text, "has qVal")
    q_cells = extract_boundary_property_indices(property_bits, q_bit)
    if len(q_cells) != q_count:
        raise ValueError("q-value property count does not match header")
    q_rows = read_qvalue_rows(
        mesh / "qval.lsb",
        element_count=q_count,
        side_count=boundary_header.side_count,
    )
    q_by_cell = scatter_qvalue_rows(
        q_rows,
        cell_count=cell_count,
        property_cells=q_cells,
    )
    return MeshContract(
        mesh=mesh,
        tree_ids=np.asarray(tree_ids),
        cell_ijk=np.asarray(cell_ijk),
        lookup=lookup,
        boundaries=boundaries,
        boundary_labels=tuple(boundary_header.labels),
        qvalues_by_cell=q_by_cell,
    )


def pressure_neighbor_indices(
    cell_indices: np.ndarray,
    cell_ijk: np.ndarray,
    lookup: dict[tuple[int, int, int], int],
    direction: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Apply pressure_eq's source-proven two-neighbor selection rule."""

    normal = np.asarray(direction, dtype=np.int64).reshape(3)
    neighbors1 = np.full(len(cell_indices), -1, dtype=np.int64)
    neighbors2 = np.full(len(cell_indices), -1, dtype=np.int64)
    for row, cell_index in enumerate(np.asarray(cell_indices, dtype=np.int64)):
        coordinate = cell_ijk[cell_index]
        neighbors1[row] = lookup.get(tuple(int(value) for value in coordinate + normal), -1)
        neighbors2[row] = lookup.get(tuple(int(value) for value in coordinate + 2 * normal), -1)
    valid = (neighbors1 >= 0) & (neighbors2 >= 0)
    return valid, neighbors1, neighbors2


def _pressure_selected_rows(
    mesh: MeshContract,
    boundary: BoundaryReconstruction,
    label: str,
    runtime_solid_cells: set[int] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    valid, neighbor1, neighbor2 = pressure_neighbor_indices(
        boundary.cell_indices,
        mesh.cell_ijk,
        mesh.lookup,
        PRESSURE_NORMALS[label],
    )
    if runtime_solid_cells:
        valid &= np.fromiter(
            (int(cell) not in runtime_solid_cells for cell in boundary.cell_indices),
            dtype=bool,
            count=len(boundary.cell_indices),
        )
    return np.flatnonzero(valid).astype(np.int64), neighbor1, neighbor2


def runtime_solid_cells(mesh: MeshContract) -> set[int]:
    """Reproduce TreElm's ``prp_solid`` marking for higher-order BC stencils.

    ``tem_build_treeHorizontalDep`` marks a boundary element solid when either
    of its required pressure-stencil neighbors is absent.  That flag is global:
    Musubi's PULL connectivity subsequently bounce-backs every link for a solid
    current element and every link whose source element is solid.
    """

    result: set[int] = set()
    for label in PORTS:
        boundary = mesh.boundaries[label]
        valid, _, _ = pressure_neighbor_indices(
            boundary.cell_indices,
            mesh.cell_ijk,
            mesh.lookup,
            PRESSURE_NORMALS[label],
        )
        result.update(int(cell) for cell in boundary.cell_indices[~valid])
    return result


def boundary_blocked_pull_mask(mesh: MeshContract) -> np.ndarray:
    """Return PULL directions whose TreElm neighbor entry is a boundary ID.

    A coordinate-adjacent tree element is not sufficient to prove a usable
    TreElm link.  ``assignBCList`` obtains a negative boundary ID from the
    neighbor table and stores its inverse direction in ``globBC%bitmask``.
    Consequently a PULL fetch in a stored bitmask direction is an implicit
    local inverse bounce even when a coordinate lookup happens to find a cell
    across a thin or diagonally touching wall.
    """

    blocked = np.zeros((len(mesh.cell_ijk), 19), dtype=bool)
    for boundary in mesh.boundaries.values():
        blocked[boundary.cell_indices, : boundary.incoming_masks.shape[1]] |= (
            boundary.incoming_masks
        )
    return blocked


def pull_fetch_pdfs_runtime(
    pdf: np.ndarray,
    mesh: MeshContract,
    target_indices: np.ndarray,
    runtime_solid: set[int],
    blocked_pull: np.ndarray | None = None,
) -> np.ndarray:
    """Mirror PULL connectivity after TreElm has set ``prp_solid``.

    The pinned ``mus_construct_connectivity`` maps a link to the current
    element's inverse PDF whenever either the current element or its source is
    solid.  TreElm boundary-neighbor entries are negative rather than ordinary
    coordinate sources and require the same local-inverse PULL mapping.
    """

    values = np.asarray(pdf, dtype=np.float64)
    targets = np.asarray(target_indices, dtype=np.int64).reshape(-1)
    blocked = (
        boundary_blocked_pull_mask(mesh)
        if blocked_pull is None
        else np.asarray(blocked_pull, dtype=bool)
    )
    if blocked.shape != (values.shape[0], values.shape[1]):
        raise ValueError("boundary-blocked PULL mask does not match PDF state")
    fetched = np.empty((len(targets), 19), dtype=np.float64)
    for row, target_value in enumerate(targets):
        target = int(target_value)
        target_is_solid = target in runtime_solid
        coordinate = mesh.cell_ijk[target]
        for direction_index, direction in enumerate(D3Q19_DIRECTIONS):
            source = mesh.lookup.get(
                tuple(int(value) for value in coordinate - direction)
            )
            if (
                target_is_solid
                or blocked[target, direction_index]
                or source is None
                or int(source) in runtime_solid
            ):
                fetched[row, direction_index] = values[
                    target, int(INVERSE_DIRECTIONS[direction_index])
                ]
            else:
                fetched[row, direction_index] = values[int(source), direction_index]
    return fetched


def _pressure_operations(
    *,
    label: str,
    boundary: BoundaryReconstruction,
    selected: np.ndarray,
    equilibrium: np.ndarray,
) -> list[tuple[str, int, int, float]]:
    operations: list[tuple[str, int, int, float]] = []
    for output_row, boundary_row in enumerate(np.asarray(selected, dtype=np.int64)):
        cell = int(boundary.cell_indices[int(boundary_row)])
        for active_dir in np.flatnonzero(boundary.incoming_masks[int(boundary_row)]):
            storage = int(INVERSE_DIRECTIONS[int(active_dir)])
            operations.append((label, cell, storage, float(equilibrium[output_row, active_dir])))
    return operations


def _wall_operations(
    pdf: np.ndarray,
    mesh: MeshContract,
    runtime_solid: set[int] | None = None,
    blocked_pull: np.ndarray | None = None,
) -> list[tuple[str, int, int, float]]:
    boundary = mesh.boundaries["wall"]
    if runtime_solid is None:
        fetched = pull_fetch_pdfs(
            pdf,
            mesh.cell_ijk,
            boundary.cell_indices,
            coordinate_lookup=mesh.lookup,
        )
    else:
        # wall_libb obtains fNgh through computeNeighBuf/FETCH.  The same
        # runtime connectivity used by pressure buffers must therefore be
        # applied here: a solid current element or solid source element maps
        # the read to the current element's inverse PDF.
        fetched = pull_fetch_pdfs_runtime(
            pdf, mesh, boundary.cell_indices, runtime_solid, blocked_pull
        )
    operations: list[tuple[str, int, int, float]] = []
    for row, cell_value in enumerate(boundary.cell_indices):
        cell = int(cell_value)
        for active_value in np.flatnonzero(boundary.incoming_masks[row]):
            active = int(active_value)
            inverse = int(INVERSE_DIRECTIONS[active])
            qvalue = float(mesh.qvalues_by_cell[cell, inverse])
            if not math.isfinite(qvalue) or qvalue <= 0.0:
                raise ValueError(f"invalid wall q-value at cell={cell}, direction={inverse}")
            f_in = float(pdf[cell, active])
            f_out = float(pdf[cell, inverse])
            f_neigh = float(fetched[row, inverse])
            if qvalue >= 0.5:
                replacement = (1.0 - 0.5 / qvalue) * f_in + (0.5 / qvalue) * f_out
            else:
                replacement = 2.0 * qvalue * f_out + (1.0 - 2.0 * qvalue) * f_neigh
            operations.append(("wall", cell, inverse, replacement))
    return operations


def replay_boundary_step(
    pdf: np.ndarray,
    mesh: MeshContract,
    *,
    dx_m: float = DX_M,
    dt_s: float = DT_S,
    density_kg_m3: float = RHO0,
    target_mass_flow_kg_s: float = TARGET_MASS_FLOW,
    outlet_pressures_pa: Mapping[str, float] = OUTLET_PRESSURES_PA,
) -> dict[str, Any]:
    """Replay all writes in geometry/runtime order from one restart state."""

    original = np.asarray(pdf, dtype=np.float64)
    state = original.copy()
    all_operations: list[tuple[str, int, int, float]] = []
    details: dict[str, Any] = {}
    selected_rows: dict[str, np.ndarray] = {}
    conversion = kg_s_per_lattice_population(density_kg_m3, dx_m, dt_s)
    target_flux_lattice = float(target_mass_flow_kg_s) / conversion
    runtime_solid = runtime_solid_cells(mesh)
    blocked_pull = boundary_blocked_pull_mask(mesh)

    for label in mesh.boundary_labels:
        boundary = mesh.boundaries[label]
        if label == "wall":
            operations = _wall_operations(
                original, mesh, runtime_solid, blocked_pull
            )
            selected_rows[label] = np.arange(len(boundary.cell_indices), dtype=np.int64)
        else:
            selected, neighbor1, neighbor2 = _pressure_selected_rows(
                mesh, boundary, label, runtime_solid
            )
            selected_rows[label] = selected
            fetched1 = pull_fetch_pdfs_runtime(
                original,
                mesh,
                neighbor1[selected],
                runtime_solid,
                blocked_pull,
            )
            fetched2 = pull_fetch_pdfs_runtime(
                original,
                mesh,
                neighbor2[selected],
                runtime_solid,
                blocked_pull,
            )
            velocity = 1.5 * velocity_from_pdf(fetched1) - 0.5 * velocity_from_pdf(fetched2)
            if label == "inlet":
                unit_eq = equilibrium_pdf(1.0, velocity)
                alpha = 0.0
                beta = 0.0
                for output_row, boundary_row in enumerate(selected):
                    cell = int(boundary.cell_indices[int(boundary_row)])
                    for active_value in np.flatnonzero(boundary.incoming_masks[int(boundary_row)]):
                        active = int(active_value)
                        storage = int(INVERSE_DIRECTIONS[active])
                        alpha += float(unit_eq[output_row, active])
                        beta -= float(state[cell, storage])
                rho = (target_flux_lattice - beta) / alpha
                eq = equilibrium_pdf(rho, velocity)
                details[label] = {"alpha": alpha, "beta": beta, "rho": rho}
            else:
                pressure_factor = float(density_kg_m3) * float(dx_m) ** 2 / float(dt_s) ** 2
                rho = float(outlet_pressures_pa[label]) / pressure_factor / CS2
                eq = equilibrium_pdf(rho, velocity)
                details[label] = {"rho": rho}
            operations = _pressure_operations(
                label=label,
                boundary=boundary,
                selected=selected,
                equilibrium=eq,
            )
        state, delta = replay_sequential_writes(state, operations)
        details.setdefault(label, {}).update(
            {
                "write_count": len(operations),
                "domain_delta_lattice": float(delta.get(label, 0.0)),
                "domain_delta_kg_s": float(delta.get(label, 0.0)) * conversion,
            }
        )
        all_operations.extend(operations)

    per_label = {
        label: float(details[label]["domain_delta_lattice"])
        for label in mesh.boundary_labels
    }
    predicted = float(sum(per_label.values()))
    return {
        "state": state,
        "operations": all_operations,
        "selected_rows": selected_rows,
        "details": details,
        "per_label_lattice": per_label,
        "per_label_kg_s_domain": {
            label: float(value) * conversion for label, value in per_label.items()
        },
        "predicted_total_lattice": predicted,
        "predicted_total_kg_s": predicted * conversion,
        "target_lattice_flux": target_flux_lattice,
        "runtime_solid_cells": sorted(runtime_solid),
    }


def window_statistics(
    samples: Sequence[Mapping[str, Any]],
    *,
    start_iteration: int,
    end_iteration: int,
    start_total_pdf: float,
    end_total_pdf: float,
) -> dict[str, Any]:
    chosen = [
        sample
        for sample in samples
        if start_iteration <= int(sample["iteration"]) <= end_iteration
    ]
    if not chosen or int(chosen[0]["iteration"]) != start_iteration or int(chosen[-1]["iteration"]) != end_iteration:
        return {"status": "UNDER-SAMPLED"}
    iterations = [int(sample["iteration"]) for sample in chosen]
    total_rate = [float(sample["predicted_total_kg_s"]) for sample in chosen]
    predicted_kg = trapezoidal_integral(iterations, total_rate) * DT_S
    observed_kg = (float(end_total_pdf) - float(start_total_pdf)) * RHO0 * DX_M**3
    inlet_mass = TARGET_MASS_FLOW * (end_iteration - start_iteration) * DT_S
    means: dict[str, float] = {}
    span = end_iteration - start_iteration
    for label in PORTS:
        domain = [float(sample["per_label_kg_s_domain"][label]) for sample in chosen]
        mean_domain = trapezoidal_integral(iterations, domain) / span
        means[label] = mean_domain if label == "inlet" else -mean_domain
    return {
        "status": "PASS",
        "sample_iterations": iterations,
        "predicted_mass_change_kg": predicted_kg,
        "observed_mass_change_kg": observed_kg,
        "closure": boundary_window_closure(predicted_kg, observed_kg, inlet_mass),
        "mean_flow_kg_s": means,
        "significant_time_averaged_backflow": significant_time_averaged_backflow(
            means["inlet"], (means["outlet_01"], means["outlet_02"], means["outlet_03"])
        ),
    }


def earliest_gate_pass(
    records: Sequence[Mapping[str, Any]],
    *,
    expected_inlet_globbc: int | None = None,
) -> Mapping[str, Any] | None:
    for record in sorted(records, key=lambda item: int(item["iteration"])):
        finite_keys = (
                "R_mass_short",
                "R_mass_long",
                "R_velocity",
                "R_pressure",
                "R_inlet",
                "R_conservation_identity",
                "boundary_window_closure",
                "minimum_pdf",
                "maximum_lattice_speed",
            )
        try:
            finite = all(math.isfinite(float(record[key])) for key in finite_keys)
        except (KeyError, TypeError, ValueError):
            finite = False
        if finite and all(
            (
                float(record["R_mass_short"]) <= MASS_GATE,
                float(record["R_mass_long"]) <= MASS_GATE,
                float(record["R_velocity"]) <= VELOCITY_GATE,
                float(record["R_pressure"]) <= PRESSURE_GATE,
                float(record["R_inlet"]) <= INLET_GATE,
                float(record["R_conservation_identity"]) <= FULL_TIMESTEP_IDENTITY_GATE,
                float(record["boundary_window_closure"]) <= BOUNDARY_WINDOW_CLOSURE_GATE,
                not bool(record["significant_time_averaged_backflow"]),
                float(record["minimum_pdf"]) > 0.0,
                float(record["maximum_lattice_speed"]) < MAXIMUM_LATTICE_SPEED,
                expected_inlet_globbc is None
                or int(record["inlet_globbc"]) == expected_inlet_globbc,
            )
        ):
            return record
    return None
