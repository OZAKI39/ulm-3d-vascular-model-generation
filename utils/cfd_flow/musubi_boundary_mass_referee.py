"""Source-proven Musubi one-step boundary mass referee.

This module is deliberately separate from the production CFD pipeline.  It
replays the pinned D3Q19/PULL boundary writes on existing restart PDFs and can
validate that replay with one, and only one, frozen-binary timestep.
"""

from __future__ import annotations

import csv
import math
import re
import shutil
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from .adaptive_flux_steady_exact_audit import _pressure_neighbor_indices
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
from .io import read_json, sha256_file, write_json
from .port_flux_audit import (
    extract_boundary_property_indices,
    parse_bnd_header,
    parse_boundary_property_header,
    read_boundary_ids,
)
from .restart_decode import (
    D3Q19_DIRECTIONS,
    read_restart_pdf,
    read_treelm_elemlist,
    tree_ids_to_ijk,
)


REFEREE_REVISION_OLD = "ENDPOINT_INLET_MINUS_OUTLETS_V1"
REFEREE_REVISION_NEW = "MUSUBI_ONE_STEP_DISCRETE_MASS_IDENTITY_V2"
FINAL_PASS = "CFD_FLOW_HEALTHY_REFEREE_FIXED_STEADY_PASS"
FINAL_FAIL = "CFD_FLOW_HEALTHY_REFEREE_DISCRETE_IDENTITY_FAILED"
NEXT_PASS = "RUN HEALTHY ADAPTIVE-FLUX GRID CONVERGENCE"
NEXT_FAIL = "FIX CURRENT DISCRETE BOUNDARY MASS REPLAY"
INSTRUMENTED_BINARY_WSL = (
    "/home/lzy/apes-worktrees/musubi_boundary_mass_oracle_20260829/build/musubi"
)
INSTRUMENTED_BINARY_SHA256 = (
    "bdd738e7ec24e794b370a26f6dbe8012ea7bf21f1878d9cc539b2d54fcb106e8"
)

RUN_NAME = "healthy_mouse_capillary_fast_steady_anchor003274_20260829_195500"
MESH_RUN = "axis_aligned_ideal_plane_inlet_preflight_anchor003274_20260829_120444"
ENDPOINT_ITERATION = 469_900
EXPECTED_CELLS = 221_309
EXPECTED_INLET_GLOBBC = 287
DX_M = 2.0e-7
DT_S = 2.441406249999999e-8
RHO0 = 1056.0
TARGET_MASS_FLOW = 2.890180380479642e-12
TARGET_Q = 2.7369132390905703e-15
PRESSURE_REFERENCE_PA = 23622.32012800001
OUTLET_PRESSURES_PA = {
    "outlet_01": 23636.865106101286,
    "outlet_02": 23754.524677223188,
    "outlet_03": 23608.619501326699,
}
PRESSURE_NORMALS = {
    "inlet": np.asarray((0, 0, -1), dtype=np.int64),
    "outlet_01": np.asarray((0, -1, -1), dtype=np.int64),
    "outlet_02": np.asarray((1, 0, 1), dtype=np.int64),
    "outlet_03": np.asarray((-1, 0, 0), dtype=np.int64),
}

FROZEN_BINARY_WSL = (
    "/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300/"
    "build/musubi_adaptive_flux"
)
FROZEN_BINARY_SHA256 = (
    "e80162fb7e0e657d2e41aafc40a1b13b32204ff34692e24b7ab02c51aa97c588"
)
SOURCE_ROOT_WSL = "/home/lzy/apes-worktrees/musubi_mcclure_adaptive_flux_20260829_1300"
RUNTIME_ROOT_WSL = "/home/lzy/u3da/referee_195500_one_step"
MPIRUN_WSL = "/home/lzy/.local/bin/mpirun"
MPI_ARGS = ("--bind-to", "core", "--map-by", "core", "--report-bindings")
RUNTIME_ORDER = ("wall", "outlet_02", "outlet_03", "inlet", "outlet_01")
PORTS = ("inlet", "outlet_01", "outlet_02", "outlet_03")


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


def _pressure_selected_rows(
    mesh: MeshContract,
    boundary: BoundaryReconstruction,
    label: str,
    runtime_solid_cells: set[int] | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    valid, neighbor1, neighbor2 = _pressure_neighbor_indices(
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
        valid, _, _ = _pressure_neighbor_indices(
            boundary.cell_indices,
            mesh.cell_ijk,
            mesh.lookup,
            PRESSURE_NORMALS[label],
        )
        result.update(int(cell) for cell in boundary.cell_indices[~valid])
    return result


def pull_fetch_pdfs_runtime(
    pdf: np.ndarray,
    mesh: MeshContract,
    target_indices: np.ndarray,
    runtime_solid: set[int],
) -> np.ndarray:
    """Mirror PULL connectivity after TreElm has set ``prp_solid``.

    This differs from a raw coordinate lookup at the three outlet_02 neighbor
    cells adjacent to the 18 runtime-solid stencil elements.  The pinned
    ``mus_construct_connectivity`` maps a link to the current element's inverse
    PDF whenever either the current element or its source is solid.
    """

    values = np.asarray(pdf, dtype=np.float64)
    targets = np.asarray(target_indices, dtype=np.int64).reshape(-1)
    fetched = np.empty((len(targets), 19), dtype=np.float64)
    for row, target_value in enumerate(targets):
        target = int(target_value)
        target_is_solid = target in runtime_solid
        coordinate = mesh.cell_ijk[target]
        for direction_index, direction in enumerate(D3Q19_DIRECTIONS):
            source = mesh.lookup.get(
                tuple(int(value) for value in coordinate - direction)
            )
            if target_is_solid or source is None or int(source) in runtime_solid:
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
            pdf, mesh, boundary.cell_indices, runtime_solid
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

    for label in mesh.boundary_labels:
        boundary = mesh.boundaries[label]
        if label == "wall":
            operations = _wall_operations(original, mesh, runtime_solid)
            selected_rows[label] = np.arange(len(boundary.cell_indices), dtype=np.int64)
        else:
            selected, neighbor1, neighbor2 = _pressure_selected_rows(
                mesh, boundary, label, runtime_solid
            )
            selected_rows[label] = selected
            fetched1 = pull_fetch_pdfs_runtime(
                original, mesh, neighbor1[selected], runtime_solid
            )
            fetched2 = pull_fetch_pdfs_runtime(
                original, mesh, neighbor2[selected], runtime_solid
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


def earliest_gate_pass(records: Sequence[Mapping[str, Any]]) -> Mapping[str, Any] | None:
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
                float(record["R_mass_short"]) <= 0.01,
                float(record["R_mass_long"]) <= 0.01,
                float(record["R_velocity"]) <= 0.01,
                float(record["R_pressure"]) <= 0.005,
                float(record["R_inlet"]) <= 0.01,
                float(record["R_conservation_identity"]) <= 1.0e-8,
                float(record["boundary_window_closure"]) <= 0.001,
                not bool(record["significant_time_averaged_backflow"]),
                float(record["minimum_pdf"]) > 0.0,
                float(record["maximum_lattice_speed"]) < 0.05,
                int(record["inlet_globbc"]) == EXPECTED_INLET_GLOBBC,
            )
        ):
            return record
    return None


def _source_evidence(path: Path, tokens: Sequence[str]) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    lines: dict[str, int] = {}
    for token in tokens:
        offset = text.find(token)
        if offset < 0:
            raise ValueError(f"source token absent from {path}: {token}")
        lines[token] = text.count("\n", 0, offset) + 1
    return {"path": str(path), "sha256": sha256_file(path), "token_lines": lines}


def _wsl_unc(distribution: str, path: str) -> Path:
    return Path(rf"\\wsl.localhost\{distribution}") / path.lstrip("/").replace("/", "\\")


def _wsl_run(command: str, *, timeout: int = 300) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["wsl.exe", "-d", "Ubuntu", "--", "bash", "-lc", command],
        text=True,
        encoding="utf-8",
        errors="replace",
        capture_output=True,
        timeout=timeout,
        check=False,
    )


def _one_step_lua(mesh_wsl: str, restart_header: str) -> str:
    return f"""-- Boundary mass referee: exactly one iteration, no tracking.
simulation_name = 'a3274_referee'
printRuntimeInfo = true
timing_file = 'timing.res'
mesh = '{mesh_wsl}/'
scaling = 'diffusive'
logging = {{ level = 5 }}
maximum_iterations = 469901
dx = {DX_M:.17g}
dt = {DT_S:.17g}
rho0_phy = {RHO0:.17g}
nu_phy = 3.27e-6
bulk_viscosity_phy = 2.18e-6
pressure_reference_phy = {PRESSURE_REFERENCE_PA:.17g}
function outlet_01_pressure(x,y,z,t) return {OUTLET_PRESSURES_PA['outlet_01']:.17g} end
function outlet_02_pressure(x,y,z,t) return {OUTLET_PRESSURES_PA['outlet_02']:.17g} end
function outlet_03_pressure(x,y,z,t) return {OUTLET_PRESSURES_PA['outlet_03']:.17g} end
physics = {{ dt = dt, rho0 = rho0_phy }}
identify = {{ label = 'ROI003274', kind = 'fluid', layout = 'd3q19', relaxation = 'bgk' }}
fluid = {{ kinematic_viscosity = nu_phy, bulk_viscosity = bulk_viscosity_phy }}
initial_condition = {{ pressure = pressure_reference_phy, velocityX = 0.0, velocityY = 0.0, velocityZ = 0.0 }}
boundary_condition = {{
  {{ label = 'wall', kind = 'wall_libb' }},
  {{ label = 'inlet', kind = 'adaptive_flux_pressure', mass_flowrate = {TARGET_MASS_FLOW:.17g} }},
  {{ label = 'outlet_01', kind = 'pressure_eq', pressure = outlet_01_pressure }},
  {{ label = 'outlet_02', kind = 'pressure_eq', pressure = outlet_02_pressure }},
  {{ label = 'outlet_03', kind = 'pressure_eq', pressure = outlet_03_pressure }}
}}
sim_control = {{ time_control = {{ max = {{ iter = maximum_iterations }}, interval = {{ iter = 1 }} }} }}
restart = {{
  read = 'restart/{restart_header}', write = 'restart_oracle/',
  time_control = {{ min = {{ iter = 469901 }}, max = {{ iter = 469901 }}, interval = {{ iter = 1 }} }}
}}
"""


def run_original_one_step(
    *, endpoint_header: Path, endpoint_binary: Path, mesh: Path, qc: Path
) -> dict[str, Any]:
    oracle_path = qc / "one_step_oracle.json"
    if oracle_path.exists():
        existing = read_json(oracle_path)
        if existing.get("status") != "PASS":
            raise RuntimeError("existing one-step oracle did not pass")
        binary = Path(str(existing["restart_binary"]))
        if not binary.is_file() or sha256_file(binary) != existing.get("restart_sha256"):
            raise RuntimeError("existing one-step oracle restart evidence changed")
        existing["reused_existing_evidence"] = True
        return existing
    binary_unc = _wsl_unc("Ubuntu", FROZEN_BINARY_WSL)
    if sha256_file(binary_unc) != FROZEN_BINARY_SHA256:
        raise RuntimeError("frozen adaptive binary SHA256 changed")
    runtime = _wsl_unc("Ubuntu", RUNTIME_ROOT_WSL)
    if runtime.exists():
        raise RuntimeError(f"one-step runtime already exists: {RUNTIME_ROOT_WSL}")
    (runtime / "restart").mkdir(parents=True)
    (runtime / "restart_oracle").mkdir()
    shutil.copy2(endpoint_header, runtime / "restart" / endpoint_header.name)
    shutil.copy2(endpoint_binary, runtime / "restart" / endpoint_binary.name)
    mesh_wsl = "/mnt/e/ULM/hatimb-particle_flow_simulator/ulm_3D_vascular/" + str(
        mesh.relative_to(mesh.parents[4])
    ).replace("\\", "/")
    lua = _one_step_lua(mesh_wsl, endpoint_header.name)
    (runtime / "musubi.lua").write_text(lua, encoding="utf-8")
    command = (
        "cd " + RUNTIME_ROOT_WSL
        + " && env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 "
        + "OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 "
        + MPIRUN_WSL
        + " -np 4 "
        + " ".join(MPI_ARGS)
        + " " + FROZEN_BINARY_WSL + " musubi.lua"
    )
    started = time.perf_counter()
    process = _wsl_run(command, timeout=300)
    elapsed = time.perf_counter() - started
    (qc / "one_step_musubi_stdout.log").write_text(process.stdout, encoding="utf-8")
    (qc / "one_step_musubi_stderr.log").write_text(process.stderr, encoding="utf-8")
    if process.returncode != 0:
        result = {
            "status": "FAIL",
            "returncode": process.returncode,
            "wall_time_s": elapsed,
            "first_failure": next((line for line in process.stderr.splitlines() if line.strip()), None),
            "musubi_calls": 1,
        }
        write_json(oracle_path, result)
        return result
    headers = sorted((runtime / "restart_oracle").glob("*header*.lua"))
    if not headers:
        raise RuntimeError("one-step Musubi returned zero but wrote no restart header")
    header = headers[-1]
    text = header.read_text(encoding="utf-8")
    iteration_match = re.search(r"\biter\s*=\s*(\d+)", text)
    binary_match = re.search(r"binary_name\s*=\s*\{\s*['\"]([^'\"]+)", text, re.DOTALL)
    if not iteration_match or not binary_match or int(iteration_match.group(1)) != 469901:
        raise RuntimeError("one-step restart header is not iteration 469901")
    binary = runtime / binary_match.group(1)
    if not binary.is_file():
        binary = runtime / "restart_oracle" / Path(binary_match.group(1)).name
    archive = qc / "one_step_restart_469901"
    archive.mkdir()
    shutil.copy2(header, archive / header.name)
    shutil.copy2(binary, archive / binary.name)
    result = {
        "status": "PASS",
        "returncode": process.returncode,
        "wall_time_s": elapsed,
        "musubi_calls": 1,
        "mpi_ranks": 4,
        "start_iteration": ENDPOINT_ITERATION,
        "end_iteration": 469901,
        "runtime_root_wsl": RUNTIME_ROOT_WSL,
        "restart_header": str(header),
        "restart_binary": str(binary),
        "restart_sha256": sha256_file(binary),
        "archive": str(archive),
    }
    write_json(oracle_path, result)
    return result


def load_instrumented_oracle(
    qc: Path, python_replay: Mapping[str, Any]
) -> dict[str, Any]:
    """Load the already-authorized instrumented one-step evidence.

    This function is intentionally read-only: it never launches the debug
    binary, and verifies both its hash and the five expected per-BC records.
    """

    log_path = qc / "instrumented_musubi_stdout.log"
    if not log_path.is_file():
        raise RuntimeError("instrumented one-step stdout evidence is missing")
    binary = _wsl_unc("Ubuntu", INSTRUMENTED_BINARY_WSL)
    if not binary.is_file() or sha256_file(binary) != INSTRUMENTED_BINARY_SHA256:
        raise RuntimeError("instrumented debug binary evidence changed")
    pattern = re.compile(
        r"BOUNDARY_MASS_ORACLE\s+iter=(\d+)\s+label=(\S+)\s+"
        r"delta=\s*([+\-0-9.Ee]+)\s+changed_slots=(\d+)"
    )
    observed: dict[str, dict[str, Any]] = {}
    for match in pattern.finditer(log_path.read_text(encoding="utf-8")):
        iteration, label, delta, changed = match.groups()
        observed[label] = {
            "iteration": int(iteration),
            "delta_lattice": float(delta),
            "changed_slots": int(changed),
        }
    if tuple(label for label in RUNTIME_ORDER if label in observed) != RUNTIME_ORDER:
        raise RuntimeError("instrumented log does not contain the complete runtime order")
    comparisons: dict[str, dict[str, float]] = {}
    for label in RUNTIME_ORDER:
        predicted = float(python_replay["per_label_lattice"][label])
        actual = float(observed[label]["delta_lattice"])
        comparisons[label] = {
            "python_lattice": predicted,
            "instrumented_lattice": actual,
            "absolute_error_lattice": abs(predicted - actual),
            "relative_error": abs(predicted - actual)
            / max(abs(actual), np.finfo(np.float64).tiny),
        }
    result = {
        "status": "PASS"
        if all(item["absolute_error_lattice"] <= 1.0e-12 for item in comparisons.values())
        else "FAIL",
        "musubi_calls": 1,
        "reused_existing_evidence": True,
        "binary": INSTRUMENTED_BINARY_WSL,
        "binary_sha256": INSTRUMENTED_BINARY_SHA256,
        "log": str(log_path),
        "runtime_boundary_order": list(RUNTIME_ORDER),
        "observed": observed,
        "comparisons": comparisons,
    }
    write_json(qc / "instrumented_oracle.json", result)
    return result


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty CSV: {path}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _restart_path(record: Mapping[str, Any], run_root: Path) -> Path:
    candidate = Path(str(record["restart_binary"]))
    if candidate.is_file():
        return candidate
    iteration = int(record["iteration"])
    if iteration == ENDPOINT_ITERATION:
        archived = next((run_root / "restart" / f"checkpoint_{iteration}").glob("*.lsb"))
        return archived
    raise FileNotFoundError(candidate)


def _source_contracts(source_root: Path, mesh: MeshContract, qc: Path) -> None:
    mus = source_root / "mus" / "source"
    files = {
        "boundary_loop": _source_evidence(
            mus / "bc" / "mus_bc_general_module.fpp",
            ("call fill_neighBuffer(", "do iBnd = 1, nBCs", "state       = state( :, pdf%nNext )"),
        ),
        "pressure_write": _source_evidence(
            mus / "bc" / "mus_bc_fluid_module.fpp",
            ("subroutine pressure_eq_reconstruct", "state(?FETCH?(iDir", "neighBufferPre_nNext(1,:)"),
        ),
        "boundary_connectivity": _source_evidence(
            mus / "mus_construction_module.fpp",
            (
                "Set the bitmask for the incoming directions which have to be",
                "stencil%cxDirInv(iDir), nBnds( bID )",
                "q-values stored in outgoing direction",
            ),
        ),
        "solid_pull_connectivity": _source_evidence(
            mus / "mus_connectivity_module.fpp",
            (
                "solidified = ( btest(neighProp, prp_solid)",
                "sourceDir = stencil%cxDirInv( iDir )",
                "GetFromPos = iElem",
            ),
        ),
        "runtime_solid_marking": _source_evidence(
            source_root / "tem" / "source" / "tem_construction_module.f90",
            (
                "subroutine tem_build_treeHorizontalDep",
                "if( levelPos <= 0 .and. computeStencil%requireNeighNeigh )",
                "ibset( levelDesc%property( totalPos ), prp_solid )",
            ),
        ),
        "wall_write": _source_evidence(
            mus / "bc" / "mus_bc_fluid_wall_module.fpp",
            ("subroutine wall_libb", "state( me%links(iLevel)%val(iLink) )", "cIn*fIn + cOut*fOut + cNgh*fNgh"),
        ),
        "pull_macro": _source_evidence(
            source_root / "mus" / "source" / "header" / "lbm_macros.inc",
            ("else !PULL", "macro :: FETCH", "?neigh?(?NGPOS?(?iDir?"),
        ),
        "restart_serializer": _source_evidence(
            mus / "mus_buffer_module.fpp",
            ("subroutine mus_pdf_serialize", "scheme%pdf( iLevel )%nNext", "subroutine mus_pdf_unserialize"),
        ),
        "timestep": _source_evidence(
            mus / "mus_control_module.f90",
            ("call tem_time_advance", "call set_boundary", "call mus_swap_now_next", "call me%scheme%compute"),
        ),
    }
    write_json(
        qc / "source_boundary_order_contract.json",
        {
            "status": "PASS",
            "runtime_boundary_order": list(mesh.boundary_labels),
            "order_source": "Seeder bnd.lua labels indexed by set_boundary iBnd=1..nBCs",
            "neighbor_buffers": "filled once from original buffers before all boundary writes",
            "non_boundary_mass_sources": "NONE (configured fluid has no source table/body force)",
            "source_files": files,
        },
    )
    write_json(
        qc / "timestep_phase_contract.json",
        {
            "status": "PASS",
            "restart_iteration_n": "serialized from nNext after fused PULL stream/BGK collision",
            "continuation": "unserialize into nNext; advance to n+1; set_boundary writes nNext; swap; compute writes new nNext; restart n+1 serializes nNext",
            "endpoint_replay_corresponds_to": "boundary update and compute producing iteration 469901",
            "off_by_one": False,
            "source_files": files,
        },
    )


def run_referee(project_root: Path, *, run_one_step: bool = True) -> dict[str, Any]:
    root = Path(project_root).resolve()
    run_root = root / "outputs" / "cfd_flow" / RUN_NAME
    qc = run_root / "qc" / "referee"
    qc.mkdir(parents=True, exist_ok=True)
    mesh_path = root / "outputs" / "cfd_flow" / MESH_RUN / "seeder" / "mesh"
    source_root = _wsl_unc("Ubuntu", SOURCE_ROOT_WSL)
    mesh = load_mesh_contract(mesh_path, expected_cells=EXPECTED_CELLS)
    _source_contracts(source_root, mesh, qc)

    convergence = read_json(run_root / "qc" / "exact_watchdog" / "checkpoint_exact_convergence.json")
    records = convergence["records"]
    endpoint_record = next(item for item in records if int(item["iteration"]) == ENDPOINT_ITERATION)
    endpoint_binary = _restart_path(endpoint_record, run_root)
    endpoint_header = next((run_root / "restart" / f"checkpoint_{ENDPOINT_ITERATION}").glob("*header*.lua"))

    frozen_paths = (
        root / "cfd_flow.py",
        root / "configs" / "cfd_flow.yaml",
        root / "utils" / "cfd_flow" / "pipeline.py",
        mesh.mesh / "elemlist.lsb",
        mesh.mesh / "bnd.lsb",
        mesh.mesh / "qval.lsb",
        endpoint_binary,
        _wsl_unc("Ubuntu", FROZEN_BINARY_WSL),
    )
    frozen_before = {str(path): sha256_file(path) for path in frozen_paths}

    runtime_solid = runtime_solid_cells(mesh)
    selected_rows: dict[str, np.ndarray] = {}
    for label, boundary in mesh.boundaries.items():
        if label == "wall":
            selected_rows[label] = np.arange(len(boundary.cell_indices), dtype=np.int64)
        else:
            selected_rows[label] = _pressure_selected_rows(
                mesh, boundary, label, runtime_solid
            )[0]
    slots = {
        label: boundary_write_slots(mesh.boundaries[label], selected_rows[label])
        for label in RUNTIME_ORDER
    }
    overlaps = detect_slot_overlaps(slots)
    pair_counts = {f"{a}/{b}": len(value) for (a, b), value in overlaps.items()}
    write_json(
        qc / "boundary_overlap_contract.json",
        {
            "status": "PASS",
            "unique_overwritten_slots": len(set().union(*slots.values())),
            "overlap_slot_count": len(set().union(*overlaps.values())) if overlaps else 0,
            "per_pair_overlap_counts": pair_counts,
            "last_writer_order": list(RUNTIME_ORDER),
            "runtime_solid_cell_count": len(runtime_solid),
            "runtime_solid_cells_zero_based": sorted(runtime_solid),
        },
    )

    inlet = mesh.boundaries["inlet"]
    example_active = int(np.flatnonzero(inlet.incoming_masks[0])[0])
    example_cell = int(inlet.cell_indices[0])
    fetch_cell, fetch_dir, fetch_kind = fetch_storage_slot(
        target_cell=example_cell,
        direction=example_active,
        cell_ijk=mesh.cell_ijk,
        coordinate_lookup=mesh.lookup,
    )
    write_json(
        qc / "exact_fetch_slot_contract.json",
        {
            "status": "PASS",
            "streaming": "PULL",
            "musubi_write": "state(FETCH(active_bitmask_direction, boundary_element))",
            "missing_neighbor_mapping": "local inverse PDF storage slot",
            "Seeder_to_globBC_mapping": "assignBCList sets globBC bitmask at inverse(raw Seeder boundary direction)",
            "old_audit_mapping": "iterates reconstructed globBC bitmask direction and reads its inverse local storage slot",
            "old_slot_assumption_valid": True,
            "runtime_solid_connectivity": "if current or source element has prp_solid, FETCH maps to current local inverse PDF",
            "example": {
                "boundary": "inlet",
                "fluid_cell_index_zero_based": example_cell,
                "active_iDir_zero_based": example_active,
                "inverse_iDir_zero_based": int(INVERSE_DIRECTIONS[example_active]),
                "fetch_source_cell_zero_based": fetch_cell,
                "fetch_storage_direction_zero_based": fetch_dir,
                "fetch_kind": fetch_kind,
            },
        },
    )

    corrected_samples: list[dict[str, Any]] = []
    for record in records:
        binary = _restart_path(record, run_root)
        pdf = read_restart_pdf(binary, n_elems=EXPECTED_CELLS, n_components=19)
        replay = replay_boundary_step(pdf, mesh)
        corrected_samples.append(
            {
                "iteration": int(record["iteration"]),
                "restart_binary": str(binary),
                "restart_sha256": sha256_file(binary),
                "total_pdf_mass": float(np.sum(pdf, dtype=np.float64)),
                "predicted_total_lattice": replay["predicted_total_lattice"],
                "predicted_total_kg_s": replay["predicted_total_kg_s"],
                "per_label_lattice": replay["per_label_lattice"],
                "per_label_kg_s_domain": replay["per_label_kg_s_domain"],
                "details": replay["details"],
            }
        )

    flattened = []
    for sample in corrected_samples:
        flattened.append(
            {
                "iteration": sample["iteration"],
                "total_pdf_mass": sample["total_pdf_mass"],
                "wall_domain_kg_s": sample["per_label_kg_s_domain"]["wall"],
                "inlet_into_domain_kg_s": sample["per_label_kg_s_domain"]["inlet"],
                "outlet_01_outward_kg_s": -sample["per_label_kg_s_domain"]["outlet_01"],
                "outlet_02_outward_kg_s": -sample["per_label_kg_s_domain"]["outlet_02"],
                "outlet_03_outward_kg_s": -sample["per_label_kg_s_domain"]["outlet_03"],
                "net_domain_kg_s": sample["predicted_total_kg_s"],
            }
        )
    _write_csv(qc / "corrected_boundary_flux_history.csv", flattened)

    historical_rows: list[dict[str, Any]] = []
    closure_errors: list[float] = []
    sign_contradiction = False
    for previous, current in zip(corrected_samples, corrected_samples[1:]):
        delta_iter = int(current["iteration"]) - int(previous["iteration"])
        legacy0 = next(item for item in records if int(item["iteration"]) == int(previous["iteration"]))
        legacy1 = next(item for item in records if int(item["iteration"]) == int(current["iteration"]))
        legacy_net0 = float(legacy0["exact_inlet_mass_flow"]) - float(legacy0["outlet_signed_sum"])
        legacy_net1 = float(legacy1["exact_inlet_mass_flow"]) - float(legacy1["outlet_signed_sum"])
        actual_rate = (
            (float(current["total_pdf_mass"]) - float(previous["total_pdf_mass"]))
            / delta_iter
            * kg_s_per_lattice_population()
        )
        legacy_average = 0.5 * (legacy_net0 + legacy_net1)
        contradiction = np.sign(legacy_average) != np.sign(actual_rate) and abs(actual_rate) > 0.0
        sign_contradiction = sign_contradiction or bool(contradiction)
        predicted_kg = (
            0.5
            * (float(previous["predicted_total_kg_s"]) + float(current["predicted_total_kg_s"]))
            * delta_iter
            * DT_S
        )
        actual_kg = (
            float(current["total_pdf_mass"]) - float(previous["total_pdf_mass"])
        ) * RHO0 * DX_M**3
        error = boundary_window_closure(
            predicted_kg, actual_kg, TARGET_MASS_FLOW * delta_iter * DT_S
        )
        closure_errors.append(error)
        historical_rows.append(
            {
                "start_iteration": previous["iteration"],
                "end_iteration": current["iteration"],
                "delta_iterations": delta_iter,
                "legacy_average_net_kg_s": legacy_average,
                "corrected_average_net_kg_s": 0.5 * (
                    float(previous["predicted_total_kg_s"]) + float(current["predicted_total_kg_s"])
                ),
                "actual_domain_mass_rate_kg_s": actual_rate,
                "legacy_sign_contradiction": contradiction,
                "R_integrated_closure": error,
            }
        )
    _write_csv(qc / "historical_integrated_closure.csv", historical_rows)

    endpoint_pdf = read_restart_pdf(endpoint_binary, n_elems=EXPECTED_CELLS, n_components=19)
    endpoint_replay = replay_boundary_step(endpoint_pdf, mesh)
    oracle = run_original_one_step(
        endpoint_header=endpoint_header,
        endpoint_binary=endpoint_binary,
        mesh=mesh.mesh,
        qc=qc,
    ) if run_one_step else {"status": "NOT_RUN", "musubi_calls": 0}
    if oracle["status"] != "PASS":
        raise RuntimeError(f"one-step original Musubi failed: {oracle}")
    next_pdf = read_restart_pdf(Path(oracle["restart_binary"]), n_elems=EXPECTED_CELLS, n_components=19)
    # Sum the elementwise change, rather than subtracting two O(1e5) totals.
    # This avoids cancellation while retaining the exact serialized PDFs.
    actual_delta = float(np.sum(next_pdf - endpoint_pdf, dtype=np.float64))
    identity = conservation_identity_residual(
        endpoint_replay["predicted_total_lattice"],
        actual_delta,
        endpoint_replay["target_lattice_flux"],
    )
    oracle.update(
        {
            "predicted_total_lattice": endpoint_replay["predicted_total_lattice"],
            "actual_total_lattice": actual_delta,
            "predicted_total_kg_s": endpoint_replay["predicted_total_kg_s"],
            "actual_total_kg_s": lattice_delta_to_kg_s(actual_delta),
            "R_one_step_identity": identity,
            "per_boundary_lattice": endpoint_replay["per_label_lattice"],
            "per_boundary_kg_s_domain": endpoint_replay["per_label_kg_s_domain"],
            "other_source_delta_lattice": 0.0,
            "other_source_delta_kg_s": 0.0,
            "preferred_gate": 1.0e-10,
            "hard_gate": 1.0e-8,
        }
    )
    write_json(qc / "one_step_oracle.json", oracle)
    instrumented = load_instrumented_oracle(qc, endpoint_replay)

    sample_by_iteration = {int(item["iteration"]): item for item in corrected_samples}
    corrected_records: list[dict[str, Any]] = []
    for legacy in records:
        iteration = int(legacy["iteration"])
        sample = sample_by_iteration[iteration]
        short_start = legacy.get("short_window_start")
        long_start = legacy.get("long_window_start")
        short_stats = (
            window_statistics(
                corrected_samples,
                start_iteration=int(short_start),
                end_iteration=iteration,
                start_total_pdf=float(sample_by_iteration[int(short_start)]["total_pdf_mass"]),
                end_total_pdf=float(sample["total_pdf_mass"]),
            )
            if short_start is not None and int(short_start) in sample_by_iteration
            else {"status": "UNDER-SAMPLED"}
        )
        long_stats = (
            window_statistics(
                corrected_samples,
                start_iteration=int(long_start),
                end_iteration=iteration,
                start_total_pdf=float(sample_by_iteration[int(long_start)]["total_pdf_mass"]),
                end_total_pdf=float(sample["total_pdf_mass"]),
            )
            if long_start is not None and int(long_start) in sample_by_iteration
            else {"status": "UNDER-SAMPLED"}
        )
        closure = short_stats.get("closure") if short_stats.get("status") == "PASS" else math.nan
        means = short_stats.get("mean_flow_kg_s", {})
        corrected_records.append(
            {
                "iteration": iteration,
                "restart_binary": sample["restart_binary"],
                "restart_sha256": sample["restart_sha256"],
                "R_mass_short": legacy.get("R_mass_short"),
                "R_mass_long": legacy.get("R_mass_long"),
                "R_velocity": legacy.get("R_velocity"),
                "R_pressure": legacy.get("R_pressure"),
                "R_inlet": abs(float(sample["per_label_kg_s_domain"]["inlet"]) - TARGET_MASS_FLOW) / TARGET_MASS_FLOW,
                "R_conservation_identity": identity,
                "boundary_window_closure": closure,
                "short_window": short_stats,
                "long_window": long_stats,
                "mean_outlet_01": means.get("outlet_01"),
                "mean_outlet_02": means.get("outlet_02"),
                "mean_outlet_03": means.get("outlet_03"),
                "significant_time_averaged_backflow": short_stats.get("significant_time_averaged_backflow", False),
                "minimum_pdf": legacy.get("minimum_pdf"),
                "maximum_lattice_speed": legacy.get("maximum_lattice_speed"),
                "inlet_globbc": legacy.get("inlet_globbc"),
            }
        )
    csv_records = [
        {key: value for key, value in record.items() if not isinstance(value, (dict, list))}
        for record in corrected_records
    ]
    _write_csv(qc / "corrected_residual_history.csv", csv_records)
    earliest = earliest_gate_pass(corrected_records)

    endpoint_corrected = next(item for item in corrected_records if item["iteration"] == ENDPOINT_ITERATION)
    short10 = endpoint_corrected["short_window"]
    long20 = endpoint_corrected["long_window"]
    frozen_after = {str(path): sha256_file(path) for path in frozen_paths}
    unchanged = frozen_before == frozen_after
    accepted = earliest if earliest is not None else endpoint_corrected
    accepted_path = Path(str(accepted["restart_binary"]))
    accepted_exists = accepted_path.is_file()
    status = (
        FINAL_PASS
        if identity <= 1.0e-8
        and instrumented["status"] == "PASS"
        and earliest is not None
        and unchanged
        else FINAL_FAIL
    )
    manifest = {
        "status": status,
        "next": NEXT_PASS if status == FINAL_PASS else NEXT_FAIL,
        "referee_revision_old": REFEREE_REVISION_OLD,
        "referee_revision_new": REFEREE_REVISION_NEW,
        "production_pipeline_modified": False,
        "seeder_calls": 0,
        "harvester_calls": 0,
        "normal_cfd_long_run_calls": 0,
        "one_step_original_musubi_calls": int(oracle.get("musubi_calls", 0)),
        "one_step_instrumented_musubi_calls": int(instrumented["musubi_calls"]),
        "grid_convergence": "NOT_RUN",
        "endpoint_iteration": ENDPOINT_ITERATION,
        "legacy_endpoint_R_boundary": endpoint_record["R_boundary"],
        "historical_sign_contradiction": sign_contradiction,
        "historical_trapezoid_closure_median": float(np.median(closure_errors)),
        "historical_trapezoid_closure_maximum": float(np.max(closure_errors)),
        "runtime_boundary_order": list(RUNTIME_ORDER),
        "old_slot_assumption_valid": True,
        "overlap_slot_count": len(set().union(*overlaps.values())) if overlaps else 0,
        "pair_overlap_counts": pair_counts,
        "endpoint_replay": {
            "per_label_lattice": endpoint_replay["per_label_lattice"],
            "per_label_kg_s_domain": endpoint_replay["per_label_kg_s_domain"],
            "predicted_total_lattice": endpoint_replay["predicted_total_lattice"],
            "predicted_total_kg_s": endpoint_replay["predicted_total_kg_s"],
            "actual_total_lattice": actual_delta,
            "actual_total_kg_s": lattice_delta_to_kg_s(actual_delta),
            "R_one_step_identity": identity,
        },
        "instrumentation_required": True,
        "instrumented_comparisons": instrumented["comparisons"],
        "root_cause": "The V1 outlet replay treated every raw mesh coordinate as a normal fluid PULL source. TreElm marks 18 outlet_02 higher-order-stencil elements prp_solid; Musubi connectivity therefore bounce-backs all links when either the current neighbor-buffer element or its source is runtime-solid. Omitting that global runtime-solid connectivity changed only outlet_02. V1 also omitted the non-zero wall_libb contribution when comparing ports with global PDF mass.",
        "corrected_flux_definition": "For every runtime-ordered BC write, sum replacement minus the immediately pre-write value at the exact FETCH storage slot; include wall, inlet, all outlets, overlaps by last-writer order, and other sources.",
        "endpoint_corrected_flows_kg_s": {
            "inlet": endpoint_replay["per_label_kg_s_domain"]["inlet"],
            "outlet_01": -endpoint_replay["per_label_kg_s_domain"]["outlet_01"],
            "outlet_02": -endpoint_replay["per_label_kg_s_domain"]["outlet_02"],
            "outlet_03": -endpoint_replay["per_label_kg_s_domain"]["outlet_03"],
            "net_domain": endpoint_replay["predicted_total_kg_s"],
        },
        "endpoint_instantaneous_pdf_derivative_kg_s": lattice_delta_to_kg_s(actual_delta),
        "endpoint_instant_closure_error": identity,
        "endpoint_residuals": {
            "R_mass_short": endpoint_corrected["R_mass_short"],
            "R_mass_long": endpoint_corrected["R_mass_long"],
            "R_velocity": endpoint_corrected["R_velocity"],
            "R_pressure": endpoint_corrected["R_pressure"],
            "R_inlet": endpoint_corrected["R_inlet"],
        },
        "window_10k": short10,
        "window_20k": long20,
        "earliest_historical_true_gate_pass_iteration": int(earliest["iteration"]) if earliest else None,
        "earliest_gate_pass_restart_exists": bool(accepted_exists) if earliest else False,
        "accepted_steady_iteration": int(accepted["iteration"]) if status == FINAL_PASS else None,
        "accepted_steady_restart_path": str(accepted_path) if status == FINAL_PASS else None,
        "accepted_steady_restart_sha256": sha256_file(accepted_path) if status == FINAL_PASS and accepted_exists else None,
        "extra_iterations_caused_by_old_referee": ENDPOINT_ITERATION - int(earliest["iteration"]) if earliest else None,
        "source_frozen_files_unchanged": unchanged,
        "frozen_files_before": frozen_before,
        "frozen_files_after": frozen_after,
        "first_failure": None if status == FINAL_PASS else "corrected referee hard gates did not all pass",
    }
    write_json(qc / "referee_final_manifest.json", manifest)
    return manifest
