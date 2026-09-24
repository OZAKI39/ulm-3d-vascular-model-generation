"""Canonical source-proven replay of one complete Musubi timestep.

The production solver and its configuration deliberately do not import this
module.  It decomposes an existing restart-to-restart transition into the
boundary-write, PULL-connectivity, collision, and source phases used by the
pinned AOS/PULL D3Q19 BGK executable.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any, Mapping

import numpy as np

from .exact_link_flux import equilibrium_pdf, velocity_from_pdf
from .io import sha256_file
from .musubi_boundary_mass_referee import (
    MeshContract,
    boundary_blocked_pull_mask,
    conservation_identity_residual,
    pull_fetch_pdfs_runtime,
    replay_boundary_step,
    runtime_solid_cells,
)
from .restart_decode import D3Q19_DIRECTIONS
from .validated_contract import FULL_TIMESTEP_IDENTITY_GATE as FULL_IDENTITY_GATE


PREFERRED_FULL_IDENTITY_GATE = 1.0e-9
BOUNDARY_FLUX_DEFINITION = "BOUNDARY_PDF_ACCOUNTING_FLUX"
PHYSICAL_FLUX_DEFINITION = "PHYSICAL_CROSS_SECTION_FLUX"
DEFERRED_PHYSICAL_FLUX = "DEFERRED_TO_GRID_CONVERGENCE_STAGE"


def full_identity_pass(
    one_step_residuals: list[float] | tuple[float, ...],
    cumulative_residual: float,
    *,
    gate: float = FULL_IDENTITY_GATE,
) -> bool:
    """Apply the unchanged hard gate to every step and the full window."""

    if not one_step_residuals:
        return False
    values = (*[float(value) for value in one_step_residuals], float(cumulative_residual))
    return all(math.isfinite(value) and value <= float(gate) for value in values)


def stable_delta(new: np.ndarray, old: np.ndarray) -> float:
    """Sum a small elementwise delta without subtracting two large totals."""

    difference = np.subtract(
        np.asarray(new, dtype=np.float64),
        np.asarray(old, dtype=np.float64),
        dtype=np.float64,
    )
    return math.fsum(float(value) for value in difference.ravel())


def stable_total(values: np.ndarray) -> float:
    """Return a reproducible compensated total for diagnostic phase tables."""

    return math.fsum(float(value) for value in np.asarray(values).ravel())


def pull_compute_targets(
    state_after_boundary: np.ndarray,
    mesh: MeshContract,
    runtime_solid: set[int] | None = None,
) -> np.ndarray:
    """Materialize the PDFs fetched for every actual single-level target.

    The frozen Base has no coarse/fine ghosts.  Source evidence and the runtime
    log establish ``nElems_solve == nFluids == 182320`` globally; halo elements
    are source storage only and are not serialized restart rows.
    """

    values = np.asarray(state_after_boundary, dtype=np.float64)
    solid = runtime_solid_cells(mesh) if runtime_solid is None else runtime_solid
    targets = np.arange(values.shape[0], dtype=np.int64)
    return pull_fetch_pdfs_runtime(values, mesh, targets, solid)


def pull_link_counts(
    mesh: MeshContract, runtime_solid: set[int] | None = None
) -> dict[str, int]:
    """Classify every solve-target PULL link by its first runtime rule."""

    solid = runtime_solid_cells(mesh) if runtime_solid is None else runtime_solid
    blocked = boundary_blocked_pull_mask(mesh)
    counts = {
        "normal_fluid_source": 0,
        "missing_source_inverse": 0,
        "current_solid_inverse": 0,
        "source_solid_inverse": 0,
        "boundary_source_inverse": 0,
    }
    for target, coordinate in enumerate(mesh.cell_ijk):
        current_solid = target in solid
        for direction_index, direction in enumerate(D3Q19_DIRECTIONS):
            source = mesh.lookup.get(
                tuple(int(value) for value in coordinate - direction)
            )
            if current_solid:
                counts["current_solid_inverse"] += 1
            elif blocked[target, direction_index]:
                counts["boundary_source_inverse"] += 1
            elif source is None:
                counts["missing_source_inverse"] += 1
            elif int(source) in solid:
                counts["source_solid_inverse"] += 1
            else:
                counts["normal_fluid_source"] += 1
    counts["total"] = sum(counts.values())
    return counts


def tau1_bgk_collision(pull_fetched_state: np.ndarray) -> np.ndarray:
    """Numerically construct the pinned omega=1 D3Q19 BGK output.

    The optimized Musubi kernel reads density and velocity from its aux field;
    that aux field was computed from the same PULL-fetched PDFs.  At omega=1,
    the ``(1-omega)*f`` term vanishes and the kernel writes equilibrium PDFs.
    """

    fetched = np.asarray(pull_fetched_state, dtype=np.float64)
    density = np.sum(fetched, axis=1, dtype=np.float64)
    velocity = velocity_from_pdf(fetched)
    return equilibrium_pdf(density, velocity)


def replay_full_timestep(
    pdf_start: np.ndarray,
    pdf_end: np.ndarray,
    mesh: MeshContract,
    *,
    dx_m: float,
    dt_s: float,
    density_kg_m3: float,
    target_mass_flow_kg_s: float,
    outlet_pressures_pa: Mapping[str, float],
) -> dict[str, Any]:
    """Replay and score all mass-relevant phases of one existing transition."""

    start = np.asarray(pdf_start, dtype=np.float64)
    end = np.asarray(pdf_end, dtype=np.float64)
    if start.shape != end.shape or start.ndim != 2 or start.shape[1] != 19:
        raise ValueError("full timestep replay requires matching D3Q19 states")

    boundary = replay_boundary_step(
        start,
        mesh,
        dx_m=dx_m,
        dt_s=dt_s,
        density_kg_m3=density_kg_m3,
        target_mass_flow_kg_s=target_mass_flow_kg_s,
        outlet_pressures_pa=outlet_pressures_pa,
    )
    after_boundary = np.asarray(boundary["state"], dtype=np.float64)
    solid = set(int(value) for value in boundary["runtime_solid_cells"])
    pulled = pull_compute_targets(after_boundary, mesh, solid)
    collided = tau1_bgk_collision(pulled)

    delta_boundary = stable_delta(after_boundary, start)
    delta_pull = stable_delta(pulled, after_boundary)
    delta_collision = stable_delta(collided, pulled)
    delta_source = 0.0
    predicted = math.fsum(
        (delta_boundary, delta_pull, delta_collision, delta_source)
    )
    actual = stable_delta(end, start)
    target = float(boundary["target_lattice_flux"])
    per_label = {
        label: float(boundary["per_label_lattice"][label])
        for label in mesh.boundary_labels
    }
    wall_raw = _wall_operations_without_runtime(start, mesh)
    wall_runtime_only = _wall_operations_runtime_only(start, mesh, solid)
    wall_corrected = _wall_operations_with_runtime(start, mesh, solid)
    runtime_solid_affected = _changed_operation_count(
        wall_raw, wall_runtime_only
    )
    boundary_blocked_affected = _changed_operation_count(
        wall_runtime_only, wall_corrected
    )

    return {
        "delta_wall": per_label["wall"],
        "delta_outlet_02": per_label["outlet_02"],
        "delta_outlet_03": per_label["outlet_03"],
        "delta_inlet": per_label["inlet"],
        "delta_outlet_01": per_label["outlet_01"],
        "delta_boundary_total": delta_boundary,
        "delta_pull_connectivity": delta_pull,
        "sum_rho_before_collision": stable_total(
            np.sum(pulled, axis=1, dtype=np.float64)
        ),
        "delta_collision": delta_collision,
        "delta_source": delta_source,
        "source_status": "SOURCE_PROVEN_NONE",
        "delta_full_predicted": predicted,
        "delta_actual": actual,
        "R_boundary_only": conservation_identity_residual(
            delta_boundary, actual, target
        ),
        "R_full_one_step_identity": conservation_identity_residual(
            predicted, actual, target
        ),
        "runtime_solid_count": len(solid),
        "affected_wall_writes": _changed_operation_count(
            wall_raw, wall_corrected
        ),
        "runtime_solid_affected_wall_writes": runtime_solid_affected,
        "boundary_blocked_affected_wall_writes": boundary_blocked_affected,
        "compute_target_count": int(start.shape[0]),
        "exact_inlet_lattice_flux": target,
        "mass_after_boundary": stable_total(after_boundary),
        "mass_after_pull": stable_total(pulled),
        "collision_constructed_numerically": True,
        "flux_definition": BOUNDARY_FLUX_DEFINITION,
        "_state_after_boundary": after_boundary,
        "_pull_fetched_state": pulled,
        "_pdf_after_collision": collided,
    }


def _wall_operations_without_runtime(
    pdf: np.ndarray, mesh: MeshContract
) -> list[tuple[str, int, int, float]]:
    from .musubi_boundary_mass_referee import _wall_operations

    return _wall_operations(pdf, mesh)


def _wall_operations_with_runtime(
    pdf: np.ndarray, mesh: MeshContract, solid: set[int]
) -> list[tuple[str, int, int, float]]:
    from .musubi_boundary_mass_referee import _wall_operations

    return _wall_operations(pdf, mesh, solid)


def _wall_operations_runtime_only(
    pdf: np.ndarray, mesh: MeshContract, solid: set[int]
) -> list[tuple[str, int, int, float]]:
    from .musubi_boundary_mass_referee import _wall_operations

    return _wall_operations(
        pdf,
        mesh,
        solid,
        np.zeros((len(mesh.cell_ijk), 19), dtype=bool),
    )


def _changed_operation_count(
    first: list[tuple[str, int, int, float]],
    second: list[tuple[str, int, int, float]],
) -> int:
    if len(first) != len(second):
        raise ValueError("wall operation cardinality changed")
    changed = 0
    for left, right in zip(first, second, strict=True):
        if left[:3] != right[:3]:
            raise ValueError("wall operation ordering changed")
        changed += float(left[3]) != float(right[3])
    return changed


def public_step_record(replay: Mapping[str, Any]) -> dict[str, Any]:
    """Drop materialized arrays before writing the replay evidence."""

    return {key: value for key, value in replay.items() if not key.startswith("_")}


def source_token_evidence(path: Path, *tokens: str) -> dict[str, Any]:
    """Bind a contract statement to exact tokens in the pinned source tree."""

    text = path.read_text(encoding="utf-8")
    token_lines: dict[str, int] = {}
    for token in tokens:
        offset = text.find(token)
        if offset < 0:
            raise ValueError(f"source token absent from {path}: {token}")
        token_lines[token] = text.count("\n", 0, offset) + 1
    return {
        "path": str(path),
        "sha256": sha256_file(path),
        "token_lines": token_lines,
    }
