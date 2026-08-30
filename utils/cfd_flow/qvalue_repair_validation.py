"""Validation gates for the research-only Seeder qVal repair."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np

from .grid_convergence import _mesh_qc, grid_specs
from .qvalue_contract_forensics import (
    sample_base_stl_ray_oracle,
    vascular_wall_qvalue_distribution,
)
from .wall_qvalue_oracle import audit_mesh_qvalues


TINY_MEDIAN_ABS_LIMIT = 0.02
TINY_P95_ABS_LIMIT = 0.05
MINIMUM_CONTINUOUS_UNIQUE = 16


def q_error_metrics(errors: np.ndarray) -> dict[str, float | int | None]:
    """Return the repair-gate statistics for signed q errors."""

    values = np.asarray(errors, dtype=np.float64)
    values = values[np.isfinite(values)]
    if not len(values):
        return {
            "count": 0,
            "bias": None,
            "rms": None,
            "median_absolute": None,
            "p95_absolute": None,
            "max_absolute": None,
        }
    absolute = np.abs(values)
    return {
        "count": int(len(values)),
        "bias": float(np.mean(values)),
        "rms": float(np.sqrt(np.mean(values * values))),
        "median_absolute": float(np.median(absolute)),
        "p95_absolute": float(np.percentile(absolute, 95.0)),
        "max_absolute": float(np.max(absolute)),
    }


def tiny_cylinder_gate(mesh_dir: Path, mesh_summary: Path) -> dict[str, Any]:
    """Apply the analytic-cylinder gate to an already-generated tiny mesh."""

    audit, errors = audit_mesh_qvalues(
        mesh_dir, mesh_summary, return_error_samples=True
    )
    metrics = q_error_metrics(errors)
    support = audit["q_value_distribution_on_declared_wall_links"]
    continuous_unique = sum(
        not math.isclose(float(value), 0.5, rel_tol=0.0, abs_tol=1.0e-10)
        and not math.isclose(float(value), 1.0, rel_tol=0.0, abs_tol=1.0e-10)
        for value in support
    )
    median = metrics["median_absolute"]
    p95 = metrics["p95_absolute"]
    bias = metrics["bias"]
    median_pass = median is not None and float(median) <= TINY_MEDIAN_ABS_LIMIT
    p95_pass = p95 is not None and float(p95) <= TINY_P95_ABS_LIMIT
    bias_pass = bias is not None and abs(float(bias)) <= TINY_MEDIAN_ABS_LIMIT
    support_pass = (
        continuous_unique >= MINIMUM_CONTINUOUS_UNIQUE
        and audit["uniform_halfway_fallback_fraction"] < 1.0
    )
    passed = (
        audit["status"].startswith("PASS")
        and support_pass
        and median_pass
        and p95_pass
        and bias_pass
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "audit": audit,
        "q_error": metrics,
        "continuous_unique_count": int(continuous_unique),
        "gates": {
            "structural": audit["status"].startswith("PASS"),
            "continuous_support": support_pass,
            "median_absolute_le_0p02": median_pass,
            "p95_absolute_le_0p05": p95_pass,
            "no_large_systematic_bias": bias_pass,
        },
    }


def topology_difference(
    old_qc: dict[str, Any], new_qc: dict[str, Any]
) -> dict[str, Any]:
    """Report topology changes without requiring identical repaired counts."""

    old_cells = int(old_qc["fluid_cell_count"])
    new_cells = int(new_qc["fluid_cell_count"])
    labels = sorted(
        set(old_qc["boundary_cell_counts"]) | set(new_qc["boundary_cell_counts"])
    )
    return {
        "old_fluid_cell_count": old_cells,
        "new_fluid_cell_count": new_cells,
        "fluid_cell_difference": new_cells - old_cells,
        "fluid_cell_relative_difference": (new_cells - old_cells) / old_cells,
        "old_connected_fluid_regions": int(old_qc["connected_fluid_regions"]),
        "new_connected_fluid_regions": int(new_qc["connected_fluid_regions"]),
        "boundary_cell_differences": {
            label: int(new_qc["boundary_cell_counts"].get(label, 0))
            - int(old_qc["boundary_cell_counts"].get(label, 0))
            for label in labels
        },
        "topology_gate_pass": int(new_qc["connected_fluid_regions"]) == 1,
        "identical_counts_not_required": True,
    }


def repaired_base_gate(
    *,
    old_mesh: Path,
    repaired_mesh: Path,
    wall_stl: Path,
    target_comparable: int = 10_000,
) -> dict[str, Any]:
    """Run read-only topology, distribution, and STL-ray gates on BASE."""

    spec = grid_specs()["base"]
    old_qc = _mesh_qc(old_mesh, spec)
    new_qc = _mesh_qc(repaired_mesh, spec)
    distribution = vascular_wall_qvalue_distribution(repaired_mesh)
    oracle = sample_base_stl_ray_oracle(
        repaired_mesh, wall_stl, target_comparable=target_comparable
    )
    difference = topology_difference(old_qc, new_qc)
    qvalue_pass = (
        distribution["status"] == "PASS"
        and oracle["status"] == "PASS_SAMPLE_SIZE"
    )
    topology_pass = new_qc["status"] == "PASS" and difference["topology_gate_pass"]
    passed = qvalue_pass and topology_pass
    return {
        "status": "PASS" if passed else "FAIL",
        "qvalue_gate_pass": qvalue_pass,
        "topology_port_gate_pass": topology_pass,
        "old_mesh_qc": old_qc,
        "repaired_mesh_qc": new_qc,
        "topology_difference": difference,
        "qvalue_distribution": distribution,
        "stl_ray_oracle": oracle,
        "first_failure": (
            None
            if passed
            else (
                "Repaired BASE has more than one connected fluid component."
                if not difference["topology_gate_pass"]
                else "Repaired BASE qVal validation failed."
            )
        ),
    }


def libb_reference_replacement(
    qvalue: float, *, f_in: float, f_out: float, f_neigh: float
) -> float:
    """Pinned Musubi ``set_bouzidi_coeff``/``wall_libb`` scalar equation."""

    q = float(qvalue)
    if not math.isfinite(q) or q <= 0.0:
        raise ValueError("qvalue must be finite and positive")
    if q >= 0.5:
        return (1.0 - 0.5 / q) * float(f_in) + (0.5 / q) * float(f_out)
    return 2.0 * q * float(f_out) + (1.0 - 2.0 * q) * float(f_neigh)
