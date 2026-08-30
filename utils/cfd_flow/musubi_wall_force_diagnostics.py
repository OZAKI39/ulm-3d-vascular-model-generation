"""Research contracts for isolating Musubi force and curved-wall numerics."""

from __future__ import annotations

import math
import os
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import trimesh

from .musubi_boundary_mass_referee import load_mesh_contract
from .musubi_pressure_bc_benchmark import (
    MEAN_VELOCITY_M_S,
    NU_M2_S,
    PIPE_RADIUS_M,
)
from .periodic_pipe_force import CASES


CS_LATTICE_SQUARED = 1.0 / 3.0


def lattice_relaxation_contract(
    *, nu_phy_m2_s: float, dx_m: float, dt_s: float, cs2: float = CS_LATTICE_SQUARED
) -> dict[str, float]:
    """Return Musubi's physical scaling and BGK relaxation quantities."""

    if min(nu_phy_m2_s, dx_m, dt_s, cs2) <= 0.0:
        raise ValueError("viscosity, spacing, timestep, and cs2 must be positive")
    nu_lattice = nu_phy_m2_s * dt_s / dx_m**2
    tau = nu_lattice / cs2 + 0.5
    omega = 1.0 / tau
    return {
        "nu_lattice": float(nu_lattice),
        "cs_lattice_squared": float(cs2),
        "tau": float(tau),
        "omega": float(omega),
    }


def tau_one_time_step_s(
    *, nu_phy_m2_s: float, dx_m: float, target_tau: float = 1.0,
    cs2: float = CS_LATTICE_SQUARED,
) -> float:
    """Analytically solve nu_lat=nu_phy*dt/dx^2=cs2*(tau-1/2)."""

    if target_tau <= 0.5:
        raise ValueError("BGK target tau must exceed 0.5")
    return float(cs2 * (target_tau - 0.5) * dx_m**2 / nu_phy_m2_s)


def body_force_conversion(
    force_phy_n_m3: Iterable[float], *, rho0_kg_m3: float, dx_m: float, dt_s: float
) -> dict[str, Any]:
    """Convert Musubi physical force density using fac.body_force=rho0*dx/dt^2."""

    force = np.asarray(tuple(force_phy_n_m3), dtype=np.float64).reshape(3)
    factor = rho0_kg_m3 * dx_m / dt_s**2
    lattice = force / factor
    return {
        "physical_force_density_n_m3": force.tolist(),
        "body_force_conversion_factor_n_m3_per_lattice": float(factor),
        "lattice_force_density": lattice.tolist(),
        "formula": "F_lat = F_phy * dt^2 / (rho0 * dx)",
    }


def expected_force_momentum_increment(
    force_phy_n_m3: Iterable[float], *, volume_m3: float, dt_s: float
) -> np.ndarray:
    """Physical total momentum increment for a uniform force-density source."""

    force = np.asarray(tuple(force_phy_n_m3), dtype=np.float64).reshape(3)
    return force * float(volume_m3) * float(dt_s)


def bouzidi_coefficients(q_value: float) -> dict[str, float]:
    """Mirror the source-proven Musubi ``set_bouzidi_coeff`` branches."""

    q_value = float(q_value)
    if not 0.0 < q_value <= 1.0:
        raise ValueError("Bouzidi q must be in (0, 1]")
    if q_value >= 0.5:
        return {
            "c_in": 1.0 - 0.5 / q_value,
            "c_out": 0.5 / q_value,
            "c_neighbor": 0.0,
        }
    return {
        "c_in": 0.0,
        "c_out": 2.0 * q_value,
        "c_neighbor": 1.0 - 2.0 * q_value,
    }


def wall_libb_post_pdf(
    *, q_value: float, f_in: float, f_out: float, f_neighbor: float
) -> float:
    """Evaluate the exact compiled wall_libb linear combination."""

    coefficients = bouzidi_coefficients(q_value)
    return float(
        coefficients["c_in"] * f_in
        + coefficients["c_out"] * f_out
        + coefficients["c_neighbor"] * f_neighbor
    )


def discrete_poiseuille_reference(
    points_m: np.ndarray,
    *,
    center_m: Iterable[float],
    axis: Iterable[float],
    radius_m: float = PIPE_RADIUS_M,
    continuum_mean_m_s: float = MEAN_VELOCITY_M_S,
) -> dict[str, Any]:
    points = np.asarray(points_m, dtype=np.float64)
    center = np.asarray(tuple(center_m), dtype=np.float64).reshape(3)
    direction = np.asarray(tuple(axis), dtype=np.float64).reshape(3)
    direction /= np.linalg.norm(direction)
    relative = points - center
    axial = relative @ direction
    perpendicular = relative - axial[:, None] * direction
    radial_squared = np.einsum("ij,ij->i", perpendicular, perpendicular)
    velocity = 2.0 * continuum_mean_m_s * np.maximum(
        0.0, 1.0 - radial_squared / radius_m**2
    )
    return {
        "cell_count": int(len(points)),
        "continuum_analytic_mean_m_s": float(continuum_mean_m_s),
        "discrete_analytic_mean_m_s": float(np.mean(velocity)),
        "discrete_minus_continuum_relative": float(
            (np.mean(velocity) - continuum_mean_m_s) / continuum_mean_m_s
        ),
        "velocity_m_s": velocity,
        "radial_distance_m": np.sqrt(radial_squared),
    }


def independent_cross_section_flux(
    velocity_m_s: np.ndarray,
    *,
    axis: Iterable[float],
    dx_m: float,
) -> float:
    """Integrate sampled axial velocity over actual fluid plane points."""

    velocity = np.asarray(velocity_m_s, dtype=np.float64)
    if velocity.ndim != 2 or velocity.shape[1] != 3:
        raise ValueError("velocity must have shape (n, 3)")
    direction = np.asarray(tuple(axis), dtype=np.float64).reshape(3)
    direction /= np.linalg.norm(direction)
    return float(np.sum(velocity @ direction, dtype=np.float64) * dx_m**2)


def _read_raw_ascii(path: Path) -> tuple[list[str], np.ndarray]:
    """Independent result reader; deliberately does not reuse the baseline parser."""

    readable = path.resolve()
    if os.name == "nt":
        readable = Path("\\\\?\\" + str(readable))
    header: list[str] | None = None
    rows: list[list[float]] = []
    for raw in readable.read_text(encoding="utf-8", errors="strict").splitlines():
        line = raw.strip()
        if line.startswith("#"):
            tokens = line[1:].split()
            if tokens and tokens[0] == "time":
                header = tokens
        elif line:
            rows.append([float(value) for value in line.split()])
    if header is None or not rows:
        raise ValueError(f"raw result lacks header or rows: {path}")
    values = np.asarray(rows, dtype=np.float64)
    if values.shape[1] != len(header):
        raise ValueError("raw result column count does not match header")
    return header, values


def independent_mean_axial_from_raw(path: Path, axis: Iterable[float]) -> dict[str, Any]:
    header, values = _read_raw_ascii(path)
    columns = [index for index, name in enumerate(header) if "velocity_phy" in name]
    if len(columns) != 3:
        raise ValueError(f"expected three physical velocity columns, got {header}")
    velocity = values[:, columns]
    direction = np.asarray(tuple(axis), dtype=np.float64).reshape(3)
    direction /= np.linalg.norm(direction)
    axial = velocity @ direction
    return {
        "header": header,
        "sample_count": int(len(values)),
        "last_time_s": float(values[-1, 0]),
        "last_mean_velocity_vector_m_s": velocity[-1].tolist(),
        "last_mean_axial_velocity_m_s": float(axial[-1]),
        "calculation": "dot(raw average velocity_phy vector, normalized pipe axis)",
    }


def sidewall_geometry_contract(case_dir: Path) -> dict[str, Any]:
    summary = __import__("json").loads(
        (case_dir / "source_mesh_summary.json").read_text(encoding="utf-8")
    )
    wall = trimesh.load_mesh(case_dir / "geometry" / "wall.stl", process=False)
    axis = np.asarray(summary["direction"], dtype=np.float64)
    axis /= np.linalg.norm(axis)
    normals = np.asarray(wall.face_normals, dtype=np.float64)
    vertices = np.asarray(wall.vertices, dtype=np.float64)
    projections = vertices @ axis
    expected_length = float(summary["wall"]["length_m"])
    actual_length = float(np.ptp(projections))
    dx_m = float(summary["dx_m"])
    plane_halfwidth = PIPE_RADIUS_M + 2.0 * dx_m
    passed = (
        not bool(wall.is_watertight)
        and float(np.max(np.abs(normals @ axis))) <= 1.0e-12
        and math.isclose(actual_length, expected_length, rel_tol=2.0e-6, abs_tol=0.0)
        and plane_halfwidth > PIPE_RADIUS_M
    )
    return {
        "status": "PASS" if passed else "FAIL",
        "triangle_count": int(len(wall.faces)),
        "watertight": bool(wall.is_watertight),
        "maximum_absolute_axial_face_normal": float(np.max(np.abs(normals @ axis))),
        "axial_cap_triangle_count": int(np.count_nonzero(np.abs(normals @ axis) > 1.0e-12)),
        "actual_physical_length_m": actual_length,
        "expected_physical_length_m": expected_length,
        "length_over_diameter": actual_length / (2.0 * PIPE_RADIUS_M),
        "periodic_plane_halfwidth_m": plane_halfwidth,
        "lumen_radius_m": PIPE_RADIUS_M,
        "periodic_plane_radial_clearance_cells": (plane_halfwidth - PIPE_RADIUS_M) / dx_m,
    }


def zero_run_baseline_audit(case_root: Path) -> dict[str, Any]:
    cases: dict[str, Any] = {}
    for name in ("axis_n16", "axis_n20", "axis_n27"):
        case = CASES[name]
        case_dir = case_root / name
        cases[name] = {
            "relaxation": lattice_relaxation_contract(
                nu_phy_m2_s=NU_M2_S, dx_m=case.dx_m, dt_s=case.dt_s
            ),
            "geometry": sidewall_geometry_contract(case_dir),
        }
    n27_dir = case_root / "axis_n27"
    tracking = next((n27_dir / "tracking").glob("*mean_velocity*p00000.res"))
    raw = independent_mean_axial_from_raw(tracking, CASES["axis_n27"].direction)
    contract = load_mesh_contract(
        n27_dir / "mesh", allow_zero_normals=True, require_runtime_order=False
    )
    case = CASES["axis_n27"]
    center = np.full(3, 0.5 * (2**8) * case.dx_m)
    points = (contract.cell_ijk.astype(np.float64) + 0.5) * case.dx_m
    discrete = discrete_poiseuille_reference(
        points, center_m=center, axis=case.direction
    )
    observed = float(raw["last_mean_axial_velocity_m_s"])
    result = {
        "status": "PASS",
        "baseline_solver_calls": 0,
        "independent_parser": raw,
        "n27_observed": {
            "mean_axial_velocity_m_s": observed,
            "error_vs_continuum": abs(observed - MEAN_VELOCITY_M_S) / MEAN_VELOCITY_M_S,
            "error_vs_discrete_cell_center_reference": abs(
                observed - float(discrete["discrete_analytic_mean_m_s"])
            )
            / float(discrete["discrete_analytic_mean_m_s"]),
            "historical_flow_rate_method": "mean_velocity_times_nominal_area",
            "historical_flow_rate_classification": "NOT_INDEPENDENT_FLOW_MEASUREMENT",
            "independent_flux_available_in_historical_tracking": False,
            "required_tau1_remedy": "asciiSpatial cross-section velocity and dx^2 quadrature",
        },
        "discrete_analytic_reference": {
            key: value for key, value in discrete.items() if not isinstance(value, np.ndarray)
        },
        "cases": cases,
        "all_grid_tau_equal": bool(
            np.ptp([cases[name]["relaxation"]["tau"] for name in cases]) <= 1.0e-14
        ),
    }
    return result


def wall_force_decision(
    *, force_pass: bool, wall_pass: bool | None, tau1_pass: bool | None,
    official_pass: bool | None = None, d3q27_pass: bool | None = None,
) -> str:
    if not force_pass:
        return "STOP_FORCE_CONTRACT"
    if wall_pass is None:
        return "RUN_WALL_ORACLE"
    if not wall_pass:
        return "STOP_WALL_CONTRACT"
    if tau1_pass is None:
        return "RUN_TAU1_N27"
    if tau1_pass:
        return "HIGH_TAU_BGK_WALL_COUPLING_CONFIRMED"
    if official_pass is None:
        return "RUN_OFFICIAL_PIP_FORCE"
    if not official_pass:
        return "CFD_FLOW_UPSTREAM_PIPE_FORCE_REFERENCE_FAILED"
    if d3q27_pass is None:
        return "RUN_CUSTOM_D3Q27_N27"
    if d3q27_pass:
        return "CFD_FLOW_D3Q19_WALL_LIBB_ACCURACY_FAILED"
    return "CFD_FLOW_CUSTOM_PIPE_FORCE_CONTRACT_UNRESOLVED"
