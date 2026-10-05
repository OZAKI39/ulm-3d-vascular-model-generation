"""Overdamped translation. Normal lubrication is a scoped rigid-interface approximation.

Acoustic forces and radial response are inputs, never modified here. All units SI.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from time import perf_counter

import numpy as np
from scipy import sparse
from scipy.linalg import lstsq
from scipy.optimize import minimize
from scipy.sparse.linalg import splu
from scipy.spatial import cKDTree

LEGACY_LABEL = "LEGACY_STOKES_TRANSLATION_BASELINE"
LUBRICATION_PHYSICS = "Stokes self-drag + normal particle-particle and particle-wall lubrication resistance"
INTERFACE_SCOPE = "tangentially immobile-interface approximation for near-contact hydrodynamic resistance"
NOT_INCLUDED = [
    "full thin-film drainage",
    "Young–Laplace contact deformation",
    "tangential lubrication",
    "shell repulsion",
    "adhesion",
    "fusion",
    "coalescence",
    "high-fidelity near-contact acoustic interaction",
]


@dataclass
class TranslationResult:
    velocity: np.ndarray
    stokes_velocity: np.ndarray
    matrix: sparse.csc_matrix
    drag: np.ndarray
    pairs: np.ndarray
    gaps: np.ndarray
    normals: np.ndarray
    wall: tuple
    diagnostics: dict = field(default_factory=dict)
    factor: object = None
    scale: float = 1.0


class TranslationBackend(ABC):
    @abstractmethod
    def solve(self, positions_m, radii_m, force_N, flow_velocity_m_s, wall=None):
        """Unconstrained translation; geometric contact is a separate operation."""


class LegacyStokesTranslation(TranslationBackend):
    name = "legacy_stokes"

    def __init__(self, viscosity_Pa_s):
        self.mu = float(viscosity_Pa_s)
        if not np.isfinite(self.mu) or self.mu <= 0:
            raise ValueError("Positive finite viscosity required")

    def solve(self, positions_m, radii_m, force_N, flow_velocity_m_s, wall=None):
        r = np.asarray(radii_m)
        if np.any(r <= 0) or not np.isfinite(r).all():
            raise ValueError("Positive finite radii required")
        d = 6 * np.pi * self.mu * r
        # Preserve the operation order of the legacy tracker exactly.
        v = np.asarray(flow_velocity_m_s) + np.asarray(force_N) / (6 * np.pi * self.mu * r[:, None])
        return TranslationResult(
            v,
            v.copy(),
            sparse.diags(np.repeat(d, 3), format="csc"),
            d,
            np.empty((0, 2), int),
            np.empty(0),
            np.empty((0, 3)),
            wall,
            {"backend": self.name, "identity": LEGACY_LABEL},
        )


def _pairs(x, cutoff):
    pairs = cKDTree(x).query_pairs(np.nextafter(cutoff, np.inf), output_type="ndarray")
    return pairs[np.lexsort((pairs[:, 1], pairs[:, 0]))] if len(pairs) else pairs.reshape(-1, 2)


def _gap_vectors(x, pairs, radius):
    delta = x[pairs[:, 1]] - x[pairs[:, 0]]
    distance = np.linalg.norm(delta, axis=1)
    if np.any(distance == 0):
        raise ValueError("COINCIDENT_PARTICLE_CENTERS")
    return distance - 2 * radius, delta / distance[:, None]


class LubricationCoupledTranslation(LegacyStokesTranslation):
    name = "lubrication_coupled"

    def __init__(self, viscosity_Pa_s, activation=0.1, floor=0.001, residual_tolerance=1e-11):
        super().__init__(viscosity_Pa_s)
        if not 0 < floor < activation:
            raise ValueError("Require 0 < floor < activation")
        self.activation, self.floor = float(activation), float(floor)
        self.residual_tolerance = float(residual_tolerance)
        self.last_failure = None

    @staticmethod
    def mono_radius(radii):
        r = np.asarray(radii, float)
        if not len(r) or np.any(r <= 0) or not np.isfinite(r).all():
            raise ValueError("Positive finite radii required")
        if np.any(np.abs(r - r[0]) > 8 * np.finfo(float).eps * r[0]):
            raise NotImplementedError("UNEQUAL_RADIUS_LUBRICATION_NOT_IMPLEMENTED")
        return float(r[0])

    def coefficient(self, gap, radius, kind="pair"):
        e = np.asarray(gap, float) / radius
        effective = np.clip(e, self.floor, self.activation)
        if kind == "pair":

            def phi(v):
                return 1 / (2 * v) - 9 / 20 * np.log(v) - 3 / 56 * v * np.log(v)

            prefactor = 3 * np.pi * self.mu * radius
        elif kind == "wall":

            def phi(v):
                return 1 / v - 1 / 5 * np.log(v) - 1 / 21 * v * np.log(v)

            prefactor = 6 * np.pi * self.mu * radius
        else:
            raise ValueError("kind must be pair or wall")
        return np.where(
            e < self.activation, prefactor * np.maximum(phi(effective) - phi(self.activation), 0), 0.0
        )

    def solve(self, positions_m, radii_m, force_N, flow_velocity_m_s, wall=None):
        started = perf_counter()
        x, r = np.asarray(positions_m, float), np.asarray(radii_m, float)
        radius = self.mono_radius(r)
        result = super().solve(x, r, force_N, flow_velocity_m_s, wall)
        n = len(r)
        pairs = _pairs(x, (2 + self.activation) * radius)
        gaps, normals = _gap_vectors(x, pairs, radius)
        a = self.coefficient(gaps, radius)
        take = a > 0
        pairs, gaps, normals, a = pairs[take], gaps[take], normals[take], a[take]
        if wall is None:
            wall = (np.full(n, np.inf), np.zeros((n, 3)), np.ones(n, bool))
        wg, wn, valid = wall
        wg, wn, valid = np.asarray(wg), np.asarray(wn), np.asarray(valid)
        if not valid.all():
            raise RuntimeError("INVALID_WALL_GAP_OR_NORMAL")
        wa = self.coefficient(wg, radius, "wall")
        wi = np.flatnonzero(wa > 0)
        rows, cols, values = [], [], []

        def block(i, j, value):
            rows.extend(np.repeat(3 * i + np.arange(3), 3))
            cols.extend(np.tile(3 * j + np.arange(3), 3))
            values.extend(value.ravel())

        for (i, j), normal, coef in zip(pairs, normals, a):
            b = coef * np.outer(normal, normal)
            block(i, i, b)
            block(j, j, b)
            block(i, j, -b)
            block(j, i, -b)
        for i in wi:
            block(i, i, wa[i] * np.outer(wn[i], wn[i]))
        correction = sparse.coo_matrix((values, (rows, cols)), shape=(3 * n, 3 * n)).tocsc()
        matrix = result.matrix + correction
        matrix.eliminate_zeros()
        rhs = (np.asarray(force_N) + result.drag[:, None] * flow_velocity_m_s).ravel()
        scale = float(np.median(result.drag))
        solve_start = perf_counter()
        factor = None
        if len(pairs) or len(wi):
            factor = splu(matrix / scale)
            velocity = factor.solve(rhs / scale).reshape(n, 3)
            solver = "scipy.sparse.linalg.splu"
        else:
            velocity = result.velocity  # bitwise far-field bypass
            solver = "exact_legacy_bypass"
        residual = float(
            np.linalg.norm(matrix @ velocity.ravel() - rhs) / max(np.linalg.norm(rhs), np.finfo(float).tiny)
        )
        diag = dict(
            backend=self.name,
            matrix_size=3 * n,
            nnz=matrix.nnz,
            solver=solver,
            relative_residual=residual,
            solve_time_s=perf_counter() - solve_start,
            failure_state=None,
            LUBRICATION_ACTIVE=int(len(pairs) + len(wi)),
            SUBFLOOR=int(
                np.count_nonzero(gaps <= self.floor * radius) + np.count_nonzero(wg <= self.floor * radius)
            ),
            pair_active=len(pairs),
            wall_active=len(wi),
            assembly_and_solve_time_s=perf_counter() - started,
        )
        if not np.isfinite(velocity).all() or residual > self.residual_tolerance:
            diag["failure_state"] = "LINEAR_SOLVER_RESIDUAL_FAILURE"
            self.last_failure = diag
            raise RuntimeError(str(diag))
        return TranslationResult(
            velocity,
            result.stokes_velocity,
            matrix,
            result.drag,
            pairs,
            gaps,
            normals,
            wall,
            diag,
            factor,
            scale,
        )

    def contact(self, result, x, radii, velocity, tolerance, horizon=None):
        """A-metric frictionless projection; finite horizon allows subfloor closure.

        At contact GU>=0. With horizon, positive gaps have GU>=-h/horizon,
        permitting approach to zero (never locking at the coefficient floor).
        """
        radius = self.mono_radius(radii)
        if horizon is None:
            pairs = _pairs(x, 2 * radius + tolerance)
        else:
            centered = velocity - np.mean(velocity, axis=0)
            reach = 2 * np.max(np.linalg.norm(centered, axis=1)) * horizon
            pairs = _pairs(x, 2 * radius + self.floor * radius + reach)
        gaps, normals = _gap_vectors(x, pairs, radius)
        wg, wn, _ = result.wall
        wall_ids = np.flatnonzero(
            wg
            <= (
                tolerance
                if horizon is None
                else np.maximum(self.floor * radius, np.linalg.norm(velocity, axis=1) * horizon)
            )
        )
        m, size = len(pairs) + len(wall_ids), velocity.size
        diag = {
            "HARD_CONTACT": int(np.count_nonzero(gaps <= tolerance) + np.count_nonzero(wg <= tolerance)),
            "contact_constraints": m,
            "contact_metric": "A",
            "contact_kkt_residual": 0.0,
            "contact_solve_time_s": 0.0,
        }
        if not m:
            return velocity, diag
        g = np.zeros((m, size))
        for k, ((i, j), normal) in enumerate(zip(pairs, normals)):
            g[k, 3 * i : 3 * i + 3], g[k, 3 * j : 3 * j + 3] = -normal, normal
        for k, i in enumerate(wall_ids, len(pairs)):
            g[k, 3 * i : 3 * i + 3] = -wn[i]
        bounds = np.zeros(m)
        if horizon is not None:
            bounds = -np.maximum(np.r_[gaps, wg[wall_ids]], 0) / horizon
            bounds[np.r_[gaps, wg[wall_ids]] <= tolerance] = 0
        vscale = max(np.max(np.abs(velocity)), np.finfo(float).tiny)
        u = velocity.ravel() / vscale
        q = g @ u - bounds / vscale
        if np.min(q) >= -1e-12:
            return velocity, diag
        started = perf_counter()
        ainvg = (
            result.factor.solve(g.T)
            if result.factor is not None
            else g.T / (np.repeat(result.drag, 3) / result.scale)[:, None]
        )
        mtx = g @ ainvg
        # Positive diagonal equilibration of the contact dual (an LCP).
        d = np.sqrt(np.diag(mtx))
        mm = mtx / d[:, None] / d[None, :]
        qq = q / d
        fit = minimize(
            lambda z: (0.5 * z @ mm @ z + qq @ z, mm @ z + qq),
            np.zeros(m),
            jac=True,
            method="L-BFGS-B",
            bounds=[(0, None)] * m,
            options={"ftol": 1e-15, "gtol": 1e-12, "maxiter": 5000, "maxls": 50},
        )
        z = np.maximum(fit.x, 0)
        # Polish the active set; do not accept an optimizer's success flag alone.
        active = z > 1e-10 * max(1.0, np.max(z))
        for _ in range(4 * m + 10):
            trial = np.zeros(m)
            if active.any():
                trial[active] = lstsq(
                    mm[np.ix_(active, active)], -qq[active], cond=1e-13, lapack_driver="gelsy"
                )[0]
            negative = active & (trial < 0)
            if negative.any():
                ratio = z[negative] / np.maximum(z[negative] - trial[negative], np.finfo(float).tiny)
                alpha = min(1.0, float(np.min(ratio)))
                z += alpha * (trial - z)
                remove = active & (z <= 1e-13 * max(1.0, np.max(z))) & (trial < 0)
                active[remove] = False
                z[~active] = 0
                continue
            z = trial
            grad = mm @ z + qq
            violation = np.where(~active, grad, np.inf)
            if np.min(violation) >= -1e-11:
                break
            active[np.argmin(violation)] = True
        u += ainvg @ (z / d)
        slack = g @ u - bounds / vscale
        kkt = max(
            float(max(0.0, -np.min(slack))),
            float(np.max(np.abs(slack[z > 1e-9 * max(1.0, np.max(z))])))
            if np.any(z > 1e-9 * max(1.0, np.max(z)))
            else 0.0,
        )
        diag.update(
            contact_kkt_residual=kkt,
            contact_solve_time_s=perf_counter() - started,
            contact_dual_optimizer=fit.message,
            contact_active_multipliers=int(np.count_nonzero(z)),
        )
        if kkt > 1e-9 or not np.isfinite(u).all():
            self.last_failure = dict(diag, failure_state="A_METRIC_CONTACT_KKT_FAILURE")
            raise RuntimeError(str(self.last_failure))
        return u.reshape(velocity.shape) * vscale, diag

    def timestep(self, x, radii, velocity, wall, remaining):
        """Predictive neighborhood plus half-current-gap guard; floor is not a target."""
        radius = self.mono_radius(radii)
        centered = velocity - np.mean(velocity, axis=0)
        reach = 2 * np.max(np.linalg.norm(centered, axis=1)) * remaining
        pairs = _pairs(x, (2 + self.activation) * radius + reach)
        gaps, normals = _gap_vectors(x, pairs, radius)
        closing = -np.einsum("ij,ij->i", velocity[pairs[:, 1]] - velocity[pairs[:, 0]], normals)
        wg, wn, valid = wall
        if not np.asarray(valid).all():
            raise RuntimeError("INVALID_WALL_GAP_OR_NORMAL")
        wc = np.einsum("ij,ij->i", velocity, wn)
        all_gaps, all_closing = np.r_[gaps, wg], np.r_[closing, wc]
        take = (all_gaps > self.floor * radius) & (all_closing > 0)
        dt = min(remaining, float(np.min(0.5 * all_gaps[take] / all_closing[take], initial=np.inf)))
        if not np.isfinite(dt) or dt <= 0:
            raise RuntimeError("NONPOSITIVE_OR_NONFINITE_TIMESTEP")
        return dt


def make_translation_backend(config):
    p = config["particles"]
    name = p.get("translation_backend", "legacy_stokes")
    mu = config["medium"]["viscosity_Pa_s"]
    if name == "legacy_stokes":
        return LegacyStokesTranslation(mu)
    if name == "lubrication_coupled":
        return LubricationCoupledTranslation(mu, **p.get("lubrication", {}))
    raise ValueError(f"Unknown translation backend: {name}")
