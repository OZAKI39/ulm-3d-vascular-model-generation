"""Display-only similarity glyphs; preserve every member and all within-group geometry."""

from __future__ import annotations

import numpy as np
import pyvista as pv
from scipy.optimize import minimize
from scipy.spatial.distance import pdist

from mevo_mb.coupling.bifurcation import BifurcationGeometry


def fit_rigid_group(geometry, positions, radii, initial_shift):
    shift = initial_shift.copy()
    for iteration in range(100):
        gap, normal, valid = geometry.wall_gap_and_normal(positions + shift, radii)
        assert valid.all()
        worst = int(np.argmin(gap))
        if gap[worst] >= 1e-8:
            return shift, iteration
        shift -= (1e-8 - gap[worst] + 1e-9) * normal[worst]

    # Near two opposing walls, successive normal corrections can converge very
    # slowly. Solve only the three rigid-translation coordinates, in mm for
    # numerical conditioning; never project or repack individual members.
    def constraints(s):
        g, n, ok = geometry.wall_gap_and_normal(positions + s / 1000, radii)
        assert ok.all()
        return (g - 1e-8) * 1000, -n

    origin = initial_shift * 1000
    fit = minimize(
        lambda s: 0.5 * np.sum((s - origin) ** 2),
        shift * 1000,
        jac=lambda s: s - origin,
        method="SLSQP",
        constraints=[dict(type="ineq", fun=lambda s: constraints(s)[0], jac=lambda s: constraints(s)[1])],
        options=dict(maxiter=500, ftol=1e-14),
    )
    gap, _, valid = geometry.wall_gap_and_normal(positions + fit.x / 1000, radii)
    if valid.all() and gap.min() >= 0:
        return fit.x / 1000, 100 + fit.nit
    raise ValueError(
        f"Rigid display group does not fit: {len(positions)} particles, min gap {gap.min()}, {fit.message}"
    )


def avoid_other_groups(geometry, positions, radii, shift, fixed, fixed_radii):
    """Translate the entire glyph to avoid falsely merging distinct groups."""
    if not len(fixed):
        return shift
    clearance = 5e-6

    def constraints(s):
        x = positions + s / 1000
        gap, normal, valid = geometry.wall_gap_and_normal(x, radii)
        assert valid.all()
        delta = x[:, None] - fixed[None]
        distance = np.linalg.norm(delta, axis=2)
        cross_gap = distance - radii[:, None] - fixed_radii[None] - clearance
        jac = delta / np.maximum(distance[..., None], 1e-30)
        return np.concatenate([(gap - 1e-8), cross_gap.ravel()]) * 1000, np.concatenate(
            [-normal, jac.reshape(-1, 3)]
        )

    initial = shift * 1000
    if constraints(initial)[0].min() >= -1e-12:
        return shift
    directions = np.array(
        [[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]], float
    )
    # Bounded local starts avoid a solver jumping to a remote feasible chamber
    # when a nearby tangential placement exists.
    for trust_mm in [0.25, 0.60, 1.2]:
        for direction in directions:
            solution = minimize(
                lambda s: 0.5 * np.sum((s - initial) ** 2),
                initial + direction * 0.05,
                jac=lambda s: s - initial,
                method="SLSQP",
                bounds=[(v - trust_mm, v + trust_mm) for v in initial],
                constraints=[
                    dict(type="ineq", fun=lambda s: constraints(s)[0], jac=lambda s: constraints(s)[1])
                ],
                options=dict(maxiter=150, ftol=1e-13),
            )
            if constraints(solution.x)[0].min() >= -1e-9:
                return solution.x / 1000
    raise ValueError(f"Cannot place {len(positions)}-bubble glyph without overlapping other groups")


def prepare(d, scale=24.0):
    geometry = BifurcationGeometry(d["config"]["geometry"])
    geometry.bind_mesh_surface(pv.read(d["surface"]))
    displayed = np.full_like(d["x"], np.nan)
    active = np.stack([s["active_mask"] for s in d["states"]])
    all_groups = []
    previous = {}
    max_shape_error = 0.0
    max_shift = 0.0
    min_gap = np.inf
    max_iterations = 0
    anchors = np.full(active.shape, -1, dtype=np.int32)
    shifts = np.zeros_like(d["x"])
    for k, groups in enumerate(d["clusters"]):
        if k % 50 == 0:
            print(f"layout t={d['t'][k]:.2f}s; scale={scale:g}", flush=True)
        current = {}
        rows = []
        fixed, fixed_radii = np.empty((0, 3)), np.empty(0)
        for g in groups:
            indices = np.array([d["id_index"][pid] for pid in g["particle_ids"]])
            real = d["x"][k, indices]
            candidates = [previous[int(pid)] for pid in g["particle_ids"] if int(pid) in previous]
            parents = {}
            for candidate in candidates:
                parents.setdefault(candidate["anchor"], []).append(candidate)
            if parents:
                parent = max(parents.values(), key=lambda values: (len(values), -values[0]["anchor"]))[0]
                anchor = parent["anchor"]
                if anchor not in indices:
                    anchor = int(indices[np.argmin(np.linalg.norm(real - d["x"][k, anchor], axis=1))])
                initial_shift = parent["shift"] * 0.99
            else:
                anchor = int(indices[np.argmin(np.linalg.norm(real - real.mean(axis=0), axis=1))])
                initial_shift = np.zeros(3)
            origin = d["x"][k, anchor]
            unshifted = origin + scale * (real - origin)
            radii = scale * d["radii"][indices]
            shift, iterations = fit_rigid_group(geometry, unshifted, radii, initial_shift)
            shift = avoid_other_groups(geometry, unshifted, radii, shift, fixed, fixed_radii)
            glyph = unshifted + shift
            fixed = np.concatenate([fixed, glyph])
            fixed_radii = np.concatenate([fixed_radii, radii])
            displayed[k, indices] = glyph
            anchors[k, indices] = anchor
            shifts[k, indices] = shift
            gap, _, valid = geometry.wall_gap_and_normal(glyph, radii)
            assert valid.all() and gap.min() >= 0
            error = float(np.max(np.abs(pdist(glyph) / scale - pdist(real)))) if len(real) > 1 else 0.0
            assert error < 1e-15
            max_shape_error = max(max_shape_error, error)
            max_shift = max(max_shift, float(np.linalg.norm(shift)))
            min_gap = min(min_gap, float(gap.min()))
            max_iterations = max(max_iterations, iterations)
            row = dict(
                particle_ids=g["particle_ids"],
                count=len(indices),
                anchor=int(d["ids"][anchor]),
                source_centroid_m=real.mean(axis=0).tolist(),
                display_centroid_m=glyph.mean(axis=0).tolist(),
                shift_m=shift.tolist(),
                scale=scale,
            )
            rows.append(row)
            for pid in g["particle_ids"]:
                current[int(pid)] = dict(anchor=anchor, shift=shift.copy())
        previous = current
        all_groups.append(rows)
    common = active[:-1] & active[1:]
    jumps = np.linalg.norm(np.diff(displayed, axis=0), axis=2)[common]
    report = dict(
        scale=scale,
        max_within_group_distance_error_m=max_shape_error,
        min_display_surface_gap_m=min_gap,
        max_rigid_boundary_shift_m=max_shift,
        max_fit_iterations=max_iterations,
        max_saved_display_step_m=float(jumps.max()),
        display_step_percentiles_m=np.percentile(jumps, [50, 90, 99, 99.9]).tolist(),
        no_members_merged=True,
        no_members_filtered=True,
        no_group_rotation=True,
        intergroup_display_clearance_m=5e-6,
        positions_are_display_glyphs=True,
        scientific_data_modified=False,
    )
    return (
        dict(
            positions_m=displayed,
            radii_m=d["radii"] * scale,
            active=active,
            anchor_indices=anchors,
            rigid_shifts_m=shifts,
            particle_id=d["ids"],
            time_s=d["t"],
        ),
        all_groups,
        report,
    )
