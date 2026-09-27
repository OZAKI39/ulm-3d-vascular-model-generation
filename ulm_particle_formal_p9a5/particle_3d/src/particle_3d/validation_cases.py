"""Deterministic validation inputs, independent analytical answers and raw records.

This module constructs checks/plots only; it never modifies the frozen field.
"""
from itertools import combinations
import numpy as np
from .field import FrozenFEMField
from .geometry import EPS

SEED = 20260920
# Fixed before measuring results. 256 eps covers a small 3x3 affine solve and
# evaluation on this well-conditioned, O(1) synthetic tetra (cond_inf = 1.94).
AFFINE_RTOL = 256 * EPS
AFFINE_ATOL = 256 * EPS


def affine_case():
    vertices = np.array([[.1, -.2, .3], [1.3, -.1, .2], [.3, 1.2, .5], [.2, .1, 1.5]])
    matrix = np.array([[2., -3., .5], [4., .25, -2.], [-1.5, 3., 1.]])
    offset = np.array([.7, -.4, 1.1])
    pressure_gradient = np.array([2.5, -1.25, 3.75])
    pressure_offset = 7.
    weights = [np.full(4, .25)]
    kinds = ["centroid"]
    for w in np.random.default_rng(SEED).dirichlet(np.ones(4), 128):
        weights.append(w); kinds.append("random_interior")
    for missing in range(4):
        w = np.ones(4) / 3; w[missing] = 0
        weights.append(w); kinds.append("face")
    for a, b in combinations(range(4), 2):
        w = np.zeros(4); w[[a, b]] = .5
        weights.append(w); kinds.append("edge")
    for w in np.eye(4):
        weights.append(w); kinds.append("vertex")
    weights = np.array(weights)
    positions = weights @ vertices
    field = FrozenFEMField(vertices, np.array([[0, 1, 2, 3]]), vertices @ matrix.T + offset,
                           pressure_offset + vertices @ pressure_gradient)
    batch = field.sample_many(positions)
    return dict(field=field, weights=weights, kinds=kinds, positions=positions, batch=batch,
                velocity=positions @ matrix.T + offset, pressure=pressure_offset + positions @ pressure_gradient,
                gradient=matrix, vorticity=np.array([5., 2., 7.]), strain=(matrix + matrix.T) / 2)


def real_node_case(field, boundaries):
    rng = np.random.default_rng(SEED + 2)
    groups = [("fixed", np.array([0, 1, 2, 10, 100, 1000, len(field.points_m) - 1])),
              ("random", rng.choice(len(field.points_m), 500, replace=False))]
    boundary_ids = []
    for role, surface in sorted(boundaries.items()):
        ids = np.asarray(surface.point_data["GlobalNodeID"], dtype=np.int64) - 1
        boundary_ids.extend(ids)
        groups.append((role, rng.choice(ids, min(200, len(ids)), replace=False)))
    interior = np.setdiff1d(np.arange(len(field.points_m)), boundary_ids)
    groups.append(("vessel_interior", rng.choice(interior, 200, replace=False)))
    rows = []
    for group, ids in groups:
        for node in ids:
            sample = field.sample(field.points_m[node])
            # Bounds derive from the selected tetra's floating geometry budget
            # times local nodal variation (not a hand-adjusted fixed threshold).
            cell = sample.tetra_id
            tau = field.geometry.weight_tolerance[cell] if cell >= 0 else 0.
            nodes = field.tetra[cell] if cell >= 0 else np.array([node])
            ubound = 4 * tau * np.max(np.ptp(field.velocity_nodes_m_s[nodes], axis=0)) + 16 * EPS * np.max(np.abs(field.velocity_nodes_m_s[nodes]))
            pbound = 4 * tau * np.ptp(field.pressure_nodes_pa[nodes]) + 16 * EPS * np.max(np.abs(field.pressure_nodes_pa[nodes]))
            rows.append(dict(group=group, node_id=int(node), tetra_id=cell, inside_lumen=sample.inside_lumen,
                             x_m=field.points_m[node, 0], y_m=field.points_m[node, 1], z_m=field.points_m[node, 2],
                             velocity_error_m_s=float(np.max(np.abs(sample.velocity_m_s - field.velocity_nodes_m_s[node]))),
                             pressure_error_pa=float(abs(sample.pressure_pa - field.pressure_nodes_pa[node])),
                             velocity_bound_m_s=float(ubound), pressure_bound_pa=float(pbound)))
    return rows


def shared_faces(field, count=8):
    """Obtain actual shared faces by sorted global IDs, independently of locator."""
    tetra = field.tetra
    faces = np.stack([tetra[:, inds] for inds in combinations(range(4), 3)], axis=1).reshape(-1, 3)
    keys = np.sort(faces, axis=1)
    order = np.lexsort(keys.T[::-1])
    ordered = keys[order]
    duplicates = np.flatnonzero(np.all(ordered[1:] == ordered[:-1], axis=1))
    chosen = np.random.default_rng(SEED + 3).choice(duplicates, count, replace=False)
    result = []
    for k in chosen:
        ids = sorted([int(order[k] // 4), int(order[k + 1] // 4)])
        nodes = ordered[k]
        xyz = field.points_m[nodes]
        center = xyz.mean(axis=0)
        normal = np.cross(xyz[1] - xyz[0], xyz[2] - xyz[0]); normal /= np.linalg.norm(normal)
        centroids = field.points_m[tetra[ids]].mean(axis=1)
        if np.dot(normal, centroids[1] - centroids[0]) < 0:
            normal = -normal
        # Keep the line inside both cells: intersect their barycentric halfspaces.
        limits = []
        for cell, sign in zip(ids, [-1, 1]):
            w = field.geometry.weights(center, cell)
            slope = field.geometry.shape_gradients_m_inv[cell] @ (sign * normal)
            inward = slope < 0
            limits.append(np.min(w[inward] / -slope[inward]))
        distance = .15 * min(limits)
        result.append(dict(ids=ids, nodes=nodes, center=center, normal=normal, distance=distance))
    return result


def shared_face_case(field, faces):
    rows, checks = [], []
    for face_id, face in enumerate(faces):
        left, right = face["ids"]
        center = face["center"]
        # Independent on-face truth is the mean of the THREE face node values.
        uref = field.velocity_nodes_m_s[face["nodes"]].mean(axis=0)
        pref = field.pressure_nodes_pa[face["nodes"]].mean()
        direct_u, direct_p = [], []
        ubounds, pbounds = [], []
        for cell in [left, right]:
            w = field.geometry.weights(center, cell)
            nodes = field.tetra[cell]
            direct_u.append(w @ field.velocity_nodes_m_s[nodes])
            direct_p.append(w @ field.pressure_nodes_pa[nodes])
            tau = field.geometry.weight_tolerance[cell]
            ubounds.append(4 * tau * np.max(np.abs(field.velocity_nodes_m_s[nodes])))
            pbounds.append(4 * tau * np.max(np.abs(field.pressure_nodes_pa[nodes])))
        checks.append(dict(face_id=face_id, left=left, right=right,
                           velocity_trace_error=float(np.max(np.abs(np.array(direct_u) - uref))),
                           pressure_trace_error=float(np.max(np.abs(np.array(direct_p) - pref))),
                           velocity_bound=float(max(ubounds)), pressure_bound=float(max(pbounds)),
                           gradient_jump_s_inv=float(np.linalg.norm(field.gradients_s_inv[right] - field.gradients_s_inv[left]))))
        for fraction in [-1., -.5, -.1, -.01, -1e-4, 0., 1e-4, .01, .1, .5, 1.]:
            distance = fraction * face["distance"]
            position = center + distance * face["normal"]
            s = field.sample(position)
            row = dict(face_id=face_id, signed_distance_m=distance, tetra_id=s.tetra_id,
                       left_id=left, right_id=right, speed_m_s=np.linalg.norm(s.velocity_m_s),
                       pressure_pa=s.pressure_pa, gradient_norm_s_inv=np.linalg.norm(s.velocity_gradient_s_inv),
                       velocity_trace_error_m_s=checks[-1]["velocity_trace_error"],
                       pressure_trace_error_pa=checks[-1]["pressure_trace_error"],
                       gradient_jump_s_inv=checks[-1]["gradient_jump_s_inv"])
            row.update({f"velocity_{axis}_m_s": s.velocity_m_s[i] for i, axis in enumerate("xyz")})
            row.update({f"gradient_{i}{j}_s_inv": s.velocity_gradient_s_inv[i, j] for i in range(3) for j in range(3)})
            rows.append(row)
    return rows, checks


def interior_positions(field, count, seed=SEED + 5):
    """Uniform in tetrahedral volume: volume-weighted cell + Dirichlet(1) weights."""
    rng = np.random.default_rng(seed)
    ids = rng.choice(len(field.tetra), count, p=field.geometry.volume_m3 / field.geometry.volume_m3.sum())
    weights = rng.dirichlet(np.ones(4), count)
    return np.einsum("ni,nij->nj", weights, field.points_m[field.tetra[ids]]), ids, weights


def classification_case(field, boundaries):
    rows = []
    def add(position, kind, expected, role="", owner=-1, offset=0.):
        sample = field.sample(position)
        rows.append(dict(x_m=position[0], y_m=position[1], z_m=position[2], kind=kind, role=role,
                         expected_inside=expected, inside_lumen=sample.inside_lumen, tetra_id=sample.tetra_id,
                         owner_tetra_id=int(owner), offset_m=offset, correct=sample.inside_lumen == expected,
                         all_physical_nan=bool(np.isnan(np.r_[sample.velocity_m_s, sample.pressure_pa,
                             sample.velocity_gradient_s_inv.ravel(), sample.vorticity_s_inv, sample.strain_rate_s_inv.ravel()]).all())))
    positions, ids, _ = interior_positions(field, 120, SEED + 4)
    for p, cell in zip(positions, ids):
        add(p, "inside_random", True, owner=cell)
    for cell in ids[:20]:
        add(field.points_m[field.tetra[cell]].mean(axis=0), "centroid", True, owner=cell)
    for role, surface in sorted(boundaries.items()):
        triangles = surface.faces.reshape(-1, 4)[:, 1:]
        selected = np.random.default_rng(SEED + 4).choice(len(triangles), min(30, len(triangles)), replace=False)
        for face_id in selected:
            xyz = surface.points[triangles[face_id]]
            p = xyz.mean(axis=0)
            cell = int(surface.cell_data["GlobalElementID"][face_id]) - 1
            normal = np.cross(xyz[1] - xyz[0], xyz[2] - xyz[0]); normal /= np.linalg.norm(normal)
            assert np.dot(normal, p - field.points_m[field.tetra[cell]].mean(axis=0)) > 0
            add(p, "boundary", True, role, cell)
            h = np.min(np.linalg.norm(xyz[[1, 2, 0]] - xyz, axis=1))
            # 1e-6 edge length lies well beyond computed floating uncertainty;
            # it is a test offset, never the sampler's inclusion tolerance.
            for ratio in [1e-6, 1e-3]:
                delta = h * ratio
                assert delta > 100 * field.geometry.candidate_padding_m
                add(p + delta * normal, "near_outside", False, role, cell, delta)
            add(p - h * 1e-3 * normal, "near_inside", True, role, cell, -h * 1e-3)
    low, high = field.points_m.min(axis=0), field.points_m.max(axis=0)
    for axis in range(3):
        for side in [-1, 1]:
            p = (low + high) / 2
            p[axis] = low[axis] - (high[axis] - low[axis]) if side < 0 else high[axis] + (high[axis] - low[axis])
            add(p, "bbox_outside", False)
    return rows


def real_flow_case(field, count=400):
    positions, ids, weights = interior_positions(field, count)
    batch = field.sample_many(positions)
    rows = []
    for k, p in enumerate(positions):
        row = dict(x_m=p[0], y_m=p[1], z_m=p[2], tetra_id=int(batch.tetra_id[k]),
                   speed_m_s=float(np.linalg.norm(batch.velocity_m_s[k])), pressure_pa=batch.pressure_pa[k],
                   vorticity_norm_s_inv=float(np.linalg.norm(batch.vorticity_s_inv[k])),
                   strain_norm_s_inv=float(np.linalg.norm(batch.strain_rate_s_inv[k])))
        for prefix, array in [("velocity", batch.velocity_m_s), ("vorticity", batch.vorticity_s_inv)]:
            unit = "m_s" if prefix == "velocity" else "s_inv"
            row.update({f"{prefix}_{axis}_{unit}": array[k, i] for i, axis in enumerate("xyz")})
        for name, array in [("gradient", batch.velocity_gradient_s_inv), ("strain", batch.strain_rate_s_inv)]:
            row.update({f"{name}_{i}{j}_s_inv": array[k, i, j] for i in range(3) for j in range(3)})
        rows.append(row)
    return rows, positions, ids, weights, batch
