"""Float64 batched postprocessing; positions and physical events are read-only.

All definitions match metrics.py/morphodynamics.py. Connected-component labels
are assigned in increasing minimum particle ID, exactly as the SciPy reference.
No physical integration, output decimation, or approximate graph is used.
"""

from time import perf_counter
import numpy as np
from scipy.sparse.csgraph import connected_components
from .metrics import DEFINITIONS
from .morphodynamics import parts, rigid_fit


def components(graph):
    """Batched breadth-first components, bounded memory via <=512 frames."""
    shape = graph.shape[:-2]
    n = graph.shape[-1]
    flat = graph.reshape(-1, n, n)
    result = np.full((len(flat), n), -1, dtype=np.int32)
    for begin in range(0, len(flat), 512):
        g = flat[begin : begin + 512]
        labels = result[begin : begin + 512]
        for c in range(n):
            rows = np.flatnonzero(np.any(labels < 0, axis=1))
            if not len(rows):
                break
            sub = g[rows]
            visited = np.zeros((len(rows), n), bool)
            visited[np.arange(len(rows)), (labels[rows] < 0).argmax(1)] = True
            while True:
                grown = visited | np.any(sub & visited[:, None, :], axis=2)
                if np.array_equal(grown, visited):
                    break
                visited = grown
            labels[rows] = np.where(visited, c, labels[rows])
    return result.reshape(*shape, n)


def stacked_shape(x):
    center = x.mean(axis=-2)
    rel = x - center[..., None, :]
    tensor = np.swapaxes(rel, -1, -2) @ rel / x.shape[-2]
    e, v = np.linalg.eigh(tensor)
    return center, rel, tensor, e[..., ::-1], v[..., :, ::-1]


def trajectory_metrics(positions, radius, times, timing=None):
    timing = {} if timing is None else timing
    start = perf_counter()
    x = np.asarray(positions)
    if x.ndim != 4 or x.shape[-1] != 3 or x.shape[2] < 1 or not np.isfinite(x).all():
        raise ValueError("Finite nonempty [time,batch,particle,xyz] positions required")
    nt, nb, n, _ = x.shape
    xx = x.reshape(-1, n, 3)
    frames = len(xx)
    timing["identity_validity_s"] = perf_counter() - start
    start = perf_counter()
    center = xx.mean(axis=1)
    rel = xx - center[:, None, :]
    rg_global = np.sqrt(np.mean(np.sum(rel * rel, axis=2), axis=1))
    timing["centroid_rg_s"] = perf_counter() - start
    start = perf_counter()
    tensor = rel.swapaxes(1, 2) @ rel / n
    eigen_global = np.linalg.eigvalsh(tensor)[:, ::-1]
    timing["global_tensor_eigen_s"] = perf_counter() - start
    start = perf_counter()
    # Same three-coordinate summation order as scipy.spatial.distance.pdist.
    distance = (xx[:, :, None, 0] - xx[:, None, :, 0]) ** 2
    distance += (xx[:, :, None, 1] - xx[:, None, :, 1]) ** 2
    distance += (xx[:, :, None, 2] - xx[:, None, :, 2]) ** 2
    np.sqrt(distance, out=distance)
    distance[:, np.arange(n), np.arange(n)] = np.inf
    gap = distance - 2 * radius
    graph = gap <= DEFINITIONS["cluster_gap_over_radius"] * radius
    labels = components(graph)
    sizes = np.zeros((frames, n), dtype=np.int64)
    np.add.at(sizes, (np.repeat(np.arange(frames), n), labels.ravel()), 1)
    count = np.sum(sizes > 0, axis=1)
    timing["pair_distances_components_s"] = perf_counter() - start
    start = perf_counter()
    centroids = np.zeros((frames, n, 3))
    tensors = np.zeros((frames, n, 3, 3))
    eigen = np.zeros((frames, n, 3))
    axes = np.zeros((frames, n, 3, 3))
    rg = np.zeros((frames, n))
    identifiable = np.zeros((frames, n), bool)
    for size in range(1, n + 1):
        ff, cc = np.nonzero(sizes == size)
        if not len(ff):
            continue
        members = np.nonzero(labels[ff] == cc[:, None])[1].reshape(-1, size)
        sx = xx[ff[:, None], members]
        ce, _, te, ee, vv = stacked_shape(sx)
        centroids[ff, cc] = ce
        tensors[ff, cc] = te
        eigen[ff, cc] = ee
        axes[ff, cc] = vv
        rg[ff, cc] = np.sqrt(np.trace(te, axis1=-2, axis2=-1))
        identifiable[ff, cc] = (ee[:, 0] > 0) & ((ee[:, 0] - ee[:, 1]) > 1e-3 * ee[:, 0])
    largest = sizes.argmax(1)
    index = np.arange(frames)
    arrays = dict(
        cluster_count=count,
        labels=labels,
        cluster_sizes=sizes,
        centroids_m=centroids,
        shape_tensors_m2=tensors,
        eigenvalues_m2=eigen,
        principal_axes=axes,
        axis_identifiable=identifiable,
        rg_m=rg,
        largest_size=sizes[index, largest],
        largest_centroid_m=centroids[index, largest],
        largest_rg_m=rg[index, largest],
        largest_eigenvalues_m2=eigen[index, largest],
        global_centroid_m=center,
        global_rg_m=rg_global,
        global_shape_tensor_m2=tensor,
        global_eigenvalues_m2=eigen_global,
        mean_nearest_neighbor_m=distance.min(2).mean(1),
        minimum_gap_m=gap.min(axis=(1, 2)),
        near_contact_pair_count=np.count_nonzero(np.triu(gap <= 0.02 * radius, 1), axis=(1, 2)),
        hard_contact_pair_count=np.count_nonzero(np.triu(gap <= 0.001 * radius, 1), axis=(1, 2)),
    )
    arrays = {k: v.reshape(nt, nb, *v.shape[1:]) for k, v in arrays.items()}
    timing["component_shape_pair_statistics_s"] = perf_counter() - start
    start = perf_counter()
    events = []
    for ti, b in np.argwhere(np.any(arrays["labels"][1:] != arrays["labels"][:-1], axis=2)):
        t = ti + 1
        old, new = arrays["labels"][t - 1, b], arrays["labels"][t, b]
        for kind, left, right in [("merge", old, new), ("split", new, old)]:
            for c in np.unique(right):
                members = np.flatnonzero(right == c)
                if len(np.unique(left[members])) > 1:
                    events.append(
                        dict(
                            kind=kind,
                            batch=int(b),
                            observed_time_s=float(times[t]),
                            bracket_s=[float(times[t - 1]), float(times[t])],
                            particle_ids=members.tolist(),
                            scope="sampled membership change; event time bounded by saved frames",
                        )
                    )
    timing["sampled_lineage_s"] = perf_counter() - start
    return arrays, events


def morphology_batch(x, labels, metrics, axis_gap=0.1, timing=None):
    """Original all-fixed-ID morphology, with the same determinant correction."""
    timing = {} if timing is None else timing
    start = perf_counter()
    center, rel, g, e, v = stacked_shape(x)
    a = x[0] - x[0].mean(axis=1)[:, None, :]
    e0, v0 = e[0], v[0]
    timing["fixed_tensor_eigen_s"] = perf_counter() - start
    start = perf_counter()
    u, singular, vt = np.linalg.svd(rel.swapaxes(-1, -2) @ a)
    s = np.broadcast_to(np.eye(3), u.shape).copy()
    s[..., 2, 2] = np.linalg.det(u @ vt)
    rotation = u @ s @ vt
    residual = rel @ rotation - a
    rmsd = np.sqrt(np.mean(np.sum(residual * residual, axis=-1), axis=-1))
    # Singular Kabsch solutions are not unique: exactly follow the scalar reference.
    degenerate = (singular[..., 0] == 0) | (singular[..., -1] <= 1e-12 * singular[..., 0])
    for t, b in np.argwhere(degenerate):
        rotation[t, b], rmsd[t, b] = rigid_fit(x[0, b], x[t, b])
    timing["kabsch_s"] = perf_counter() - start
    start = perf_counter()
    with np.errstate(divide="ignore", invalid="ignore"):
        valid = (e[..., 0] > 0) & ((e[..., 0] - e[..., 1]) / e[..., 0] >= axis_gap)
        valid &= (e0[..., 0] - e0[..., 1]) / e0[..., 0] >= axis_gap
    angle = np.degrees(np.arccos(np.clip(np.abs(np.sum(v[..., :, 0] * v0[..., :, 0], axis=-1)), 0, 1)))
    angle = np.where(valid, angle, np.nan)
    chosen = metrics["cluster_sizes"].argmax(-1)
    tt, bb = np.indices(chosen.shape)
    de = metrics["eigenvalues_m2"][tt, bb, chosen]
    result = dict(
        fixed_centroid_m=center,
        fixed_tensor_m2=g,
        fixed_eigenvalues_m2=e,
        fixed_rg_m=np.sqrt(e.sum(-1)),
        fixed_aspect_ratio=np.sqrt(e[..., 0] / np.maximum(e[..., -1], 1e-300)),
        fixed_principal_axes=v,
        principal_axis_identifiable=valid,
        principal_axis_angle_deg=angle,
        kabsch_rotation_deg=np.degrees(
            np.arccos(np.clip((np.trace(rotation, axis1=-2, axis2=-1) - 1) / 2, -1, 1))
        ),
        aligned_internal_rmsd_m=rmsd,
        rotation_matrix=rotation,
        dynamic_centroid_m=metrics["centroids_m"][tt, bb, chosen],
        dynamic_tensor_m2=metrics["shape_tensors_m2"][tt, bb, chosen],
        dynamic_eigenvalues_m2=de,
        dynamic_rg_m=np.sqrt(de.sum(-1)),
        dynamic_member_mask=labels == chosen[..., None],
    )
    timing["aspect_axis_dynamic_s"] = perf_counter() - start
    timing["degenerate_kabsch_frames"] = int(degenerate.sum())
    return result


def topology_fast(x0, frame_times, frame_labels, edge_times, edge_ids, end_time, persistence_s=0.002):
    """Same accepted-step timeline/dwell; skip set construction only when labels equal."""
    graph = np.linalg.norm(x0[:, None] - x0[None, :], axis=-1) <= 2.1 * 1.78e-6
    lab = connected_components(graph, directed=False)[1]
    timeline_t, timeline_labels, events = [0.0], [lab.copy()], []
    cursor = 0
    for frame, t in enumerate(frame_times):
        while cursor < len(edge_times) and edge_times[cursor] <= t + 1e-14:
            et = float(edge_times[cursor])
            stop = int(np.searchsorted(edge_times, et, side="right"))
            ids = edge_ids[cursor:stop]
            assert np.all(ids[:, 0] == 0)
            old = lab
            # Preserve sequential last-write semantics even for duplicate pairs.
            for _, i, j, on in ids:
                graph[i, j] = graph[j, i] = bool(on)
            lab = connected_components(graph, directed=False)[1]
            if not np.array_equal(old, lab):
                for kind, left, right in [("merge", old, lab), ("split", lab, old)]:
                    for c in np.unique(right):
                        members = np.flatnonzero(right == c)
                        others = np.unique(left[members])
                        if len(others) > 1:
                            events.append(
                                dict(
                                    kind=kind,
                                    time_s=et,
                                    particle_ids=members.tolist(),
                                    groups=[members[left[members] == other].tolist() for other in others],
                                )
                            )
            timeline_t.append(et)
            timeline_labels.append(lab.copy())
            cursor = stop
        assert np.array_equal(lab, frame_labels[frame]) or parts(lab) == parts(frame_labels[frame]), (
            "Topology/frame mismatch",
            frame,
        )
    assert cursor == len(edge_times)
    times = np.asarray(timeline_t)
    labels = np.asarray(timeline_labels)
    for event in events:
        idx = int(np.searchsorted(times, event["time_s"], side="right")) - 1
        ids = event["particle_ids"]
        invalid = None
        for k in range(idx, len(times)):
            current = labels[k]
            if event["kind"] == "merge":
                valid = len(np.unique(current[ids])) == 1
            else:
                a = [set(current[g]) for g in event["groups"]]
                valid = all(not a[i].intersection(a[j]) for i in range(len(a)) for j in range(i + 1, len(a)))
            if not valid:
                invalid = float(times[k])
                break
        observed = (invalid if invalid is not None else end_time) - event["time_s"]
        event["observed_duration_s"] = float(observed)
        event["persistence_class"] = (
            "persistent"
            if observed >= persistence_s - 1e-14
            else "transient"
            if invalid is not None
            else "right_censored"
        )
    return (
        dict(
            events=events,
            raw_edge_changes=len(edge_times),
            verified_frames=len(frame_times),
            persistence_s=persistence_s,
            counts={
                f"{kind}_{p}": sum(e["kind"] == kind and e["persistence_class"] == p for e in events)
                for kind in ["merge", "split"]
                for p in ["persistent", "transient", "right_censored"]
            },
            scope="Discrete accepted-step graph changes; persistence never replaces raw events; not physical bubble breakup/fusion",
        ),
        times,
        labels,
    )


def process(x, radius, times, edge_ids, edge_times):
    started = perf_counter()
    timing = {}
    metrics, lineage = trajectory_metrics(x, radius, times, timing)
    morph = morphology_batch(x, metrics["labels"], metrics, timing=timing)
    start = perf_counter()
    top = []
    for i in range(x.shape[1]):
        take = (edge_ids[:, 0] == i) & (edge_times <= times[-1] + 1e-14)
        ei = edge_ids[take].copy()
        ei[:, 0] = 0
        tp, _, _ = topology_fast(x[0, i], times, metrics["labels"][:, i], edge_times[take], ei, times[-1])
        top.append(tp)
    timing["accepted_step_topology_s"] = perf_counter() - start
    timing["total_s"] = perf_counter() - started
    return metrics, lineage, morph, top, timing
