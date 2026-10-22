"""Topology dwell, maturity, and separation of rigid motion from fixed-ID deformation."""

import numpy as np
from scipy.sparse.csgraph import connected_components


def parts(labels):
    return {frozenset(np.flatnonzero(labels == k)) for k in np.unique(labels)}


def topology(x0, frame_times, frame_labels, edge_times, edge_ids, end_time, persistence_s=0.002):
    radius = 1.78e-6
    graph = np.linalg.norm(x0[:, None] - x0[None, :], axis=-1) <= 2.1 * radius
    lab = connected_components(graph, directed=False)[1]
    timeline_t = [0.0]
    timeline_labels = [lab.copy()]
    events = []
    cursor = 0
    for frame, t in enumerate(frame_times):
        while cursor < len(edge_times) and edge_times[cursor] <= t + 1e-14:
            et = float(edge_times[cursor])
            stop = int(np.searchsorted(edge_times, et, side="right"))
            old = lab
            for batch, i, j, on in edge_ids[cursor:stop]:
                assert batch == 0
                graph[i, j] = graph[j, i] = bool(on)
            lab = connected_components(graph, directed=False)[1]
            for kind, left, right in [("merge", old, lab), ("split", lab, old)]:
                for c in np.unique(right):
                    ids = np.flatnonzero(right == c)
                    others = np.unique(left[ids])
                    if len(others) > 1:
                        events.append(
                            dict(
                                kind=kind,
                                time_s=et,
                                particle_ids=ids.tolist(),
                                groups=[ids[left[ids] == other].tolist() for other in others],
                            )
                        )
            timeline_t.append(et)
            timeline_labels.append(lab.copy())
            cursor = stop
        assert parts(lab) == parts(frame_labels[frame]), ("Topology/frame mismatch", frame)
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


def mature_index(data, history_times, history_labels, criterion, radius):
    times = data["time_s"]
    labels = data["labels"]
    rg = data["largest_rg_m"]
    n = labels.shape[1]
    need = int(np.ceil(criterion["minimum_fraction"] * n))
    window = criterion["stable_window_s"]
    for k, t in enumerate(times):
        if t < window:
            continue
        sizes = np.bincount(labels[k])
        c = int(sizes.argmax())
        ids = np.flatnonzero(labels[k] == c)
        if len(ids) < need:
            continue
        begin = int(np.searchsorted(times, t - window - 1e-14))
        a = max(0, int(np.searchsorted(history_times, t - window, side="right")) - 1)
        b = int(np.searchsorted(history_times, t + 1e-14, side="right"))
        if any(set(np.flatnonzero(lab == lab[ids[0]])) != set(ids) for lab in history_labels[a:b]):
            continue
        values = rg[begin : k + 1]
        if (values.max() - values.min()) / max(float(values.mean()), 1e-300) > criterion["rg_relative_range"]:
            continue
        if data["quality_samples"][k, :, 3].max() > radius * 0.01:
            continue
        return k, ids
    return None, None


def rigid_fit(reference, current):
    a = reference - reference.mean(0)
    b = current - current.mean(0)
    u, _, vt = np.linalg.svd(b.T @ a)
    s = np.eye(3)
    s[2, 2] = np.linalg.det(u @ vt)
    rotation = u @ s @ vt
    residual = b @ rotation - a
    return rotation, float(np.sqrt(np.mean(np.sum(residual * residual, axis=1))))


def shape(x):
    center = x.mean(0)
    a = x - center
    g = a.T @ a / len(x)
    e, v = np.linalg.eigh(g)
    e = e[::-1]
    v = v[:, ::-1]
    return center, g, e, v


def morphology_series(positions, initial_ids, labels, axis_gap=0.1):
    reference = positions[0, initial_ids]
    _, g0, e0, v0 = shape(reference)
    records = []
    for x, lab in zip(positions, labels):
        fixed = x[initial_ids]
        center, g, e, v = shape(fixed)
        rotation, rmsd = rigid_fit(reference, fixed)
        chosen = np.bincount(lab[initial_ids]).argmax()
        dynamic_ids = np.flatnonzero(lab == chosen)
        dc, dg, de, dv = shape(x[dynamic_ids])
        dynamic_mask = np.isin(np.arange(len(x)), dynamic_ids)
        valid = e[0] > 0 and (e[0] - e[1]) / e[0] >= axis_gap and (e0[0] - e0[1]) / e0[0] >= axis_gap
        records.append(
            dict(
                fixed_centroid_m=center,
                fixed_tensor_m2=g,
                fixed_eigenvalues_m2=e,
                fixed_rg_m=np.sqrt(e.sum()),
                fixed_aspect_ratio=np.sqrt(e[0] / max(e[-1], 1e-300)),
                fixed_principal_axes=v,
                principal_axis_identifiable=valid,
                principal_axis_angle_deg=np.degrees(np.arccos(np.clip(abs(v[:, 0] @ v0[:, 0]), 0, 1)))
                if valid
                else np.nan,
                kabsch_rotation_deg=np.degrees(np.arccos(np.clip((np.trace(rotation) - 1) / 2, -1, 1))),
                aligned_internal_rmsd_m=rmsd,
                rotation_matrix=rotation,
                dynamic_centroid_m=dc,
                dynamic_tensor_m2=dg,
                dynamic_eigenvalues_m2=de,
                dynamic_rg_m=np.sqrt(de.sum()),
                dynamic_member_mask=dynamic_mask,
            )
        )
    return {k: np.asarray([r[k] for r in records]) for k in records[0]}


def contrast(reference, current, ids):
    a = reference[ids]
    b = current[ids]
    ca, ga, ea, _ = shape(a)
    cb, gb, eb, _ = shape(b)
    rotation, rms = rigid_fit(a, b)
    return dict(
        centroid_distance_m=float(np.linalg.norm(cb - ca)),
        centroid_delta_m=(cb - ca).tolist(),
        rg_relative=float(abs(np.sqrt(np.trace(gb) / np.trace(ga)) - 1)),
        eigenvalues_relative=float(np.linalg.norm(eb - ea) / np.linalg.norm(ea)),
        tensor_relative=float(np.linalg.norm(gb - ga) / np.linalg.norm(ga)),
        aligned_internal_rmsd_m=rms,
        rigid_rotation_deg=float(np.degrees(np.arccos(np.clip((np.trace(rotation) - 1) / 2, -1, 1)))),
    )
