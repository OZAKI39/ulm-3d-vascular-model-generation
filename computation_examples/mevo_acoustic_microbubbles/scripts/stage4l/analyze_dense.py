"""Same-ID dense-cluster metrics on identical physical-time sampling grids."""

import argparse
import json
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist

from scripts.stage4l.dense_probe import BASE, write, sha


def edges(x, cutoff):
    return {tuple(p) for p in cKDTree(x).query_pairs(cutoff, output_type="ndarray")}


def largest(x, cutoff):
    pairs = np.array(sorted(edges(x, cutoff)), int).reshape(-1, 2)
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(x), len(x))).tocsr()
    _, labels = connected_components(graph, directed=False)
    return np.flatnonzero(labels == np.argmax(np.bincount(labels)))


def lifetime(graphs, times):
    live = {edge: times[0] for edge in graphs[0]}
    initial = set(live)
    rows = []
    for k in range(1, len(times)):
        for edge in set(live) - graphs[k]:
            rows.append(
                dict(
                    pair=list(edge),
                    duration_s=float(times[k] - live.pop(edge)),
                    left_censored=edge in initial,
                    right_censored=False,
                )
            )
            initial.discard(edge)
        for edge in graphs[k] - set(live):
            live[edge] = times[k]
    for edge, start in live.items():
        rows.append(
            dict(
                pair=list(edge),
                duration_s=float(times[-1] - start),
                left_censored=edge in initial,
                right_censored=True,
            )
        )
    return rows


def analyze(run, cadence=1e-5, member=None, tag=""):
    summary = json.loads((run / "SUMMARY.json").read_text())
    with h5py.File(run / "trajectories.h5", "r") as h:
        t = h["time_s"][:]
        x = h["positions_m"][:]
        active = h["active"][0]
        r = float(h["radii_m"][0])
        ids = h["particle_id"][:]
    if member is None:
        member = np.flatnonzero(active)[largest(x[0, active], 2.1 * r)]
    if len(member) < 2:
        raise ValueError("No multi-particle cluster at the selected checkpoint")
    if len(t) < 2:
        return dict(status="NO_ACCEPTED_STEPS")
    samples = np.arange(t[0], t[-1] + cadence * 0.01, cadence)
    # Sampling for statistics only; raw accepted positions are retained unchanged.
    sampled = np.empty((len(samples), len(member), 3))
    for j, index in enumerate(member):
        for a in range(3):
            sampled[:, j, a] = np.interp(samples, t, x[:, index, a])
    centered = sampled - sampled.mean(axis=1, keepdims=True)
    reference = centered[0]
    rms = np.sqrt(np.mean(np.sum((centered - reference) ** 2, axis=2), axis=1))
    rigid = []
    for frame in centered:
        u, _, vt = np.linalg.svd(reference.T @ frame)
        fix = np.diag([1.0, 1.0, np.linalg.det(u @ vt)])
        q = u @ fix @ vt
        rigid.append(float(np.sqrt(np.mean(np.sum((reference @ q - frame) ** 2, axis=1)))))
    rg = np.sqrt(np.mean(np.sum(centered**2, axis=2), axis=1))
    graph = [edges(frame, 2.1 * r) for frame in sampled]
    contact_graph = [edges(frame, 2 * r + 1e-11) for frame in sampled]
    turnover = [len(a.symmetric_difference(b)) / max(1, len(a | b)) for a, b in zip(graph[:-1], graph[1:])]
    swapped = [len(a.symmetric_difference(b)) for a, b in zip(graph[:-1], graph[1:])]
    edge_life = lifetime(graph, samples)
    contact_life = lifetime(contact_graph, samples)
    for row in edge_life + contact_life:
        row["pair"] = [int(ids[member[k]]) for k in row["pair"]]
    distances = np.stack([np.quantile(pdist(frame), [0, 0.1, 0.25, 0.5, 0.75, 0.9, 1]) for frame in sampled])
    cluster_sizes = [len(largest(frame, 2.1 * r)) for frame in sampled]
    key = f"{tag}{cadence * 1e6:g}us"
    np.savez_compressed(
        run / f"cluster_metrics_{key}.npz",
        time_s=samples,
        particle_ids=ids[member],
        centered_rms_m=rms,
        rigid_fit_residual_m=rigid,
        radius_of_gyration_m=rg,
        pair_distance_quantiles_m=distances,
        largest_component_size=cluster_sizes,
        turnover_fraction=turnover,
    )
    write(
        run / f"edge_lifetimes_{key}.json",
        dict(
            near_contact=edge_life,
            hard_contact=contact_life,
            semantics="Observed, censored lifetimes within this probe; not uncensored material lifetimes",
        ),
    )
    result = dict(
        status=summary["status"],
        sampling_cadence_s=cadence,
        initial_cluster_size=len(member),
        same_ID_particle_ids=ids[member].tolist(),
        same_ID_centroid_removed_RMS_end_m=float(rms[-1]),
        same_ID_centroid_removed_RMS_mean_m=float(np.mean(rms)),
        rigid_fit_residual_end_m=rigid[-1],
        rigid_fit_residual_mean_m=float(np.mean(rigid)),
        neighbor_turnover_mean=float(np.mean(turnover)),
        edge_changes_per_s=float(sum(swapped) / (samples[-1] - samples[0])),
        near_contact_observed_edge_lifetime_mean_s=float(np.mean([v["duration_s"] for v in edge_life])),
        hard_contact_observed_edge_lifetime_mean_s=float(np.mean([v["duration_s"] for v in contact_life]))
        if contact_life
        else None,
        radius_of_gyration_initial_m=float(rg[0]),
        radius_of_gyration_final_m=float(rg[-1]),
        gyration_relative_change=float(rg[-1] / rg[0] - 1),
        largest_component_final_size=cluster_sizes[-1],
        pair_distance_quantiles_initial_m=distances[0].tolist(),
        pair_distance_quantiles_final_m=distances[-1].tolist(),
        interpolation="Linear interpolation of raw accepted states for statistics only; same physical cadence in both arms",
    )
    write(run / f"CLUSTER_METRICS_{key}.json", result)
    return result


def compare(legacy, lub, out):
    a = json.loads((legacy / "SUMMARY.json").read_text())
    b = json.loads((lub / "SUMMARY.json").read_text())
    metrics = {}
    for cadence in [1e-5, 2e-6]:
        aa, bb = analyze(legacy, cadence), analyze(lub, cadence)
        keys = [
            "same_ID_centroid_removed_RMS_end_m",
            "rigid_fit_residual_end_m",
            "neighbor_turnover_mean",
            "edge_changes_per_s",
            "near_contact_observed_edge_lifetime_mean_s",
        ]
        metrics[f"{cadence * 1e6:g}us"] = dict(
            legacy=aa,
            lubrication=bb,
            ratios_lub_over_legacy={k: bb[k] / aa[k] if aa.get(k) else None for k in keys if k in bb},
        )
    with (
        np.load(legacy / "INITIAL_FORCES_AND_VELOCITIES.npz") as aa,
        np.load(lub / "INITIAL_FORCES_AND_VELOCITIES.npz") as bb,
    ):
        identical = {
            k: bool(np.array_equal(aa[k], bb[k]))
            for k in [
                "primary_force_N",
                "secondary_force_N",
                "flow_velocity_m_s",
                "local_pressure_complex_Pa",
                "volume_response_complex_m3",
            ]
        }
    with np.load(legacy / "INPUT_SNAPSHOT.npz") as aa, np.load(lub / "INPUT_SNAPSHOT.npz") as bb:
        snapshots_identical = {k: bool(np.array_equal(aa[k], bb[k])) for k in aa.files}
        initial = aa["positions_m"][aa["active"]]
        active_indices = np.flatnonzero(aa["active"])
        radius = float(aa["radii_m"][0])
    pairs = np.array(sorted(edges(initial, 2.1 * radius)), int).reshape(-1, 2)
    graph = coo_matrix(
        (np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(len(initial), len(initial))
    ).tocsr()
    _, labels = connected_components(graph, directed=False)
    clusters = sorted(
        [
            active_indices[labels == label]
            for label in np.unique(labels)
            if np.count_nonzero(labels == label) >= 3
        ],
        key=len,
        reverse=True,
    )
    all_clusters = []
    for k, member in enumerate(clusters):
        tag = f"initial_cluster{k}_"
        all_clusters.append(
            dict(legacy=analyze(legacy, 1e-5, member, tag), lubrication=analyze(lub, 1e-5, member, tag))
        )
    steps = [json.loads(line) for line in (lub / "steps.jsonl").read_text().splitlines()]
    weights = np.array([row["dt_s"] for row in steps])
    diagnostics = {}
    for k in [
        "LUBRICATION_ACTIVE",
        "SUBFLOOR",
        "HARD_CONTACT",
        "relative_residual",
        "contact_kkt_residual",
        "secondary_force_max_N",
    ]:
        vals = np.array([row.get(k, 0) for row in steps])
        diagnostics[k] = dict(
            max=float(vals.max(initial=0)), time_weighted_mean=float(np.average(vals, weights=weights))
        )
    result = dict(
        legacy=str(legacy),
        lubrication=str(lub),
        case=a["case"],
        legacy_summary=a,
        lubrication_summary=b,
        metrics=metrics,
        all_initial_clusters_size_at_least_3=all_clusters,
        initial_snapshot_arrays_identical=snapshots_identical,
        initial_inputs_and_forces_bitwise_identical=identical,
        diagnostics=diagnostics,
        analysis_sha256=sha(Path(__file__)),
        gate_status="REQUIRES_SCIENTIFIC_REVIEW",
        overall_stage4_acceptance="NOT_EVALUATED",
    )
    write(out, result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy", type=Path, required=True)
    parser.add_argument("--lubrication", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.output.is_relative_to(BASE):
        raise ValueError("Output must be under lubrication_upgrade")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    compare(args.legacy, args.lubrication, args.output)


if __name__ == "__main__":
    main()
