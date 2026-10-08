"""Separate prescribed transport, fixed-ID deformation, and graph connectivity."""

import argparse
from pathlib import Path
import h5py
import numpy as np
from mevo_mb.fast_sim.morphodynamics import morphology_series, topology
from .common import BASE, ROOT, js, write, sha


def read(path):
    with h5py.File(path) as h:
        data = {
            key: h[key][:]
            for key in [
                "time_s",
                "positions_m",
                "particle_id",
                "radii_m",
                "active",
                "wall_gap_m",
                "wall_stats",
            ]
        }
        data["flow"] = {key: h["flow"][key][:] for key in h["flow"]}
        data["metrics"] = {key: h["metrics"][key][:] for key in h["metrics"]}
        data["edge_t"] = h["edge_events/time_s"][:]
        data["edge_ids"] = h["edge_events/batch_i_j_linked"][:]
    data["x"] = data["positions_m"][:, 0]
    return data


def series(d):
    x, time = d["x"], d["time_s"]
    n = x.shape[1]
    labels = d["metrics"]["labels"][:, 0]
    morph = morphology_series(x, np.arange(n), labels)
    displacement = morph["fixed_centroid_m"] - morph["fixed_centroid_m"][0]
    distances = np.linalg.norm(x[:, :, None] - x[:, None, :], axis=-1)
    mask = np.eye(n, dtype=bool)
    distances[:, mask] = np.inf
    nearest = distances.min(-1)
    nearest_id = distances.argmin(-1)
    neighbors = distances <= 2.1 * float(d["radii_m"][0])
    changed = np.count_nonzero(neighbors != neighbors[0], axis=(1, 2)) // 2
    u = d["flow"]["velocity_m_s"][:, 0]
    grad = d["flow"]["gradient_s_inv"][:, 0]
    strain = 0.5 * (grad + grad.swapaxes(-1, -2))
    du = np.linalg.norm(u[:, :, None] - u[:, None, :], axis=-1).max(axis=(1, 2))
    out = dict(
        time_s=time,
        **morph,
        centroid_displacement_m=displacement,
        centroid_speed_m_s=np.linalg.norm(np.gradient(morph["fixed_centroid_m"], time, axis=0), axis=1),
        streamwise_displacement_m=displacement[:, 0],
        transverse_displacement_m=np.linalg.norm(displacement[:, 1:], axis=1),
        nearest_neighbor_m=nearest,
        nearest_neighbor_id=nearest_id,
        nearest_neighbor_identity_changes=np.count_nonzero(nearest_id != nearest_id[0], axis=1),
        pair_neighbor_changes_from_initial=changed,
        cluster_count=d["metrics"]["cluster_count"][:, 0],
        largest_size=d["metrics"]["largest_size"][:, 0],
        largest_rg_m=d["metrics"]["largest_rg_m"][:, 0],
        mean_flow_m_s=u.mean(1),
        maximum_pair_flow_difference_m_s=du,
        flow_gradient_s_inv=grad,
        flow_strain_s_inv=strain,
        mean_flow_gradient_s_inv=grad.mean(1),
        mean_flow_strain_s_inv=strain.mean(1),
        mean_flow_strain_frobenius_s_inv=np.linalg.norm(strain, axis=(-1, -2)).mean(1),
    )
    graph, _, _ = topology(x[0], time, labels, d["edge_t"], d["edge_ids"], time[-1], 0.002)
    summary = dict(
        duration_s=float(time[-1] - time[0]),
        fixed_particle_count=n,
        centroid_delta_m=displacement[-1],
        centroid_distance_m=float(np.linalg.norm(displacement[-1])),
        streamwise_delta_m=float(displacement[-1, 0]),
        transverse_delta_m=float(out["transverse_displacement_m"][-1]),
        mean_path_speed_m_s=float(
            np.sum(np.linalg.norm(np.diff(morph["fixed_centroid_m"], axis=0), axis=1)) / (time[-1] - time[0])
        ),
        initial_rg_m=float(morph["fixed_rg_m"][0]),
        final_rg_m=float(morph["fixed_rg_m"][-1]),
        rg_ratio=float(morph["fixed_rg_m"][-1] / morph["fixed_rg_m"][0]),
        initial_eigenvalues_m2=morph["fixed_eigenvalues_m2"][0],
        final_eigenvalues_m2=morph["fixed_eigenvalues_m2"][-1],
        initial_tensor_m2=morph["fixed_tensor_m2"][0],
        final_tensor_m2=morph["fixed_tensor_m2"][-1],
        initial_aspect_ratio=float(morph["fixed_aspect_ratio"][0]),
        final_aspect_ratio=float(morph["fixed_aspect_ratio"][-1]),
        final_aligned_internal_rmsd_m=float(morph["aligned_internal_rmsd_m"][-1]),
        initial_cluster_count=int(out["cluster_count"][0]),
        final_cluster_count=int(out["cluster_count"][-1]),
        initial_largest_size=int(out["largest_size"][0]),
        final_largest_size=int(out["largest_size"][-1]),
        initial_nearest_neighbor_mean_m=float(nearest[0].mean()),
        final_nearest_neighbor_mean_m=float(nearest[-1].mean()),
        final_neighbor_pair_changes=int(changed[-1]),
        final_nearest_neighbor_id_changes=int(out["nearest_neighbor_identity_changes"][-1]),
        mean_flow_at_start_m_s=u[0].mean(0),
        mean_flow_at_end_m_s=u[-1].mean(0),
        maximum_pair_flow_difference_range_m_s=[float(du.min()), float(du.max())],
        mean_flow_gradient_at_start_s_inv=grad[0].mean(0),
        mean_flow_gradient_at_end_s_inv=grad[-1].mean(0),
        mean_strain_norm_range_s_inv=[
            float(out["mean_flow_strain_frobenius_s_inv"].min()),
            float(out["mean_flow_strain_frobenius_s_inv"].max()),
        ],
        topology_counts=graph["counts"],
        active_all=bool(d["active"].all()),
        minimum_saved_pair_gap_m=float(d["metrics"]["minimum_gap_m"].min()),
        minimum_saved_wall_gap_m=float(d["wall_gap_m"].min()),
    )
    return out, summary, graph


def subset(d, start, stop=None):
    indices = np.arange(len(d["time_s"]))[slice(start, stop)]
    shift = float(d["time_s"][indices[0]])
    out = {
        key: (
            value[indices]
            if isinstance(value, np.ndarray) and key not in ["particle_id", "radii_m", "edge_t", "edge_ids"]
            else value
        )
        for key, value in d.items()
    }
    out["time_s"] = d["time_s"][indices] - shift
    for group in ["metrics", "flow"]:
        out[group] = {key: value[indices] for key, value in d[group].items()}
    edges = (d["edge_t"] > shift + 1e-14) & (d["edge_t"] <= d["time_s"][indices[-1]] + 1e-14)
    out["edge_t"] = d["edge_t"][edges] - shift
    out["edge_ids"] = d["edge_ids"][edges]
    return out


def main(label, f4_reference):
    run = BASE / label
    manifest = js(run / "RESULT.json")
    assert manifest["status"] == "COMPLETE"
    d = read(run / "main/trajectory.h5")
    b = read(run / "branch_B/trajectory.h5")
    idx = manifest["switch_index"]
    datasets = {"main": d, "ON": subset(d, 0, idx + 1), "A": subset(d, idx), "B": b}
    np.testing.assert_array_equal(d["x"][idx], b["x"][0])
    rows = {}
    for name, data in datasets.items():
        arrays, rows[name], graph = series(data)
        np.savez_compressed(run / (name + "_ANALYSIS.npz"), **arrays)
        write(run / (name + "_TOPOLOGY.json"), graph)
    reference = Path(f4_reference)
    with h5py.File(reference / "trajectory.h5") as h:
        np.testing.assert_array_equal(h["positions_m"][0, 0], d["x"][0])
        frame = int(np.argmin(abs(h["time_s"][:] - 0.3)))
        assert abs(h["time_s"][frame] - 0.3) < 1e-14
        noflow = dict(
            time_s=float(h["time_s"][frame]),
            largest_size=int(h["metrics/largest_size"][frame, 0]),
            cluster_count=int(h["metrics/cluster_count"][frame, 0]),
            rg_m=float(h["metrics/global_rg_m"][frame, 0]),
        )
        np.savez_compressed(
            run / "F4_ON_NO_FLOW_REFERENCE.npz",
            time_s=h["time_s"][: frame + 1],
            positions_m=h["positions_m"][: frame + 1],
            largest_size=h["metrics/largest_size"][: frame + 1, 0],
            cluster_count=h["metrics/cluster_count"][: frame + 1, 0],
            rg_m=h["metrics/global_rg_m"][: frame + 1, 0],
        )
    refcfg = js(reference / "CONFIG.json")
    cfg = js(run / "main/CONFIG.json")
    for key in ["amplitude_Pa", "phase_rad", "enabled"]:
        assert refcfg["events"][0][key] == cfg["events"][0][key]
    comparison = dict(
        status="COMPLETE",
        switch_s=0.3,
        same_state_bitwise=True,
        A=rows["A"],
        B=rows["B"],
        ON_with_flow=rows["ON"],
        ON_without_flow_F4=noflow,
        F4_reference=dict(
            path=str(reference), sha256=sha(reference / "trajectory.h5"), same_initial_and_PZT=True
        ),
        interpretations=dict(
            A="Uninterrupted main suffix: all acoustics zero, prescribed flow retained",
            B="Same physical state, all acoustics and background flow zero; residual numerical contact relaxation may remain",
            topology="Distance-based connected components (gap <= 0.1R), persistent after 2 ms; split/merge does not mean bubble fragmentation/fusion",
            strain="Gradient of prescribed CFD P1 flow, distinct from measured bubble deformation; no wall hydrodynamic lubrication",
        ),
    )
    write(BASE / "ON_OFF_FLOW_COMPARISON.json", comparison)
    write(
        BASE / "CLUSTER_TRANSPORT_METRICS.json",
        dict(
            status="COMPLETE",
            runs=rows,
            full_series={name: str((run / (name + "_ANALYSIS.npz")).relative_to(ROOT)) for name in datasets},
            fixed_ID_set=list(range(32)),
            flow_gradient="Per bubble piecewise-constant P1 gradients; discontinuous between tetrahedra",
            definitions="Rg and eigenvalues use all 32 fixed IDs; dynamic clusters and per-cluster tensors retained separately in HDF5 metrics",
            frames=[0, idx, len(d["time_s"]) - 1],
            run_path=str(run),
            main_sha256=sha(run / "main/trajectory.h5"),
        ),
    )
    print(
        {
            name: {
                key: row[key]
                for key in [
                    "centroid_distance_m",
                    "rg_ratio",
                    "final_largest_size",
                    "final_cluster_count",
                    "topology_counts",
                ]
            }
            for name, row in rows.items()
        }
    )


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--f4-reference", required=True)
    a = ap.parse_args()
    main(a.label, a.f4_reference)
