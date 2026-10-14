"""Fixed-ID analysis; model contrasts are distinct from numerical sensitivity."""

import numpy as np
import h5py
from scipy.sparse.csgraph import connected_components
from mevo_mb.fast_sim.morphodynamics import rigid_fit, topology
from .common import BASE, OLD, CONFIG, read, write


def load(path):
    with h5py.File(path / "trajectory.h5") as h:
        x = h["positions_m"][:]
        t = h["time_s"][:]
        ids = h["particle_id"][:]
    if x.ndim == 4:
        x = x[:, 0]
    np.testing.assert_array_equal(ids, np.arange(32))
    return t, x


def metrics(t, x):
    center = x.mean(1)
    q = x - center[:, None]
    tensor = np.einsum("tni,tnj->tij", q, q) / 32
    eig, axes = np.linalg.eigh(tensor)
    identifiable = (eig[:, -1] - eig[:, -2]) > 1e-3 * eig[:, -1]
    return dict(
        time_s=t,
        centroid_m=center,
        tensor_m2=tensor,
        rg_m=np.sqrt(eig.sum(1)),
        aspect=np.sqrt(eig[:, -1] / eig[:, 0]),
        directional_xy=np.sqrt(tensor[:, 0, 0] / tensor[:, 1, 1]),
        principal_axis=axes[:, :, -1],
        axis_identifiable=identifiable,
        angle_to_x_deg=np.where(
            identifiable, np.degrees(np.arccos(np.clip(abs(axes[:, 0, -1]), 0, 1))), np.nan
        ),
        centered_positions_m=q,
        kabsch_from_initial_m=np.array([rigid_fit(x[0], row)[1] for row in x]),
    )


def compare(a, b, kind):
    ta, xa = load(a)
    tb, xb = load(b)
    ia = []
    ib = []
    for i, t in enumerate(ta):
        j = np.searchsorted(tb, t - 1e-13)
        if j < len(tb) and abs(tb[j] - t) < 1e-13:
            ia.append(i)
            ib.append(j)
    xa, xb = xa[ia], xb[ib]
    times = ta[ia]
    ma, mb = metrics(times, xa), metrics(times, xb)
    arrays = dict(
        time_s=times,
        fixed_id_m=np.linalg.norm(xa - xb, axis=-1).max(1),
        COM_m=np.linalg.norm(ma["centroid_m"] - mb["centroid_m"], axis=1),
        centered_fixed_id_m=np.linalg.norm(
            ma["centered_positions_m"] - mb["centered_positions_m"], axis=-1
        ).max(1),
        kabsch_m=np.array([rigid_fit(x, y)[1] for x, y in zip(xa, xb)]),
        Rg_m=abs(ma["rg_m"] - mb["rg_m"]),
        aspect=abs(ma["aspect"] - mb["aspect"]),
        directional_xy=abs(ma["directional_xy"] - mb["directional_xy"]),
        tensor_m2=np.linalg.norm(ma["tensor_m2"] - mb["tensor_m2"], axis=(-2, -1)),
    )
    summary = {key: float(arr.max(initial=0)) for key, arr in arrays.items() if key != "time_s"}
    summary.update(
        common_end_s=float(times[-1]),
        common_frames=len(times),
        interpretation=kind,
        complete_2ms=bool(times[-1] >= 0.002 - 1e-13),
        fixed_id_pass=bool(times[-1] >= 0.002 - 1e-13 and summary["fixed_id_m"] <= 5e-8),
    )
    if kind.startswith("Different"):
        summary["fixed_id_pass"] = None  # cross-model contrast is not a numerical gate
    return summary, arrays


def topology_record(t, x):
    r = read(CONFIG)["radius_m"]
    ii, jj = np.triu_indices(32, 1)
    edge = np.linalg.norm(x[:, ii] - x[:, jj], axis=-1) <= 2.1 * r
    ts = []
    ids = []
    labs = []
    for k in range(len(x)):
        adj = np.eye(32, dtype=bool)
        adj[ii, jj] = edge[k]
        adj[jj, ii] = edge[k]
        labs.append(connected_components(adj, directed=False)[1])
        if k:
            changed = np.flatnonzero(edge[k] != edge[k - 1])
            ts.extend([t[k]] * len(changed))
            ids.extend([[0, int(ii[j]), int(jj[j]), int(edge[k, j])] for j in changed])
    ledger, _, _ = topology(
        x[0], t, np.array(labs), np.array(ts), np.array(ids, dtype=int).reshape(-1, 4), t[-1], 0.002
    )
    ledger["components_initial_final"] = [len(np.unique(labs[k])) for k in [0, -1]]
    ledger["max_components"] = max(len(np.unique(lab)) for lab in labs)
    ledger["initial_network_persists_entire_2ms"] = bool(
        t[-1] >= 0.002 - 1e-13 and ledger["max_components"] == 1
    )
    return ledger, np.array(labs), np.array(ts), np.array(ids, dtype=int).reshape(-1, 4)


def main():
    analysis = BASE / "analysis"
    analysis.mkdir(exist_ok=True)
    runs = {}
    cost = []
    for campaign in sorted(BASE.glob("*/CAMPAIGN.json")):
        c = read(campaign)
        cost.append(dict(attempt=campaign.parent.name, **c))
        for path in sorted((campaign.parent / "runs").glob("*/RESULT.json")):
            record = read(path)
            key = campaign.parent.name + "/" + path.parent.name
            t, x = load(path.parent)
            m = metrics(t, x)
            topo, labels, ets, eids = topology_record(t, x)
            np.savez_compressed(
                analysis / (key.replace("/", "__") + ".npz"),
                **m,
                labels=labels,
                edge_times_s=ets,
                edge_ids=eids,
            )
            write(analysis / (key.replace("/", "__") + "_TOPOLOGY.json"), topo)
            record.update(
                rg_initial_final_m=m["rg_m"][[0, -1]],
                aspect_initial_final=m["aspect"][[0, -1]],
                xy_initial_final=m["directional_xy"][[0, -1]],
                tensor_initial_final_m2=m["tensor_m2"][[0, -1]],
                axis_identifiable_initial_final=m["axis_identifiable"][[0, -1]],
                angle_to_x_initial_final_deg=m["angle_to_x_deg"][[0, -1]],
                topology=topo,
            )
            runs[key] = record
    write(analysis / "RUN_ANALYSIS.json", runs)
    write(
        analysis / "COSTS.json",
        dict(
            attempts=cost,
            total_process_s=sum(c["wall_s"] for c in cost),
            integration_s=sum(r["integration_s"] for r in runs.values()),
        ),
    )
    comparisons = {}
    validation = BASE / "validation_v1/runs"
    if validation.exists():
        for tag in ["STRONG", "OFF", "WEAK"]:
            a = validation / (tag + "_DT1")
            for suffix, kind in [("DT05", "Time-step sensitivity"), ("FINE", "Field-cache sensitivity")]:
                key = tag + "_" + suffix
                summary, arrays = compare(a, validation / (tag + "_" + suffix), kind)
                comparisons[key] = summary
                np.savez_compressed(analysis / (key + ".npz"), **arrays)
            old = {"STRONG": "OPTIMIZED_REPLAY", "OFF": "OFF", "WEAK": "SHAPE_REPLAY"}[tag]
            key = tag + "_M0_VS_B2"
            summary, arrays = compare(
                OLD / "campaign_v1/runs" / old,
                a,
                "Different contact/mobility approximations; NOT pure numerical error",
            )
            comparisons[key] = summary
            np.savez_compressed(analysis / (key + ".npz"), **arrays)
    write(analysis / "COMPARISONS.json", comparisons)
    historical = {}
    for name in ["OFF", "BASELINE", "OPTIMIZED_REPLAY", "SHAPE_REPLAY"]:
        path = OLD / "campaign_v1/runs" / name
        t, x = load(path)
        m = metrics(t, x)
        _, labels, _, _ = topology_record(t, x)
        with h5py.File(path / "trajectory.h5") as h:
            ledger, _, _ = topology(
                x[0], t, labels, h["edge_event_time_s"][:], h["edge_event_ids"][:], t[-1], 0.002
            )
        historical[name] = dict(
            topology=ledger,
            rg_initial_final_m=m["rg_m"][[0, -1]],
            xy_initial_final=m["directional_xy"][[0, -1]],
            aspect_initial_final=m["aspect"][[0, -1]],
            result=read(path / "RESULT.json"),
        )
    write(analysis / "HISTORICAL_M0.json", historical)
    matched = {}
    for model, root, names in [
        ("M0", OLD / "campaign_v1/runs", ["OPTIMIZED_REPLAY", "OPTIMIZED_HALF_DT", "OPTIMIZED_FINE_CACHE"]),
        ("B2", BASE / "validation_v1/runs", ["STRONG_DT1", "STRONG_DT05", "STRONG_FINE"]),
    ]:
        if not (root / names[0]).exists():
            continue
        points = []
        for name in names:
            t, x = load(root / name)
            j = np.flatnonzero(abs(t - 1e-4) < 1e-13)
            assert len(j) == 1
            points.append(x[j[0]])
        matched[model] = dict(
            time_s=1e-4,
            dt_same_ID_m=float(np.linalg.norm(points[0] - points[1], axis=1).max()),
            cache_same_ID_m=float(np.linalg.norm(points[0] - points[2], axis=1).max()),
        )
    write(
        analysis / "MATCHED_0p1MS.json",
        dict(data=matched, scope="Exactly matching saved time; not 2ms admission"),
    )


if __name__ == "__main__":
    main()
