"""Fixed-state reduction ablation, graph replay, and model-equivalence checks."""

import argparse
import numpy as np
from scipy.sparse.csgraph import connected_components
from mevo_mb.fast_sim.f3_engine import F3Simulator
from mevo_mb.fast_sim.models import Basis, ResponseSection
from mevo_mb.fast_sim.nearfield import F2Physics, pair_reference
from .common import BASE, F2, PROTOCOL, write, js, state, identity
from .snapshot import evaluate


def main(label):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    p = F2Physics()
    x = state()[0]
    table = ResponseSection(BASE / "inputs/response.h5", p)
    basis = Basis(BASE / "inputs/basis.h5", p.frequency_Hz)
    event = js(F2 / "mature_v1/FIELD_SELECTION.json")["baseline"]
    snapshots = {}
    result = {}
    for mode in ["fast_atomic", "stable_pairs", "stable_components", "stable_reduction"]:
        sim = F3Simulator(x, p, table, basis, execution_mode=mode)
        raw = [evaluate(sim, event) for _ in range(PROTOCOL["fixed_state_repeats"])]
        snapshots[mode] = {k: np.asarray([r[k] for r in raw]) for k in raw[0]}
        ref = pair_reference(x, raw[0]["harmonics_m3"], raw[0]["pressure_complex_Pa"], p, sim.chi_zero)
        np.testing.assert_allclose(raw[0]["secondary_N"], ref, rtol=1e-11, atol=1e-24)
        graph = np.linalg.norm(x[:, None] - x[None, :], axis=-1) < (2 + p.lubrication_activation) * p.radius_m
        lab = connected_components(graph, directed=False)[1]
        for r in raw:
            assert np.array_equal(
                r["component_roots"][:, None] == r["component_roots"][None, :], lab[:, None] == lab[None, :]
            )
        result[mode] = {
            k: dict(
                bitwise_equal=bool(all(np.array_equal(v[0], a) for a in v[1:])),
                max_absolute_difference=float(np.max(abs(v - v[0]))),
            )
            for k, v in snapshots[mode].items()
            if k not in ["component_counts", "component_roots"]
        }
        if mode == "stable_reduction":
            assert all(r["bitwise_equal"] for r in result[mode].values())
        np.savez_compressed(out / (mode + "_fixed_state.npz"), **snapshots[mode])
    np.testing.assert_allclose(
        snapshots["fast_atomic"]["total_N"][0],
        snapshots["stable_reduction"]["total_N"][0],
        rtol=1e-11,
        atol=1e-24,
    )
    sim = F3Simulator(x, p, table, basis)
    graph_run = sim.run(0.0005, 1e-5, [event], 0.0001)
    sim.reset()
    direct = sim.run(0.0005, 1e-5, [event], 0.0001, graphs=False)
    np.testing.assert_array_equal(graph_run["positions_m"], direct["positions_m"])
    np.testing.assert_array_equal(graph_run["edge_event_ids"], direct["edge_event_ids"])
    np.testing.assert_array_equal(graph_run["edge_event_time_s"], direct["edge_event_time_s"])
    np.savez_compressed(
        out / "graph_direct.npz",
        graph_positions_m=graph_run["positions_m"],
        direct_positions_m=direct["positions_m"],
        time_s=direct["time_s"],
    )
    cloud = np.array([[0, 80e-6, 40e-6], [10e-6, 80e-6, 40e-6]])
    off = dict(event, amplitude_Pa=[0] * 4)
    s = F3Simulator(cloud, p, table, basis)
    zero = s.run(5e-5, 1e-5, [off], 1e-5)
    assert np.array_equal(zero["positions_m"][-1], cloud[None])
    s.reset()
    timed = s.run(5e-5, 1e-5, [off, dict(event, time_s=17e-6)], 1e-5)
    assert 17e-6 in timed["time_s"] and timed["summary"]["controls"][1]["time_s"] == 17e-6
    write(
        out / "RESULT.json",
        dict(
            status="PASS",
            identity=identity(),
            fixed_state=result,
            checks=[
                "F2 NumPy pair formula agreement",
                "Canonical component membership versus CPU graph",
                "Stable repeat bitwise equality at fixed state",
                "CUDA graph and direct 50-step states/events identical",
                "Dispersed PZT OFF stationary",
                "17 us control switch exact",
            ],
            attribution_scope="Reduction ablations isolate fixed-state floating sum variation. Full trajectory repeats follow separately. Integer event append order is sorted; contact/mobility have fixed per-particle ID loops.",
            graph_direct_steps=50,
            stable_vs_atomic_initial_force_relative=float(
                np.linalg.norm(
                    snapshots["fast_atomic"]["total_N"][0] - snapshots["stable_reduction"]["total_N"][0]
                )
                / np.linalg.norm(snapshots["stable_reduction"]["total_N"][0])
            ),
        ),
    )
    print("PASS: F3 fixed-state, reduction, graph and clock checks", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    main(a.label)
