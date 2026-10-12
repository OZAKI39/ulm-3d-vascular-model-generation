"""Four preregistered independent seeds; prepare once, retain every failure."""

import argparse
from time import perf_counter
import numpy as np
import h5py
from scipy.spatial.distance import pdist
from mevo_mb.fast_sim.morphodynamics import rigid_fit
from .common import BASE, ROOT, identity, objects, case, js, write, sha
from .pipeline import execute
from .postprocess_reference import process as reference
from mevo_mb.fast_sim.f10d_metrics import process
from .validate import compare_results


def initial(seed, center):
    rng = np.random.default_rng(seed)
    accepted = []
    attempts = 0
    while len(accepted) < 32 and attempts < 100000:
        point = rng.uniform(-12e-6, 12e-6, 3)
        attempts += 1
        if np.linalg.norm(point) > 12e-6:
            continue
        if accepted and np.min(np.linalg.norm(np.array(accepted) - point, axis=1)) < 2.05 * 1.78e-6:
            continue
        accepted.append(point)
    if len(accepted) != 32:
        raise RuntimeError("Preregistered placement did not complete")
    return np.array(accepted) + center, attempts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    proto = js(ROOT / "configs/fast_pzt/f10d/protocol.json")["diversity"]
    center = np.array(case("O_1000")["environment"]["positions_m"]).mean(0)
    raw = []
    prepared = []
    prep = []
    comparisons = []
    seeds = []
    start_all = perf_counter()
    start = perf_counter()
    obj = objects()
    static = perf_counter() - start
    for seed in proto["seeds"]:
        x, attempts = initial(seed, center)
        raw.append(x)
        family = f"F10D_seed_{seed}"
        record = dict(
            sample_id=family + "_preparation",
            positions_m=x,
            events=case("O_1000")["events"],
            split_group=family,
            initial_state_family=family,
            sample_class="EXPLORATORY_NUMERICAL_SAMPLE",
            initial_state_source=dict(
                seed=seed,
                method=proto["method"],
                attempts=attempts,
                minimum_gap_m=float(pdist(x).min() - 2 * 1.78e-6),
            ),
        )
        seeds.append(record)
        path = out / f"prepare_{seed}"
        timing = execute(path, [record], obj, duration=proto["preparation_s"])
        summary = js(path / "SUMMARY.json")
        prep.append(
            dict(
                seed=seed,
                timing=timing,
                status=summary["status"],
                failures=summary["failures"],
                reason=summary["reason"],
            )
        )
        with h5py.File(path / "trajectory.h5") as h:
            last = h["positions_m"][-1, 0]
            largest = int(h["metrics/largest_size"][-1, 0])
        prep[-1]["largest_size"] = largest
        if timing["valid_sample_count"] == 1:
            prepared.append(
                dict(
                    record,
                    positions_m=last,
                    initial_state_source=dict(
                        seed=seed,
                        preparation_path=str(path.relative_to(BASE)),
                        preparation_sha256=sha(path / "trajectory.h5"),
                        original_generated_state=record["initial_state_source"],
                        preparation_duration_s=proto["preparation_s"],
                        largest_size=largest,
                    ),
                )
            )
    preparation_elapsed = perf_counter() - start_all
    write(
        out / "INITIAL_STATES.json",
        dict(
            records=seeds,
            method=proto["method"],
            scope="Four independent pseudorandom initial arrangements under one numerical model and physical environment; no independent experimental validation",
        ),
    )
    np.savez_compressed(out / "INITIAL_STATES.npz", positions_m=np.array(raw), seeds=proto["seeds"])
    pairwise = [
        dict(seeds=[proto["seeds"][i], proto["seeds"][j]], aligned_rmsd_m=rigid_fit(raw[i], raw[j])[1])
        for i in range(4)
        for j in range(i + 1, 4)
    ]
    assert all(r["aligned_rmsd_m"] > 1e-6 for r in pairwise)
    rows = []
    start = perf_counter()
    for control in proto["controls"]:
        if not prepared:
            break
        records = [
            dict(r, sample_id=r["initial_state_family"] + "_" + control, events=case(control)["events"])
            for r in prepared
        ]
        path = out / control
        timing = execute(path, records, obj, duration=proto["duration_s"])
        summary = js(path / "SUMMARY.json")
        rows.append(
            dict(
                control=control,
                sample_ids=[r["sample_id"] for r in records],
                timing=timing,
                status=summary["status"],
                reason=summary["reason"],
                failures=summary["failures"],
                per_sample_quality=summary["per_sample_quality"],
            )
        )
    integration_output_elapsed = perf_counter() - start
    # Recompute both postprocessors on actual new states, including accepted failed prefixes.
    for path in sorted(out.iterdir()):
        if not path.is_dir() or not (path / "trajectory.h5").exists():
            continue
        with h5py.File(path / "trajectory.h5") as h:
            x = h["positions_m"][:]
            times = h["time_s"][:]
            ei = h["edge_events/batch_i_j_linked"][:]
            et = h["edge_events/time_s"][:]
        old = reference(x, 1.78e-6, times, ei, et)
        new = process(x, 1.78e-6, times, ei, et)
        comparisons.append(dict(path=path.name, **compare_results(old, new)))
    valid = sum(row["timing"]["valid_sample_count"] for row in rows)
    # Each seed contributes multiple correlated control branches, never multiple independent seeds.
    unique_valid = set()
    for row in rows:
        if row["timing"]["valid_sample_count"]:
            unique_valid.update(r["initial_state_family"] for r in prepared)
    production_s = preparation_elapsed + integration_output_elapsed
    status = "PASS" if all(c["status"] == "PASS" for c in comparisons) else "FAIL"
    write(
        out / "RESULTS.json",
        dict(
            status=status,
            preparation=prep,
            rows=rows,
            comparisons=comparisons,
            pairwise_initial_shape=pairwise,
            requested_seeds=4,
            prepared_valid=len(prepared),
            valid_control_trajectories=valid,
            valid_initial_families=len(unique_valid),
            static_prepare_s=static,
            preparation_total_s=preparation_elapsed,
            control_trajectories_s=integration_output_elapsed,
            total_generation_s=production_s,
            valid_trajectories_per_gpu_hour=3600 * valid / production_s,
            unique_initial_families_per_gpu_hour=3600 * len(unique_valid) / production_s,
            scope="Small-sample extrapolation including preparation; 3 controls share each seed. Offline validation time reported separately, not production integration. No independent validation set or model calibration.",
            training_ready=False,
        ),
    )
    assert status == "PASS"
    print(
        dict(
            prepared=len(prepared),
            valid_control_trajectories=valid,
            valid_initial_families=len(unique_valid),
            total_generation_s=production_s,
        ),
        flush=True,
    )


if __name__ == "__main__":
    main()
