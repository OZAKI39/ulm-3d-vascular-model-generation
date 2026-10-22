"""Real OFF/S/O integrations, exact family routing and mixed-input order recovery."""

import argparse
from time import perf_counter
import numpy as np
from mevo_mb.fast_sim.f10c_families import control_families
from .common import BASE, identity, objects, case, write, js
from .pipeline import benchmark_records, execute, compare_hdf


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    records = []
    for control, size in zip(["OFF", "S_1000", "O_1000"], [5, 7, 8]):
        records.extend(benchmark_records(size, control, case(control)))
    order = np.random.default_rng(110).permutation(len(records))
    records = [records[i] for i in order]
    write(
        out / "INPUTS.json",
        dict(
            records=records,
            scope="BENCHMARK_ONLY",
            control_budget="Actual archived F10A EXPLORATORY_HIGH_DRIVE source commands; 1 kPa refers to initial local bubble pressure, not source amplitude cap",
        ),
    )
    start = perf_counter()
    obj = objects()
    static = perf_counter() - start
    start = perf_counter()
    groups = control_families(records, max_batch=16)
    group_s = perf_counter() - start
    lookup = {r["sample_id"]: r for r in records}
    paths = {}
    members = {}
    rows = []
    start = perf_counter()
    for i, f in enumerate(groups):
        rr = [lookup[s] for s in f["sample_ids"]]
        path = out / f"family_{i:02d}"
        timing = execute(path, rr, obj)
        summary = js(path / "SUMMARY.json")
        for j, sample in enumerate(f["sample_ids"]):
            paths[sample] = path
            members[sample] = j
        rows.append(
            dict(
                control=rr[0]["initial_state_source"]["control"],
                control_family_id=f["control_family_sha256"],
                sample_ids=f["sample_ids"],
                batch_size=len(rr),
                capacity=16,
                fill_fraction=len(rr) / 16,
                completed=len(rr) if summary["status"] == "COMPLETE" else 0,
                valid=timing["valid_sample_count"],
                aborted=len(rr) - timing["valid_sample_count"],
                timing=timing,
            )
        )
    grouped_execution = perf_counter() - start
    serial = []
    comparisons = []
    for i, r in enumerate(records):
        path = out / f"serial_{i:03d}"
        timing = execute(path, [r], obj)
        serial.append(timing)
        comparison = compare_hdf(path, paths[r["sample_id"]], [members[r["sample_id"]]])
        comparisons.append(dict(sample_id=r["sample_id"], **comparison))
    restored = [
        dict(
            sample_id=r["sample_id"],
            path=str(paths[r["sample_id"]].relative_to(BASE)),
            member=members[r["sample_id"]],
        )
        for r in records
    ]
    assert [r["sample_id"] for r in restored] == [r["sample_id"] for r in records]
    valid = sum(r["valid"] for r in rows)
    serial_s = sum(t["complete_run_wall_s"] for t in serial)
    status = "PASS" if all(c["status"] == "PASS" for c in comparisons) else "FAIL"
    write(
        out / "RESULTS.json",
        dict(
            status=status,
            rows=rows,
            serial=serial,
            comparisons=comparisons,
            original_order_mapping=restored,
            grouping_s=group_s,
            static_prepare_s=static,
            grouped_execution_s=grouped_execution,
            grouped_total_with_static_s=group_s + static + grouped_execution,
            serial_total_s=serial_s,
            actual_grouping_speedup_vs_serial=serial_s / (group_s + grouped_execution),
            valid_sample_count=valid,
            resident_valid_samples_per_gpu_hour=3600 * valid / (group_s + grouped_execution),
            total_valid_samples_per_gpu_hour=3600 * valid / (group_s + static + grouped_execution),
            independent_initial_families=0,
            training_ready=False,
        ),
    )
    assert status == "PASS"
    print(dict(status=status, valid=valid, grouped_s=grouped_execution, serial_s=serial_s), flush=True)


if __name__ == "__main__":
    main()
