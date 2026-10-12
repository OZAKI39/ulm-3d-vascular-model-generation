"""Gated B=1/4/16/(64), real serial and simultaneous full-particle trajectories."""

import argparse
import gc
from time import perf_counter
import h5py
import numpy as np
import warp as wp
from .common import BASE, ROOT, case, objects, identity, write, js
from .runner import run


def compare(serial, batch, member):
    differences = []
    tested = []
    with h5py.File(serial / "trajectory.h5") as a, h5py.File(batch / "trajectory.h5") as b:
        for group in ["trajectory", "metrics", "morphology"]:
            for key in a[group]:
                aa = a[group + "/" + key][:, 0]
                bb = b[group + "/" + key][:, member]
                if not np.array_equal(aa, bb, equal_nan=True):
                    differences.append(group + "/" + key)
                tested.append(group + "/" + key)
        for key in ["time_s", "physical_dt_s"]:
            if not np.array_equal(a[key][:], b[key][:]):
                differences.append(key)
        for key in ["label", "time_s"]:
            if not np.array_equal(a["boundary/" + key][0], b["boundary/" + key][member]):
                differences.append("boundary/" + key)
        bi = b["edge_events/batch_i_j_linked"][:]
        select = bi[:, 0] == member
        bb = bi[select].copy()
        bb[:, 0] = 0
        for key, aa, bb in [
            ("edge_ids", a["edge_events/batch_i_j_linked"][:], bb),
            ("edge_times", a["edge_events/time_s"][:], b["edge_events/time_s"][:][select]),
        ]:
            if not np.array_equal(aa, bb):
                differences.append(key)
    return dict(
        member=member,
        status="BATCH_NUMERICAL_FAIL" if differences else "BATCH_NUMERICAL_PASS",
        bitwise=not differences,
        differences=differences,
        tested_arrays=len(tested) + 6,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--max-batch", type=int, default=16)
    ap.add_argument("--only64", action="store_true")
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    started = perf_counter()
    obj = objects()
    static = perf_counter() - started
    cfg = case("O_1000")
    x = np.array(cfg["environment"]["positions_m"])
    events = cfg["events"]
    variants = np.stack([x + np.array([0, dy, 0]) for dy in [0, 2e-7, -2e-7, 4e-7]])
    write(
        out / "INPUTS.json",
        dict(
            positions_m=variants,
            scope="BENCHMARK_ONLY repeated legal translations",
            split_group="F10A_translated_original_Y32",
            training_ready=False,
        ),
    )
    rows = []
    if a.only64:
        prior = js(BASE / "batch_v1/RESULTS.json")["rows"][-1]
        assert (
            prior["B"] == 16 and prior["status"] == "BATCH_NUMERICAL_PASS" and prior["resident_speedup"] >= 2
        )
    for size in [64] if a.only64 else [1, 4, 16, 64]:
        if size > a.max_batch:
            break
        if (
            size == 64
            and not a.only64
            and not (rows[-1]["status"] == "BATCH_NUMERICAL_PASS" and rows[-1]["resident_speedup"] >= 2)
        ):
            break
        serial = []
        serial_paths = []
        for i in range(size):
            path = out / f"B{size}_serial_{i:03d}"
            t = perf_counter()
            timing = run(path, variants[i % 4], events, obj, original=True)
            timing["outer_wall_s"] = perf_counter() - t
            serial.append(timing)
            serial_paths.append(path)
            gc.collect()
        batch = out / f"B{size}_batch"
        t = perf_counter()
        timing = run(batch, np.stack([variants[i % 4] for i in range(size)]), events, obj)
        timing["outer_wall_s"] = perf_counter() - t
        checks = [compare(path, batch, i) for i, path in enumerate(serial_paths)]
        status = (
            "BATCH_NUMERICAL_PASS"
            if all(c["bitwise"] for c in checks) and timing["valid_sample_count"] == size
            else "BATCH_NUMERICAL_FAIL"
        )
        serial_s = sum(r["outer_wall_s"] for r in serial)
        # GPU memory snapshot after all buffers have been released is not a peak.
        # Actual allocator high-water mark is collected separately by job telemetry.
        row = dict(
            B=size,
            status=status,
            comparison=checks,
            serial=serial,
            batch=timing,
            serial_total_s=serial_s,
            batch_total_s=timing["outer_wall_s"],
            resident_speedup=serial_s / timing["outer_wall_s"],
            valid_samples_per_gpu_hour=3600 * timing["valid_sample_count"] / timing["outer_wall_s"],
            cold_static_plus_batch_s=static + timing["outer_wall_s"],
            cold_valid_samples_per_gpu_hour=3600
            * timing["valid_sample_count"]
            / (static + timing["outer_wall_s"]),
            shared_static_prepare_s=static,
            peak_gpu_memory_MiB=None,
            warp_pool_peak_MiB=wp.get_mempool_used_mem_high("cuda:0") / 2**20,
            warp_pool_current_MiB=wp.get_mempool_used_mem_current("cuda:0") / 2**20,
        )
        rows.append(row)
        write(
            out / "RESULTS.json",
            dict(
                rows=rows,
                static_prepare_s=static,
                scope="Actual sequential integrations versus simultaneous GPU batch, identical full32 output; disk compiler cache warm; cold denotes new CUDA context/static input preparation, not forced compiler-cache purge.",
                training_ready=False,
            ),
        )
        print(
            dict(
                B=size,
                status=status,
                serial_s=serial_s,
                batch_s=timing["outer_wall_s"],
                speedup=row["resident_speedup"],
            ),
            flush=True,
        )
        if status != "BATCH_NUMERICAL_PASS":
            break
    # The unmodified F6 serial baseline must also reproduce archived F10A positions.
    with (
        h5py.File(out / ("B64_serial_000/trajectory.h5" if a.only64 else "B1_serial_000/trajectory.h5")) as h,
        h5py.File(ROOT / "results/fast_pzt/f10a/extended_10_v1/O_1000/trajectory.h5") as ref,
    ):
        np.testing.assert_array_equal(h["positions_m"][:], ref["positions_m"][:])
    write(
        out / "COMPLETE.json",
        dict(
            status="COMPLETE",
            F10A_positions_bitwise=True,
            tested_sizes=[row["B"] for row in rows],
            training_ready=False,
        ),
    )


if __name__ == "__main__":
    main()
