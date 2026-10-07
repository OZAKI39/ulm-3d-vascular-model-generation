"""One completed 1-second simple scene, then eight identical-seed serial/batch samples."""

import argparse
import json
import numpy as np
import h5py
from .common import BASE, write
from .experiment import run_case


def main(label, performance):
    out = BASE / label
    out.mkdir(parents=True, exist_ok=False)
    perf = json.loads((BASE / performance / "RESULT.json").read_text())
    chosen = perf["selected"]
    cfg = dict(
        count=32,
        seed=20261007,
        cloud_radius_m=20e-6,
        duration_s=0.05,
        dt_s=chosen["dt_s"],
        sample_interval_s=0.0005,
        physics={"contact_iterations": chosen["contact_iterations"]},
        events=[{"time_s": 0, "pressure_Pa": 1000}],
    )
    prediction = perf["predicted_1s_N32_wall_s"]
    if prediction > 90:
        write(
            out / "RESULT.json",
            dict(
                complete=False,
                status="LONG_NOT_RUN",
                predicted_N32_1s_wall_s=prediction,
                reason="Measured short-run rate exceeds the preregistered short long-run budget; prediction is not completion",
            ),
        )
        return
    longrun = run_case(out / "N32_1s", dict(cfg, duration_s=1.0, sample_interval_s=0.001))
    seeds = list(range(20261007, 20261015))
    rows = []
    for i, seed in enumerate(seeds):
        rows.append(run_case(out / f"serial_{i}", dict(cfg, seed=seed)))
    batched = run_case(out / "batch8", dict(cfg, seeds=seeds))
    with h5py.File(out / "batch8/trajectory.h5") as h:
        bx = h["positions_m"][:]
        bt = h["time_s"][:]
        batch_metrics = {k: h["metrics/" + k][:] for k in ["cluster_count", "largest_size", "global_rg_m"]}
    differences = []
    for i in range(8):
        with h5py.File(out / f"serial_{i}/trajectory.h5") as h:
            np.testing.assert_array_equal(h["positions_m"][0, 0], bx[0, i])
            np.testing.assert_array_equal(h["time_s"][:], bt)
            sx = h["positions_m"][:, 0]
            differences.append(
                dict(
                    seed=seeds[i],
                    max_same_ID_position_difference_m=float(np.linalg.norm(sx - bx[:, i], axis=-1).max()),
                    final_largest_size_difference=int(
                        abs(h["metrics/largest_size"][-1, 0] - batch_metrics["largest_size"][-1, i])
                    ),
                    final_global_rg_relative_difference=float(
                        abs(h["metrics/global_rg_m"][-1, 0] - batch_metrics["global_rg_m"][-1, i])
                        / h["metrics/global_rg_m"][-1, 0]
                    ),
                )
            )
    serial_wall = sum(r["dataset_wall_s"] for r in rows)
    batch_wall = batched["dataset_wall_s"]
    serial_valid = sum(r["numerical_quality"] == "PASS" for r in rows)
    batch_valid = sum(status == "PASS" for status in batched["sample_numerical_quality"])
    write(
        out / "RESULT.json",
        dict(
            complete=True,
            long_run=longrun,
            serial=rows,
            batched=batched,
            batch_seed_list=seeds,
            serial_dataset_wall_s=serial_wall,
            batch_dataset_wall_s=batch_wall,
            serial_valid_samples=serial_valid,
            batch_certified_valid_samples=batch_valid,
            raw_sample_throughput_speedup=serial_wall / batch_wall,
            serial_valid_samples_per_GPU_hour=serial_valid / serial_wall * 3600,
            batch_valid_samples_per_GPU_hour=batch_valid / batch_wall * 3600,
            paired_differences=differences,
            scope="Dataset throughput includes actual sampled metrics and HDF5; setup/first compile reported separately",
        ),
    )


if __name__ == "__main__":
    a = argparse.ArgumentParser()
    a.add_argument("--label", required=True)
    a.add_argument("--performance", required=True)
    q = a.parse_args()
    main(q.label, q.performance)
