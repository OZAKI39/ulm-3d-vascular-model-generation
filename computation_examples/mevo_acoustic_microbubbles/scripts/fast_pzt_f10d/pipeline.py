"""Opt-in data pipeline around the existing F10C runner; no solver copy."""

from time import perf_counter
import json
import h5py
import numpy as np
from mevo_mb.fast_sim.f10d_metrics import process
from mevo_mb.fast_sim.f10c_families import control_families
from scripts.fast_pzt_f10c.runner import run
from .common import identity, write
from .validate import arrays_equal


def execute(out, records, obj, mode="optimized", duration=0.01):
    families = control_families(records, max_batch=64)
    if len(families) != 1:
        raise ValueError("A single execution requires one synchronous control family")
    f = families[0]
    sources = [r.get("initial_state_source", {}) for r in records]
    initial_families = [r["initial_state_family"] for r in records]
    benchmark = all(r["sample_class"] == "BENCHMARK_ONLY" for r in records)
    metadata = dict(
        control_family_id=f["control_family_sha256"],
        initial_state_family=initial_families,
        split_groups=f["split_groups"],
        split_group=f["split_groups"][0]
        if len(set(f["split_groups"])) == 1
        else "PER_SAMPLE_SEE_split_groups",
        labels=["BENCHMARK_ONLY" if benchmark else "EXPLORATORY_NUMERICAL_SAMPLE", "PHYSICS_UNCALIBRATED"],
        independent_scientific_initial_states=0 if benchmark else len(set(initial_families)),
        initial_state_sources=sources,
        postprocessing_mode=mode,
        original_physics_mode="F6_M0_F10C_full_particle_stride1",
    )
    source = identity()
    started = perf_counter()
    timing = run(
        out,
        f["positions_m"],
        f["events"],
        obj,
        sample_ids=f["sample_ids"],
        postprocess=process if mode == "optimized" else None,
        dataset_metadata=metadata,
        identity_provider=lambda: source,
        dataset_options=None,  # Measured byte shuffle was slower and larger; retain gzip1.
        duration=duration,
    )
    timing["complete_run_wall_s"] = perf_counter() - started
    timing["serialization_and_tail_s"] = timing["complete_run_wall_s"] - timing["resident_body_s"]
    timing["valid_samples_per_gpu_hour"] = 3600 * timing["valid_sample_count"] / timing["complete_run_wall_s"]
    timing["output_bytes"] = (out / "trajectory.h5").stat().st_size
    write(out / "PIPELINE_TIMING.json", timing)
    return timing


def compare_hdf(a, b, mapping=None):
    checks = []
    with h5py.File(a / "trajectory.h5") as old, h5py.File(b / "trajectory.h5") as new:
        ma = json.loads(old.attrs["metadata_json"])
        mb = json.loads(new.attrs["metadata_json"])
        mapping = list(range(old["positions_m"].shape[1])) if mapping is None else mapping
        for group in ["trajectory", "metrics", "morphology"]:
            assert set(old[group]) == set(new[group])
            for key in old[group]:
                aa = old[group + "/" + key][:]
                bb = new[group + "/" + key][:, mapping]
                if group == "trajectory":
                    equal = np.array_equal(aa, bb, equal_nan=True)
                    checks.append(
                        dict(key=group + "/" + key, status="PASS" if equal else "FAIL", bitwise=equal)
                    )
                else:
                    checks.append(arrays_equal(aa, bb, group + "/" + key))
        for key in ["time_s", "physical_dt_s", "control_index", "particle_id", "radii_m"]:
            good = np.array_equal(old[key][:], new[key][:])
            checks.append(dict(key=key, status="PASS" if good else "FAIL", bitwise=good))
        for key in ["label", "time_s"]:
            good = np.array_equal(old["boundary/" + key][:], new["boundary/" + key][:][mapping])
            checks.append(dict(key="boundary/" + key, status="PASS" if good else "FAIL", bitwise=good))
        for i, j in enumerate(mapping):
            aa = old["edge_events/batch_i_j_linked"][:]
            bb = new["edge_events/batch_i_j_linked"][:]
            sa = aa[:, 0] == i
            sb = bb[:, 0] == j
            good = np.array_equal(aa[sa, 1:], bb[sb, 1:]) and np.array_equal(
                old["edge_events/time_s"][:][sa], new["edge_events/time_s"][:][sb]
            )
            checks.append(dict(key=f"edge_events/{i}", status="PASS" if good else "FAIL", bitwise=good))
            for key in ["quality"]:
                assert ma["per_sample_quality"][i][key] == mb["per_sample_quality"][j][key]
            assert ma["per_sample_quality"][i]["status"] == mb["per_sample_quality"][j]["status"]
            assert ma["topology"][i] == mb["topology"][j]
        assert ma["effective_control_events"] == mb["effective_control_events"]
    return dict(
        status="PASS" if all(r["status"] == "PASS" for r in checks) else "FAIL",
        checks=checks,
        quality_controls_topology_exact=True,
        compared_members=len(mapping),
    )


def benchmark_records(size, control, cfg):
    x = np.array(cfg["environment"]["positions_m"])
    variants = [x + np.array([0, dy, 0]) for dy in [0, 2e-7, -2e-7, 4e-7]]
    return [
        dict(
            sample_id=f"benchmark_{control}_{i:03d}",
            positions_m=variants[i % 4],
            events=cfg["events"],
            split_group="F10A_translated_original_Y32",
            initial_state_family="F10A_translated_original_Y32",
            sample_class="BENCHMARK_ONLY",
            initial_state_source=dict(control=control, translation_m=[0, [0, 2e-7, -2e-7, 4e-7][i % 4], 0]),
        )
        for i in range(size)
    ]
