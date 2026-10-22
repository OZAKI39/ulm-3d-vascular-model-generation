"""Six bounded original-initial-state replays; stop at each first failed substep."""

import argparse
from dataclasses import replace
from time import perf_counter
import h5py
import numpy as np
from mevo_mb.fast_sim.f10f_smallsteps import SmallStepsBatch
from mevo_mb.fast_sim.f10d_metrics import process
from scripts.fast_pzt_f10c.runner import run
from scripts.fast_pzt_f10e.prefix import compare
from .common import BASE, HISTORY, ROOT, objects, identity, js, sha, write


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    assert js(BASE / "checks_v1/RESULTS.json")["status"] == "PASS"
    source = identity()
    write(out / "IDENTITY.json", source)
    protocol = js(ROOT / "configs/fast_pzt/f10f/protocol.json")
    selection = js(HISTORY / "SELECTION.json")
    start = perf_counter()
    obj = objects()
    static = perf_counter() - start
    rows = []
    for case in ["branch", "preparation"]:
        chosen = selection[case]
        parent = ROOT / chosen["relative_path"]
        assert sha(parent / "trajectory.h5") == chosen["sha256"]
        cfg = js(parent / "CONFIG.json")
        summary = js(parent / "SUMMARY.json")
        for mode, setting in protocol["modes"].items():
            name = chosen["case_id"] + "_" + mode
            path = out / name
            pp = replace(obj[0], contact_iterations=setting["contact_iterations"])
            timing = run(
                path,
                np.array(cfg["initial_positions_m"]),
                cfg["events"],
                (pp, *obj[1:]),
                sample_ids=summary["sample_id"],
                duration=chosen["target_end_s"],
                physical_dt_s=setting["dt_s"],
                simulator_class=SmallStepsBatch,
                substep_protocol=True,
                export_observer=True,
                diagnostic=True,
                numerical_mode=mode,
                postprocess=process,
                identity_provider=lambda: source,
                dataset_metadata=dict(
                    initial_state_family=summary["initial_state_family"],
                    split_groups=summary["split_groups"],
                    split_group=summary["split_group"],
                    labels=["EXPLORATORY", "PHYSICS_UNCALIBRATED"],
                    independent_scientific_initial_states=1,
                    initial_state_sources=[dict(parent_sha256=chosen["sha256"], original_initial_state=True)],
                    original_physics_mode=mode,
                    postprocessing_mode="F10D_optimized",
                ),
            )
            ob = js(path / "OBSERVER.json")
            result = js(path / "SUMMARY.json")
            row = dict(
                case=chosen["case_id"],
                mode=mode,
                status=result["status"],
                setting=setting,
                first_rejection=ob["first_rejection_step"],
                timing=timing,
                attempted_steps=ob["attempted_steps"],
                accepted_steps=ob["accepted_steps"],
                target_end_s=chosen["target_end_s"],
            )
            if mode == "B0":
                old = HISTORY / "contact_v1" / (chosen["case_id"] + "_C0")
                check = compare(old / "trajectory.h5", path / "trajectory.h5")
                assert check["status"] == "PASS", check
                assert (
                    abs(
                        js(old / "OBSERVER.json")["first_rejection_step"]["attempted_end_s"]
                        - ob["first_rejection_step"]["attempted_end_s"]
                    )
                    < 1e-14
                )
                write(path / "B0_REGRESSION.json", check)
            with h5py.File(path / "trajectory.h5") as h:
                row["accepted_saved_end_s"] = float(h["time_s"][-1])
            rows.append(row)
            write(
                out / "RESULTS.json",
                dict(
                    rows=rows,
                    static_prepare_s=static,
                    process_body_s=perf_counter() - start,
                    training_ready=False,
                ),
            )
            print({k: row[k] for k in ["case", "mode", "status", "attempted_steps"]}, flush=True)


if __name__ == "__main__":
    main()
