import argparse
from dataclasses import asdict
from time import perf_counter
import h5py
import numpy as np
from mevo_mb.fast_sim.f10b_contact import ContactDiagnosticPhysics
from mevo_mb.fast_sim.f10e_contact import ContactAuditBatch
from mevo_mb.fast_sim.f10d_metrics import process
from scripts.fast_pzt_f10c.runner import run
from .common import BASE, ROOT, identity, objects, write, js, sha
from .prefix import compare


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--phase", choices=["baseline", "alternatives", "extension"], required=True)
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    source = identity()
    write(out / "IDENTITY.json", source)
    selection = js(BASE / "SELECTION.json")
    proto = js(ROOT / "configs/fast_pzt/f10e/protocol.json")
    started = perf_counter()
    obj = objects()
    static = perf_counter() - started
    jobs = []
    if a.phase == "baseline":
        jobs = [
            (selection["preparation"], proto["alternatives"][0]),
            (selection["branch"], proto["alternatives"][0]),
        ]
    elif a.phase == "alternatives":
        jobs = [(selection["branch"], mode) for mode in proto["alternatives"][1:]]
    else:
        decision = js(BASE / "EXTENSION_DECISION.json")
        assert decision["admitted"]
        mode = next(m for m in proto["alternatives"] if m["label"] == decision["selected_mode"])
        jobs = [(selection["branch"], mode)]
    rows = []
    for chosen, mode in jobs:
        parent = ROOT / chosen["relative_path"]
        cfg = js(parent / "CONFIG.json")
        summary = js(parent / "SUMMARY.json")
        assert sha(parent / "trajectory.h5") == chosen["sha256"]
        target = 0.01 if a.phase == "extension" else chosen["target_end_s"]
        pp = ContactDiagnosticPhysics(
            **(asdict(obj[0]) | dict(contact_iterations=mode["contact_iterations"]))
        )
        current_obj = (pp, *obj[1:])
        name = chosen["case_id"] + "_" + mode["label"]
        path = out / name
        t = perf_counter()
        timing = run(
            path,
            np.array(cfg["initial_positions_m"]),
            cfg["events"],
            current_obj,
            sample_ids=summary["sample_id"],
            duration=target,
            physical_dt_s=mode["dt_s"],
            simulator_class=ContactAuditBatch,
            export_observer=True,
            audit_start_s=chosen["audit_start_s"] if a.phase != "extension" else 1.0,
            numerical_mode=mode["mode"],
            postprocess=process,
            identity_provider=lambda: source,
            dataset_metadata=dict(
                initial_state_family=summary["initial_state_family"],
                split_groups=summary["split_groups"],
                split_group=summary["split_group"],
                labels=["EXPLORATORY", "PHYSICS_UNCALIBRATED"],
                independent_scientific_initial_states=1,
                initial_state_sources=[dict(parent_sha256=chosen["sha256"], original_initial_state=True)],
                original_physics_mode=mode["mode"],
                postprocessing_mode="F10D_optimized",
            ),
        )
        timing["complete_run_wall_s"] = perf_counter() - t
        write(path / "PIPELINE_TIMING.json", timing)
        result = js(path / "SUMMARY.json")
        observer = js(path / "OBSERVER.json")
        row = dict(
            case=name,
            mode=mode,
            status=result["status"],
            first_rejection=observer["first_rejection_step"],
            timing=timing,
            target_end_s=target,
            parent_sha256=chosen["sha256"],
            original_initial_positions_sha256=sha(parent / "CONFIG.json"),
        )
        if mode["label"] == "C0":
            check = compare(parent / "trajectory.h5", path / "trajectory.h5")
            write(path / "DIAGNOSTIC_PREFIX_REGRESSION.json", check)
            assert check["status"] == "PASS", check
            row["diagnostic_on_off_prefix"] = "PASS"
        with h5py.File(path / "trajectory.h5") as h:
            row["accepted_saved_end_s"] = float(h["time_s"][-1])
        rows.append(row)
        write(
            out / "RESULTS.json",
            dict(
                rows=rows,
                static_prepare_s=static,
                process_body_s=perf_counter() - started,
                training_ready=False,
            ),
        )
        print(
            dict(
                case=name,
                status=result["status"],
                first_rejection=row["first_rejection"],
                wall_s=timing["complete_run_wall_s"],
            ),
            flush=True,
        )


if __name__ == "__main__":
    main()
