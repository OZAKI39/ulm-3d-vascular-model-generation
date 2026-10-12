"""Conditional original-initial-state trajectory experiments; original kernels reused."""

import argparse
from dataclasses import replace
from time import perf_counter
import h5py
import json
import numpy as np
from mevo_mb.fast_sim.f10h_lubrication import LubricationDiagnostic
from mevo_mb.fast_sim.f10d_metrics import process
from scripts.fast_pzt_f9.common import objects, clean
from scripts.fast_pzt_f10c.runner import run
from .common import BASE, ROOT, js, write, identity, sha


def record(out, chosen, obj, mode, duration, source, simulator=LubricationDiagnostic, **extra):
    parent = ROOT / chosen["relative_path"]
    assert sha(parent / "trajectory.h5") == chosen["sha256"]
    cfg = js(parent / "CONFIG.json")
    summ = js(parent / "SUMMARY.json")
    dt = 5e-6 if mode.endswith("DT5") else 1e-5
    start = perf_counter()
    timing = run(
        out,
        np.array(cfg["initial_positions_m"]),
        cfg["events"],
        obj,
        duration=duration,
        physical_dt_s=dt,
        graphs=False,
        simulator_class=simulator,
        numerical_mode=mode,
        substep_protocol=True,
        export_observer=True,
        sample_ids=summ["sample_id"],
        postprocess=process,
        identity_provider=lambda: source,
        dataset_metadata=dict(
            initial_state_family=summ["initial_state_family"],
            split_groups=summ["split_groups"],
            split_group=summ["split_group"],
            labels=["EXPLORATORY", "PHYSICS_UNCALIBRATED", "DIAGNOSTIC"],
            independent_scientific_initial_states=0,
            initial_state_sources=[dict(parent_sha256=chosen["sha256"], original_initial_state=True)],
            original_physics_mode=mode,
            postprocessing_mode="F10D_optimized",
        ),
        **extra,
    )
    ob = js(out / "OBSERVER.json")
    result = js(out / "SUMMARY.json")
    meta = dict(
        numerical_scheme_id=mode,
        lubrication_iterations=4 if "LUB4" in mode else 2,
        lubrication_residual="Recorded per physical step; offline same-operator audit",
        contact_model="original Jacobi PBD12",
        contact_stiffness=None,
        contact_geometry_mode="original simultaneous sphere-distance PBD",
        contact_outer_iterations=12,
        physical_dt=dt,
        contact_quality=result["numerical_quality"],
        morphology_admission="PENDING",
        timestep_convergence="PENDING",
        execution_commit=source["execution_commit"],
        input_hashes=source["inputs"],
        completed_requested_window=result["status"] == "COMPLETE",
        complete_trajectory=duration == 0.01 and result["status"] == "COMPLETE",
        failure_reason=result["reason"],
        split_group=summ["split_group"],
        training_ready=False,
        observer_mode="Full per-step diagnostic; CUDA span includes CPU readbacks, not pure kernel time",
    )
    if mode.startswith("B2_GEOM"):
        meta.update(
            lubrication_iterations=None,
            lubrication_residual="exact same-discrete-operator Cholesky",
            contact_model="geometric overdamped implicit penalty",
            contact_stiffness=extra["kc_N_m"],
            contact_geometry_mode="current final-position distances and normals; initial-step resistance fixed",
            contact_outer_iterations="per-step in SOFT_STEPS",
        )
    write(out / "NUMERICAL_IDENTITY.json", meta)
    with h5py.File(out / "trajectory.h5", "a") as h:
        h.attrs["f10h_numerical_identity_json"] = json.dumps(clean(meta))
        h.require_group("numerics").attrs["identity_json"] = json.dumps(clean(meta))
        if not mode.startswith("B2_GEOM"):
            h["bubbles"].attrs["metadata_json"] = json.dumps(
                replace(obj[0], mobility_iterations=meta["lubrication_iterations"]).metadata()
            )
    result.update(f10h_numerical_identity=meta)
    write(out / "SUMMARY.json", result)
    return dict(
        case=chosen["case_id"],
        mode=mode,
        target_s=duration,
        status=result["status"],
        accepted_steps=ob["accepted_steps"],
        attempted_steps=ob["attempted_steps"],
        first_rejection=ob["first_rejection_step"],
        timing=timing,
        record_wall_s=perf_counter() - start,
        path=str(out.relative_to(ROOT)),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--phase", choices=["short", "dt5", "long"], default="short")
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    source = identity()
    write(out / "IDENTITY.json", source)
    obj = objects()
    selection = js(ROOT / "results/fast_pzt/f10e/SELECTION.json")
    cases = [selection["branch"], selection["preparation"]]
    if a.phase != "short":
        prior = js(BASE / "short_v1/RESULTS.json")["rows"]
        cases = [
            c for c in cases if any(r["case"] == c["case_id"] and r["status"] == "COMPLETE" for r in prior)
        ]
        assert cases, "No passed LUB4_DT10 window"
    if a.phase == "long":
        cases = [c for c in cases if c["case_id"] == "O_1000_1102"]
    rows = []
    for c in cases:
        mode = "M0_LUB4_PBD12_DT5" if a.phase == "dt5" else "M0_LUB4_PBD12_DT10"
        rows.append(
            record(
                out / (c["case_id"] + "_" + mode),
                c,
                obj,
                mode,
                0.01 if a.phase == "long" else c["target_end_s"],
                source,
            )
        )
        write(out / "RESULTS.json", dict(rows=rows, training_ready=False))
        print(
            c["case_id"],
            mode,
            rows[-1]["status"],
            rows[-1]["accepted_steps"],
            rows[-1]["first_rejection"],
            flush=True,
        )


if __name__ == "__main__":
    main()
