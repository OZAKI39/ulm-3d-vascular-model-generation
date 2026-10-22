"""Independent initial-only predictions and one output-matched full M0 timer.

No trajectory HDF5 is opened here. Accuracy comparison is a later offline stage.
Every timed sample uses a fresh process and constructs its own static providers.
"""

import argparse
import json
import os
import subprocess
import sys
from time import perf_counter
import h5py
import numpy as np
import warp as wp
from mevo_mb.fast_sim.f10b_surrogate import CentroidSurrogate, admit
from mevo_mb.fast_sim.f6_engine import F6Simulator
from .common import ROOT, BASE, case, js, write, identity, objects, clean


def child(out, mode, key):
    started = perf_counter()
    proto = js(ROOT / "configs/fast_pzt/f10b/protocol.json")
    cfg = case(key)
    ev = cfg["events"][0]
    x = np.array(cfg["environment"]["positions_m"])
    mat = js(BASE / "checks_v1/MATURITY.json")
    obj = objects()
    p, b, r, g, f = obj
    static_s = perf_counter() - started
    start = perf_counter()
    admission = admit(
        x,
        ev,
        obj,
        proto["admission"],
        mobility=cfg["mobility"],
        duration_s=0.01,
        maturity=mat["mature"],
        validated_initial=np.array_equal(x, np.array(case("OFF")["environment"]["positions_m"])),
        validated_control=key in proto["surrogate"]["controls"],
    )
    assert admission["status"] == "ADMITTED_FOR_SHADOW_TEST", admission
    dt = proto["surrogate"]["dt_s"]
    sample = proto["surrogate"]["sample_s"]
    duration = 0.01
    if mode == "FULL":
        sim = F6Simulator(
            x, p, r, b, flow=f, vessel=g, mobility_mode="f5_local", gradient_mode="consistent_trilinear"
        )
        sim.control(ev)
        model_prepare = perf_counter() - start
        start = perf_counter()
        sim.graph(dt, round(sample / dt))
        graph = perf_counter() - start
        sim.state()  # initial error check
        positions = [x.copy()]
        values = [x.mean(0)]
        integration = 0.0
        sync = 0.0
        for _ in range(round(duration / sample)):
            start = perf_counter()
            sim.advance(dt, round(sample / dt))
            wp.synchronize_device(sim.device)
            integration += perf_counter() - start
            start = perf_counter()
            xx = sim.state()[0].copy()
            values.append(xx.mean(0))
            positions.append(xx)
            sync += perf_counter() - start
        q = sim.quality()
        assert (
            q["error_bits"] == 0
            and q["max_accepted_penetration_m"] <= 0.01 * p.radius_m
            and q["max_contact_pass_correction_m"] <= 0.1 * p.radius_m
        )
        assert sim.active.numpy().all() and q["wall"]["cumulative_correction_m"] == 0
        data = dict(
            time_s=np.arange(len(values)) * sample,
            centroid_m=np.array(values),
            label="FULL_PARTICLE_REFERENCE",
            quality=q,
            timing=dict(graph_prepare_s=graph, integration_s=integration, sync_readback_s=sync),
        )
    else:
        sim = CentroidSurrogate(x, p, r, b, f, ev, mode, proto["admission"], admission)
        model_prepare = perf_counter() - start
        data = sim.rollout(duration, dt, sample)
    timing = data["timing"] | dict(static_input_prepare_s=static_s, model_admission_prepare_s=model_prepare)
    timing["preparation_s"] = static_s + model_prepare + timing["graph_prepare_s"]
    start = perf_counter()
    out.mkdir(parents=True, exist_ok=False)
    metadata = dict(
        schema="f10b_centroid_v1",
        label=data["label"],
        evidence=["MODEL_EXPLORATORY", "PHYSICS_UNCALIBRATED"],
        training_ready=False,
        control=ev,
        medium=r.metadata,
        mode=mode,
        dt_s=dt,
        duration_s=duration,
        admission=admission,
        initial_template_label="FROZEN_TEMPLATE_ASSUMPTION",
        initial_template_offsets_m=x - x.mean(0),
        quality=data["quality"],
        execution_commit=identity()["execution_commit"],
    )
    with h5py.File(out / "centroid.h5", "w") as h:
        h["time_s"] = data["time_s"]
        h["centroid_m"] = data["centroid_m"]
        h.attrs["metadata_json"] = json.dumps(clean(metadata))
    write(
        out / "CONFIG.json",
        dict(
            initial_positions_m=x,
            event=ev,
            model=mode,
            protocol=proto,
            source_config=cfg,
            admission=admission,
        ),
    )
    timing["output_s"] = perf_counter() - start
    timing["body_s"] = perf_counter() - started
    write(out / "TIMING.json", timing)
    # This optional validation output is excluded equally from the matched timer.
    if mode == "FULL":
        np.savez_compressed(
            out / "reference_validation_positions.npz", positions_m=np.array(positions), time_s=data["time_s"]
        )
    print(json.dumps(dict(mode=mode, control=key, timing=timing, quality=clean(data["quality"]))), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--mode")
    ap.add_argument("--key")
    a = ap.parse_args()
    out = BASE / a.label
    if a.mode:
        child(out / (a.mode + "_" + a.key), a.mode, a.key)
        return
    assert js(BASE / "checks_v1/RESULTS.json")["status"] == "PASS"
    out.mkdir(parents=True, exist_ok=False)
    write(out / "IDENTITY.json", identity())
    hardware = subprocess.check_output(
        ["nvidia-smi", "--query-gpu=name,driver_version,utilization.gpu,memory.used", "--format=csv"],
        text=True,
    )
    rows = []
    # Single full-reference integration. Six predictions each have independent state.
    jobs = [("FULL", "O_1000")] + [(m, k) for m in ["G0", "G1"] for k in ["OFF", "S_1000", "O_1000"]]
    for mode, key in jobs:
        begin = perf_counter()
        with open(out / (mode + "_" + key + ".log"), "w") as log:
            subprocess.run(
                [
                    sys.executable,
                    "-m",
                    "scripts.fast_pzt_f10b.rollout",
                    "--label",
                    a.label,
                    "--mode",
                    mode,
                    "--key",
                    key,
                ],
                stdout=log,
                stderr=subprocess.STDOUT,
                check=True,
                env=os.environ.copy(),
            )
        row = dict(
            mode=mode,
            control=key,
            subprocess_wall_s=perf_counter() - begin,
            **js(out / (mode + "_" + key) / "TIMING.json"),
        )
        rows.append(row)
        write(
            out / "TIMING.json",
            dict(
                hardware=hardware,
                rows=rows,
                output_scope="201 centroid/time records with controls, medium, quality and identity; full model validation positions saved separately after body timer",
                replicates=1,
                scope="fresh process per sample; static preparation included; no CUDA event/kernel-only claim",
            ),
        )
        print(row, flush=True)
    write(
        out / "COMPLETE.json",
        dict(
            status="COMPLETE",
            independent_predictions=6,
            new_full_reference_runs=1,
            future_reference_reads=False,
        ),
    )


if __name__ == "__main__":
    main()
