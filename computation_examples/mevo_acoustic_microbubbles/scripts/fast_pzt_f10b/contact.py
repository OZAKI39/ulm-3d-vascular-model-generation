"""Four registered short tests, stopped at the first actual rejected step."""

import argparse
from dataclasses import asdict
import time
import numpy as np
import warp as wp
from mevo_mb.fast_sim.f10b_contact import ContactDiagnosticPhysics, ContactTrace, distances
from .common import ROOT, BASE, case, js, write, objects, identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--cases", nargs="+", default=["A0", "A1", "A2", "A3"])
    a = ap.parse_args()
    out = BASE / a.label
    out.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    proto = js(ROOT / "configs/fast_pzt/f10b/protocol.json")
    cfg = case("O_3000", "pilot_v3")
    x = np.array(cfg["environment"]["positions_m"])
    write(out / "IDENTITY.json", identity())
    write(out / "CONFIG.json", dict(original=cfg, diagnostic=proto["contact"]))
    p, b, r, g, f = objects()
    results = {}
    for label, dt, iters in proto["contact"]["cases"]:
        if label not in a.cases:
            continue
        pp = ContactDiagnosticPhysics(**(asdict(p) | dict(contact_iterations=iters)))
        sim = ContactTrace(
            x,
            pp,
            r,
            b,
            flow=f,
            vessel=g,
            flow_enabled=True,
            acoustic_enabled=True,
            mobility_mode="f5_local",
            gradient_mode="consistent_trilinear",
        )
        sim.control(cfg["events"][0])
        traces = []
        failure = None
        for step in range(round(proto["contact"]["max_s"] / dt)):
            sim.step(dt)
            wp.synchronize()
            t = sim.trace
            q = sim.quality()
            proposal = t["proposal_positions_m"]
            t.update(
                step=step,
                start_s=step * dt,
                attempted_end_s=(step + 1) * dt,
                quality=q,
                secondary_force_N=t["total_force_N"] - t["primary_force_N"],
                prediction=distances(t["predicted_positions_m"], pp.radius_m),
                final=distances(proposal, pp.radius_m),
            )
            t["balanced_velocity_m_s"] = (t["predicted_positions_m"] - t["initial_positions_m"]) / dt
            t["max_free_displacement_m"] = float(np.linalg.norm(t["free_velocity_m_s"] * dt, axis=1).max())
            failed = (
                q["error_bits"] != 0
                or t["final"]["max_overlap_m"] > 0.01 * pp.radius_m
                or q["max_contact_pass_correction_m"] > 0.1 * pp.radius_m
                or q["wall"]["max_pass_correction_m"] > 0.1 * pp.radius_m
                or not np.isfinite(proposal).all()
            )
            t["accepted"] = not failed
            traces.append(t)
            if failed:
                failure = (step + 1) * dt
                break
        result = dict(
            label=label,
            dt_s=dt,
            contact_rounds=iters,
            first_rejection_s=failure,
            accepted_duration_s=(len(traces) - (failure is not None)) * dt,
            attempted_steps=len(traces),
            classification="NUMERICAL_CONTACT_FAILURE" if failure is not None else "SHORT_NUMERICAL_PASS",
            physics="PHYSICAL_MODEL_UNCALIBRATED",
            long_trajectory_admitted=False,
            traces=traces,
        )
        write(out / (label + ".json"), result)
        results[label] = {k: v for k, v in result.items() if k != "traces"}
        print(results[label], flush=True)
    write(
        out / "RESULTS.json",
        dict(
            cases=results,
            process_body_s=time.perf_counter() - started,
            scope="per-round synchronized diagnostic, not performance measurement",
        ),
    )


if __name__ == "__main__":
    main()
