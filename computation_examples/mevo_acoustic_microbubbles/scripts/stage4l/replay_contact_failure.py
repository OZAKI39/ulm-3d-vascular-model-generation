"""Reproduce one rejected Stage 4L step and preserve the contact QP as a fixture."""

import argparse
from pathlib import Path
import json

import numpy as np
from scipy.sparse import save_npz

from scripts.stage4l.dense_probe import BASE, write


def main():
    from mevo_mb.acoustics.basis import ComplexAcousticField
    from mevo_mb.bubbles.ordered_contact_sweep import OrderedContactSweep
    from mevo_mb.coupling.integrators import GlobalSynchronousIntegrator
    from mevo_mb.flow.fields import VTKFlowField
    import mevo_mb.bubbles.tracker as tracking

    p = argparse.ArgumentParser()
    p.add_argument("run", type=Path)
    args = p.parse_args()
    c = json.loads((args.run / "RUN_IDENTITY.json").read_text())["config"]
    field = ComplexAcousticField.load(c["effective_acoustic_field"]["path"])
    f = c["flow"]
    flow = VTKFlowField(
        f["snapshot"],
        velocity_array=f["velocity_array"],
        coordinate_scale_to_m=f["coordinate_scale_to_m"],
        velocity_scale_to_m_s=f["velocity_scale_to_m_s"],
        spatial_method=f["spatial_method"],
    )
    tracker = GlobalSynchronousIntegrator(c, field, flow)
    sweep = OrderedContactSweep(BASE / "unit_tests/ordered_contact_cache")
    tracking.project_contacts = sweep.project
    saved = np.load(args.run / "particles/last_accepted_state.npz")
    x, r, t = saved["positions_m"], saved["radii_m"], float(saved["time_s"])
    out = BASE / "unit_tests/contact_solver_regression" / args.run.name
    out.mkdir(parents=True, exist_ok=False)
    original = tracker.translation.contact
    counter = 0

    def contact(result, x, radii, velocity, tolerance, horizon=None):
        nonlocal counter
        counter += 1
        try:
            return original(result, x, radii, velocity, tolerance, horizon)
        except Exception:
            save_npz(out / "matrix.npz", result.matrix)
            np.savez_compressed(
                out / "contact_input.npz",
                positions_m=x,
                radii_m=radii,
                velocity_m_s=velocity,
                wall_gap_m=result.wall[0],
                wall_normal=result.wall[1],
                wall_valid=result.wall[2],
                drag=result.drag,
                scale=result.scale,
                tolerance_m=tolerance,
                horizon_s=np.nan if horizon is None else horizon,
            )
            write(
                out / "REPRODUCTION.json",
                dict(
                    status="REPRODUCED",
                    source_run=str(args.run),
                    last_accepted_time_s=t,
                    call=counter,
                    failure=tracker.translation.last_failure,
                ),
            )
            raise

    tracker.translation.contact = contact
    state = tracker.evaluate(x, r, t)
    nominal = c["particles"]["dt_s"]
    target = (int(np.floor(t / nominal)) + 1) * nominal
    dt = tracker.select_timestep(x, r, state, target - t)
    x, _, _, _, _ = tracker.step(x, r, t, dt, state)
    tracker.evaluate(x, r, t + dt)
    raise RuntimeError("Original failed step did not reproduce; do not use a guessed fixture")


if __name__ == "__main__":
    main()
