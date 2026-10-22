"""Whole-step stiff squeezing and wall approach vs independent scalar ODE."""

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scripts.stage4l.dense_probe import write
from mevo_mb.bubbles.motion_r2 import MotionR2, SemiImplicit
from mevo_mb.bubbles.motion_sweep_r2 import pair_sweep
from mevo_mb.bubbles.constrained_motion import SubdivideInterval

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage5/motion_solver_r2"
R = 1.78e-6
MU = 0.001
DRAG = 6 * np.pi * MU * R


class Geometry:
    def __init__(self, wall):
        self.wall = wall

    def wall_gap_and_normal(self, x, r):
        return (
            (-x[:, 2] - r, np.tile([0.0, 0.0, 1.0], (len(x), 1)), np.ones(len(x), bool))
            if self.wall
            else (np.full(len(x), np.inf), np.zeros_like(x), np.ones(len(x), bool))
        )

    def sweep(self, x, y, r, tol):
        pair = pair_sweep(x, y, r)
        wall = min(self.wall_gap_and_normal(x, r)[0].min(), self.wall_gap_and_normal(y, r)[0].min())
        return dict(
            pass_gate=min(pair, wall) >= -tol,
            minimum_pair_swept_gap_m=pair,
            minimum_wall_swept_gap_bound_m=wall,
            projection_max_m=0.0,
        )


def experiment(task, output=None):
    wall, factor, grid = task
    v = 0.001
    duration = 0.001
    geometry = Geometry(wall)
    model = MotionR2(MU, components=False)
    x = np.array([[0.0, 0.0, -1.05 * R]]) if wall else np.array([[0.0, 0.0, 0.0], [2.05 * R, 0.0, 0.0]])
    flow = np.array([[0.0, 0.0, v]]) if wall else np.array([[v, 0.0, 0.0], [-v, 0.0, 0.0]])
    r = np.full(len(x), R)
    ids = np.arange(len(x))

    def forcing(x, r, t):
        return dict(
            primary_force_N=np.zeros_like(x), secondary_force_N=np.zeros_like(x), flow_velocity_m_s=flow
        )

    step = SemiImplicit(
        model,
        geometry,
        forcing,
        geometry.sweep,
        nonlinear_strategy="quasi_newton",
        atol_R=0.002 * factor,
        rtol=0.01 * factor,
        gap_rtol=0.02 * factor,
    )

    def rhs(t, h):
        c = model.coefficient(h, R, "wall" if wall else "pair")[0]
        return [-v / (1 + c / DRAG)] if wall else [-2 * v / (1 + 2 * c / DRAG)]

    def event(t, h):
        return h[0]

    event.terminal = True
    event.direction = -1
    reference = solve_ivp(
        rhs, [0, duration], [0.05 * R], events=event, dense_output=True, rtol=2e-12, atol=1e-20, max_step=1e-6
    )

    def ref(t):
        return float(reference.sol(t)[0]) if t <= reference.t[-1] else 0.0

    large, _, _, large_diag = step.advance(x, r, ids, 0.0, duration, "single_large_diagnostic")
    large_gap = -large[0, 2] - R if wall else np.linalg.norm(large[1] - large[0]) - 2 * R
    times = [0.0]
    gaps = [0.05 * R]
    errors = []
    dt = 1e-4 * grid
    t = 0.0
    rejected = []
    accepted = []
    for trial in range(10000):
        if t >= duration - 1e-15:
            break
        dt = min(dt, duration - t, 1e-4 * grid)
        try:
            states, d = step.trial(x, r, ids, t, dt)
        except SubdivideInterval as exc:
            rejected.append(dict(time_s=t, dt_s=dt, reason=exc.reason, diagnostics=exc.diagnostics))
            dt *= 0.5
            continue
        for i, (y, u, cache, sub) in enumerate(states):
            at = t + dt * (i + 1) / 2
            gap = -y[0, 2] - R if wall else np.linalg.norm(y[1] - y[0]) - 2 * R
            times.append(at)
            gaps.append(gap)
            errors.append(abs(gap - ref(at)))
            accepted.append(
                dict(time_s=at, dt_s=dt / 2, kkt=sub["kkt"], nonlinear_kkt=sub["qp"]["nonlinear_kkt"])
            )
        x = states[-1][0]
        model.cache = states[-1][2]
        t += dt
        dt *= min(2.0, max(0.5, 0.85 / max(d["error_ratio"], 1e-8) ** 0.5))
    out = Path(output) if output is not None else BASE / "diagnosis/manufactured_v2"
    out.mkdir(exist_ok=True)
    name = f"{'wall' if wall else 'pair'}_tol{factor:g}_grid{grid:g}"
    np.savez_compressed(
        out / f"{name}.npz", time_s=times, gap_m=gaps, reference_gap_m=[ref(at) for at in times]
    )
    summary = dict(
        case="wall" if wall else "pair",
        tolerance_factor=factor,
        max_dt_s=1e-4 * grid,
        completed_s=t,
        raw_large_step_gap_m=float(large_gap),
        raw_large_step_error_m=float(abs(large_gap - ref(duration))),
        raw_large_step_safe=large_diag["safety"]["pass_gate"],
        max_gap_error_m=max(errors),
        endpoint_error_m=errors[-1],
        accepted_steps=len(accepted),
        min_gap_m=min(gaps),
        min_dt_s=min(q["dt_s"] for q in accepted),
        kkt_max=max(q["kkt"] for q in accepted),
        nonlinear_kkt_max=max(q["nonlinear_kkt"] for q in accepted),
        projection_max_m=0.0,
        accepted=accepted,
        rejected=rejected,
    )
    write(out / f"{name}.json", summary)
    return {k: v for k, v in summary.items() if k not in ["accepted", "rejected"]}


if __name__ == "__main__":
    with ProcessPoolExecutor(max_workers=4) as pool:
        rows = list(
            pool.map(
                experiment,
                [
                    (wall, f, g)
                    for wall in [False, True]
                    for f, g in [(1.0, 1.0), (0.25, 1.0), (0.0625, 1.0), (1.0, 0.5), (1.0, 0.25)]
                ],
            )
        )
    write(BASE / "diagnosis/manufactured_v2/SUMMARY.json", rows)
    for row in rows:
        print(row, flush=True)
