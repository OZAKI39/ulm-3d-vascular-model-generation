"""Bounded control proposals followed by actual, independently repeated trajectories."""

import time
import numpy as np
from scipy.optimize import minimize
from .common import BASE, CONFIG, event, read, write, sha


class SearchBudget(Exception):
    pass


def run_navigation(c):
    plan = read(CONFIG.parent / "navigation_admission.json")
    assert plan["strong_group_position_eligible"] and plan["strong_basic_pass"]
    # Admission is backed by immutable full trajectories, not optimizer estimates.
    for name, digest in plan["trajectory_sha256"].items():
        assert sha(BASE / "validation_v1/runs" / name / "trajectory.h5") == digest
    basis = c.bases["coarse"]
    pb, gb = basis.query_basis(c.x0)
    zeta = 6 * np.pi * c.physics.viscosity_Pa_s * c.physics.radius_m
    strong = event("OPTIMIZED_REPLAY", 0.002)
    initial = np.r_[np.array(strong["amplitude_Pa"]) / 1500.0, strong["phase_rad"][1:]]
    bounds = [(0.0, 1.0)] * 4 + [(-np.pi, np.pi)] * 3
    for axis, target in c.p["XY_targets_m"].items():
        start = time.perf_counter()
        records = []
        best = [np.inf, None]

        def objective(q):
            if len(records) >= c.p["max_search_objective_evaluations"]:
                raise SearchBudget()
            control = dict(strong, amplitude_Pa=(q[:4] * 1500).tolist(), phase_rad=[0.0, *q[4:].tolist()])
            record = dict(parameters=q.copy(), control=control)
            try:
                u = basis.vector(control["amplitude_Pa"], control["phase_rad"], control["enabled"])
                pressure = u @ pb
                grad = np.einsum("i,inj->nj", u, gb)
                volumes = c.response.query(pressure)
                primary = -0.5 * np.real(volumes[:, 0, None] * grad.conj())
                # Exact instantaneous group balance under this resistance/force model;
                # frozen-state 2ms estimate is only a candidate proposal, not a trajectory.
                estimate = 0.002 * primary.mean(0) / zeta
                value = float(np.sum((estimate - np.array(target)) ** 2) / (2e-7) ** 2)
                record.update(
                    status="VALID_INITIAL_RESPONSE",
                    surrogate_endpoint_m=estimate,
                    initial_pressure_min_max_Pa=[float(abs(pressure).min()), float(abs(pressure).max())],
                    objective=value,
                )
            except ValueError as ex:
                value = 1e6
                record.update(status="REJECTED_RESPONSE", reason=str(ex), objective=value)
            records.append(record)
            if value < best[0]:
                best[:] = [value, control]
            return value

        try:
            result = minimize(
                objective,
                initial,
                method="L-BFGS-B",
                bounds=bounds,
                options=dict(maxiter=3, maxfun=32, eps=1e-4, ftol=1e-8),
            )
            stop = str(result.message)
        except SearchBudget:
            stop = "32 objective evaluations reached; best legal proposal retained"
        write(
            c.out / (axis + "_SEARCH.json"),
            dict(
                target_m=target,
                records=records,
                best_control=best[1],
                stopping_reason=stop,
                wall_s=time.perf_counter() - start,
                model="Frozen-state instantaneous COM force balance; proposal only, no trajectory claim",
            ),
        )
        if best[1] is None:
            continue
        r = c.run("NAV_" + axis + "_DT1", "OPTIMIZED_REPLAY", control=best[1])
        if r["status"] != "PASS":
            continue
        # Quality and actual morphology cost are evaluated in layered postprocessing.
        r = c.run("NAV_" + axis + "_DT05", "OPTIMIZED_REPLAY", control=best[1], dt=0.5e-6)
        if r["status"] == "PASS":
            c.run("NAV_" + axis + "_FINE", "OPTIMIZED_REPLAY", control=best[1], cache="fine")
