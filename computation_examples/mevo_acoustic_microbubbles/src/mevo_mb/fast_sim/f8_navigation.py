"""Bounded all-bubble RPY-aware feedback. Controls only; no imposed motion."""

import time
import numpy as np
from scipy.optimize import minimize
from .f8_predictor import MatchedPredictor


class MatchedFeedback:
    def __init__(self, basis, response, physics, target, flow_function, policy):
        self.target = np.asarray(target, float)
        self.flow = flow_function
        self.policy = policy
        self.predictor = MatchedPredictor(basis, response, physics, mode=policy["prediction_mobility"])
        self.previous = np.r_[np.ones(4), np.zeros(3)]
        self.records = []

    def update(self, positions, t):
        begin = time.perf_counter()
        x = np.asarray(positions)
        center = x.mean(axis=0)
        self.predictor.prepare(x, self.flow(x))
        desired = self.policy["gain_s_inv"] * (self.target - center)
        gate_task = self.policy.get("task_type") == "flow_guided_transverse"
        if gate_task:
            uflow = self.flow(x).mean(axis=0)
            remaining = max(self.policy["interval_s"], (self.target[0] - center[0]) / max(uflow[0], 1e-12))
            desired = (self.target - center) / remaining
            desired[0] = uflow[0]
        dx = x - center
        rg = np.sqrt(np.mean(np.sum(dx * dx, axis=1)))
        cap = self.policy["amplitude_cap_Pa"]
        scale = self.policy.get("velocity_scale_m_s", 1e-5)
        best = [np.inf, self.previous.copy(), None]
        evaluations = 0
        invalid = 0
        floor = self.policy.get("minimum_control_pressure_Pa", 0.0)
        samples = self.predictor.samples.numpy()
        samples = samples[..., 0] + 1j * samples[..., 1]

        def pressure_constraint(z):
            return (np.abs(samples @ (z[:4] * cap * np.exp(1j * np.r_[0.0, z[4:]]))).min() - floor) / 100.0

        def evaluate(z):
            nonlocal evaluations, invalid
            evaluations += 1
            try:
                v = self.predictor.evaluate(z[:4] * cap, np.r_[0.0, z[4:]])
                mean = v.mean(axis=0)
                relative = v - mean
                rate = float(np.mean(np.sum(dx * relative, axis=1)) / rg)
                pressure = np.linalg.norm(self.predictor.pressure.numpy(), axis=1)
                cost = float(
                    (
                        np.sum(((mean - desired)[1:] / scale) ** 2)
                        if gate_task
                        else np.sum(((mean - desired) / scale) ** 2)
                    )
                    + self.policy.get("expansion_penalty", 0.1) * (max(0.0, rate) / scale) ** 2
                    + 0.1 * (max(0.0, 100 - float(pressure.min())) / 100) ** 2
                )
                row = dict(
                    predicted_mean_velocity_m_s=mean,
                    predicted_rg_rate_m_s=rate,
                    predicted_relative_rms_m_s=np.sqrt(np.mean(np.sum(relative**2, axis=1))),
                    min_pressure_Pa=pressure.min(),
                    max_pressure_Pa=pressure.max(),
                )
                if cost < best[0] and pressure.min() >= floor - 1e-6:
                    best[:] = [cost, z.copy(), row]
                return cost
            except ValueError:
                invalid += 1
                return 1e12

        initial = evaluate(self.previous)
        if floor > 0:
            fit = minimize(
                evaluate,
                self.previous,
                method="SLSQP",
                bounds=[(0.0, 1.0)] * 4 + [(-np.pi, np.pi)] * 3,
                constraints=[dict(type="ineq", fun=pressure_constraint)],
                options=dict(maxiter=self.policy["optimizer_maxiter"], ftol=1e-9),
            )
        else:
            fit = minimize(
                evaluate,
                self.previous,
                method="L-BFGS-B",
                bounds=[(0.0, 1.0)] * 4 + [(-np.pi, np.pi)] * 3,
                options=dict(
                    maxiter=self.policy["optimizer_maxiter"],
                    maxfun=self.policy["optimizer_maxfun"],
                    ftol=1e-10,
                    gtol=1e-8,
                ),
            )
        feasibility_repair = None
        if best[2] is None and self.policy.get("analytic_pressure_feasibility_repair", False):
            # Linear source superposition gives an exact amplitude-only projection
            # onto the existing lower-pressure constraint. This is an actual source
            # command within its original cap, NOT an acoustic-force multiplier.
            zfix = fit.x.copy()
            minimum = float(np.abs(samples @ (zfix[:4] * cap * np.exp(1j * np.r_[0.0, zfix[4:]]))).min())
            if minimum > 0:
                factor = floor * (1 + 1e-10) / minimum
                if factor >= 1 and np.max(zfix[:4] * factor) <= 1:
                    zfix[:4] *= factor
                    evaluate(zfix)  # Recheck table, pair closure, full RPY and pressure.
                    feasibility_repair = dict(
                        source_amplitude_factor=factor,
                        before_min_pressure_Pa=minimum,
                        after_min_pressure_Pa=None if best[2] is None else best[2]["min_pressure_Pa"],
                        optimizer_iterations_added=0,
                    )
        if best[2] is None:
            self.failed_update = dict(
                time_s=float(t),
                positions_m=x.copy(),
                previous=self.previous.copy(),
                optimizer_message=str(fit.message),
                optimizer_success=bool(fit.success),
                optimizer_iterations=int(fit.nit),
                evaluations=evaluations,
                invalid_candidates=invalid,
                feasibility_repair=feasibility_repair,
                minimum_pressure_constraint=pressure_constraint(fit.x),
                wall_s=time.perf_counter() - begin,
            )
            raise RuntimeError("No legal RPY control candidate; no silent Stokes fallback")
        z = best[1]
        self.previous = z.copy()
        self.records.append(
            dict(
                time_s=float(t),
                centroid_m=center.tolist(),
                target_error_m=float(np.linalg.norm(self.target - center)),
                desired_velocity_m_s=desired.tolist(),
                objective=best[0],
                feasibility_repair=feasibility_repair,
                initial_objective=initial,
                optimizer_success=bool(fit.success),
                optimizer_iterations=int(fit.nit),
                optimizer_message=str(fit.message),
                evaluations=evaluations,
                invalid_candidates=invalid,
                wall_s=time.perf_counter() - begin,
                **best[2],
            )
        )
        return dict(
            time_s=float(t),
            amplitude_Pa=(z[:4] * cap).tolist(),
            phase_rad=np.r_[0.0, z[4:]].tolist(),
            enabled=[True] * 4,
            secondary_enabled=True,
            flow_enabled=self.policy["flow_enabled"],
        )
