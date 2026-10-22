"""Bounded centroid feedback using actual consistent pressure/primary response.

No prescribed particle velocity: the optimizer sets only four source controls;
all motion follows the existing force/mobility/contact equations.
"""

import time
import numpy as np
from scipy.optimize import minimize
from .consistent_field import query_numpy
from .f6_engine import F6Simulator


class BoundedFeedback:
    def __init__(self, basis, response, physics, target, flow_function, policy):
        self.basis, self.response, self.p = basis, response, physics
        self.target = np.array(target)
        self.flow = flow_function
        self.policy = policy
        self.previous = np.r_[np.ones(4), np.zeros(3)]
        self.records = []

    def update(self, positions, t):
        begin = time.perf_counter()
        x = np.mean(positions, axis=0, keepdims=True)
        samples = [
            query_numpy(v, self.basis.origin, self.basis.spacing, x) for v in self.basis.db.pressure_basis
        ]
        P = np.array([v[0][0] for v in samples])
        G = np.stack([v[1][0] for v in samples], axis=-1)
        uf = np.asarray(self.flow(x))[0]
        desired = self.policy["gain_s_inv"] * (self.target - x[0])
        cap = self.policy["amplitude_cap_Pa"]

        def objective(z):
            u = z[:4] * cap * np.exp(1j * np.r_[0.0, z[4:]])
            pressure = P @ u
            if 0 < abs(pressure) < self.response.knots[0]:
                return 1e8
            try:
                v = self.response.query(np.array([pressure]))[0, 0]
            except ValueError:
                return 1e8
            acoustic = (
                -0.5 * np.real(v * (G @ u).conj()) / (6 * np.pi * self.p.viscosity_Pa_s * self.p.radius_m)
            )
            return float(
                np.sum((uf + acoustic - desired) ** 2) / (10e-6) ** 2
                + 0.1 * max(0.0, 100 - abs(pressure)) ** 2 / 100**2
            )

        fit = minimize(
            objective,
            self.previous,
            method="L-BFGS-B",
            bounds=[(0.0, 1.0)] * 4 + [(-np.pi, np.pi)] * 3,
            options={
                "maxiter": self.policy["optimizer_maxiter"],
                "maxfun": self.policy["optimizer_maxfun"],
                "ftol": 1e-10,
                "gtol": 1e-8,
            },
        )
        z = fit.x if np.isfinite(fit.fun) and fit.fun <= objective(self.previous) else self.previous
        self.previous = z.copy()
        event = dict(
            time_s=float(t),
            amplitude_Pa=(z[:4] * cap).tolist(),
            phase_rad=np.r_[0.0, z[4:]].tolist(),
            enabled=[True] * 4,
            secondary_enabled=True,
            flow_enabled=self.policy["flow_enabled"],
        )
        self.records.append(
            dict(
                time_s=float(t),
                centroid_m=x[0].tolist(),
                target_error_m=float(np.linalg.norm(x[0] - self.target)),
                objective=objective(z),
                optimizer_success=bool(fit.success),
                evaluations=int(fit.nfev),
                wall_s=time.perf_counter() - begin,
            )
        )
        return event


class FeedbackSimulator(F6Simulator):
    def __init__(self, *args, controller, **kwargs):
        self.controller = controller
        self.feedback_events = []
        self.feedback_active = False
        super().__init__(*args, **kwargs)

    def state(self):
        x = super().state()
        if self.feedback_active:
            t = float(self.time_origin.numpy()[0] + self.step_count.numpy()[0] * self._current_dt)
            due = self.feedback_next * self.controller.policy["interval_s"]
            if t >= due - 1e-13 and t < self.feedback_end - 1e-13:
                if abs(t - due) > 1e-11:
                    raise RuntimeError("Feedback event not reached exactly; choose aligned save cadence")
                event = self.controller.update(x[0], due)
                self.feedback_events.append(self.control(event))
                self.feedback_next += 1
        return x

    def run(self, duration_s, dt_s, events, sample_interval_s=1e-4, *args, **kwargs):
        interval = self.controller.policy["interval_s"]
        if not np.isclose(
            interval / sample_interval_s, round(interval / sample_interval_s), atol=1e-12, rtol=0
        ):
            raise ValueError("Feedback interval must align with saved-state cadence")
        if kwargs.get("adaptive", False):
            raise ValueError("Feedback/adaptive combination not admitted; independently validated modes")
        initial = self.controller.update(self.initial[0], 0.0)
        self.feedback_events = [initial]
        self.feedback_next = 1
        self.feedback_end = duration_s
        self.feedback_active = True
        try:
            result = super().run(duration_s, dt_s, [initial], sample_interval_s, *args, **kwargs)
        finally:
            self.feedback_active = False
        self.feedback_events[0] = result["summary"]["controls"][0]
        result["summary"]["controls"] = self.feedback_events
        result["summary"]["feedback_policy"] = self.controller.policy
        result["summary"]["feedback_optimization"] = self.controller.records
        ci = (
            np.searchsorted(
                [e["time_s"] for e in self.feedback_events], result["time_s"] + 1e-13, side="right"
            )
            - 1
        )
        result["control_index"] = ci
        for key in ["flow_enabled", "acoustic_enabled", "secondary_enabled"]:
            result["flow"][key] = np.array([self.feedback_events[i][key] for i in ci])
        return result
