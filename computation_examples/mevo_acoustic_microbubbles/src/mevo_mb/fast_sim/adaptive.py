"""Optional error-controlled step doubling with complete diagnostic rollback.

GPU integration stays in existing kernels; this conservative host controller is
for bounded reference/adaptive work, not claimed to outperform fixed CUDA graphs.
Only accepted half steps advance physical time. Trial work/rejections are explicit.
"""

import numpy as np
import warp as wp
from .engine import FastSimulator

STATE_ARRAYS = (
    "x",
    "error",
    "adjacency",
    "event_count",
    "step_count",
    "time_origin",
    "diag",
    "near_diag",
    "free_path",
    "active",
    "pending",
    "fractions",
    "exit_label",
    "exit_time",
    "wall_stats",
    "cache_x",
    "cache_d",
    "wall_gap",
)


class StepDoubling:
    def __init__(self, sim, *, absolute_error_m=2e-9, relative_error=1e-4, minimum_step_s=1.25e-6):
        if min(absolute_error_m, relative_error, minimum_step_s) <= 0:
            raise ValueError("Positive adaptive tolerances required")
        self.s = sim
        self.atol = absolute_error_m
        self.rtol = relative_error
        self.minimum = minimum_step_s
        arrays = [getattr(sim, k) for k in STATE_ARRAYS if hasattr(sim, k)]
        arrays += [sim.flow_query.cache, sim.flow_query.stats, sim.flow_query.valid, sim.flow_query.failed_x]
        self.state = [(a, wp.empty_like(a)) for a in arrays]
        self.time = 0.0
        self.steps = []
        self.attempts = []
        self.trial_steps = 0

    def checkpoint(self):
        for source, backup in self.state:
            wp.copy(backup, source)

    def restore(self):
        for source, backup in self.state:
            wp.copy(source, backup)

    def trial(self, dt, count, graphs):
        self.s.time_origin.assign(np.array([self.time]))
        self.s.step_count.zero_()
        FastSimulator.advance(self.s, dt, count, graphs)
        self.trial_steps += count
        return self.s.x.numpy().reshape(self.s.batch, self.s.n, 3), int(self.s.error.numpy()[0])

    def advance(self, nominal, count, graphs):
        end = self.time + nominal * count
        h = min(2 * nominal, end - self.time)
        while self.time < end - 1e-15:
            h = min(h, end - self.time)
            self.checkpoint()
            coarse, err1 = self.trial(h, 1, graphs)
            self.restore()
            fine, err2 = self.trial(h / 2, 2, graphs)
            error = float(np.max(np.linalg.norm(fine - coarse, axis=-1)))
            rel = fine - fine.mean(axis=1, keepdims=True)
            scale = float(np.sqrt(np.mean(np.sum(rel**2, axis=-1))))
            budget = self.atol + self.rtol * scale
            q = self.s.quality()
            good = (
                not (err1 or err2)
                and np.isfinite(error)
                and error <= budget
                and q["max_accepted_penetration_m"] <= 0.01 * self.s.p.radius_m
                and q["max_contact_pass_correction_m"] <= 0.1 * self.s.p.radius_m
            )
            self.attempts.append((self.time, h, error, budget, int(good), err1, err2))
            if good:
                self.steps.extend([h / 2, h / 2])
                self.time += h
                if error < budget * 0.25:
                    h = min(2 * h, 2 * nominal)
            else:
                self.restore()
                if h / 4 < self.minimum - 1e-15:
                    raise RuntimeError("ADAPTIVE_MINIMUM_STEP: tolerance/domain/quality failure retained")
                h /= 2
