"""Observe a fixed-ID centroid gate crossing at accepted physical steps on GPU.

The linear intersection is a diagnostic derived from two true accepted states;
it never changes a position or inserts a trajectory sample.
"""

import warp as wp
from . import kernels as k
from .engine import FastSimulator
from .f6_engine import F6Simulator
from .f7_navigation import RecordedFeedbackSimulator


@wp.kernel
def observe_gate(
    before: wp.array(dtype=wp.vec3d),
    after: wp.array(dtype=wp.vec3d),
    active: wp.array(dtype=int),
    error: wp.array(dtype=int),
    step: wp.array(dtype=int),
    origin: wp.array(dtype=wp.float64),
    dt: wp.float64,
    gate: wp.float64,
    n: int,
    flag: wp.array(dtype=int),
    interval: wp.array(dtype=wp.float64),
    left: wp.array(dtype=wp.vec3d),
    right: wp.array(dtype=wp.vec3d),
):
    if flag[0] != 0 or error[0] != 0:
        return
    a = wp.float64(0.0)
    b = wp.float64(0.0)
    for i in range(n):
        if active[i] != 1:
            flag[0] = -1  # No surviving-ID recentering after any event.
            return
        a += before[i][0] / wp.float64(n)
        b += after[i][0] / wp.float64(n)
    if a < gate and b >= gate:
        alpha = (gate - a) / (b - a)
        interval[0] = origin[0] + wp.float64(step[0]) * dt
        interval[1] = interval[0] + dt
        interval[2] = alpha
        flag[0] = 1
        for i in range(n):
            left[i] = before[i]
            right[i] = after[i]


class GateMixin:
    def __init__(self, *args, x_gate=None, **kwargs):
        self.gate_ready = False
        self.x_gate = x_gate
        super().__init__(*args, **kwargs)
        if x_gate is not None:
            if self.batch != 1:
                raise ValueError("F8 fixed-ID gate observer is admitted for one 32-ID sample")
            self.gate_before = wp.zeros(self.n, dtype=wp.vec3d, device=self.device)
            self.gate_left = wp.zeros_like(self.gate_before)
            self.gate_right = wp.zeros_like(self.gate_before)
            self.gate_flag = wp.zeros(1, dtype=int, device=self.device)
            self.gate_interval = wp.zeros(3, dtype=wp.float64, device=self.device)
            wp.load_module(module=observe_gate.module, device=self.device)
            self.gate_ready = True

    def reset(self):
        super().reset()
        if self.gate_ready:
            self.gate_flag.zero_()
            self.gate_interval.zero_()

    def launch(self, kernel, dim, args):
        watching = self.gate_ready and kernel is k.accept
        if watching:
            wp.copy(self.gate_before, self.x)
        result = super().launch(kernel, dim, args)
        if watching:
            FastSimulator.launch(
                self,
                observe_gate,
                1,
                [
                    self.gate_before,
                    self.x,
                    self.active,
                    self.error,
                    self.step_count,
                    self.time_origin,
                    self._current_dt,
                    self.x_gate,
                    self.n,
                    self.gate_flag,
                    self.gate_interval,
                    self.gate_left,
                    self.gate_right,
                ],
            )
        return result

    def gate_record(self):
        if not self.gate_ready:
            return None
        flag = int(self.gate_flag.numpy()[0])
        interval = self.gate_interval.numpy()
        row = dict(
            status="CROSSED" if flag == 1 else "ID_EVENT_CENSORED" if flag == -1 else "RIGHT_CENSORED",
            x_gate_m=self.x_gate,
            fixed_id_count=self.n,
            method="Linear intersection of the fixed-ID centroid between consecutive accepted physical steps; derived intersection, not a saved sample frame",
        )
        if flag == 1:
            left, right = self.gate_left.numpy(), self.gate_right.numpy()
            fraction = float(interval[2])
            positions = left + fraction * (right - left)
            row.update(
                bracket_s=interval[:2],
                bracket_width_s=float(interval[1] - interval[0]),
                fraction=fraction,
                time_s=float(interval[0] + fraction * (interval[1] - interval[0])),
                accepted_before_m=left,
                accepted_after_m=right,
                derived_positions_m=positions,
                derived_centroid_m=positions.mean(0),
            )
            assert abs(positions.mean(0)[0] - self.x_gate) < 1e-14
        return row


class GateSimulator(GateMixin, F6Simulator):
    pass


class GateFeedbackSimulator(GateMixin, RecordedFeedbackSimulator):
    pass
