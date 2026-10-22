"""Opt-in F2 engine extension; F1 execution stays independently selectable."""

from time import perf_counter
import numpy as np
import warp as wp
from .engine import FastSimulator
from . import kernels as k
from . import f2_kernels as f2


class F2Simulator(FastSimulator):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.contact_mobility = wp.empty(self.total, dtype=wp.mat33d, device=self.device)
        self.near_diag = wp.empty((self.total, 3), dtype=wp.float64, device=self.device)
        self.free_path = wp.zeros(self.total, dtype=wp.float64, device=self.device)
        self.chi_zero = (
            self.response.query(np.array([self.response.knots[0] + 0j]))[0, 0] / self.response.knots[0]
        )
        self._contact_prepared = False
        self._current_dt = 0.0
        self._snapshots = None
        self.reset()
        start = perf_counter()
        wp.load_module(module=f2, device=self.device)
        wp.synchronize_device(self.device)
        self.compile_s += perf_counter() - start

    def reset(self):
        super().reset()
        if hasattr(self, "near_diag"):
            d = np.zeros((self.total, 3))
            d[:, [0, 2]] = np.inf
            self.near_diag.assign(d)
            self.free_path.zero_()

    def pairs(self):
        p = self.p
        if p.force_model == "f1_smooth_near_field":
            return super().pairs()
        if p.secondary:
            self.launch(
                f2.pair_force,
                (self.batch, self.n, self.n),
                [
                    self.x,
                    self.harmonics,
                    self.pressure,
                    self.force,
                    self.n,
                    p.radius_m,
                    2 * np.pi * p.frequency_Hz,
                    p.density_kg_m3,
                    p.transition_over_sum,
                    p.soft_core_over_sum,
                    p.feedback_end_over_sum,
                    wp.vec2d(self.chi_zero.real, self.chi_zero.imag),
                    float(self.response.knots[0]),
                    float(self.response.knots[-1]),
                    self.near_diag,
                    self.error,
                ],
            )

    def step(self, dt):
        self._contact_prepared = False
        self._current_dt = dt
        super().step(dt)

    def launch(self, kernel, dim, args):
        if kernel is k.contact and self.p.contact_model == "f2_mobility_aware_contact":
            if not self._contact_prepared:
                super().launch(
                    f2.mobility,
                    self.total,
                    [
                        args[0],
                        self.contact_mobility,
                        self.n,
                        self.p.radius_m,
                        self.p.lubrication_activation,
                        self.p.lubrication_floor,
                    ],
                )
                self._contact_prepared = True
            return super().launch(f2.contact, dim, [*args, self.contact_mobility])
        super().launch(kernel, dim, args)
        if kernel is k.free_velocity:
            super().launch(f2.free_path, self.total, [self.free, self.free_path, self._current_dt])

    def state(self):
        try:
            x = super().state()
        except RuntimeError as exc:
            raise RuntimeError(
                str(exc) + "; F2 bits: closure_singular=64 coupled_pressure_domain=128 mobility_metric=256"
            ) from exc
        if getattr(self, "_snapshots", None) is not None:
            d = self.diag.numpy()
            self._snapshots.append(d.reshape(self.batch, self.n, -1).copy())
        return x

    def quality(self, sample=None):
        q = super().quality(sample)
        d = self.near_diag.numpy()
        path = self.free_path.numpy()
        if sample is not None:
            d = d[sample * self.n : (sample + 1) * self.n]
            path = path[sample * self.n : (sample + 1) * self.n]
        for key, col, op in [
            ("coupled_effective_pressure_min_Pa", 0, np.min),
            ("coupled_effective_pressure_max_Pa", 1, np.max),
            ("closure_min_abs_determinant", 2, np.min),
        ]:
            value = float(op(d[:, col]))
            q[key] = value if np.isfinite(value) else None
        q["cumulative_free_Stokes_displacement_m"] = float(path.sum())
        q["correction_over_free_Stokes_path"] = q["cumulative_correction_m"] / max(float(path.sum()), 1e-300)
        return q

    def run(self, *args, **kwargs):
        self._snapshots = []
        result = super().run(*args, **kwargs)
        if len(self._snapshots) != len(result["time_s"]):
            raise RuntimeError("Quality/position sampling misalignment")
        q = np.array(self._snapshots)
        q[0] = 0
        result["quality_samples"] = q
        self._snapshots = None
        return result
