"""F3 GPU execution modes for the frozen F2 model; no changes to physics or integration."""

from dataclasses import asdict
from time import perf_counter
import numpy as np
import warp as wp
from .f2_engine import F2Simulator
from .nearfield import F2Physics
from . import kernels as k
from . import f3_kernels as f3


class F3Simulator(F2Simulator):
    def __init__(self, *args, execution_mode="stable_reduction", **kwargs):
        if execution_mode not in ["fast_atomic", "stable_reduction", "stable_pairs", "stable_components"]:
            raise ValueError("Unknown numerical execution mode")
        self.execution_mode = execution_mode
        super().__init__(*args, **kwargs)
        if asdict(self.p) != asdict(F2Physics()):
            raise ValueError("F3 experiment physics is frozen to selected F2 model")
        self.pair_buffer = wp.zeros((self.total, self.n), dtype=wp.vec3d, device=self.device)
        start = perf_counter()
        wp.load_module(module=f3, device=self.device)
        wp.synchronize_device(self.device)
        self.compile_s += perf_counter() - start

    def pairs(self):
        if self.execution_mode in ["fast_atomic", "stable_components"]:
            return super().pairs()
        p = self.p
        self.launch(
            f3.pair_buffer,
            (self.batch, self.n, self.n),
            [
                self.x,
                self.harmonics,
                self.pressure,
                self.pair_buffer,
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
        self.launch(f3.pair_sum, self.total, [self.pair_buffer, self.force, self.n])

    def launch(self, kernel, dim, args):
        if kernel is k.component_balance and self.execution_mode in ["stable_reduction", "stable_components"]:
            return super().launch(f3.component_balance_ordered, dim, [*args, self.n])
        return super().launch(kernel, dim, args)

    def run(self, *args, **kwargs):
        result = super().run(*args, **kwargs)
        result["summary"].update(
            execution_mode=self.execution_mode,
            numerical_semantics="Frozen F2 discrete model; ordered pair/member sums only in stable mode; contact and mobility kernels unchanged",
            repeatability_scope="Same GPU, compiler, inputs, block layout and software; no cross-device bitwise promise",
        )
        return result
