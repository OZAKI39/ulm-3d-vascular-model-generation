"""Explicitly opt-in allocation limit for bounded synthetic fixed-state benchmarks.

Not a production integrator, not B independent trajectories. F10J equations,
PCG256, true/recursive residual guards, Armijo and quality budgets are inherited.
"""

from .f10j_gpu import OptimizedB2


class ScalingB2(OptimizedB2):
    max_fixed_state_bubbles = 1024
