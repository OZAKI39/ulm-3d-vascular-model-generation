"""Explicit identity for the retained Stage 4 synchronous reference integrator."""
from mevo_mb.coupling.bifurcation_experiment import BifurcationExperimentTracker


class GlobalSynchronousIntegrator(BifurcationExperimentTracker):
    """P0: inherits the original evaluator, RK step, contact solver and run loop unchanged.

    The historical absorbing classification boundary is preserved for regression.
    It is not reinterpreted as proven physical terminality. The Stage 4I exact
    interior-wall certificate remains an explicitly installed execution helper.
    """
