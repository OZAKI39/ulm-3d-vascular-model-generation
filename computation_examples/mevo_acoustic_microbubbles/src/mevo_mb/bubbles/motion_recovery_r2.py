"""Bounded S0 recovery at detected discontinuities of the frozen WALL normal.

S1 retains its original endpoint-equilibrium requirement. S0 steps instead solve
the original start-state QP, undergo their own step-doubling and full sweeps, and
are explicitly labelled. No failed S1 candidate is accepted or relabelled.
"""

RECOVERY_CAP_S = 1e-7
BURST_STEPS = 8
MAX_BURSTS = 1024


def selected_method(controller, default="S1"):
    return "S0" if controller.get("recovery_steps_left", 0) else default


def recovery_ceiling(tolerance_factor):
    # Actual refinement of explicit recovery as well as implicit macro steps.
    return RECOVERY_CAP_S * (0.5 if tolerance_factor < 1 else 1.0)


def begin_recovery(controller, exc, failed_dt, ceiling, tolerance_factor=1.0):
    if (
        getattr(exc, "reason", None) != "NONLINEAR_NOT_CONVERGED"
        or exc.diagnostics.get("wall_normal_jump", 0.0) <= 0.01
        or controller.get("recovery_steps_left", 0)
    ):
        return False
    count = controller.get("recovery_bursts", 0)
    if count >= MAX_BURSTS:
        raise RuntimeError("NONSMOOTH_RECOVERY_BUDGET")
    controller.update(
        recovery_bursts=count + 1,
        recovery_steps_left=BURST_STEPS,
        recovery_failed_dt=failed_dt,
        next_dt=min(recovery_ceiling(tolerance_factor), ceiling),
    )
    return True


def accepted_recovery(controller, method, dt, proposed_dt, max_dt):
    if method == "S0" and controller.get("recovery_steps_left", 0):
        controller["recovery_steps_left"] -= 1
        controller["recovery_accepted_macro_steps"] = controller.get("recovery_accepted_macro_steps", 0) + 1
        controller["recovery_duration_s"] = controller.get("recovery_duration_s", 0.0) + dt
        if not controller["recovery_steps_left"]:
            return min(max_dt, max(proposed_dt, controller["recovery_failed_dt"]))
    return min(max_dt, proposed_dt)
