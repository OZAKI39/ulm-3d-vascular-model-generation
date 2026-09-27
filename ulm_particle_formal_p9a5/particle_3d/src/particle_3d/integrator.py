"""Explicit Euler position update with caller-supplied dt; no default timestep."""
from .microbubble import MicrobubbleState, equilibrium_state, positive_scalar, vector3


def euler_position(position_m, velocity_m_s, dt_s):
    dt = positive_scalar(dt_s, "dt_s")
    return vector3(vector3(position_m, "position_m") + dt * vector3(velocity_m_s, "velocity_m_s"), "new_position_m")


def advance_single_microbubble(state: MicrobubbleState, fem_field, dt_s):
    """Return a fresh state; never modify input.

    Position uses ONLY old-position velocity: x_new=x_old+dt*u(x_old).
    Returned V/Omega are refreshed at x_new so they describe that position.
    This endpoint refresh is not an extra position integration stage.
    Outside start/end raises FieldUnavailableError. Segment classification is
    a separate validation diagnostic; this API applies no boundary response.
    """
    dt = positive_scalar(dt_s, "dt_s")
    old = equilibrium_state(state, fem_field.sample(state.position_m))
    position = euler_position(old.position_m, old.velocity_m_s, dt)
    provisional = MicrobubbleState(old.particle_id, position, old.radius_m,
                                   old.velocity_m_s, old.angular_velocity_s_inv)
    return equilibrium_state(provisional, fem_field.sample(position))
