"""One rigid RBC: old-position Euler translation and old-Omega quaternion step."""
from .microbubble import FieldUnavailableError, positive_scalar
from .integrator import euler_position
from .rbc import RBCState
from .rbc_orientation import angular_velocity, advance_orientation, short_axis


def equilibrium_rbc(state, sample):
    if not sample.inside_lumen or sample.tetra_id < 0:
        raise FieldUnavailableError("RBC center has no frozen fluid field")
    omega = angular_velocity(state.short_axis, state.geometry.jeffery_lambda,
                             sample.velocity_gradient_s_inv, sample.vorticity_s_inv)
    return RBCState(state.particle_id, state.position_m, state.quaternion_wxyz,
                    sample.velocity_m_s, omega, state.geometry)


def advance_single_rbc(state: RBCState, fem_field, dt_s):
    """Return independent state. q is wxyz body→world; world delta left-multiplies.

Both position and orientation use old-point fields. Endpoint resampling only
refreshes V/Omega for the returned position/q; it adds no integration stage.
No wall response. Separate validation code classifies segments before stepping.
"""
    dt = positive_scalar(dt_s, "dt_s")
    old = equilibrium_rbc(state, fem_field.sample(state.position_m))
    position = euler_position(old.position_m, old.velocity_m_s, dt)
    q = advance_orientation(old.quaternion_wxyz, old.angular_velocity_s_inv, dt)
    provisional = RBCState(old.particle_id, position, q, old.velocity_m_s, old.angular_velocity_s_inv, old.geometry)
    return equilibrium_rbc(provisional, fem_field.sample(position))
