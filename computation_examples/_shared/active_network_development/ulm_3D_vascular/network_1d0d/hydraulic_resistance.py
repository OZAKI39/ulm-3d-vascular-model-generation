"""All hydraulic inputs and outputs are SI; conversion occurs at the boundary."""
import numpy as np

MU_PA_S = 0.00345312
RHO_KG_M3 = 1056.0
ROI_TARGET_Q_M3_S = 1.551359160440232e-14


def um_to_m(value):
    return np.asarray(value, dtype=float) * 1e-6


def linear_radius_resistance(length_m, r0_m, r1_m, mu=MU_PA_S):
    """Exact integral, avoiding cancellation when r0 approaches r1.

    Factoring (r1**3-r0**3)/(r1-r0) gives this continuous expression.
    Zero length is allowed for subsegment integrals, never for solver edges.
    """
    length, r0, r1, viscosity = np.broadcast_arrays(length_m, r0_m, r1_m, mu)
    if (not all(np.all(np.isfinite(x)) for x in (length, r0, r1, viscosity))
            or np.any(length < 0) or np.any(r0 <= 0) or np.any(r1 <= 0)
            or np.any(viscosity <= 0)):
        raise ValueError('Lengths must be nonnegative; radii and viscosity positive finite SI')
    result = 8 * viscosity * length / np.pi * (r0*r0 + r0*r1 + r1*r1) / (3*r0**3*r1**3)
    if not np.all(np.isfinite(result)):
        raise ValueError('Nonfinite hydraulic resistance')
    return result
