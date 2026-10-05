from dataclasses import replace

import numpy as np
import pytest

from mevo_mb.bubbles.forces import primary_bjerknes
from mevo_mb.bubbles.translation import QuasiSteadyStokesTranslation, primary_cycle_average
from mevo_mb.validation.doinikov_dayton import DoinikovDaytonParameters, DoinikovDaytonSingleBubbleReference


def test_cycle_integral_keeps_all_and_only_matching_pressure_harmonics():
    t = np.linspace(0, 4e-6, 4097)
    w = 2 * np.pi * 1e6
    v1, v2 = (2 + 3j) * 1e-18, (-1 + 0.5j) * 1e-18
    g1 = np.array([2 + 1j, -3 + 0.4j, 0]) * 1e8
    g2 = np.array([0.0, 2.0 - 1j, 1.0]) * 1e8
    v = 3e-17 + np.real(v1 * np.exp(1j * w * t) + v2 * np.exp(2j * w * t))
    g = np.real(np.exp(1j * w * t[:, None]) * g1)
    np.testing.assert_allclose(
        primary_cycle_average(t, v, g), primary_bjerknes(v1, g1), rtol=1e-12, atol=1e-24
    )
    g += np.real(np.exp(2j * w * t[:, None]) * g2)
    expected = primary_bjerknes(v1, g1) + primary_bjerknes(v2, g2)
    np.testing.assert_allclose(primary_cycle_average(t, v, g), expected, rtol=1e-12, atol=1e-24)


def test_stokes_drag_units_and_reciprocal_scaling():
    drag = QuasiSteadyStokesTranslation(0.001)
    radii = np.array([1e-6, 2e-6])
    force = np.array([[6 * np.pi * 1e-9, 0, 0], [6 * np.pi * 1e-9, 0, 0]])
    np.testing.assert_allclose(drag.velocity(force, radii)[:, 0], [1.0, 0.5])


def test_reference_equilibrium_and_analytic_derivatives():
    p = DoinikovDaytonParameters(pressure_Pa=0)
    model = DoinikovDaytonSingleBubbleReference(p)
    assert abs(model.G(model.r0, 0, 0, 0, 0)) < 1e-10
    assert np.linalg.norm(model.rhs(0, [model.r0, 0, 0, 0])) < 1e-3
    vals = [model.r0 * 1.01, 0.1, 1e-8, 0.01, 0.1e-6]
    analytic = model.derivatives_G(*vals)
    for i, step in enumerate([1e-12, 1e-5, 1e-10, 1e-4, 1e-12]):
        a, b = vals.copy(), vals.copy()
        a[i] += step
        b[i] -= step
        fd = (model.G(*a) - model.G(*b)) / (2 * step)
        np.testing.assert_allclose(analytic[i], fd, rtol=1e-5, atol=1e-7)
    # Constant-radius, unforced added mass + Stokes relaxation from Eq32.
    v = 0.01
    acceleration = model.rhs(0, [model.r0, 0, 0, v])[3]
    mass = model.mb + 2 * np.pi / 3 * p.density_liquid_kg_m3 * p.outer_radius_m**3
    np.testing.assert_allclose(
        acceleration, -6 * np.pi * p.viscosity_liquid_Pa_s * p.outer_radius_m * v / mass
    )


@pytest.mark.literature
@pytest.mark.slow
def test_reference_weak_field_matches_published_equation63():
    p = replace(DoinikovDaytonParameters(), pressure_Pa=5000.0)
    model = DoinikovDaytonSingleBubbleReference(p)
    z = model.solve()
    t = z["time_s"]
    force = -4 * np.pi / 3 * z["outer_radius_m"] ** 3 * model.gradient(z["displacement_m"], t)
    q = t >= 10 / p.frequency_Hz
    actual = np.trapezoid(force[q], t[q]) / (t[q][-1] - t[q][0])
    assert abs(actual / model.linear_force() - 1) < 0.03
