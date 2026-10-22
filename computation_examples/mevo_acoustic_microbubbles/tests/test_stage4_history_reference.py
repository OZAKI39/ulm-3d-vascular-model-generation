"""Numerical reference checks; no production dependency on the history module."""
import numpy as np
import importlib.util
from pathlib import Path


def load(path,name):
    spec = importlib.util.spec_from_file_location(name,Path(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


periodic = load("validation/history_branch_robustness/periodic.py","branch_history_reference")
periodic_translation,slow_history_step = periodic.periodic_translation,periodic.slow_history_step
solve = load("validation/history_force_sensitivity/model.py","existing_history_reference").solve


def test_constant_radius_periodic_reference_matches_analytic_impedance():
    radius,frequency,mu,rho = 1.78e-6,2.5e6,.001,1000.
    volume = 4*np.pi*radius**3/3
    mass = rho*volume/2
    omega = 2*np.pi*frequency
    for history in [False,True]:
        result = periodic_translation(np.full(512,radius),frequency,history=history,n=65)
        impedance = 6*np.pi*mu*radius+1j*omega*mass
        if history:
            impedance += 6*np.pi*mu*radius**2/np.sqrt(mu/rho)*np.sqrt(1j*omega)
        phasors = np.array([-volume,-1j*volume])/impedance
        expected = np.real(np.exp(2j*np.pi*result["phase"][:,None])*phasors)
        np.testing.assert_allclose(result["velocity_per_gradient"],expected,rtol=2e-10,atol=1e-24)
        assert result["relative_linear_residual"] < 1e-12


def test_slow_history_matches_independent_existing_constant_radius_reference():
    radius,dt = 1.78e-6,1e-7
    time = np.arange(201)*dt
    force = 1e-11*(1-np.cos(2*np.pi*5e4*time))
    for history in [False,True]:
        old = solve(time,np.full(len(time),radius),force,history=history)
        slip = np.zeros((len(time),1,3))
        increments = np.zeros((len(time)-1,1,3))
        for j in range(1,len(time)):
            slip[j] = slow_history_step(np.full((1,3),force[j]),slip[j-1],slip[max(j-2,0)],
                                        increments,j,dt,radius,history=history)
            increments[j-1] = slip[j]-slip[j-1]
        np.testing.assert_allclose(slip[:,0,0],old["velocity_m_s"],rtol=1e-11,atol=1e-16)


def test_variable_radius_periodic_reference_has_correct_quasistatic_mean():
    """Independent closed form exercises transformed-time weights and phase sign.

    With vanishing added mass and history OFF, v=-2 R(t)^2 grad(p)/(9 mu).
    For R=R0[1+eps*cos(wt+phi)], both mean gradient coefficients are analytic.
    This toy limit is a numerical test, not a physical low-density simulation.
    """
    radius,mu,epsilon,phase = 1.78e-6,.001,.10,.47
    theta = 2*np.pi*np.arange(1024)/1024
    # Raw APECSS clock is shifted by -T/4 inside the periodic reference.
    raw_radius = radius*(1+epsilon*np.cos(theta+phase+np.pi/2))
    result = periodic_translation(raw_radius,2.5e6,history=False,n=129,density=1e-9,viscosity=mu)
    expected = -2*radius**2*epsilon/(9*mu)*np.array([np.cos(phase),np.sin(phase)])
    np.testing.assert_allclose(result["mean_velocity_per_gradient"],expected,rtol=2e-7,atol=1e-24)
