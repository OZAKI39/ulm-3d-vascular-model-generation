"""Validation-only periodic variable-radius translation in transformed time.

No production import is permitted. Garbin (2009), Eq.1 and Cattaneo & Supponen
arXiv:2506.16413v2, Eqs.11,14 specify conservative added momentum and the Abel
kernel in theta=int R^-2 dt. For an infinitely periodic input the Abel operator
has Fourier multiplier sqrt(pi)*sqrt(i*omega_theta). The zero mode is zero.
This is an adiabatic local periodic closure, not a full moving microcycle solve.
"""
import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.interpolate import CubicSpline, PchipInterpolator


def periodic_translation(radius_cycle_m,frequency_Hz,*,history=True,n=129,density=1000.,viscosity=.001):
    if n < 17 or n % 2 == 0:
        raise ValueError("Odd Fourier size >=17 avoids an ambiguous Nyquist derivative")
    r = np.asarray(radius_cycle_m,float)
    if np.any(r<=0) or not np.isfinite(r).all():
        raise ValueError("Positive finite prescribed radii required")
    cycle = CubicSpline(np.linspace(0,1,len(r)+1),np.r_[r,r[0]],bc_type="periodic")
    # APECSS convention -P*sin(wt); shift by -T/4 for the real-cosine pressure reference.
    dense_phase = np.linspace(0,1,8193)
    dense_r = cycle(np.mod(dense_phase-.25,1.))
    period = 1/frequency_Hz
    theta = cumulative_trapezoid(1/dense_r**2,dense_phase*period,initial=0.)
    theta_period = theta[-1]
    uniform = np.arange(n)*theta_period/n
    phase = PchipInterpolator(theta,dense_phase)(uniform)
    radius = cycle(np.mod(phase-.25,1.))
    mass = .5*density*4*np.pi*radius**3/3
    omega = 2*np.pi*np.fft.fftfreq(n,d=theta_period/n)
    basis = np.fft.fft(np.eye(n),axis=0)
    derivative = np.fft.ifft((1j*omega)[:,None]*basis,axis=0).real
    half = np.fft.ifft(np.sqrt(1j*omega)[:,None]*basis,axis=0).real
    operator = derivative*(mass/radius)[None,:]/radius[:,None]**2 + 6*np.pi*viscosity*np.eye(n)
    if history:
        operator += 6*np.pi*viscosity/np.sqrt(viscosity/density)*half
    volume = 4*np.pi*radius**3/3
    force_basis = volume[:,None]*np.column_stack([-np.cos(2*np.pi*phase),np.sin(2*np.pi*phase)])
    q = np.linalg.solve(operator,force_basis)
    velocity = q/radius[:,None]
    mean = (velocity*radius[:,None]**2).sum(axis=0)/(radius**2).sum()
    residual = np.linalg.norm(operator@q-force_basis)/np.linalg.norm(force_basis)
    return dict(mean_velocity_per_gradient=mean,relative_linear_residual=float(residual),
                mean_force_per_gradient=(force_basis*radius[:,None]**2).sum(axis=0)/(radius**2).sum(),
                theta=uniform,phase=phase,radius_m=radius,velocity_per_gradient=velocity,
                force_per_gradient=force_basis,operator=operator)


def slow_history_step(force,slip,slip_previous,increments,step,dt,radius,*,history=True,density=1000.,viscosity=.001):
    """BDF2 added mass and exact piecewise-linear Abel weights for slow slip.

    Radius is equilibrium R0 here. Fast radial/history effects belong to the
    local periodic closure, so adding this slow-envelope kernel is a two-scale
    approximation. No finite-memory cutoff or empirical multiplier is used.
    """
    mass = .5*density*4*np.pi*radius**3/3
    drag = 6*np.pi*viscosity*radius
    coefficient = 6*np.pi*viscosity*radius**2/np.sqrt(np.pi*viscosity/density)
    alpha = 1. if step==1 else 1.5
    rhs = mass*(slip if step==1 else 2*slip-.5*slip_previous)/dt
    past,last = np.zeros_like(slip),0.
    if history:
        lag = np.arange(step,0,-1,dtype=float)
        weights = 2/np.sqrt(dt)/(np.sqrt(lag)+np.sqrt(lag-1))
        past = np.tensordot(weights[:-1],increments[:step-1],axes=(0,0))
        last = weights[-1]
    return (force+rhs-coefficient*past+coefficient*last*slip)/(alpha*mass/dt+drag+coefficient*last)
