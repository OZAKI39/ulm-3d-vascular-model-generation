"""Optional equivalent-sphere fragment drag, not irregular-clot hydrodynamics.

Schiller–Naumann correlation is documented by the OpenFOAM Foundation:
https://cpp.openfoam.org/v12/SchillerNaumann_8C_source.html
The low-Re limit is Stokes drag. No torque approximation is added in this version.
"""
import numpy as np
from ..fragment_fluid import relaxation_half_step


def drag_coefficient(mass, volume, slip, rho, mu):
    if mass <= 0 or volume <= 0 or rho <= 0 or mu <= 0:
        raise ValueError('Positive fragment mass, volume, density and viscosity required')
    diameter = (6*volume/np.pi)**(1/3)
    speed = float(np.linalg.norm(slip))
    Re = rho*diameter*speed/mu
    if Re < 1000:
        beta = 3*np.pi*mu*diameter*(1+.15*Re**.687)
        regime = 'STOKES_LIMIT_WITH_SCHILLER_NAUMANN_CORRECTION' if Re < .1 else 'SCHILLER_NAUMANN'
    else:
        beta = .5*rho*(np.pi*diameter**2/4)*.44*speed
        regime = 'CONSTANT_CD_0.44_HIGH_RE_SURROGATE'
    return beta, diameter, Re, regime


def fragment_drag_half_step(cloud, positions, velocity, mass, components, category, provider, time_s, dt, config):
    if dt <= 0:
        raise ValueError('Positive transport increment required')
    result = velocity.copy()
    debris = (category >= 2) & ~cloud.fixed
    uf = provider.sample(positions, time_s).velocity
    result = relaxation_half_step(result, uf, debris, dt, config['transport']['tau_h_s'])
    rows = []
    for comp in components:
        ids = np.asarray(comp['particle_ids'], int)
        if comp['attached_to_base'] or not np.all(category[ids] == 1):
            continue
        if np.any(cloud.fixed[ids]):
            raise AssertionError('Fixed point in a detached resolved component')
        M, V = float(mass[ids].sum()), float(cloud.volume[ids].sum())
        com = np.average(positions[ids], axis=0, weights=mass[ids])
        vc = np.average(velocity[ids], axis=0, weights=mass[ids])
        uc = np.average(uf[ids], axis=0, weights=cloud.volume[ids])
        slip = uc-vc
        beta, diameter, Re, regime = drag_coefficient(M, V, slip, config['pipe']['density_kg_m3'], config['pipe']['viscosity_Pa_s'])
        delta_v = slip*(-np.expm1(-beta*dt/M))
        result[ids] += delta_v
        net_impulse = np.sum(mass[ids, None]*(result[ids]-velocity[ids]), axis=0)
        if not np.allclose(net_impulse, M*delta_v, rtol=1e-12, atol=1e-25):
            raise AssertionError('Fragment drag distribution lost momentum')
        rows.append(dict(fragment_id=comp['fragment_id'], particle_count=len(ids), volume_m3=V, mass_kg=M,
                         center_of_mass_m=com.tolist(), center_velocity_m_s=vc.tolist(), sampled_fluid_velocity_m_s=uc.tolist(),
                         equivalent_diameter_m=diameter, projected_area_m2=np.pi*diameter**2/4,
                         Reynolds_number=Re, drag_regime=regime, linear_drag_coefficient_kg_s=beta,
                         instantaneous_drag_N=(beta*slip).tolist(), integrated_drag_impulse_kg_m_s=net_impulse.tolist(),
                         sample_method='Volume-weighted velocity at all component particle positions',
                         torque_model='NOT_INCLUDED', integration='Exact frozen-coefficient COM relaxation over this half-step'))
    return result, rows
