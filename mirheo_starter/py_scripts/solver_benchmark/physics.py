"""One SI definition, native frozen Mirheo mapping, independent continuum target."""
import math
from pathlib import Path
import numpy as np
import yaml
from py_scripts.fluid_physics.common import PROJECT_ROOT, sha256_file
from py_scripts.fluid_physics.analysis import periodic_shape, block_statistics, t95
from py_scripts.fluid_comparison.models import frozen_targets, candidates


def load_config(path):
    path = Path(path)
    if not path.is_absolute(): path = PROJECT_ROOT / path
    c = yaml.safe_load(path.read_text())
    validate_config(c)
    c['_config_path'] = str(path)
    c['_config_sha256'] = sha256_file(path)
    return c


def validate_config(c):
    if c['campaign_id']!='solver_benchmark_20260909':raise ValueError('FIXED_CAMPAIGN_NO_BUDGET_RESET')
    b=c['budget']
    if not 0<b['cpu_limit_s']<=1200 or b['task_limit_s']!=600 or b['concurrency']!=1 or b['gpu_requested_s']!=600:
        raise ValueError('BUDGET_OR_CONCURRENCY_OUT_OF_SCOPE')
    if b['automatic_retry'] is not False:raise ValueError('AUTOMATIC_RETRY_NOT_AUTHORIZED')
    if c['physics']['acceleration_si']<=0:raise ValueError('INVALID_ACCELERATION')
    if c['physics']['periodic']!=[True,True,True]:raise ValueError('PERIODIC_CASE_REQUIRED')
    if c['sampling']['bins']%2:raise ValueError('HALF_BOX_BINS_REQUIRED')
    if any(x not in [1,2,4] for x in c['lbm']['ranks']):raise ValueError('RANK_SCAN_OUT_OF_SCOPE')


def definition(c):
    p = PROJECT_ROOT / c['frozen_model_config']
    old = yaml.safe_load(p.read_text())
    case, mapping, units = frozen_targets(old)
    candidate = next(x for x in candidates(old, mapping) if x['id'] == 'sdpd_target_linear')
    domain_star = old['sampling']['domain_star']
    box = units.to_si(domain_star, 'length').tolist()
    a = c['physics']['acceleration_si']
    rho, nu = case['density_kg_m3'], case['kinematic_viscosity_m2_s']
    force = units.to_star(units.to_si(candidate['m_star'], 'particle_mass') * a, 'force')
    L = box[1]
    p = dict(rho_si=rho, nu_si=nu, mu_si=rho*nu, temperature_K=case['temperature_K'],
             box_si=box, acceleration_si=a, force_density_si=rho*a,
             particle_force_star=force, particle_force_si=units.to_si(force, 'force'),
             domain_star=domain_star, bins=c['sampling']['bins'],
             bin_width_si=L/c['sampling']['bins'], relaxation_time_si=L**2/(4*math.pi**2*nu),
             expected_max_velocity_si=a*L**2/(32*nu), expected_half_mean_si=a*L**2/(48*nu),
             background_temperature_scope='LBM material background only; thermal fluctuations NOT_APPLICABLE',
             force_sign='y > Ly/2: +x; y <= Ly/2: -x; all directions periodic',
             frozen_config_sha256=sha256_file(PROJECT_ROOT/c['frozen_model_config']),
             mapping=mapping, candidate=candidate)
    p['expected_Re_box'] = p['expected_max_velocity_si']*L/nu
    p['expected_SDPD_Mach'] = p['expected_max_velocity_si']/units.to_si(candidate['sound_speed_star'], 'velocity')
    return p, units


def lattice_definition(p, dx, tau=1):
    N = round(p['box_si'][1]/dx)
    if N%2 or N%p['bins'] or any(not math.isclose(L, N*dx, rel_tol=1e-12) for L in p['box_si']):
        raise ValueError('PERIODIC_GRID_OR_BINS_MISMATCH')
    dt = (tau-.5)*dx*dx/(3*p['nu_si'])
    if dt <= 0: raise ValueError('INVALID_TAU')
    return dict(N=N, dx_si=dx, dt_si=dt, tau=tau, acceleration_lbm=p['acceleration_si']*dt*dt/dx,
                expected_Mach=p['expected_max_velocity_si']*dt/dx*math.sqrt(3),
                lattice_nodes=N**3, node_origin_si=dx/2, periodic_unique_nodes=N)


def target(p):
    width=p['bin_width_si']; y=(np.arange(p['bins'])+.5)*width
    return y, p['acceleration_si']/p['nu_si']*periodic_shape(y,p['box_si'][1],width)


def profile_metrics(p, measured):
    y, truth=target(p); v=np.asarray(measured,float)
    if v.shape!=truth.shape or not np.isfinite(v).all(): raise ValueError('INVALID_PROFILE')
    # Scale shape before fitting: SI quadratics are ~1e-12; avoid rank loss.
    slope, offset=np.linalg.lstsq(np.column_stack([truth,np.ones(len(truth))]),v,rcond=None)[0]
    fit=slope*truth+offset
    h=len(v)//2
    means=[float(v[:h].mean()),float(v[h:].mean())]
    wanted=[float(truth[:h].mean()),float(truth[h:].mean())]
    area=p['box_si'][1]*p['box_si'][2]/2
    return dict(profile_relative_l2=float(np.linalg.norm(v-truth)/np.linalg.norm(truth)),
                half_mean_velocity_si=means, half_flow_si=[x*area for x in means],
                full_signed_flow_si=float(v.mean()*2*area),
                half_flow_relative_error=max(abs(means[i]/wanted[i]-1) for i in range(2)),
                apparent_nu_si=float(p['nu_si']/slope) if slope>0 else None,
                apparent_nu_relative_error=float(abs(1/slope-1)) if slope>0 else None,
                fitted_profile_si=fit.tolist(), fitted_offset_si=float(offset),
                target_profile_si=truth.tolist(), measured_profile_si=v.tolist(), y_si=y.tolist())


def qualified_ratio(a,b):
    if not a.get('qualified') or not b.get('qualified'): return None
    x,y=a.get('time_to_qualified_solution_s'),b.get('time_to_qualified_solution_s')
    return x/y if x is not None and y is not None and x>0 and y>0 else None


def valid_blocks(t,x,p,c):
    out=block_statistics(t,x,min_duration=c['criteria']['block_relaxation_multiples']*p['relaxation_time_si'],min_blocks=c['criteria']['min_blocks'])
    # Old helper deliberately retains descriptive CI; this comparison only publishes valid CI.
    if out['status']!='SUFFICIENT': out['ci95_halfwidth']=None
    # The shared helper's historical keys say star; this adapter supplied SI seconds.
    return {k.replace('_star','_si'):v for k,v in out.items()}
