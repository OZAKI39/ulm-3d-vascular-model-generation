"""Frozen-integrity mechanical equilibrium interface.

The legacy integrator and SciPy reference reuse the existing force kernel.
ADR is deliberately unavailable until Kilic & Madenci (2010) methods can be
read. No damping/mass formula is reconstructed from secondary sources.
"""
from dataclasses import dataclass
import copy
import numpy as np
from scipy.optimize import least_squares
from ..geometry import make_cloud
from ..mechanics import timestep_limit
from ..fragment_mechanics import prepare_supported_shape
from ..fragment_topology import surface_particles
from .fast_mechanics import evaluate_supported
from .flow import ManufacturedStreamingProvider,ConsistentFragmentFluid
from .resolved_flow import ResolvedStreamingFieldProvider
from .energy import potential,kinetic
from .erosion import tensile_strain,SettlingMonitor
from .erosion_runner import graph_state


class NeedPaperInput(RuntimeError):
    """Full primary method is required before any ADR implementation."""


@dataclass
class ForceState:
    internal_force_N: np.ndarray
    external_force_N: np.ndarray
    epsilon_plus: np.ndarray
    stored_energy_J: float
    minimum_J: float
    surface_traction_Pa: np.ndarray
    loaded_surface: np.ndarray

    @property
    def residual(self):return self.internal_force_N+self.external_force_N


@dataclass
class EquilibriumResult:
    x: np.ndarray
    v: np.ndarray | None
    method: str
    converged: bool
    history: list
    checks: list
    summary: dict


def make_fixture(config):
    cloud=make_cloud(config['clot'])
    cloud.fixed=cloud.X[:,2] < config['clot']['origin_m'][2]+config['clot']['fixed_base_thickness_m']-cloud.spacing*1e-10
    return cloud


class FrozenForceEvaluator:
    """Same NOSB mechanics and position-dependent traction as erosion v1."""
    def __init__(self,cloud,integrity,config,provider=None):
        self.cloud=cloud;self.g=np.array(integrity,float,copy=True);self.g.setflags(write=False)
        self.c=copy.deepcopy(config)
        if config['streaming']['oscillatory_amplitude']!=0 or config['streaming']['mean_fraction']!=1:
            raise ValueError('Only static mean loading is supported')
        self.shape=prepare_supported_shape(cloud,self.g,config['safety'])
        attached,*_=graph_state(cloud,self.g,cloud.X,np.zeros_like(cloud.X),self.shape[3],config)
        if not attached.all():raise ValueError('Current equilibrium fixture requires wall-attached topology')
        self.attached=attached
        self.exposed,_=surface_particles(cloud,self.g,config['surface'])
        if provider is None:
            if config['streaming']['provider']=='manufactured':provider=ManufacturedStreamingProvider(config)
            elif config['streaming']['provider']=='resolved_field':
                s=config['streaming'];provider=ResolvedStreamingFieldProvider(s['field_path'],s['metadata_path'],config['pipe']['viscosity_Pa_s'],cloud.X)
            else:raise ValueError('Unsupported field provider')
        if isinstance(provider,ManufacturedStreamingProvider):
            if provider.c['streaming']['oscillatory_amplitude']!=0 or provider.c['streaming']['mean_fraction']!=1:
                raise ValueError('Injected manufactured provider must be mean-only')
        elif not isinstance(provider,ResolvedStreamingFieldProvider):raise ValueError('Verified static field provider required')
        self.fluid=ConsistentFragmentFluid(config,provider)

    def __call__(self,x):
        check_positions(self.cloud,x,self.c)
        m=evaluate_supported(self.cloud,x,self.g,self.c['material'],self.c['safety'],self.shape)
        force,normals,traction=self.fluid.surface_force(self.cloud,x,self.g,self.exposed,self.attached,0.)
        eps=tensile_strain(np.linalg.norm(x[self.cloud.pairs[:,1]]-x[self.cloud.pairs[:,0]],axis=1)/self.cloud.length)
        minimum=float(m[2][self.shape[3]>0].min())
        loaded=self.exposed & self.attached & ~self.cloud.fixed & (np.linalg.norm(normals,axis=1)>0)
        return ForceState(m[0],force,eps,sum(potential(self.cloud,m)),minimum,traction,loaded)


def check_positions(cloud,x,c):
    if not np.isfinite(x).all():raise FloatingPointError('Nonfinite positions')
    if np.linalg.norm(x-cloud.X,axis=1).max()>c['safety']['maximum_displacement_m']:
        raise FloatingPointError('Displacement guard exceeded')
    if np.linalg.norm(x[:,1:3],axis=1).max()>c['pipe']['radius_m']*c['safety']['maximum_radius_ratio']:
        raise FloatingPointError('Particle left pipe cross-section')
    if not np.array_equal(x[cloud.fixed],cloud.X[cloud.fixed]):raise AssertionError('Fixed coordinates changed')


def diagnostics(cloud,x,state,c,velocity=None):
    free=~cloud.fixed;R=state.residual[free];external=state.external_force_N[free]
    norm=float(np.linalg.norm(R));ext=float(np.linalg.norm(external));ec=c['mechanical_equilibrium']
    traction=np.linalg.norm(state.surface_traction_Pa[state.loaded_surface],axis=1)
    if velocity is None:K=None;ratio=None
    else:
        K=kinetic(cloud.volume*c['clot']['density_kg_m3'],velocity)
        ratio=K/max(abs(state.stored_energy_J),ec['energy_floor_J'])
    return dict(force_residual_L2_N=norm,force_residual_RMS_N=norm/np.sqrt(int(free.sum())),
        force_residual_max_N=float(np.linalg.norm(R,axis=1).max()),external_force_L2_N=ext,
        relative_force_residual=norm/max(ext,ec['force_floor_N']),force_floor_N=ec['force_floor_N'],
        physical_kinetic_energy_J=K,kinetic_to_stored_energy_ratio=ratio,adr_pseudo_kinetic_measure=None,
        maximum_epsilon_plus=float(state.epsilon_plus.max()),mean_epsilon_plus=float(state.epsilon_plus.mean()),
        stored_energy_J=float(state.stored_energy_J),minimum_J=state.minimum_J,
        maximum_displacement_m=float(np.linalg.norm(x-cloud.X,axis=1).max()),
        fixed_max_error_m=float(np.max(abs(x[cloud.fixed]-cloud.X[cloud.fixed]))),
        surface_traction_mean_Pa=float(traction.mean()) if len(traction) else 0.,
        surface_traction_max_Pa=float(traction.max()) if len(traction) else 0.,particle_count=len(x),
        erosion_exposure_time_s=0.,maximum_damage=0.)


def validate_options(cloud,x0,v0,integrity,c,method,verification_only):
    ec=c['mechanical_equilibrium']
    if method=='adaptive_dynamic_relaxation':
        raise NeedPaperInput('NEED_PAPER_INPUT: Kilic & Madenci (2010), DOI 10.1016/j.tafmec.2010.08.001; ADR formulas not implemented')
    if method not in ['legacy_fixed_damping','scipy_reference']:raise ValueError('Unknown equilibrium method')
    if not cloud.fixed.any() or cloud.fixed.all():raise ValueError('Both fixed and free particles required')
    if x0.shape!=cloud.X.shape or v0.shape!=cloud.X.shape:raise ValueError('Invalid state shape')
    if np.asarray(integrity).shape!=(len(cloud.pairs),) or np.any(np.asarray(integrity)!=1.):
        raise ValueError('This first independent benchmark permits D=0 only; erosion integration remains gated')
    for key in ['relative_force_residual_tolerance','force_floor_N','driver_drift_tolerance','driver_absolute_floor','kinetic_to_stored_energy_tolerance','energy_floor_J']:
        if not np.isfinite(ec[key]) or ec[key]<=0:raise ValueError(f'Positive finite {key} required')
    for key in ['check_interval','minimum_consecutive_converged_checks','maximum_iterations']:
        if isinstance(ec[key],bool) or ec[key]<1 or int(ec[key])!=ec[key]:raise ValueError(f'Positive integer {key} required')
    if method=='scipy_reference':
        if not verification_only:raise ValueError('SciPy reference is restricted to explicit verification_only=True')
        if 3*int((~cloud.fixed).sum())>ec['scipy_reference_max_free_dofs']:raise ValueError('SciPy reference is restricted to small cases')
    check_positions(cloud,x0,c)


def equilibrate(cloud,x0,v0,integrity,config,force_callback=None,method=None,verification_only=False):
    """Never update integrity, topology, D or the erosion exposure clock.

    legacy_fixed_damping integrates the entire prescribed window even if
    equilibrium was reached earlier, for a controlled window sensitivity test.
    SciPy's algebraic root has no velocity or mechanical/pseudo-time history.
    """
    c=config;ec=c['mechanical_equilibrium'];method=method or ec['method']
    validate_options(cloud,x0,v0,integrity,c,method,verification_only)
    x=np.array(x0,copy=True);v=np.array(v0,copy=True);g=np.array(integrity,copy=True);g.setflags(write=False)
    evaluate=force_callback or FrozenForceEvaluator(cloud,g,c)
    if method=='scipy_reference':return _scipy(cloud,x,evaluate,c)
    dt=c['simulation']['dt_s'];window=ec['legacy_window_s'];n=round(window/dt)
    if n<1 or not np.isclose(n*dt,window,rtol=0,atol=1e-14):raise ValueError('dt must divide fixed window')
    if dt>timestep_limit(cloud,c):raise ValueError('Unsafe mechanical time step')
    if n>ec['maximum_iterations']:raise ValueError('Fixed window exceeds predeclared maximum_iterations')
    mass=cloud.volume*c['clot']['density_kg_m3'];decay=np.exp(-c['simulation']['damping_s_inv']*dt/2)
    state=evaluate(x);acc=state.residual/mass[:,None];acc[cloud.fixed]=0.;v[cloud.fixed]=0.
    history=[];checks=[];consecutive=0;first_pass=None
    monitor=SettlingMonitor(g>0,ec['driver_absolute_floor']);monitor.add(0.,0.,state.epsilon_plus,kinetic(mass,v))
    initial=diagnostics(cloud,x,state,c,v)
    initial.update(equilibrium_iteration=0,mechanical_time_s=0.,equilibrium_pseudo_time=None,phase='NO_DAMAGE_WINDOW',driver_drift=None)
    history.append(initial)
    for iteration in range(1,n+1):
        v*=decay;v+=.5*dt*acc;x+=dt*v;x[cloud.fixed]=cloud.X[cloud.fixed];v[cloud.fixed]=0.
        state=evaluate(x);acc=state.residual/mass[:,None];acc[cloud.fixed]=0.
        if np.linalg.norm(acc,axis=1).max()>c['safety']['maximum_acceleration_m_s2']:raise FloatingPointError('Acceleration guard exceeded')
        v+=.5*dt*acc;v*=decay;v[cloud.fixed]=0.
        row=diagnostics(cloud,x,state,c,v)
        row.update(equilibrium_iteration=iteration,mechanical_time_s=iteration*dt,equilibrium_pseudo_time=None,phase='NO_DAMAGE_WINDOW',driver_drift=None)
        monitor.add(iteration*dt,iteration*dt,state.epsilon_plus,row['physical_kinetic_energy_J'])
        if iteration%ec['check_interval']==0:
            drift=monitor.summary(ec['driver_drift_tolerance'])['driver_relative_drift']
            okay=(row['relative_force_residual']<=ec['relative_force_residual_tolerance'] and
                  row['kinetic_to_stored_energy_ratio']<=ec['kinetic_to_stored_energy_tolerance'] and
                  drift<=ec['driver_drift_tolerance'])
            consecutive=consecutive+1 if okay else 0
            row['driver_drift']=drift
            checks.append(dict(iteration=iteration,mechanical_time_s=iteration*dt,driver_drift=drift,
                relative_force_residual=row['relative_force_residual'],kinetic_ratio=row['kinetic_to_stored_energy_ratio'],
                all_criteria_passed=bool(okay),consecutive_passes=consecutive))
            if consecutive>=ec['minimum_consecutive_converged_checks'] and first_pass is None:first_pass=iteration
            monitor=SettlingMonitor(g>0,ec['driver_absolute_floor'])
            monitor.add(iteration*dt,iteration*dt,state.epsilon_plus,row['physical_kinetic_energy_J'])
        history.append(row)
    final_checked=bool(checks and checks[-1]['iteration']==n)
    passed=bool(final_checked and consecutive>=ec['minimum_consecutive_converged_checks'])
    summary=dict(status='EQUILIBRIUM_CRITERIA_PASSED' if passed else 'INCOMPLETE_MECHANICAL_EQUILIBRIUM',
        method=method,converged=passed,iterations=n,mechanical_time_s=n*dt,equilibrium_pseudo_time=None,
        first_converged_iteration=first_pass,erosion_exposure_time_s=0.,damage_advanced=False,
        fixed_topology=True,physical_density_unchanged=True,final=history[-1],
        convergence_scope='Force, physical-mass kinetic ratio and driver drift; no independent ADR verification')
    return EquilibriumResult(x,v,method,passed,history,checks,summary)


def _scipy(cloud,x,evaluate,c):
    ec=c['mechanical_equilibrium'];free=~cloud.fixed;h=cloud.spacing;origin=cloud.X[free].copy()
    initial=evaluate(x);force_scale=max(np.linalg.norm(initial.external_force_N[free]),ec['force_floor_N'])
    history=[];evaluations=0
    def unpack(q):
        xx=cloud.X.copy();xx[free]=origin+h*q.reshape((-1,3));return xx
    def fun(q):
        nonlocal evaluations
        evaluations+=1;xx=unpack(q);state=evaluate(xx)
        return state.residual[free].ravel()/force_scale
    def capture(q):
        xx=unpack(q);state=evaluate(xx);row=diagnostics(cloud,xx,state,c)
        row.update(reference_iteration=len(history),residual_evaluations=evaluations,
                   mechanical_time_s=None,equilibrium_pseudo_time=None,phase='SCIPY_REFERENCE',driver_drift=None)
        history.append(row)
    def callback(intermediate_result):capture(intermediate_result.x)
    q0=((x[free]-origin)/h).ravel();capture(q0)
    solution=least_squares(fun,q0,callback=callback,**ec['scipy'])
    final_x=unpack(solution.x);capture(solution.x);final=history[-1]
    passed=bool(solution.success and final['relative_force_residual']<=ec['relative_force_residual_tolerance'])
    summary=dict(status='SCIPY_REFERENCE_FORCE_BALANCE_PASSED' if passed else 'SCIPY_REFERENCE_NOT_AVAILABLE',method='scipy_reference',
        converged=passed,scipy_success=bool(solution.success),scipy_status=int(solution.status),scipy_termination_reason=solution.message,
        scipy_nfev=int(solution.nfev),actual_residual_evaluations=evaluations,accepted_iterations=max(0,len(history)-2),
        scipy_cost=float(solution.cost),scipy_optimality=float(solution.optimality),
        final_jacobian_condition=float(np.linalg.cond(solution.jac)),force_scale_N=float(force_scale),
        coordinate_scale_m=h,free_coordinate_count=len(q0),erosion_exposure_time_s=0.,damage_advanced=False,
        fixed_topology=True,physical_density_unchanged=True,final=final,
        convergence_scope='Small algebraic equilibrium reference only; no physical velocity, time or ADR iterations',
        mechanical_time_s=None,equilibrium_pseudo_time=None)
    return EquilibriumResult(final_x,None,'scipy_reference',passed,history,[],summary)
