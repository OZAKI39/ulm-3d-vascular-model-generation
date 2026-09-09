"""Pinned local call-chain evidence and fixed-snapshot formulas, never an integrator."""
import time
from pathlib import Path
import numpy as np
from py_scripts.fluid_physics.common import PROJECT_ROOT,sha256_file,read_json
from .calculations import ou_discrete_temperature_ratio,shared_reference_covariance
from py_scripts.fluid_physics.analysis import t95


def source_contract():
    entries=[
      ('src/mirheo/core/mirheo.cpp','void Mirheo::run','run sets dt, calls setup, then sim.run; the state object persists.'),
      ('src/mirheo/core/simulation.cpp','void Simulation::init','Rebuilds RunData, cell lists, interaction bindings and plugins; does not rerun UniformIC.'),
      ('src/mirheo/core/simulation.cpp','void Simulation::run','Loop begins at existing currentStep/currentTime; scheduler runs once per step; final cell-list rebuild can change summation order.'),
      ('src/mirheo/core/simulation.cpp','static void buildDependencies','Intermediate density accumulation/gather precedes final forces; exactly one integration task follows force accumulation.'),
      ('src/mirheo/core/integrators/vv.cu','p.u +=','Executed update is v += F dt/m; x += v_new dt. No second kick in this pure-liquid call chain.'),
      ('src/mirheo/core/integrators/integration_kernel.h','real4 vel =','Kernel reads existing velocity and old position, applies the transform, writes xyz and preserves packed IDs.'),
      ('src/mirheo/core/integrators/forcing_terms/periodic_poiseuille.h','operator()','Periodic force uses the y half-domain sign and the selected x component; parameter is force on a particle.'),
      ('src/mirheo/core/pvs/particle_vector.cpp','void ParticleVector::setVelocities_vector','Setter copies xyz, preserving the fourth packed ID word, with no half-step conversion or temperature rescaling.'),
      ('src/mirheo/core/utils/common.cpp','const std::string forces','The reserved force channel is named __forces, not forces. Required by the new diagnostic saver; first failed attempt is preserved.'),
      ('src/mirheo/bindings/particle_vectors.cpp','getPerParticlePositions','Getters download owned local arrays; coordinates convert local to global, velocities expose xyz only.'),
      ('src/mirheo/bindings/cuda_array_interface.cpp','toCudaArrayInterface','Shape is logical buffer size, stride[0] is C++ sizeof(T); fourth real4 component is not a physical coordinate.'),
      ('src/mirheo/core/initial_conditions/helpers.cpp','static long genSeed','Position RNG seed is compute rank + std::hash(PV name); each unit cell gets floor/ceil(n) random positions. Uniform is not a lattice.'),
      ('src/mirheo/core/initial_conditions/helpers.cpp','p.u.x =','Uniform starts at zero velocity; the frozen Python worker sets Gaussian velocities and removes COM once.'),
      ('src/mirheo/core/interactions/pairwise/kernels/density.h','operator()','Kernel number density sums W including self; SDPD separately multiplies by particle mass for EOS.'),
      ('src/mirheo/core/interactions/pairwise/kernels/sdpd.h','const real pi =','EOS receives m*d; dissipative pair coefficient is -5*mu*(di^-2+dj^-2)*Wprime/r.'),
      ('src/mirheo/core/interactions/pairwise/kernels/sdpd.h','computeFRfact','Noise amplitude is sqrt(2*5*mu*kBT/dt), multiplied once by dt in the velocity update.'),
      ('src/mirheo/core/interactions/utils/step_random_gen.cpp','real StepRandomGen::generate','Stored Mersenne Twister state advances when currentTime changes; repeated calls at the same time reuse the sample.'),
      ('src/mirheo/core/utils/cuda_rng.h','inline __D__ real mean0var1( real seed, int','Logistic pair hash sorts IDs in the SDPD caller, then produces intended zero-mean/unit-variance bounded arcsine noise; device variance not empirically validated here.'),
      ('src/mirheo/core/interactions/pairwise/kernels/stress_wrapper.h','0.5','Stress wraps the same pair force; half contribution per ordered particle pair, not a second thermostat.'),
      ('src/mirheo/plugins/particle_channel_saver.cpp','void ParticleChannelSaverPlugin::beforeIntegration','Copies density and positions at the same pre-integration phase into persistent channels.'),
      ('src/mirheo/core/pvs/data_manager.h','enum class PersistenceMode','Active saved channels follow particle redistribution/reordering. Saved positions have no shift flag, so use periodic minimum-image distances in this single periodic domain.'),
      ('src/mirheo/plugins/stats.cu','void SimulationStats::afterIntegration','Owned-particle kinetic energy and momentum sampled after integration; serialized time belongs to the following step boundary.'),
      ('src/mirheo/plugins/virial_pressure.cu','void VirialPressurePlugin::afterIntegration','Sums trace(stress)/3 without V or kinetic term. Stress was computed from pre-integration coordinates/velocities.'),
    ]
    rows=[]
    for relative,needle,meaning in entries:
        path=PROJECT_ROOT/'vendor/Mirheo'/relative;lines=path.read_text().splitlines()
        hits=[i+1 for i,line in enumerate(lines) if needle in line]
        if not hits:raise ValueError('SOURCE_ANCHOR_MISSING '+relative+' '+needle)
        rows.append({'path':str(path),'sha256':sha256_file(path),'line':hits[0],'anchor':needle,'meaning':meaning})
    return {'status':'LOCAL_CALL_CHAIN_REVIEWED','sources':rows,
      'integrator_semantics':'One kick followed by drift. For a conservative force it has the algebra of symplectic Euler / staggered leapfrog, but this API does not reconstruct centered integer-time velocities; velocity-dependent stochastic forces use stored pre-kick velocities.',
      'native_bug_status':'INSUFFICIENT_EVIDENCE: method semantics and finite-step bias are concerns, not proof of a native library defect.',
      'repeated_run_status':'RULED_OUT_WITHIN_TESTED_SCOPE: no time, particle state or SDPD generator reset in this no-object path. Cell-list order may change floating-point trajectories; exact continuous-vs-chunk equivalence has not been experimentally tested.',
      'pressure_phase':'Virial is pre-kick pair stress, kinetic momentum flux is post-kick; adjacent time phases differ by dt. The mean estimator is discrete-time mechanical pressure, not independently calibrated thermodynamic absolute pressure.',
      'unavailable_old_observations':['initial coordinates before first force','particle velocity snapshots','force snapshots','random generator checkpoints','per-component kinetic variances'],
      'no_native_library_changes':True}


def snapshot_audit(path,spec,*,initial_state=False):
    """All-pairs fixed-state force/friction/structure audit of a real saved frame."""
    start=time.monotonic();x=np.load(path,allow_pickle=False).astype(float)
    observed=None if initial_state else np.load(path.with_name(path.name.replace('_positions.npy','_kernel_number.npy')),allow_pickle=False)
    N=len(x);box=np.array(spec['domain_star']);rc=spec['rc_star'];a=spec['candidate'];mass=spec['m_star'];norm=21/(2*np.pi*rc**3)
    density=np.zeros(N);neighbors=np.zeros(N,dtype=int);nearest=np.zeros(N);friction=np.zeros((N,3,3));force=np.zeros((N,3));virial=0.
    for begin in range(0,N,128):
        dx=x[begin:begin+128,None,:]-x[None,:,:];dx-=np.rint(dx/box)*box
        r=np.linalg.norm(dx,axis=2);q=r/rc;inside=q<1
        density[begin:begin+len(r)]=np.sum(np.where(inside,norm*(1-q)**4*(1+4*q),0),axis=1)
        nonself=r>1e-12;neighbors[begin:begin+len(r)]=np.sum(inside&nonself,axis=1)
        nearest[begin:begin+len(r)]=np.min(np.where(nonself,r,np.inf),axis=1)
    p=a['sound_speed_star']**2*(mass*density-a['rho_0_star'])
    for begin in range(0,N,128):
        dx=x[begin:begin+128,None,:]-x[None,:,:];dx-=np.rint(dx/box)*box
        r=np.linalg.norm(dx,axis=2);q=r/rc;active=(q<1)&(r*r>=1e-6)
        invr=np.divide(1.,r,out=np.zeros_like(r),where=active);er=dx*invr[...,None]
        derivative=np.where(active,-20*norm*q*(1-q)**3/rc,0.)
        gamma=-5*a['viscosity_mu_star']*(density[begin:begin+len(r),None]**-2+density[None,:]**-2)*derivative*invr
        fC=-(p[begin:begin+len(r),None]/density[begin:begin+len(r),None]**2+p[None,:]/density[None,:]**2)*derivative
        force[begin:begin+len(r)]=np.sum(fC[...,None]*er,axis=1)
        virial+=float(np.sum(fC*r))/(6*np.prod(box))
        for j in range(3):
            for k in range(3):friction[begin:begin+len(r),j,k]=np.sum(gamma*er[:,:,j]*er[:,:,k],axis=1)/mass
    eigen=np.linalg.eigvalsh(friction);mean_mode=float(np.sum(np.trace(friction,axis1=1,axis2=2))/(3*(N-1)))
    bound=float(2*np.max(np.trace(friction,axis1=1,axis2=2)));dt=spec['dt_star']
    return {'snapshot':str(path),'snapshot_sha256':sha256_file(path),'phase':read_json(path.parent/'initial_state.json') if initial_state else read_json(path.with_name(path.name.replace('_positions.npy','.json'))),
      'category':'DERIVED_FROM_REAL_POSITION_SNAPSHOT','N':N,'density_relative_error_max':None if observed is None else float(np.max(abs(density-observed)/observed)),
      'density_comparison_reason':'Initial state is before first native density calculation; CPU-derived density only.' if initial_state else 'Saved pre-kick positions and densities, same persistent particle order and phase.',
      'kernel_mean':float(density.mean()),'kernel_cv':float(density.std()/density.mean()),'global_number_density':N/np.prod(box),
      'neighbors_mean':float(neighbors.mean()),'neighbors_range':[int(neighbors.min()),int(neighbors.max())],
      'nearest_distance_quantiles_star':np.quantile(nearest,[0,.1,.5,.9,1]).tolist(),
      'conservative_force_RMS_star':float(np.sqrt(np.mean(np.sum(force**2,axis=1)))),
      'conservative_virial_star':virial,'model_pressure_at_local_density_mean_star':float(p.mean()),
      'conservative_potential_per_particle_star':float(np.mean(a['sound_speed_star']**2*(mass*np.log(density)+a['rho_0_star']/density))),
      'potential_definition':'For fixed Linear EOS parameters, U=sum cs^2*(m*ln(d)+rho0/d), up to an additive constant; only differences within a fixed model are interpretable.',
      'friction_local_eigenvalue_quantiles_star':np.quantile(eigen,[0,.1,.5,.9,1]).tolist(),
      'friction_mean_mode_rate_star':mean_mode,'friction_full_operator_eigenvalue_upper_bound_star':bound,
      'dt_times_mean_mode_rate':dt*mean_mode,'dt_times_upper_bound':dt*bound,
      'fixed_OU_mean_mode_temperature_ratio':ou_discrete_temperature_ratio(mean_mode,dt) if dt*mean_mode<2 else None,
      'fixed_OU_scope':'Analytic fixed-friction OU mode only; positions, collective modes and actual native random correlations are not frozen in the fluid. This is a bias scale estimate, not a fluid-temperature correction.',
      'CPU_wall_s':time.monotonic()-start}


def pressure_audit(old_assessments,corrected_assessments,corrected_tasks,case,u):
    a=next(x for x in corrected_assessments if x['method']=='SDPD');old=next(x for x in old_assessments if x['method']=='SDPD')
    points=a['eos']['points'];points=sorted(points,key=lambda p:p['rho_star'])
    # One common random-generator seed is used in the separate initialized runs.
    # Covariance below is an explicitly conditional independent-run approximation.
    terms=[]
    for p in points:
        task=next(t for t in corrected_tasks if t['task_id']==p['task_id']);s=task['pressure'];means=np.array(s['block_means'])*u.scales['pressure']
        terms.append({'task_id':p['task_id'],'mean_pa':p['mean_pressure_pa'],'ci_halfwidth_pa':p['ci95_halfwidth_star']*u.scales['pressure'] if p['ci95_halfwidth_star'] is not None else None,
          'mean_variance_pa2':float(np.var(means,ddof=1)/len(means)) if len(means)>=2 else None,'blocks':len(means),
          'temperature_mean_star':task['temperature']['mean'],'temperature_status':task['temperature_status'],
          'actual_density_kg_m3':task['density_si'],'requested_n_star':task['parameters']['task']['n_star'],'actual_n_star':task['actual_n_star']})
    differences=[];cov=None
    if len(terms)==3 and all(t['mean_variance_pa2'] is not None for t in terms):
        low,ref,high=terms;cov=shared_reference_covariance([t['mean_variance_pa2'] for t in terms])
        for endpoint,sign in [(low,-1),(high,1)]:
            gauge=min(case['outlet_gauge_pressures_pa'].values()) if sign==-1 else max(case['outlet_gauge_pressures_pa'].values())
            va,vb=endpoint['mean_variance_pa2'],ref['mean_variance_pa2'];df=(va+vb)**2/(va*va/(endpoint['blocks']-1)+vb*vb/(ref['blocks']-1))
            hw=t95(int(df))*np.sqrt(va+vb);delta=endpoint['mean_pa']-ref['mean_pa']
            differences.append({'endpoint':endpoint['task_id'],'pressure_difference_pa':delta,'required_gauge_pa':gauge,
              'nominal_margin_pa':sign*(delta-gauge),'independent_run_Welch_95_halfwidth_pa':float(hw),
              'conservative_endpoint_plus_reference_halfwidth_pa':endpoint['ci_halfwidth_pa']+ref['ci_halfwidth_pa'],
              'nominal_margin_exceeds_conditional_Welch_uncertainty':bool(sign*(delta-gauge)>=hw)})
    baseline=next(t for t in corrected_tasks if t['task_id']=='sdpd_equilibrium');s=baseline['parameters'];vol=float(np.prod(s['domain_star']));cs=s['candidate']['sound_speed_star']
    raw=np.genfromtxt(Path(baseline['directory'])/'pressure/pv.csv',delimiter=',',names=True)
    csv_half=float(np.max(.5*10.**(np.floor(np.log10(abs(raw['pressure'])))-6))/vol*u.scales['pressure'])
    return {'old_EOS':old['eos'],'corrected_EOS':a['eos'],'points':terms,'endpoint_margins':differences,
      'shared_reference_difference_covariance_pa2':cov.tolist() if cov is not None else None,
      'covariance_assumption':'Diagonal run covariance assumes independent runs. The same native seed may induce common-random-number covariance; it was not measured. Existing endpoint interval bounds do not require this independence but are conservative, with no joint 95% claim.',
      'background_at_global_density_pa':u.to_si(cs**2*(s['m_star']*baseline['actual_n_star']-s['candidate']['rho_0_star']),'pressure'),
      'CSV_virial_rounding_max_half_quantum_pa':csv_half,'single_float_pressure_resolution_scale_pa':float(np.finfo(np.float32).eps*abs(terms[1]['mean_pa'])),
      'precision_scope':'Formatting and one-float resolution estimates, not a full GPU reduction-error bound. Raw pressure uncertainty is compared separately.',
      'formal_interpretation':'No qualified isothermal EOS: density points still have temperature/relaxation problems. No extrapolation or offset manipulation is used.',
      'mechanical_definition':'sum ordered-pair half stress trace/(3V) + m sum|v-COM|^2/(3V); add kinetic flux once, never add EOS pressure on top of the virial.'}
