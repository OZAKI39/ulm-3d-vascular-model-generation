"""CPU-only model contracts and independent small-snapshot checks."""
from pathlib import Path
import math
import numpy as np
from py_scripts.fluid_physics.common import PROJECT_ROOT, read_json, sha256_file
from py_scripts.fluid_physics.units import Units, consistency
from py_scripts.fluid_physics.calibration import verify_package


def frozen_targets(c):
    p=Path(c['historical_result']);verify_package(p)
    if sha256_file(p/'package_sha256.json')!=c['historical_manifest_sha256']:
        raise ValueError('HISTORICAL_PACKAGE_IDENTITY_CONFLICT')
    case=read_json(p/'physical_case.json');mapping=read_json(p/'unit_mapping.json')
    u=Units(mapping['L0'],mapping['M0'],mapping['t0'])
    check=consistency(u,n_star=mapping['n_star'],m_star=mapping['m_star'],kBT_star=mapping['kBT_star'],
                      rho_si=case['density_kg_m3'],temperature_K=case['temperature_K'],
                      nu_star=mapping['required_nu_star'],nu_si=case['kinematic_viscosity_m2_s'])
    if not check['all_consistent']:raise ValueError('LOCKED_UNIT_CONFLICT')
    if not math.isclose(mapping['required_mu_star'],mapping['required_nu_star']*mapping['rho_star'],rel_tol=1e-13):
        raise ValueError('DYNAMIC_KINEMATIC_VISCOSITY_CONFLICT')
    if not math.isclose(u.to_star(case['dynamic_viscosity_pa_s'],'dynamic_viscosity'),mapping['required_mu_star'],rel_tol=1e-13):
        raise ValueError('LOCKED_MU_CONFLICT')
    return case,mapping,u


def validate_model(candidate):
    if candidate['method'] not in ('DPD','SDPD'):raise ValueError('UNKNOWN_METHOD')
    for key in ['m_star','n_star','kBT_star','rc_star','dt_star','nu_prior_star']:
        if not math.isfinite(candidate[key]) or candidate[key]<=0:raise ValueError('INVALID '+key)
    if candidate['method']=='SDPD':
        den=candidate.get('density_interaction')
        if not den:raise ValueError('SDPD_REQUIRES_DENSITY')
        if den['kernel']!='WendlandC2' or den['rc_star']!=candidate['rc_star']:
            raise ValueError('DENSITY_KERNEL_OR_CUTOFF_MISMATCH')
        if candidate['EOS']!='Linear' or candidate['viscosity_mu_star']<=0 or candidate['sound_speed_star']<=0:
            raise ValueError('INVALID_SDPD_EOS_OR_VISCOSITY')
        if candidate.get('also_apply_dpd',False):raise ValueError('DUPLICATED_THERMOSTAT')
    return candidate


def candidates(c, mapping):
    shared={k:mapping[k] for k in ['m_star','n_star','kBT_star']};shared['rc_star']=c['space']['rc_star']
    rows=[]
    for raw in c['candidates']:
        a={**shared,**raw}
        if a['method']=='SDPD':
            a['viscosity_mu_star']=mapping['required_mu_star']
            a['nu_prior_star']=mapping['required_nu_star']
        validate_model(a);rows.append(a)
    return rows


def linear_eos(mass_density,sound_speed,rho_0):
    return sound_speed**2*(np.asarray(mass_density,float)-rho_0)


def snapshot_density(positions,domain,rc):
    """Periodic Wendland C2 sum including W(0); CPU audit only, never dynamics."""
    x=np.asarray(positions,float);box=np.asarray(domain,float)
    if x.ndim!=2 or x.shape[1]!=3 or not np.isfinite(x).all() or np.any(box<=2*rc):
        raise ValueError('INVALID_SMALL_PERIODIC_SNAPSHOT')
    d=np.empty(len(x));norm=21/(2*np.pi*rc**3)
    for start in range(0,len(x),96):
        dr=x[start:start+96,None,:]-x[None,:,:];dr-=box*np.rint(dr/box)
        q=np.linalg.norm(dr,axis=2)/rc;z=np.maximum(0,1-q)
        d[start:start+96]=np.sum(norm*z**4*(1+4*q),axis=1)
    return d


def native_contract():
    vendor=PROJECT_ROOT/'vendor/Mirheo'
    specs={
      'tests/sdpd/rest.py':['kind="Density"','kind="SDPD"'],
      'tests/sdpd/double_poiseuille.py':['viscosity=10.0'],
      'tests/sdpd/rigid.py':['createParticleChannelSaver'],
      'src/mirheo/core/interactions/pairwise/kernels/sdpd.h':['pressure_(di * dst.m)','viscosity * zeta_','computeFRfact'],
      'src/mirheo/core/interactions/pairwise/kernels/density.h':['number density','densityKernel_(rij, invrc_)'],
      'src/mirheo/core/interactions/pairwise/kernels/density_kernels.h':['WendlandC2DensityKernel','21.0 / (2.0 * M_PI)'],
      'src/mirheo/core/interactions/pairwise/kernels/type_traits.h':['needSelfInteraction<WendlandC2DensityKernel>'],
      'src/mirheo/core/interactions/pairwise/drivers.h':['if (needSelfInteraction<Interaction>::value)'],
      'src/mirheo/core/interactions/pairwise/kernels/pressure_EOS.h':['cSq_ * (rho - rho0_)','p0_ * (r7 - 1._r)'],
      'src/mirheo/core/interactions/pairwise/kernels/parameters.h':['SDPDParams'],
      'src/mirheo/core/interactions/pairwise/factory.cpp':['makePairwiseSDPDInteraction'],
      'src/mirheo/core/interactions/pairwise/density.cu':['Stage::Intermediate'],
      'src/mirheo/core/interactions/pairwise/sdpd.cu':['PairwiseStressWrapper','channel_names::densities'],
      'src/mirheo/bindings/interactions.cpp':['same density kernel'],
      'src/mirheo/bindings/data_manager.cpp':['__getitem__'],
      'src/mirheo/bindings/cuda_array_interface.cpp':['__cuda_array_interface__'],
      'src/mirheo/bindings/particle_vectors.cpp':['downloadFromDevice'],
      'src/mirheo/plugins/particle_channel_saver.cpp':['beforeIntegration','PersistenceMode::Active'],
      'src/mirheo/plugins/stats.cu':['afterIntegration','serializeAndSend'],
      'src/mirheo/plugins/virial_pressure.cu':['afterIntegration','savedTime_','(s.xx + s.yy + s.zz) / 3.0'],
      'src/mirheo/core/interactions/pairwise/kernels/stress_wrapper.h':['0.5_r * dr.x * f.x'],
      'src/mirheo/core/simulation.cpp':['gatherInteractionIntermediate','currentTime +=','pluginsBeforeIntegration'],
      'src/mirheo/core/integrators/forcing_terms/periodic_poiseuille.h':['return ef + original','magnitude_']}
    evidence=[]
    for path,needles in specs.items():
        p=vendor/path;lines=p.read_text().splitlines();matches={n:[i+1 for i,line in enumerate(lines) if n in line] for n in needles}
        if not all(matches.values()):raise ValueError('NATIVE_SOURCE_CONTRACT_NOT_FOUND '+path+' '+str(matches))
        evidence.append({'path':str(p),'sha256':sha256_file(p),'evidence_lines':matches})
    return {'status':'SOURCE_VERIFIED_RUNTIME_PROBE_REQUIRED','source_files':evidence,
      'registration':'Register separate Density and SDPD interactions and attach both to the same liquid PV pair; no DPD force on SDPD pairs.',
      'density_kernel':'WendlandC2: W=21/(2*pi*rc^3)*(1-q)^4*(1+4*q), 0<=q<=1; matching rc in Density and SDPD.',
      'self_contribution':'Included once by needSelfInteraction<WendlandC2DensityKernel> and computeSelfInteractions.',
      'density_stage':'Intermediate sum is cleared/recomputed each force evaluation, accumulated and gathered before SDPD final forces.',
      'density_definitions':{'global_number':'N/V','global_mass':'N*m/V','local_number':'sum_j W(rij), includes self','local_EOS_mass':'m_i * d_i; not forced equal to N*m/V'},
      'parameter_units':{'viscosity':'dynamic viscosity M/(L*T), input exact mu_star','kBT':'energy M*L^2/T^2','sound_speed':'L/T','rho_0':'mass density M/L^3; EOS parameter, not automatically initialization density'},
      'EOS':{'Linear':'p_i=c_s^2*(m_i*d_i-rho_0)','QuasiIncompressible':'p_i=p0*((m_i*d_i/rho_r)^7-1)'},
      'stress':'SDPD supports stress=True and stress_period through PairwiseStressWrapper; 0.5*r_ab*F_ab for each particle.',
      'random_force':'sqrt(2*5*mu*kBT/dt)*sqrt(-(d_i^-2+d_j^-2)*Wprime/r)*unit-variance pair noise; no extra dt factor or thermostat.',
      'mechanical_pressure':'[sum_i trace(native pair stress_i)/3 + sum_i m*|v_i-u_COM|^2/3]/V in equilibrium. Native stress includes conservative EOS, dissipative and random forces, no kinetic term.',
      'thermal_accounting':'EOS p_i is a force-law field. It is not added again to native virial pressure. Kinetic momentum flux is independently included once; mechanical pressure need not equal prescribed EOS at coarse resolution.',
      'fields':'PV coordinates and velocities via synchronized getters; local per_particle channel CUDA array interface. Save density and corresponding pre-integration positions with native ParticleChannelSaver; copy to CPU using existing CUDA runtime, no new packages.',
      'time_labels':'Native Stats and virial share afterIntegration hook. Stats serialized label advances by dt; virial saves pre-increment time. Check actual CSV offset separately per run. Saved density/positions at beforeIntegration are labelled (completed_steps-1)*dt; end-of-chunk velocities at completed_steps*dt.',
      'remaining_limitations':['No independent bulk-viscosity tuning/validation','No wall kernel-truncation validation','No membrane inside/outside coupling validation','No open pressure-reservoir implementation']}
