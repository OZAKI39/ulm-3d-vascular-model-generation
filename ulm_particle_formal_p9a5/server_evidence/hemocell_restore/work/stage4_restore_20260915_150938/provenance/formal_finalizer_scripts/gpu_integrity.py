"""Independent saved-field/CPU TestB comparisons; no solver state writes."""
from pathlib import Path
import gzip,json,math
import numpy as np
from remote_common import *
from verify_run import reconstruct_flux

def load_field(path):
    path=Path(path)
    if path.exists():return fields(path)
    raw=gzip.decompress(Path(str(path)+'.gz').read_bytes())
    n=int.from_bytes(raw[:8],'little');a=np.frombuffer(raw,FIELD_DTYPE,offset=8)
    assert len(a)==n and len(raw)==8+40*n
    return a

def compare_short(root,run,reference=None):
    root=Path(root);run=Path(run);c=contract(root);t=read(root/'contracts/GPU_NUMERICAL_COMPARISON_CONTRACT.json')
    ref=Path(reference) if reference else root/'reference/test_B'
    a=load_field(field_path(run,5000));b=load_field(field_path(ref,5000))
    assert len(a)==len(b)==182694 and np.array_equal(a['index'],b['index']),'TestB field domain mismatch'
    rho_err=float(np.linalg.norm(a['rho']-b['rho'])/np.linalg.norm(b['rho']))
    vel_err=float(np.linalg.norm(a['u']-b['u'])/np.linalg.norm(b['u']))
    ah=history(run);bh=history(ref);sa=ah[ah['iteration']==5000][0];sb=bh[bh['iteration']==5000][0]
    errors={}
    for n in ['rho_min','rho_max','rho_mean']:errors[n]=float(abs(sa[n]-sb[n]));assert errors[n]<=t['fields']['rho_min_max_abs_difference_max'],n
    for n in ['total_mass','control_volume_mass']:errors[n]=float(abs(sa[n]-sb[n])/abs(sb[n]));assert errors[n]<=t['mass']['relative_difference_max'],n
    for n in ['u_max_m_s','mean_speed_m_s']:errors[n]=float(abs(sa[n]-sb[n])/abs(sb[n]));assert errors[n]<=t['fields']['velocity_max_normalized_abs_difference_max'],n
    qnames=[f'{p}_{k}dx'+('_dx2' if v else '') for p in PORTS for k in OFFSETS for v in [0,1]]
    flux_err=max(float(abs(sa[n]-sb[n])/c['Qtarget_m3_s']) for n in qnames)
    assert rho_err<=t['fields']['rho_normalized_l2_max'] and vel_err<=t['fields']['velocity_normalized_l2_max']
    assert flux_err<=t['flux']['abs_Q_gpu_minus_cpu_over_Qtarget_max']
    return dict(status='PASS',iteration=5000,reference=str(ref),reference_field_sha256=sha(field_path(ref,5000)) if field_path(ref,5000).exists() else sha(Path(str(field_path(ref,5000))+'.gz')),rho_normalized_l2=rho_err,velocity_normalized_l2=vel_err,all24_flux_error_over_Qtarget=flux_err,scalar_errors=errors,formal_Qin=float(sa['Qin_4dx']),TestB_Qin=float(sb['Qin_4dx']),field_comparison='ALL_182694_PHYSICAL_FLUID_NODES',tolerances='Existing GPU numerical comparison contract; unchanged')

def inspect_snapshot(root,run,step):
    root=Path(root);run=Path(run);c=contract(root)
    data=history(run);scalar=data[data['iteration']==step]
    assert len(scalar)==1;scalar=scalar[0]
    f=fields(field_path(run,step));sf=fields(field_path(run,step,True))
    with np.load(run/'contracts/control_volume_nodes.npz') as mask:
        ids=mask['flat_index'];cv=ids[mask['control_volume_mask']!=0]
    assert np.array_equal(f['index'],ids),'GPU physical domain/ownership mismatch'
    assert len(sf)==len(np.unique(sf['index'])),'GPU sampled field alias'
    scale=c['dx_effective_m']/c['dt_s'];unit=c['rho_kg_m3']*c['dx_effective_m']**3
    mass=float(f['rho'].sum()*unit);cv_mass=float(f['rho'][np.searchsorted(ids,cv)].sum()*unit)
    assert np.isclose(mass,scalar['total_mass'],rtol=2e-11,atol=1e-25)
    assert np.isclose(cv_mass,scalar['control_volume_mass'],rtol=2e-11,atol=1e-25)
    speed=np.linalg.norm(f['u'],axis=1)*scale
    values={'rho_min':float(f['rho'].min()),'rho_max':float(f['rho'].max()),'rho_mean':float(f['rho'].mean()),'u_max_m_s':float(speed.max()),'mean_speed_m_s':float(speed.mean())}
    for k,v in values.items():assert np.isclose(v,scalar[k],rtol=2e-11,atol=1e-18),k
    quad=np.loadtxt(run/'contracts/multiplane_quadrature.tsv')
    shape=np.array(read(run/'contracts/geometry_reuse_contract.json')['lattice_shape'])
    q,m=reconstruct_flux(quad,sf,shape,c['dx_effective_m'],c['dt_s'],c['rho_kg_m3'])
    eq=np.array([scalar[f'{p}_{k}dx'+('_dx2' if v else '')] for p in PORTS for k in OFFSETS for v in [0,1]]);eq[:6]*=-1
    em=np.array([scalar[f'mass_outward_g{i}'] for i in range(24)])
    qe=float(np.max(np.abs(q-eq))/c['Qtarget_m3_s']);me=float(np.max(np.abs(m-em))/(c['rho_kg_m3']*c['Qtarget_m3_s']))
    assert qe<=1e-9 and me<=1e-9,'GPU independent quadrature mismatch'
    return dict(status='PASS',step=step,full_field_nodes=len(f),sample_nodes=len(sf),field_sha256=sha(field_path(run,step)),sample_sha256=sha(field_path(run,step,True)),total_mass=mass,control_volume_mass=cv_mass,volume_flux_error_over_Qtarget=qe,mass_flux_error_over_rho_Qtarget=me,statistics=values)
