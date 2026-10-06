"""Read saved states; no additional equilibrium or force campaign."""
import csv,json,time
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from run_audit import OUT,ROOT,CFG,KINDS,H,L,DELTA,MAT,load,save,csvsave,primary,quality,field_for,norm,require
from pd_clot.pfpd_candidate.nonlocal_boundary.residual_decomposition import balances
from pd_clot.pfpd_candidate.nonlocal_boundary.diagnostics import load_checks

def c1_profile(kind):
    old=ROOT/'verification/baqp_near_boundary_kinematic_energy_audit_v1'
    prefix='AXIAL' if kind=='nonaffine_axial' else 'SHEAR'
    with np.load(old/f'{prefix}_NODE_DIAGNOSTICS.npz') as a:
        node={k:a[k] for k in a.files};psi_scale=float(a['Psi_exact'].max())
    with np.load(old/'BOND_GEOMETRY.npz') as a:geom={k:a[k] for k in a.files}
    with np.load(old/f'{prefix}_BOND_DIAGNOSTICS.npz') as a:bond={k:a[k] for k in a.files}
    # Read the accepted depth field by its documented array key.
    depthkey=next(k for k in geom if k in ('depth_over_delta','bond_depth_over_delta','midpoint_depth_over_delta'))
    rows=[]
    for sample,d,depth in [('Fbar',node,node['depth_over_delta']),('Ftilde',bond,geom[depthkey])]:
        bins=np.rint(depth*DELTA/H*2).astype(int)
        for k in np.unique(bins):
            mask=bins==k
            rows.append(dict(case=kind,backend='C1_TARGET',sample=sample,depth_h=k/2,depth_over_delta=k/2*H/DELTA,
                F_error_max_scaled=float(d['actual_F_error'][mask].max()/.004),Psi_plus_error_max_scaled=float(abs(d['Psi_actual_bias'][mask]).max()/psi_scale)))
    return rows

def run():
    started=time.monotonic();g,f,b=primary();allprofiles=[];recovery=[];sums=load('PD2_TARGET_SUMMARIES.json')
    with np.load(OUT/'PD2_NODEWISE_BOUNDARY_DECOMPOSITION.npz') as a:
        affine_roundoff={}
        for k in KINDS[:2]:
            external=a[k+'__f_external_NL'];conversion=a[k+'__f_BC_localized']
            affine_roundoff[k]=dict(external_minus_conversion_max_N=float(norm(external-conversion).max()),
                external_minus_conversion_scaled=float(norm(external-conversion).max()/sums[k]['scales']['F_ref_N']),
                actual_external_load_checks=load_checks(g,f.X,external,field_for(k)),
                interpretation='Analytic affine q=0. Measured oracle roundoff was retained in conversion; the common NL load path includes -qV of roundoff magnitude. Actual saved external load is reported without retrospectively changing it.')
        save('PD2_AFFINE_LOAD_ROUNDOFF_ACCOUNTING.json',affine_roundoff)
        for k in KINDS:
            summary,rows,_,psi=quality(g,f,b,a[k+'__target'],k)
            csvsave(k.upper()+'_TARGET_FIELD_PROFILE.csv',rows)
            allprofiles += [dict(case=r['case'],backend='PD2_TARGET',sample=r['sample'],depth_h=r['depth_h'],depth_over_delta=r['depth_over_delta'],F_error_max_scaled=r['F_error_max_scaled'],Psi_plus_error_max_scaled=r['Psi_plus_error_max_scaled']) for r in rows]
        # General nonaffine force / reference-moment / work identities.
        acc={};work=[]
        for k in KINDS[2:]:
            field=field_for(k);sc=sums[k]['scales'];arrays={key:a[k+'__'+key] for key in ['f_FD','f_BC_localized','f_body_local','f_model']}
            continuum=field.continuum();kw={}
            if k=='nonaffine_shear':kw=dict(variations=[(n,u) for n,u,_ in field.variations(f.X)],moment_source=continuum['reference_moment_source_Nm'])
            force,moment,w=balances(f.X,arrays,continuum,sc['F_ref_N'],L,**kw)
            acc[k]=dict(force=force,moment=moment,continuum=continuum)
            work += [dict(case=k,**r) for r in w]
            require(max(force['closure_scaled'],moment['closure_scaled'],max(r['closure_scaled'] for r in w))<=1e-10,'Nonaffine accounting closure failed')
        save('PD2_NONAFFINE_ACCOUNTING.json',acc);csvsave('PD2_NONAFFINE_VIRTUAL_WORK.csv',work)
    for k in CFG['solve_order']:
        a=load('PD2_'+k.upper()+'.json');require(a['PASS'],'Unaccepted equilibrium result')
        recovery.append(dict(case=k,force_residual=a['force_residual'],Fbar_error_scaled=a['fields']['Fbar']['F_error_max_scaled'],Ftilde_error_scaled=a['fields']['Ftilde']['F_error_max_scaled'],Psi_plus_error_scaled=max(v['Psi_plus_error_max_scaled'] for v in a['fields'].values()),displacement_shift_m=a['displacement_shift_max_m'],energy_target_J=a['energy_target_J'],energy_solved_J=a['energy_solved_J'],external_work_target_J=a['external_work_target_J'],external_work_solved_J=a['external_work_solved_J'],PASS=a['PASS']))
        with (OUT/(k.upper()+'_FIELD_PROFILE.csv')).open() as file:
            allprofiles += [dict(case=k,backend='PD2_SOLVED',sample=r['sample'],depth_h=float(r['depth_h']),depth_over_delta=float(r['depth_over_delta']),F_error_max_scaled=float(r['F_error_max_scaled']),Psi_plus_error_max_scaled=float(r['Psi_plus_error_max_scaled'])) for r in csv.DictReader(file)]
    with ProcessPoolExecutor(8,mp_context=mp.get_context('spawn')) as pool:
        for rows in pool.map(c1_profile,KINDS[2:]):allprofiles+=rows
    csvsave('C1_PD2_FIELD_QUALITY.csv',allprofiles);csvsave('PD2_EQUILIBRIUM_RECOVERY.csv',recovery)
    with np.load(OUT/'NONAFFINE_SHEAR_RECOVERED_STATE.npz') as n,np.load(OUT/'CLASSICAL_SHEAR_RECOVERED_STATE.npz') as c,np.load(OUT/'PD2_NODEWISE_BOUNDARY_DECOMPOSITION.npz') as a:
        model=a['nonaffine_shear__f_model'];free=~n['fixed'];force_diff=c['force']-n['force']
        comparison=dict(NL_to_classical_displacement_max_m=float(norm(c['x']-n['x']).max()),NL_to_classical_Fbar_max=float(norm(c['Fbar']-n['Fbar']).max()),
            external_difference_equals_model_scaled=float(norm(c['external']-n['external']-model).max()/sums['nonaffine_shear']['scales']['F_ref_N']),
            force_response_closure_relative_to_model=float(np.linalg.norm((force_diff+model)[free])/np.linalg.norm(model[free])),
            NL_recovery_over_classical_shift=float(norm(n['x']-n['target']).max()/norm(c['x']-c['target']).max()),
            classical_target_residual_on_solver_scale=float(np.linalg.norm(model[free])/load('PD2_CLASSICAL_SHEAR.json')['solver_force_scale_N']))
        save('PD2_CLASSICAL_VS_NONLOCAL_SHEAR.json',comparison)
        z=f.X[:,2];levels=np.unique(z);profiles=[]
        psi_n=np.load(OUT/'NONAFFINE_SHEAR_FIELD_DIAGNOSTICS.npz');psi_c=np.load(OUT/'CLASSICAL_SHEAR_FIELD_DIAGNOSTICS.npz')
        for j,Z in enumerate(levels):
            m=z==Z
            profiles.append(dict(level=j,Z_um=Z*1e6,NL_displacement_error_max_m=float(norm(n['x'][m]-n['target'][m]).max()),classical_displacement_error_max_m=float(norm(c['x'][m]-c['target'][m]).max()),NL_Fbar_error_max=float(psi_n['Fbar__F_error'][m].max()),classical_Fbar_error_max=float(psi_c['Fbar__F_error'][m].max()),NL_Psi_error_max_Pa=float(abs(psi_n['Fbar__Psi_plus_error'][m]).max()),classical_Psi_error_max_Pa=float(abs(psi_c['Fbar__Psi_plus_error'][m]).max())))
        csvsave('PD2_CLASSICAL_SHEAR_SHIFT_PROFILE.csv',profiles)
    save('POSTPROCESS_COMPLETE.json',dict(PASS=True,wall_seconds=time.monotonic()-started))

if __name__=='__main__':run()
