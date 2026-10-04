"""C5 activation diagnostics only; no historical or production writers used."""
from pathlib import Path
from unittest.mock import patch
from datetime import datetime, timezone
import json, time, hashlib, resource
import numpy as np
from scipy.optimize import brentq
from scripts import spectral_split_audit as spectral
from scripts import spatial_audit_v2 as v2
from pd_clot.pfpd_candidate import branch_stability_tracking as tracker
from pd_clot.pfpd_candidate import coupling_stability_audit as ca
from pd_clot.pfpd_candidate import adr_adapter_at1 as adr
from pd_clot.regularization.equilibrium import check_positions

ROOT=spectral.ROOT
OLD=spectral.OUT
OUT=ROOT/'verification/pfpd_at1_c5_activation_stationarity_audit_v1'
RESULTS=ROOT/'results/pfpd_at1_c5_activation_stationarity_audit_v1'
Problem=spectral.SpectralProblem
read=ca.read_json
save=ca.save_json
SCAN=[.00430,.00435,.00438,.00440,.00442,.00444,.00446,.00448,.00450,.00452,.00455]
GUARD=dict(safety=dict(maximum_displacement_m=.01,maximum_radius_ratio=1.),pipe=dict(radius_m=1.))


def declare():
    cfg=dict(task='PFPD_AT1_C5_ACTIVATION_STATIONARITY_AUDIT_V1',created_utc=datetime.now(timezone.utc).isoformat(),
        force_tolerance=1e-9,phase_tolerance=1e-8,q_scaled_active_tolerance=ca.audit_config()['active_q_scaled_tolerance'],
        s_tolerance=ca.audit_config()['active_s_tolerance'],q_scale='D_b=max(Psi_plus_b,Yc); q_active_tol_b=1e-7*D_b',
        scan=SCAN,activation_bracket_width=1e-7,maximum_bisections=12,
        activation_point='After bracket: SciPy brentq on equilibrated min(q), xtol=1e-13; obtain a weak-lower diagnostic state.',
        active_selection_A3='A1 union A2 plus bonds sharing at least one endpoint and q_scaled<=1e-3',
        active_selection_A3_near_threshold=1e-3,maximum_active_set_updates=10,
        root_scan_offsets=[1e-6,2e-6,5e-6,1e-5,2e-5,5e-5],continuation_step=1e-6,maximum_additional_continuation_states=5,
        reduced_solver=dict(method='trf',tr_solver='exact',ftol=1e-13,xtol=1e-13,gtol=1e-13,max_nfev=120,x_scale=1.,loss='linear'),
        diagnostic_Jx_step=3e-6,diagnostic_Jx_validation_steps=[1e-5,3e-6,1e-6],
        singular_value_relative_residual_tolerance=1e-6,svds_maxiter=500,
        classification_ill_conditioned_threshold=1e8,
        reproduction_tolerances=dict(energy_relative=1e-10,force_absolute=1e-12,KKT_absolute=1e-12,q_absolute=1e-12,Y_relative=1e-10),
        physical_model_changes=False,production_changes=False,stability_assessed=False,
        original_solver_config=str(OLD/'provenance/PREDECLARED_CONFIG.json'),
        references=dict(least_squares='https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html',
                        svds='https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.svds.html'))
    save(OUT/'provenance/PREDECLARED_CONFIG.json',cfg)


def q_fields(p,x,s):
    Y=p.state(x,s).W0.copy()
    q=-2*(1-s)*Y+2*p.yc
    D=np.maximum(Y,p.yc)
    return Y,q,D,q/D


def classify(s,lower,qs,sc):
    cfg=read(OUT/'provenance/PREDECLARED_CONFIG.json')
    near=s-lower<=cfg['s_tolerance'];tol=cfg['q_scaled_active_tolerance']
    return dict(STRONG_LOWER=near&(qs>tol),WEAK_LOWER=near&(abs(qs)<=tol),
                VIOLATED_LOWER=near&(qs<-tol),INTERIOR=(~near)&(s<sc-cfg['s_tolerance']))


def record_fields(p,x,s):
    Y,q,D,qs=q_fields(p,x,s);groups=classify(s,p.lower,qs,p.sc)
    arg=s-qs;distance=np.minimum(abs(arg-p.lower),abs(arg-p.sc))
    clip=np.where(arg<p.lower,-1,np.where(arg>p.sc,1,0))
    row=dict(strain=p.strain,metrics=p.metrics(x,s),min_q=float(q.min()),
        q_percentile_levels=[0,1,5,25,50,75,95,99,100],q_percentiles=np.percentile(q,[0,1,5,25,50,75,95,99,100]).tolist(),
        min_q_scaled=float(qs.min()),negative_q_count=int((q<0).sum()),weak_count=int((abs(qs)<=1e-7).sum()),
        max_Y_over_Yc=float((Y/p.yc).max()),counts={k:int(v.sum()) for k,v in groups.items()},
        IDs={k:np.flatnonzero(v).tolist() for k,v in groups.items()},
        negative_q_IDs=np.flatnonzero(q<0).tolist(),
        kink_counts=dict(lower_clipped=int((clip==-1).sum()),interior_map=int((clip==0).sum()),upper_clipped=int((clip==1).sum()),
                         exactly_kink=int((distance==0).sum()),near_kink=int((distance<=1e-7).sum())))
    arrays=dict(x=x,s=s,Y=Y,q=q,q_scaled=qs,q_scale=D,clip_branch=clip,kink_distance=distance,
                bond_midpoints=p.family.X[p.family.pairs].mean(axis=1),bond_orientation=p.family.xi/p.family.length[:,None],pairs=p.family.pairs)
    return row,arrays


def reproduce():
    rows=[]
    for label in ['C5_W_004','C5_A_001','C5_W_005','C5_P_005']:
        before=read(OLD/f'branch_tracking/{label}/STATIONARY.json');a=ca.read_npz(OLD/f'branch_tracking/{label}/STATE.npz')
        p=Problem(5,before['strain'],a['lower']);row,arr=record_fields(p,a['x'],a['s'])
        oldqs=np.array(read(OLD/f'active_sets/{label}.json')['q_scaled']);s=a['s']
        # Reconstruct archived Y from q_scaled=-2(1-s)Y/max(Y,Yc)+2Yc/max(Y,Yc).
        above=oldqs<2*s
        oldY=np.where(above,2*p.yc/(oldqs+2*(1-s)),p.yc*(2-oldqs)/(2*(1-s)))
        oldq=oldqs*np.maximum(oldY,p.yc)
        errors=dict(energy_relative=abs(row['metrics']['total_J']/before['metrics']['total_J']-1),
            force_absolute=abs(row['metrics']['force_residual']-before['metrics']['force_residual']),
            KKT_absolute=abs(row['metrics']['projected_phase_KKT']-before['metrics']['projected_phase_KKT']),
            q_absolute=float(abs(arr['q']-oldq).max()),Y_relative=float(abs(arr['Y']-oldY).max()/arr['Y'].max()))
        passed=all(v<=read(OUT/'provenance/PREDECLARED_CONFIG.json')['reproduction_tolerances'][k] for k,v in errors.items())
        row.update(label=label,errors=errors,passed=passed,original_accepted=before['accepted'],original_solver=before['solver'])
        save(OUT/f'reproduction/{label}.json',row);np.savez_compressed(OUT/f'reproduction/{label}.npz',**arr,lower=a['lower'])
        rows.append(row)
    result=dict(passed=all(r['passed'] for r in rows),states=rows,
        Y_reference='Reconstructed algebraically from archived q_scaled and s; not a new independent constitutive oracle.')
    save(OUT/'reproduction/ALL.json',result)
    assert result['passed'],'ACTIVATION_AUDIT_REPRODUCTION_FAILED'
    print('REPRODUCTION PASSED',flush=True)


def intact(strain,label):
    folder=OUT/f'intact/{label}'
    if (folder/'SUMMARY.json').exists():
        return read(folder/'SUMMARY.json'),ca.read_npz(folder/'STATE.npz')
    p=Problem(5,strain);start=time.time();s=np.zeros(p.ns)
    answer=adr.equilibrate_tensile(p.cloud,p.family,p.material,p.top,v2.affine(p),s,strain,p.config)
    x=answer.x;row,arrays=record_fields(p,x,s);state=p.state(x,s)
    check_positions(adr.boundary_view(p.cloud,p.top,strain),x,GUARD)
    gates=dict(ADR=bool(answer.converged),force=row['metrics']['force_residual']<=1e-9,
               fixed=bool(np.array_equal(x[p.cloud.fixed],p.boundary[p.cloud.fixed])),positive_J=bool(np.all(np.linalg.det(state.Ftilde)>0)),
               phase_frozen=bool(np.count_nonzero(s)==0),safety=True)
    row.update(label=label,intact_accepted=all(gates.values()),gates=gates,wall_seconds=time.time()-start,
               peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    save(folder/'SUMMARY.json',row);save(folder/'ADR_SUMMARY.json',answer.summary)
    np.savez_compressed(folder/'STATE.npz',**arrays)
    assert row['intact_accepted'],'INTACT_EQUILIBRIUM_FAILED'
    print(f'INTACT {label} eps={strain:.12g} minq={row["min_q"]:.9g} qneg={row["negative_q_count"]} weak={row["weak_count"]}',flush=True)
    return row,arrays


def scan():
    assert read(OUT/'reproduction/ALL.json')['passed']
    rows=[]
    for i,e in enumerate(SCAN):
        r,_=intact(e,f'scan_{i:02d}');rows.append(r);save(OUT/'activation/SCAN.json',rows)
    crossing=next(((a,b) for a,b in zip(rows[:-1],rows[1:]) if a['min_q']>=0 and b['min_q']<0),None)
    if crossing is None:raise RuntimeError('Activation not bracketed; inspect before finite extension')
    low,high=crossing;history=[]
    for i in range(12):
        if high['strain']-low['strain']<=1e-7:break
        r,_=intact((low['strain']+high['strain'])/2,f'bisect_{i:02d}')
        old=[low['strain'],high['strain']]
        if r['min_q']>=0:low=r
        else:high=r
        history.append(dict(before=old,after=[low['strain'],high['strain']],min_q=r['min_q']))
    result=dict(interval=[low['strain'],high['strain']],width=high['strain']-low['strain'],lower=low,upper=high,
                history=history,kind='LOWER_BOUND_KKT_ACTIVATION_BRACKET_NOT_STABILITY')
    save(OUT/'activation/BRACKET.json',result)
    count=[0]
    def fun(e):
        if e==low['strain']:return low['min_q']
        if e==high['strain']:return high['min_q']
        count[0]+=1
        return intact(e,f'weak_locator_{count[0]:02d}')[0]['min_q']
    e=brentq(fun,low['strain'],high['strain'],xtol=1e-13)
    weak,arrays=intact(e,'activation_point')
    save(OUT/'activation/POINT.json',weak)
    spatial=[]
    for name,r in [('last_strong',low),('activation',weak),('first_q_negative',high)]:
        a=arrays if r is weak else ca.read_npz(OUT/f'intact/{r["label"]}/STATE.npz')
        ids=np.flatnonzero(a['q_scaled']<=1e-7)
        spatial.append(dict(role=name,label=r['label'],strain=r['strain'],candidate_count=len(ids),fraction=len(ids)/len(a['q']),
             candidates=[dict(ID=int(b),midpoint_m=a['bond_midpoints'][b].tolist(),orientation=a['bond_orientation'][b].tolist(),
                              Psi_plus=float(a['Y'][b]),q=float(a['q'][b])) for b in ids]))
    save(OUT/'activation/SPATIAL.json',spatial)
    print('ACTIVATION BRACKET',result['interval'],'POINT',e,'WEAK',weak['weak_count'],flush=True)


if __name__=='__main__':
    declare();reproduce();scan()
