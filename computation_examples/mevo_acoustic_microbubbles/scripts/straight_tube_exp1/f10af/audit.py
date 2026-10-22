"""Independent accepted-step balance, geometry, IDs and compact truthful data.
Far-state audit uses A=I only after checking every pair; near states use original
historical resistance assembly. This changes no integration step or threshold.
"""
import time,json,concurrent.futures,multiprocessing
import numpy as np,h5py
from mevo_mb.fast_sim.f10g_contact import resistance
from scripts.straight_tube_exp1.f10aa.analyze import independent
from .common import *
from .metrics import spatial,rotation

def one(r):
    start=time.perf_counter();path=BASE/r['trajectory_path'];p=policy();zeta=6*np.pi*.001*R;maxratio=0.;maxcomp=0.;minreaction=0.;maxoverlap=0.;mingap=np.inf;maxwall=0.;continuity=0.;checked=0;near=0
    dest=BASE/'dataset'/(r['id']+'.h5');dest.parent.mkdir(exist_ok=True)
    with h5py.File(path) as h,h5py.File(dest,'w') as d:
        tx=h['time_s'][:];x=h['positions_m'];n=x.shape[1];ii,jj=np.triu_indices(n,1);f=h['force_samples'];assert np.array_equal(h['particle_id'][:],np.arange(n));assert len(tx)==r['accepted_steps']+1;assert abs(tx[-1]-r['completed_s'])<1e-13
        for lo in range(0,len(tx)-1,128):
            hi=min(lo+128,len(tx)-1);q=x[lo:hi+1];v=f['velocity_m_s'][lo:hi];force=f['total_force_N'][lo:hi];lam=f['pair_contact_reaction_N'][lo:hi];before=q[:-1,jj]-q[:-1,ii];after=q[1:,jj]-q[1:,ii];nb=np.linalg.norm(before,axis=2);na=np.linalg.norm(after,axis=2);gap=na-2*R;scale=np.maximum(np.maximum(np.linalg.norm(force,axis=2).max(1),np.abs(lam).max(1)),1e-18);residual=zeta*v-force
            for k in range(hi-lo):
                if (nb[k]-2*R<.1*R).any():residual[k]=(zeta*resistance(q[k],R)['A']@v[k].ravel()).reshape(n,3)-force[k];near+=1
                if np.any(lam[k]):
                    pair=lam[k,:,None]*after[k]/na[k,:,None];np.add.at(residual[k],ii,pair);np.add.at(residual[k],jj,-pair)
            ratio=np.linalg.norm(residual,axis=2).max(1)/(p['force_residual_atol_N']+p['force_residual_rtol']*scale);maxratio=max(maxratio,float(ratio.max()));maxcomp=max(maxcomp,float((abs(lam*gap)/scale[:,None]).max()));minreaction=min(minreaction,float(lam.min()));mingap=min(mingap,float(gap.min()));maxoverlap=max(maxoverlap,float(-gap.min()));maxwall=max(maxwall,float((np.linalg.norm(q[:,:,1:],axis=2)+R-.0008).max()),float((abs(q[:,:,0])+R-.01625).max()));continuity=max(continuity,float(np.linalg.norm(np.diff(q,axis=0)-v*np.diff(tx[lo:hi+1])[:,None,None],axis=2).max()));assert np.all(f['gate_pass'][lo:hi]);checked+=hi-lo
        ix=np.unique(np.rint(np.linspace(0,len(tx)-1,min(401,len(tx)))).astype(int));xx=x[ix];d['time_s']=tx[ix];d['positions_m']=xx;d['particle_id']=h['particle_id'][:];d['raw_frame_index']=ix
        d['step_metrics']=h['step_metrics'][:];d['step_metrics'].attrs['columns']=h['step_metrics'].attrs['columns'];d['step_failed_gates']=h['step_failed_gates'][:]
        fi=np.minimum(ix,len(tx)-2);metrics=[]
        for j,k in enumerate(ix):
            vel=f['velocity_m_s'][fi[j]];s=spatial(xx[j]);s.pop('spectrum_S');s.pop('spectrum_freq_m1');s.pop('spectrum_theta_rad');s.update(time_s=tx[k],rotation=rotation(xx[j],vel),velocity_time_s=float(f['time_s'][fi[j]]));metrics.append(s)
        for key in ['time_s','pressure_Pa','primary_force_N','secondary_force_N','total_force_N','velocity_m_s','contact_force_N']:
            data=np.asarray([f[key][j] for j in fi]);d.create_dataset('force_samples/'+key,data=data,compression='gzip')
        for k,v in h.attrs.items():d.attrs[k]=v
        d.attrs.update(raw_relative_path=r['trajectory_path'],raw_sha256=sha(path),sampling='Exact accepted nodes, no interpolation; endpoint velocity from last accepted-step start',task_json=json.dumps(r['task']))
    write(BASE/'dataset'/(r['id']+'_metrics.json'),metrics)
    check=dict(PASS=bool(maxratio<=1 and maxcomp<=p['actual_geometry_complementarity_m'] and minreaction>=-p['dual_negative_N'] and maxoverlap<=p['hard_geometry_m'] and maxwall<=1e-12),force_gate_fraction_max=maxratio,actual_complementarity_m=maxcomp,min_reaction_N=minreaction,overlap_m=maxoverlap,min_gap_m=mingap,wall_penetration_m=maxwall,step_identity_residual_m=continuity,checked_steps=checked,near_contact_steps=near,trajectory_status=r['status'],scope='Accepted prefix; failed proposals never admitted' if r['status']!='PASS' else 'Complete trajectory')
    reference=None
    if r['id'] in ['pair_contact','pair_r0_X_PAIR']:reference=independent(path.parent);assert check['PASS']==reference['PASS'] and abs(check['force_gate_fraction_max']-reference['force_gate_fraction_max'])<1e-9
    return dict(id=r['id'],audit=check,original_audit_crosscheck=reference,wall_s=time.perf_counter()-start,raw_sha256=sha(path),raw_bytes=path.stat().st_size,dataset_sha256=sha(dest))
def main():
    begin=time.perf_counter();verify();rows=[read(p) for p in sorted((BASE/'hf').glob('*/SUMMARY.json'))];known={z['id']:z for z in read(BASE/'AUDIT.json')['rows']} if (BASE/'AUDIT.json').exists() else {};todo=[r for r in rows if r.get('trajectory_path') and r['id'] not in known]
    with concurrent.futures.ProcessPoolExecutor(max_workers=3,mp_context=multiprocessing.get_context('spawn')) as ex:
        for z in ex.map(one,todo):known[z['id']]=z;print('AUDIT',z['id'],z['audit']['PASS'],flush=True)
    write(BASE/'AUDIT.json',dict(rows=list(known.values()),checked_steps=sum(z['audit']['checked_steps'] for z in known.values()),all_accepted_prefixes_pass=all(z['audit']['PASS'] for z in known.values()),wall_s=time.perf_counter()-begin,identity=identity()))
    write(BASE/'ALL_HF_RUNS.json',rows);write(BASE/'ALL_LF_RUNS.json',[read(p) for p in sorted((BASE/'lf').glob('*/SUMMARY.json'))])
if __name__=='__main__':main()
