"""Independent complete-window metrics, all-step KKT audit and timing summaries."""
import time
import numpy as np,h5py
from scipy.spatial.distance import pdist,squareform
from mevo_mb.fast_sim.f10g_contact import resistance
from scripts.straight_tube_exp1.f10z_scaling.analyze import summarize,compare,load
from scripts.straight_tube_exp1.f10z_scaling.metrics import edge_events
from .common import *

def independent(folder):
    p=policy();r,t,x=load(folder);n=len(x[0]);ii,jj=np.triu_indices(n,1)
    max_force=0.;max_ratio=0.;max_comp=0.;min_reaction=0.;max_reaction=0.;max_gaps=0.
    with h5py.File(folder/'trajectory.h5') as h:
        f=h['force_samples'];lam=f['pair_contact_reaction_N'][:];forces=f['total_force_N'][:];v=f['velocity_m_s'][:]
        assert np.all(abs(f['time_s'][:]-t[:-1])<1e-13) and len(lam)==len(x)-1
        assert np.all(f['gate_pass'][:]);zeta=6*np.pi*.001*R
        for k in range(len(lam)):
            # Historical resistance assembly; independently scatter exact final normals.
            A=resistance(x[k],R)['A'];d=x[k+1,jj]-x[k+1,ii];norm=np.linalg.norm(d,axis=1);normal=d/norm[:,None];gap=norm-2*R
            force=np.zeros((n,3));np.add.at(force,ii,-lam[k,:,None]*normal);np.add.at(force,jj,lam[k,:,None]*normal)
            residual=np.linalg.norm((zeta*A@v[k].ravel()).reshape(n,3)-forces[k]-force,axis=1).max()
            scale=max(np.linalg.norm(forces[k],axis=1).max(),abs(lam[k]).max(initial=0),1e-18)
            max_force=max(max_force,residual);max_ratio=max(max_ratio,residual/(p['force_residual_atol_N']+p['force_residual_rtol']*scale))
            max_comp=max(max_comp,abs(lam[k]*gap).max(initial=0)/scale)
            min_reaction=min(min_reaction,lam[k].min(initial=0));max_reaction=max(max_reaction,lam[k].max(initial=0));max_gaps=max(max_gaps,-gap.min(initial=0))
        check=dict(force_residual_max_N=max_force,force_gate_fraction_max=max_ratio,actual_complementarity_m=max_comp,
            min_reaction_N=min_reaction,max_reaction_N=max_reaction,overlap_m=max_gaps,checked_steps=len(lam),all_saved_gates=True,trajectory_sha256=r['trajectory_sha256'])
        check['trajectory_status']=r['status']
        check['PASS']=bool(max_ratio<=1 and max_comp<=p['actual_geometry_complementarity_m'] and min_reaction>=-p['dual_negative_N'] and max_gaps<=p['hard_geometry_m'])
    return check

def extra_compare(a,b):
    r,t,x=load(a);s,u,y=load(b);j=np.searchsorted(u,t-1e-13);y=y[j]
    nearest=0.;ec=[];active=[];forces=[];reactions=[];pairs=[];field=[]
    for xx,yy in zip(x,y):
        da=squareform(pdist(xx));db=squareform(pdist(yy));np.fill_diagonal(da,np.inf);np.fill_diagonal(db,np.inf)
        nearest=max(nearest,abs(da.min(1)-db.min(1)).max())
        iu=np.triu_indices(len(xx),1);ec.append(np.count_nonzero((da[iu]<=2*R+1e-9)!=(db[iu]<=2*R+1e-9)))
    with h5py.File(a/'trajectory.h5') as ha,h5py.File(b/'trajectory.h5') as hb:
        fa=ha['force_samples'];fb=hb['force_samples'];tb=fb['time_s'][:];ia=np.searchsorted(tb,fa['time_s'][:]-1e-13)
        for key in ['contact_force_N','pair_contact_reaction_N','total_force_N','pressure_Pa']:
            va=fa[key][:];vb=fb[key][:][ia];forces.append(float(abs(va-vb).max(initial=0)))
        la=fa['pair_contact_reaction_N'][:];lb=fb['pair_contact_reaction_N'][:][ia]
        active=int(np.count_nonzero((la>0)!=(lb>0)))
    return dict(nearest_neighbor_max_difference_m=float(nearest),actual_contact_geometry_mismatch_pair_frames=int(sum(ec)),
        active_reaction_mismatch_pair_steps=active,net_contact_force_max_difference_N=forces[0],pair_reaction_max_difference_N=forces[1],
        total_force_max_difference_N=forces[2],pressure_max_difference_Pa=forces[3],scope='All accepted fixed IDs/frames; force comparisons at matching start times; dt refinements compare different intervals, not identical reactions')

def main():
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--core-only',action='store_true');args=ap.parse_args()
    verify();out=BASE/'analysis';out.mkdir(exist_ok=True);start=time.perf_counter();runs={};audits=read(out/'AUDITS.json') if (out/'AUDITS.json').exists() else {};metrics={}
    files=sorted(BASE.glob('full/*/runs/*/RESULT.json'))+(sorted(BASE.glob('extension/*/runs/*/RESULT.json'))+sorted(BASE.glob('threaded/*/runs/*/RESULT.json'))+sorted(BASE.glob('parallel/w*/runs/*/RESULT.json')) if not args.core_only else [])
    for pth in files:
        key=str(pth.parent.relative_to(BASE));r,m,tx=summarize(pth.parent);runs[key]=r
        # Full audit for every trajectory, including original CPU and negative candidates.
        if key not in audits or audits[key].get('trajectory_sha256')!=r['trajectory_sha256']:audits[key]=independent(pth.parent)
        np.savez_compressed(out/(key.replace('/','__')+'.npz'),**{k:v for k,v in m.items() if k not in ['centered_positions_m','edges']})
        write(out/'AUDITS.json',audits);write(out/'RUNS.json',runs)
        print(key,audits[key]['PASS'],flush=True)
    extra={};performance={}
    for n in [32,64,128]:
        duration=read(CONFIG)['full_windows_s'][str(n)];name=f'N{n}_PAIR_{round(duration*1e6)}us_DT1'
        a=BASE/'full/B0/runs'/name
        performance[str(n)]={}
        ref=read(a/'RESULT.json')
        for backend in ['B0','B1','PACKED','B3_CuPy','B2_t2','B2_t4']:
            b=BASE/f'{"threaded" if backend.startswith("B2_") else "full"}/{backend}/runs'/name
            if not b.exists():continue
            r=read(b/'RESULT.json');performance[str(n)][backend]=dict(wall_s=r['wall_s'],integration_s=r['integration_s'],initialization_s=r['initialization_s'],save_s=r['wall_s']-r['integration_s'],step=r['step_latency'],cpu_s=r['cpu_s'],speedup=ref['wall_s']/r['wall_s'],physical_s_per_wall_s=r['completed_s']/r['wall_s'],status=r['status'])
            if backend!='B0':extra[f'{n}_{backend}']=extra_compare(a,b)
        for variant in ['DT05','FINE']:
            extra[f'{n}_{variant}']=extra_compare(BASE/'full/B1/runs'/name,BASE/'full/B1/runs'/name.replace('DT1',variant))
    write(out/'EXTRA_COMPARISONS.json',extra);write(out/'PERFORMANCE.json',performance)
    write(out/'COMPLETE.json',dict(status='CORE_COMPLETE' if args.core_only else 'COMPLETE',all_independent_steps_pass=all(v['PASS'] for v in audits.values()),checked_steps=sum(v['checked_steps'] for v in audits.values()),wall_s=time.perf_counter()-start,identity=identity()))
if __name__=='__main__':main()
