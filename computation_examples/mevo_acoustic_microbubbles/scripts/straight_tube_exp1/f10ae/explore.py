"""Bounded absolute-extension diagnostics and genuine two-stage recovery search."""
import time,numpy as np,h5py
from scipy.stats import qmc
from scipy.linalg import cho_factor,cho_solve
from mevo_mb.fast_sim.accelerated_contact import resistance
from scripts.straight_tube_exp1.f10ab.run import AuditEngine
from scripts.straight_tube_exp1.f10s.common import context
from .common import *
from .run import task,batch,schedule
def static_scan():
    c=context()
    with h5py.File(BASE/'INITIAL_STATES.h5') as h:x=h['S2_r0'][:]
    y=x-x.mean(0);lam,vec=np.linalg.eigh(y.T@y/len(x));lam=lam[::-1];vec=vec[:,::-1]
    a=cho_factor(resistance(x,R)['A']);zeta=6*np.pi*c['physics'].viscosity_Pa_s*R
    engine=AuditEngine(x,c['physics'],c['response'],c['basis'],c['cavity'],policy=policy())
    points=qmc.scale(qmc.Sobol(4,scramble=True,seed=461911).random_base2(7),[0,0,0,0],[1100,360,1100,360]);rows=[]
    write(BASE/'EXTENSION_PLAN.json',dict(points=points,columns=['X_common_Pa_equivalent','X_relative_deg','Y_common_Pa_equivalent','Y_relative_deg'],frequency_Hz=2500000.,static_count=128,full_trajectory_budget=8,selection='Largest initial positive principal RMS rate after original resistance; real 2 ms trajectories decide, not static force alone',identity=identity()))
    dst=BASE/'static';dst.mkdir(exist_ok=False)
    for i,(ax,px,ay,py) in enumerate(points):
        t=task(f'extension_{i:03}',amps=[ay,ay,ax,ax],phase=[0,np.deg2rad(py),0,np.deg2rad(px)])
        try:
            engine.control(schedule(t).segments[0].event());f=engine.fields(x);v=cho_solve(a,f['total_force_N'].ravel()/zeta).reshape(x.shape);vrel=v-v.mean(0)
            qdot=(y.T@vrel+vrel.T@y)/len(x);ldot=np.diag(vec.T@qdot@vec)/(2*np.sqrt(lam))
            p,g=c['basis'].query_basis(x);u=c['basis'].vector(t['amplitudes'],t['phases']);gradient=np.einsum('i,inj->nj',u,g)
            np.savez_compressed(dst/f'{i:03}.npz',positions_m=x,gradient_Pa_m=gradient,free_velocity_m_s=v,**f)
            rows.append(dict(id=t['id'],task=t,valid=True,axis_rates_m_s=ldot,max_rate_m_s=float(max(ldot)),COM_rate_m_s=v.mean(0),pressure_min_Pa=float(abs(f['pressure_Pa']).min()),pressure_max_Pa=float(abs(f['pressure_Pa']).max()),centered_force_N=f['total_force_N']-f['total_force_N'].mean(0)))
        except Exception as e:rows.append(dict(id=t['id'],task=t,valid=False,reason=repr(e)))
    chosen=sorted([r for r in rows if r['valid']],key=lambda r:r['max_rate_m_s'],reverse=True)[:8]
    write(BASE/'STATIC_EXTENSION.json',dict(rows=rows,selected_ids=[r['id'] for r in chosen]))
    return [r['task'] for r in chosen]
def main():
    verify();start=time.perf_counter();plans=static_scan();ex=batch(plans);write(BASE/'EXTENSION_RUNS.json',ex)
    baseline=read(BASE/'BASELINE_RUNS.json');a=next(r for r in baseline if r['id']=='baseline_X180')
    first=dict(amplitudes=a['task']['amplitudes'],phases=a['task']['phases'],duration=.002)
    pts=qmc.scale(qmc.Sobol(4,scramble=True,seed=461912).random_base2(5),[0,0,0,0],[1100,360,700,360]);plans=[]
    for i,(ax,px,ay,py) in enumerate(pts):
        if i<16:ay=0.
        second=dict(amplitudes=[float(ay),float(ay),float(ax),float(ax)],phases=[0,float(np.deg2rad(py)),0,float(np.deg2rad(px))],duration=.002)
        plans.append(task(f'recovery_{i:02}',duration=.004,segments=[first,second]))
    write(BASE/'RECOVERY_PLAN.json',dict(plans=plans,first_stage=a['id'],discovery_loss='Sum of Rg/axes/covariance relative errors plus Kabsch/initial Rg; no numerical-failure endpoint',controls=['baseline_X180','schedule_C'],frequency_Hz=2.5e6,switch_model='QUASISTEADY_RADIAL_RESPONSE_AT_SWITCH',budget=32,identity=identity()))
    rr=batch(plans);write(BASE/'RECOVERY_RUNS.json',rr)
    search=read(BASE/'SEARCH_SOBOL.json')['rows']+read(BASE/'SEARCH_AX.json')['rows']
    uniq={r['id']:r for r in search+baseline+ex if r['status']=='PASS' and r['task']['duration']==.002 and not r['task'].get('diagnostic_low')}
    rows=list(uniq.values());feas=[r for r in rows if r['targets']['COM_TARGET_PASS']]
    t1=min(feas,key=lambda r:r['targets']['T1_loss']);t2=min([r for r in feas if r['connected_at_end']],key=lambda r:r['targets']['T2_loss'])
    t3=max(rows,key=lambda r:max(r['targets']['axis_delta_m']))
    valid=[r for r in rr if r['status']=='PASS']
    def loss(r):
        e=r['recovery']['final_errors'];return e['Rg']+e['axes']+e['covariance']+e['Kabsch_m']/r['initial']['Rg_3d_m']
    t4=min(valid,key=loss) if valid else None
    selected={k:dict(source=r['id'],task=r['task'],discovery_targets=r['targets'],discovery_recovery=r.get('recovery')) for k,r in [('T1',t1),('T2',t2),('T3',t3),('T4',t4)] if r}
    # Freeze before touching validation or test trajectories.
    write(BASE/'SELECTED_FROZEN.json',dict(selected=selected,criterion='Frozen discovery actions; no tuning on validation/test',identity=identity(),wall_s=time.perf_counter()-start))
    print('EXPLORATION_COMPLETE',{k:v['source'] for k,v in selected.items()},flush=True)
if __name__=='__main__':main()
