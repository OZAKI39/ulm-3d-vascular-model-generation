"""Full-physics F10AA Campaign.run inherited verbatim; extended caches only."""
import time,traceback,concurrent.futures,multiprocessing
import numpy as np,h5py
from scripts.straight_tube_exp1.f10aa.trajectory import Campaign as Original
from scripts.straight_tube_exp1.f10ab.run import AuditEngine
from scripts.straight_tube_exp1.f10ae.schedules import schedule
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from .common import *
from .field import response_context
from .metrics import endpoint
class Campaign(Original):
    def __init__(self,out,t):
        verify();self.started=time.perf_counter();self.out=out;out.mkdir(parents=True,exist_ok=False);cfg=read(CONFIG)
        self.cfg=dict(max_trajectories=1,per_run_budget_s=cfg['budget']['per_HF_wall_s'],remote_budget_s=cfg['budget']['per_HF_wall_s']);self.backend='B1';self.state=t['state'];self.p=policy();f=t.get('frequency',2500000.);self.c=response_context(f);self.ident=identity();self.ident['execution_commit']=self.ident['source_commit'];self.records=[]
        folder=BASE/'fields'/str(int(f));assert read(folder/'FIELD_GATE.json')['PASS'],'Extended field not admitted'
        self.bases={k:CachedWaterBasis(folder/(k+'.h5')) for k in ['coarse','fine']};self.c['basis']=self.bases['coarse']
        if t.get('historical_cache'):self.bases={k:CachedWaterBasis(OLD/'cache'/(k+'.h5')) for k in ['coarse','fine']}
        with h5py.File(BASE/'INITIAL_STATES.h5') as h:x=h[self.state][:]
        if t.get('offset_m'):x=x+np.asarray(t['offset_m'])
        self.states={len(x):x};self.schedules={'CONTROL':schedule(t)}
    def engine(self,n,cache='coarse'):
        e=AuditEngine(self.states[n],self.c['physics'],self.c['response'],self.bases[cache],self.c['cavity'],policy=self.p);e.gate_records=[];self.last_engine=e;return e
def task(id,action='X_PAIR',state='N64_wide_r0',duration=.002,**kw):
    a=read(BASE/'ACTIONS.json')[action]
    t=dict(id=id,action=action,state=state,amplitudes=a['amplitudes'],phases=a['phases'],frequency=a['frequency'],duration=duration);t.update(kw);return t
def perform(t):
    start=time.perf_counter();out=BASE/'hf'/t['id'];r=dict(id=t['id'],task=t,status='FAIL',final=None)
    if (out/'SUMMARY.json').exists():
        prior=read(out/'SUMMARY.json');assert prior['task']==clean(t)
        if 'trajectory_path' in prior:assert sha(BASE/prior['trajectory_path'])==prior['trajectory_sha256']
        return prior
    try:
        c=Campaign(out,t);n=next(iter(c.states));raw=c.run(n,'CONTROL',t['duration'],t.get('variant','DT1'));r.update(raw);folder=out/'runs'/raw['name'];path=folder/'trajectory.h5'
        write(folder/'F10AF_INPUT.json',dict(task=t,identity=identity(),protocol_sha256=sha(BASE/'PREREGISTRATION.json'),initial_sha256=sha(BASE/'INITIAL_STATES.h5'),response=c.c['response'].metadata,field=c.bases[raw['cache']].metadata))
        write(folder/'SOLVER_GATES.json',c.last_engine.gate_records)
        with h5py.File(path) as h:x0=h['positions_m'][0];x=h['positions_m'][-1];vel=h['force_samples/velocity_m_s'][-1] if 'velocity_m_s' in h['force_samples'] else None
        r.update(trajectory_path=str(path.relative_to(BASE)),initial=endpoint(x0,x0),final=endpoint(x0,x,vel) if raw['status']=='PASS' else None,scope='FULL_F10AA_PHYSICS')
    except Exception as e:r.update(reason=repr(e),traceback=traceback.format_exc())
    r['task_wall_s']=time.perf_counter()-start;write(out/'SUMMARY.json',r);print('HF',t['id'],r['status'],r.get('completed_s'),flush=True);return clean(r)
def batch(tasks,workers=3):
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as pool:return list(pool.map(perform,tasks))
