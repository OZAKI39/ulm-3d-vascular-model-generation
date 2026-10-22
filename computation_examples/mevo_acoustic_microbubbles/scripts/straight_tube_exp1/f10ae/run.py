"""F10AA integration reused unchanged; no new mechanics or acoustic kernels."""
import time,traceback,concurrent.futures,multiprocessing
import numpy as np,h5py
from scipy.spatial.distance import cdist
from scipy.sparse.csgraph import connected_components
from scripts.straight_tube_exp1.f10aa.trajectory import Campaign as Original
from scripts.straight_tube_exp1.f10ab.run import AuditEngine
from scripts.straight_tube_exp1.f10ab.metrics import measure,series
from scripts.straight_tube_exp1.f10ac_qd.group_diagnostics import partition
from scripts.straight_tube_exp1.f10ad_multifrequency.interface import context,FrequencySegment
from mevo_mb.control.scheduled_pzt import Schedule
from mevo_mb.acoustics.water_disc_proxy import CachedWaterBasis
from .common import *
from .metrics import target_metrics,recovery
from .schedules import schedule
class Campaign(Original):
    def __init__(self,out,t):
        verify();self.started=time.perf_counter();self.out=out;out.mkdir(parents=True,exist_ok=False)
        budget=read(CONFIG)['search']['per_trajectory_wall_s']
        self.cfg=dict(max_trajectories=1,per_run_budget_s=budget,remote_budget_s=budget)
        self.backend='B1';self.state=t.get('state','S2_r0');self.p=policy();f=t.get('frequency',2.5e6);self.c=context(f,t.get('shell',1.))
        admission=BASE/('LOW_PRESSURE_ADMISSION.json' if (BASE/'LOW_PRESSURE_ADMISSION.json').exists() else 'LOW_PRESSURE_FINAL.json')
        admitted=admission.exists() and read(admission)['production_extension_admitted']
        if f==2.5e6 and t.get('shell',1.)==1. and (t.get('joined_diagnostic') or admitted and not t.get('historical_response')):
            from .response import ExtendedResponse
            self.c['response']=ExtendedResponse(self.c['response'],self.c['physics'])
        if t.get('diagnostic_low'):
            from scripts.straight_tube_exp1.f10ad_multifrequency.interface import FixedResponse
            self.c['response']=FixedResponse(BASE/'low_pressure/diagnostic_response.h5',self.c['physics'])
        self.ident=identity();self.ident['execution_commit']=self.ident['source_commit'];self.records=[]
        # F10AD shell variants at 2.5 MHz use AD coarse/fine as a matched pair.
        old=f==2.5e6 and t.get('shell',1.)==1.
        fine=OLD/'cache/fine.h5' if old else AD/str(round(f))/'fine.h5'
        self.bases={'coarse':self.c['basis'],'fine':CachedWaterBasis(fine)}
        self.input_library=BASE/'INITIAL_STATES.h5'
        with h5py.File(self.input_library) as h:x=h[self.state][:]
        if t.get('offset_m'):x=x+np.asarray(t['offset_m'])
        self.states={len(x):x};self.schedules={'CONTROL':schedule(t)}
    def engine(self,n,cache='coarse'):
        e=AuditEngine(self.states[n],self.c['physics'],self.c['response'],self.bases[cache],self.c['cavity'],policy=self.p)
        e.gate_records=[];self.last_engine=e;return e
def perform(t):
    start=time.perf_counter();out=BASE/'runs'/t['id'];r=dict(id=t['id'],task=t,status='FAIL',scientific_status='INITIALIZATION_FAILURE',final=None,targets=None)
    if (out/'SUMMARY.json').exists():
        prior=read(out/'SUMMARY.json')
        if prior['task']!=clean(t):raise RuntimeError('Run ID collision with different task')
        if prior.get('trajectory_path') and sha(BASE/prior['trajectory_path'])!=prior['trajectory_sha256']:raise RuntimeError('Immutable trajectory digest mismatch')
        return prior
    try:
        c=Campaign(out,t);raw=c.run(64,'CONTROL',t['duration'],t.get('variant','DT1'));r.update(raw)
        folder=out/'runs'/raw['name'];path=folder/'trajectory.h5'
        write(folder/'F10AE_INPUT.json',dict(task=t,identity=identity(),protocol_sha256=sha(BASE/'PREREGISTRATION.json'),initial_sha256=sha(c.input_library),field=c.bases[raw['cache']].metadata,response=c.c['response'].metadata,switch_model=read(CONFIG)['switch_model']))
        write(folder/'SOLVER_GATES.json',c.last_engine.gate_records)
        post=time.perf_counter()
        with h5py.File(path) as h:x=h['positions_m'][:];tt=h['time_s'][:]
        complete=raw['status']=='PASS';m=series(tt,x,10)
        all_components=np.array([connected_components(cdist(z,z)<=2.1*R,directed=False,return_labels=False) for z in x])
        m['every_step_components']=all_components
        np.savez_compressed(out/'METRICS.npz',**{k:v for k,v in m.items() if k!='aligned_positions_m'})
        initial=measure(x[0],x[0]);endpoint=measure(x[0],x[-1])
        for z in [initial,endpoint]:z.pop('aligned_positions_m',None)
        status='VALID_COMPLETE_TRAJECTORY' if complete else ('OUT_OF_MODEL_DOMAIN' if any('APECSS' in s or 'cache' in s or 'Outside' in s for s in raw['reasons']) else 'NUMERICAL_OR_BUDGET_FAILURE')
        r.update(trajectory_path=str(path.relative_to(BASE)),initial=initial,final=endpoint if complete else None,prefix_endpoint=endpoint if not complete else None,scientific_status=status,targets=target_metrics(x[0],x[-1]) if complete else None,connected_at_end=complete and endpoint['proximity_components']==1,connected_every_step=complete and bool(np.all(all_components==1)),postprocessing_s=time.perf_counter()-post)
        if complete:r['component_diagnostic']=partition(x[0],x[-1],2.1)
        if complete and t.get('segments'):
            ts=t['segments'][0]['duration'];idx=int(np.argmin(abs(tt-ts)))
            assert abs(tt[idx]-ts)<1e-13
            r['recovery']=recovery(x[0],x[idx],x[-1]);r['switch_index']=idx
    except Exception as e:
        r.update(reason=repr(e),traceback=traceback.format_exc(),final=None,targets=None)
    r['task_wall_s']=time.perf_counter()-start;write(out/'SUMMARY.json',r)
    print('F10AE_RESULT',t['id'],r['status'],r.get('targets'),flush=True)
    return clean(r)
def batch(tasks,workers=4):
    if not tasks:return []
    with concurrent.futures.ProcessPoolExecutor(max_workers=workers,mp_context=multiprocessing.get_context('spawn')) as ex:return list(ex.map(perform,tasks))
def task(id,state='S2_r0',amps=(0,0,700,700),phase=(0,0,0,np.pi),duration=.002,**kw):
    return dict(id=id,state=state,amplitudes=list(amps),phases=list(phase),duration=duration,**kw)
