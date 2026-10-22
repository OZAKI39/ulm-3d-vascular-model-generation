"""Historical regressions, segmented-clock controls, pilot and low-force dynamics."""
import time,h5py,numpy as np
from scripts.straight_tube_exp1.f10z_scaling.analyze import compare
from .common import *
from .run import task,batch
def compare_rows(a,b,root_a=BASE,root_b=BASE):
    if a['status']!='PASS' or b['status']!='PASS':return dict(admitted=False,reason='Incomplete trajectory')
    z=compare((root_a/a['trajectory_path']).parent,(root_b/b['trajectory_path']).parent)
    z.update(reference=a['id'],comparison=b['id']);return z
def main():
    verify();start=time.perf_counter();plans=[]
    elites=read(AC/'ELITE_VALIDATION_PREREGISTRATION.json')['selected']
    for k,r in enumerate(elites):t=r['task'].copy();t['id']=f'baseline_elite{k}';t['historical_response']=True;plans.append(t)
    plans += [task('baseline_OFF',amps=[0]*4),task('baseline_Xsingle',amps=[0,0,700,0]),task('baseline_X180')]
    for v in ['DT05','FINE']:plans.append(task('baseline_X180_'+v,variant=v))
    a=dict(amplitudes=[0,0,700,700],phases=[0,0,0,np.pi],duration=.002)
    off=dict(amplitudes=[0]*4,phases=[0]*4,duration=.002)
    plans += [task('schedule_A',duration=.004),task('schedule_B',duration=.004,segments=[a,a]),task('schedule_C',duration=.004,segments=[a,off])]
    low=read(BASE/'LOW_PRESSURE_ANALYSIS.json')
    if all(low['gates'][k] for k in ['LOW_PRESSURE_APECSS_CONVERGED','LOW_PRESSURE_FORCE_VALIDATED','LOW_PRESSURE_INTERPOLATION_VALIDATED']):
        plans += [task('low_dynamic_'+v,amps=[0,0,.1,.1],phase=[0]*4,variant=v,joined_diagnostic=True) for v in ['DT1','DT05','FINE']]
        t=elites[2]['task'].copy();t.update(id='joined_historical_elite2',joined_diagnostic=True);plans.append(t)
        plans.append(task('joined_OFF',amps=[0]*4,joined_diagnostic=True))
        for old in read(AC/'ALL_RUNS.json'):
            if old['status']!='PASS':
                t=old['task'].copy();t.update(id='joined_retry_'+old['id'],joined_diagnostic=True);plans.append(t)
    write(BASE/'BASELINE_PLAN.json',dict(plans=plans,identity=identity(),workers=4))
    rows=batch(plans);write(BASE/'BASELINE_RUNS.json',rows);by={r['id']:r for r in rows};reg=[]
    for k,old in enumerate(elites):
        new=by[f'baseline_elite{k}'];z=compare_rows(old,new,AC,BASE)
        with h5py.File(AC/old['trajectory_path']) as h,h5py.File(BASE/new['trajectory_path']) as j:z['bitwise_positions']=bool(np.array_equal(h['positions_m'][:],j['positions_m'][:]))
        z['frequency']=old['task'].get('frequency',2.5e6);reg.append(z)
    sched=compare_rows(by['schedule_A'],by['schedule_B'])
    with h5py.File(BASE/by['schedule_C']['trajectory_path']) as h:
        x=h['positions_m'][:];t=h['time_s'][:];i=np.argmin(abs(t-.002));sched['OFF_stage_displacement_m']=float(abs(x[i:]-x[i]).max());sched['clock_strictly_increasing']=bool(np.all(np.diff(t)>0));sched['switch_count']=int(np.count_nonzero(abs(t-.002)<1e-15));sched['actual_end_s']=float(t[-1])
    checks=[compare_rows(by['baseline_X180'],by['baseline_X180_'+v]) for v in ['DT05','FINE']]
    lows=[compare_rows(by['low_dynamic_DT1'],by['low_dynamic_'+v]) for v in ['DT05','FINE']] if 'low_dynamic_DT1' in by else []
    low['gates']['HISTORICAL_2P5MHZ_REGRESSION_PASS']=all(z['admitted'] and z['bitwise_positions'] for z in reg if z['frequency']==2.5e6)
    joined_reg=compare_rows(by['baseline_elite2'],by['joined_historical_elite2']) if 'joined_historical_elite2' in by else dict(admitted=False)
    low['gates']['JOINED_HISTORICAL_TRAJECTORY_PASS']=joined_reg['admitted']
    low['gates']['LOW_PRESSURE_DYNAMIC_TEST_PASS']=len(lows)==2 and all(z['admitted'] for z in lows)
    if 'joined_historical_elite2' in by:
        from .response import ExtendedResponse
        from scripts.straight_tube_exp1.f10s.common import context
        c=context();o=c['response'];j=ExtendedResponse(o,c['physics']);p=np.geomspace(.1,40000,501).astype(complex)*np.exp(.37j)
        low['gates']['HISTORICAL_QUERY_BITWISE_PASS']=bool(np.array_equal(j.query(p),o.query(p)))
        low['gates']['JOINED_OFF_EXACT_ZERO']=by['joined_OFF']['status']=='PASS' and np.linalg.norm(by['joined_OFF']['displacement_m'])==0
        low['joined_metadata']=j.metadata
    if 'joined_historical_elite2' in by:
        from scripts.straight_tube_exp1.f10ab.run import AuditEngine
        from mevo_mb.fast_sim.nearfield import pair_reference
        from .run import schedule
        with h5py.File(BASE/'INITIAL_STATES.h5') as h:x=h['S2_r0'][:]
        forces=[]
        for amplitude in [.01,.03,.1]:
            t=task('force_diagnostic',amps=[0,0,amplitude,amplitude],phase=[0]*4)
            e=AuditEngine(x,c['physics'],j,c['basis'],c['cavity'],policy=policy());e.control(schedule(t).segments[0].event());gpu=e.fields(x)
            pp,gg=c['basis'].query_basis(x);u=c['basis'].vector(t['amplitudes'],t['phases']);pressure=u@pp;grad=np.einsum('i,inj->nj',u,gg);v=j.query(pressure)
            f1=-.5*np.real(v[:,0,None]*np.conj(grad));chi=j.query(np.array([complex(j.knots[0])]))[0,0]/j.knots[0];f2=pair_reference(x,v,pressure,c['physics'],chi)
            errors=[float(np.linalg.norm(ref-gpu[key])/max(np.linalg.norm(ref),1e-80)) for ref,key in [(f1,'primary_force_N'),(f2,'secondary_force_N')]]
            forces.append(dict(source_amplitude=amplitude,pressure_min_Pa=float(abs(pressure).min()),pressure_max_Pa=float(abs(pressure).max()),F1_rms_N=float(np.sqrt(np.mean(f1*f1))),F2_rms_N=float(np.sqrt(np.mean(f2*f2))),cpu_gpu_relative_errors=errors,PASS=max(errors)<1e-8))
            np.savez_compressed(BASE/'low_pressure'/f'force_{amplitude}.npz',positions_m=x,pressure_Pa=pressure,gradient_Pa_m=grad,F1_N=f1,F2_N=f2)
        low['actual_force_checks']=forces;low['gates']['LOW_PRESSURE_FORCE_CPU_GPU_PASS']=all(z['PASS'] for z in forces)
    low['production_extension_admitted']=all(v is True for v in low['gates'].values())
    low['decision']='Admit only explicit F10AE opt-in 0.001 Pa extension; all historical intervals, files and queries preserved; no frequency extrapolation' if low['production_extension_admitted'] else 'Retain historical fail-closed domain: incomplete low-pressure admission'
    write(BASE/'LOW_PRESSURE_FINAL.json',low)
    write(BASE/'BASELINE_VALIDATION.json',dict(regressions=reg,schedule=sched,numerical_checks=checks,low_dynamics=lows,joined_regression=joined_reg,wall_s=time.perf_counter()-start,identity=identity()))
    print('BASELINE_COMPLETE',reg,sched,flush=True)
if __name__=='__main__':main()
