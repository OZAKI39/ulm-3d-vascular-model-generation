"""Frozen actions: independent seeds, final once-only test and matched physics."""
import time,numpy as np,h5py
from .common import *
from .run import batch
from .baseline import compare_rows
def main():
    verify();start=time.perf_counter();selected=read(BASE/'SELECTED_FROZEN.json')['selected'];plans=[]
    for target,s in selected.items():
        for state in ['S2_r0']+[f'S2_validation_{i}' for i in range(3)]:
            for v in ['DT1','DT05','FINE']:
                t=s['task'].copy();t.update(id=f'validate_{target}_{state}_{v}',state=state,variant=v);plans.append(t)
    write(BASE/'VALIDATION_PLAN.json',dict(plans=plans,frozen_sha256=sha(BASE/'SELECTED_FROZEN.json'),identity=identity()))
    rows=batch(plans);write(BASE/'VALIDATION_RUNS.json',rows);lookup={r['id']:r for r in rows};checks=[];groups=[]
    for target in selected:
        rep=[]
        for state in ['S2_r0']+[f'S2_validation_{i}' for i in range(3)]:
            a=lookup[f'validate_{target}_{state}_DT1'];cc=[];axis_noise=np.zeros(3);shape_noise=0.;kabsch_noise=0.
            for v in ['DT05','FINE']:
                b=lookup[f'validate_{target}_{state}_{v}'];z=compare_rows(a,b);checks.append(z);cc.append(z)
                if a['status']=='PASS' and b['status']=='PASS':
                    axis_noise=np.maximum(axis_noise,abs(np.array(a['targets']['axis_delta_m'])-b['targets']['axis_delta_m']));shape_noise=max(shape_noise,abs(a['targets']['aspect_gain']-b['targets']['aspect_gain']));kabsch_noise=max(kabsch_noise,abs(a['targets']['kabsch_m']-b['targets']['kabsch_m']))
            numeric=all(z['admitted'] for z in cc)
            target_pass=False
            if a['status']=='PASS':
                if target in ['T1','T2']:target_pass=a['targets'][target+'_DISCOVERY_PASS'] and abs(a['targets']['aspect_gain'])>max(1e-4,5*shape_noise)
                elif target=='T3':target_pass=bool(np.any(np.array(a['targets']['axis_delta_m'])>np.maximum(1e-9,5*axis_noise)))
                elif target=='T4':target_pass=all(a['recovery']['raw_recovery_pass'].values())
            rep.append(dict(state=state,id=a['id'],numerical_pass=numeric,target_pass=target_pass,confirmed=numeric and target_pass,axis_noise_m=axis_noise,aspect_noise=shape_noise,kabsch_noise_m=kabsch_noise,targets=a.get('targets'),recovery=a.get('recovery'),connected_every_step=a.get('connected_every_step'),connected_at_end=a.get('connected_at_end')))
        groups.append(dict(target=target,replicates=rep,discovery_and_validation_confirmed=all(r['confirmed'] for r in rep)))
    write(BASE/'VALIDATION.json',dict(groups=groups,comparisons=checks,wall_s=time.perf_counter()-start,identity=identity()))
    # Final tests run once with frozen controls even if validation fails: no reselection.
    tests=[]
    for target,s in selected.items():
        for i in range(3):
            t=s['task'].copy();t.update(id=f'final_{target}_{i}',state=f'S2_test_{i}',variant='DT1');tests.append(t)
    write(BASE/'FINAL_TEST_PLAN.json',dict(plans=tests,frozen_sha256=sha(BASE/'SELECTED_FROZEN.json'),only_once=True,identity=identity()))
    final=batch(tests);write(BASE/'FINAL_TEST_RUNS.json',final)
    # Same action, matched discrete-frequency field and radial response; no interpolation.
    stress=[]
    for target in ['T1','T2']:
        for f in [2400000.,2500000.,2600000.]:
            t=selected[target]['task'].copy();t.update(id=f'freq_{target}_{int(f)}',frequency=f);stress.append(t)
        for shell in [.8,1.2]:
            t=selected[target]['task'].copy();t.update(id=f'shell_{target}_{shell:.1f}',shell=shell);stress.append(t)
        for offset in [-5e-6,5e-6]:
            t=selected[target]['task'].copy();t.update(id=f'position_{target}_{offset:.2g}',offset_m=[0,offset,0]);stress.append(t)
    # Explicit Y near-node diagnostic even if the selected T1 uses only X sources.
    from .run import task
    for offset in [-5e-6,0.,5e-6]:stress.append(task(f'Ynode_{offset:.2g}',amps=[700,700,0,0],phase=[0,np.pi/2,0,0],offset_m=[0,offset,0]))
    write(BASE/'ROBUSTNESS_PLAN.json',dict(plans=stress,frequency_microperturbation='NOT_RUN',reason='Only previously admitted 2.4/2.5/2.6 MHz tables',identity=identity()))
    sr=batch(stress);write(BASE/'ROBUSTNESS_RUNS.json',sr);ref=[]
    for r in sr:
        if r['status']=='PASS':
            for v in ['DT05','FINE']:t=r['task'].copy();t.update(id=t['id']+'_'+v,variant=v);ref.append(t)
    refined=batch(ref);write(BASE/'ROBUSTNESS_REFINED.json',refined);comparisons=[]
    for r in sr:
        for b in refined:
            if b['id'] in [r['id']+'_DT05',r['id']+'_FINE']:comparisons.append(compare_rows(r,b))
    write(BASE/'ROBUSTNESS_CONVERGENCE.json',comparisons)
    write(BASE/'VALIDATION_COMPLETE.json',dict(wall_s=time.perf_counter()-start,identity=identity()))
if __name__=='__main__':main()
