"""Official Ax/BoTorch ask/tell versus SciPy Sobol, same bounded physics budget."""
import time, numpy as np
from scipy.stats import qmc
from .common import *
from .run import task,batch
def action(identifier,p):
    return task(identifier,amps=[0,0,p['amplitude'],p['amplitude']],phase=[0,0,0,np.deg2rad(p['phase_deg'])])
def parameters(t):return dict(amplitude=float(t['amplitudes'][2]),phase_deg=float(np.rad2deg(t['phases'][3])))
def tell(client,index,r):
    if r['status']!='PASS':
        client.mark_trial_failed(trial_index=index)
        return 'FAILED_NO_ENDPOINT_OBSERVATION'
    z=r['targets'];client.complete_trial(trial_index=index,raw_data=dict(loss=(z['T1_loss'],0.),COM_nm=(z['COM_m']*1e9,0.)))
    return 'COMPLETE'
def main():
    verify();start=time.perf_counter();cfg=read(CONFIG)['search']
    assert (BASE/'LOW_PRESSURE_FINAL.json').exists(),'Low pressure decision must precede search'
    import torch
    torch.set_num_threads(1)
    from ax.api.client import Client
    from ax.api.configs import RangeParameterConfig
    n=cfg['comparison_evaluations'];shared=cfg['shared_sobol_initial']
    points=qmc.scale(qmc.Sobol(2,scramble=True,seed=cfg['seed']).random_base2(int(np.log2(n))),[200,150],[1100,240])
    sobol=[action(f'sobol_{i:02}',dict(amplitude=float(a),phase_deg=float(p))) for i,(a,p) in enumerate(points)]
    write(BASE/'SEARCH_PLAN.json',dict(identity=identity(),sobol=sobol,ax_shared_ids=[t['id'] for t in sobol[:shared]],per_method_observations=n,physical_calls=2*n-shared,selection='Feasible minimum T1 loss; optimizer advantage from first target hit and total wall, single seed only',constraints=['COM_nm <= 20'],low_pressure_decision_sha256=sha(BASE/'LOW_PRESSURE_FINAL.json')))
    common_start=time.perf_counter();common=batch(sobol[:shared]);common_s=time.perf_counter()-common_start
    write(BASE/'SEARCH_COMMON.json',dict(rows=common,wall_s=common_s))
    sub=time.perf_counter();extra=batch(sobol[shared:]);sobol_s=common_s+time.perf_counter()-sub
    sr=common+extra;write(BASE/'SEARCH_SOBOL.json',dict(rows=sr,wall_s=sobol_s,optimizer_s=0.,shared_reused=0,new_physical_calls=n))
    client=Client(random_seed=cfg['seed'])
    client.configure_experiment(name='F10AE_T1',parameters=[RangeParameterConfig(name='amplitude',parameter_type='float',bounds=(200.,1100.)),RangeParameterConfig(name='phase_deg',parameter_type='float',bounds=(150.,240.))])
    client.configure_optimization(objective='-loss',outcome_constraints=['COM_nm <= 20'])
    client.configure_generation_strategy(initialization_budget=shared,initialize_with_center=False,use_existing_trials_for_initialization=True)
    records=[];ar=list(common);optimizer_s=0.;sub=time.perf_counter()
    for r in common:
        idx=client.attach_trial(parameters=parameters(r['task']));status=tell(client,idx,r);records.append(dict(trial=idx,id=r['id'],reused=True,status=status))
    while len(ar)<n:
        if time.perf_counter()-sub>cfg['optimizer_wall_budget_s']:break
        tick=time.perf_counter();suggestions=client.get_next_trials(max_trials=min(4,n-len(ar)));ask_s=time.perf_counter()-tick;optimizer_s+=ask_s
        if not suggestions:break
        plans=[action(f'ax_{idx:02}',p) for idx,p in suggestions.items()]
        rs=batch(plans)
        tick=time.perf_counter()
        for (idx,p),r in zip(suggestions.items(),rs):records.append(dict(trial=idx,parameters=p,id=r['id'],reused=False,status=tell(client,idx,r),ask_batch_s=ask_s))
        optimizer_s+=time.perf_counter()-tick;ar+=rs
        client.save_to_json_file(filepath=str(BASE/'AX_EXPERIMENT.json'))
        write(BASE/'AX_ASK_TELL.json',records)
        write(BASE/'SEARCH_AX_PROGRESS.json',dict(rows=ar,optimizer_s=optimizer_s,wall_s=common_s+time.perf_counter()-sub))
    ax_s=common_s+time.perf_counter()-sub
    write(BASE/'SEARCH_AX.json',dict(rows=ar,wall_s=ax_s,optimizer_s=optimizer_s,shared_reused=shared,new_physical_calls=len(ar)-shared))
    def summary(name,rows,wall,opt):
        valid=[r for r in rows if r['status']=='PASS'];feasible=[r for r in valid if r['targets']['COM_TARGET_PASS']];hits=[i+1 for i,r in enumerate(rows) if r['status']=='PASS' and r['targets']['T1_DISCOVERY_PASS']]
        return dict(method=name,observations=len(rows),complete=len(valid),feasible=len(feasible),raw_target_hits=len(hits),first_target_observation=min(hits) if hits else None,best=min(feasible,key=lambda r:r['targets']['T1_loss'])['id'] if feasible else None,best_loss=min([r['targets']['T1_loss'] for r in feasible],default=None),wall_s=wall,optimizer_s=opt)
    summaries=[summary('Sobol',sr,sobol_s,0),summary('Ax_BoTorch',ar,ax_s,optimizer_s)]
    write(BASE/'OPTIMIZER_COMPARISON.json',dict(methods=summaries,fair_observation_budget=len(sr)==len(ar)==n,single_optimizer_seed=True,wall_s=time.perf_counter()-start,identity=identity()))
    print('SEARCH_COMPLETE',summaries,flush=True)
if __name__=='__main__':main()
