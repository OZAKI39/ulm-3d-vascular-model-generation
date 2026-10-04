"""Reproduce A/B, then run the two predeclared stationary schedules and brackets."""
from concurrent.futures import ProcessPoolExecutor
import argparse
import hashlib
import numpy as np
from pd_clot.pfpd_candidate.coupling_stability_audit import (OUT as AUDIT_OUT,OLD,AuditProblem,
    physics_config,load_audit_state,audit_config,read_json,read_npz,save_json)
from pd_clot.pfpd_candidate.branch_stability_tracking import (
    OUT,settings,displacement_predictor,solve_stationary_branch_state,audit_stationary_state)
from scripts.run_coupling_stability_audit import worker_budget


def save_solution(label,p,result,parent=None,initial_kind=None):
    folder=OUT/'branch_states'/label;folder.mkdir(parents=True,exist_ok=False)
    np.savez_compressed(folder/'STATE.npz',**{k:result[k] for k in ['x','s','z','z_initial','lower']})
    summary={k:v for k,v in result.items() if k not in ['x','s','z','z_initial','lower','history']}
    summary.update(label=label,strain=p.strain,previous_accepted_state=parent,initial_kind=initial_kind,
        phase_predictor_is_previous_phase=bool(np.array_equal(result['z_initial'][p.nx:],p.lower)) if initial_kind=='continuation_predictor' else None,
        lower_sha256=hashlib.sha256(p.lower.tobytes()).hexdigest())
    save_json(folder/'STATIONARITY.json',summary);save_json(folder/'EVALUATIONS.json',result['history'])
    return summary


def reproduction(task):
    key,start=task;p,d=load_audit_state(key);c=settings();label=f'reproduce_{key}_{start}'
    if start=='saved':x,s=d['x'],d['s']
    else:
        state=audit_config()['states'][key];prev_strain=physics_config()['schedules']['loading_unloading'][state['load_step']-1]
        x=displacement_predictor(p,d['previous']['x'],prev_strain);s=d['previous']['s_pf_bond']
    result=solve_stationary_branch_state(p,x,s);save_solution(label,p,result,d['path'],start)
    if not result['accepted']:return dict(state=key,start=start,passed=False,gates=result['gates'],reason=result['reason'])
    curvature=audit_stationary_state(label,p,result);old=p.metrics(d['x'],d['s'])
    old_lambda=read_json(AUDIT_OUT/f'hessian/{key}/eps_3e-06/SUMMARY.json')['minimum_reduced_eigenvalue']
    errors=dict(x_absolute_m=float(abs(result['x']-d['x']).max()),s_absolute=float(abs(result['s']-d['s']).max()),
        energy_relative=abs(result['metrics']['total_J']-old['total_J'])/abs(old['total_J']),
        force_absolute=abs(result['metrics']['force_residual']-old['force_residual']),
        KKT_absolute=abs(result['metrics']['projected_phase_KKT']-old['projected_phase_KKT']),
        lambda_absolute=abs(curvature['reduced_eigenvalues_min'][1]-old_lambda))
    passed=all(v<=c['reproduction_tolerances'][k] for k,v in errors.items())
    return dict(state=key,start=start,passed=passed,errors=errors,metrics=result['metrics'],classification=curvature['classification'])


def continuity(previous,current):
    old=read_npz(OUT/'branch_states'/previous['state']/'STATE.npz');new=read_npz(OUT/'branch_states'/current['state']/'STATE.npz')
    p=AuditProblem(physics_config(),current['strain'],new['lower']);dx=current['strain']-previous['strain']
    pred=displacement_predictor(p,old['x'],previous['strain']);change={}
    for group,ids in previous['active']['bond_IDs'].items():
        before=set(ids);after=set(current['active']['bond_IDs'][group]);change[group]=dict(entered=sorted(after-before),left=sorted(before-after))
    switched=set()
    for r in change.values():switched.update(r['entered']);switched.update(r['left'])
    event=dict(previous=previous['state'],current=current['state'],strain_increment=dx,
        x_increment_L2_m=float(np.linalg.norm(new['x']-old['x'])),phase_increment_L2=float(np.linalg.norm(new['s']-old['s'])),
        maximum_phase_increment=float(abs(new['s']-old['s']).max()),minimum_phase_increment=float((new['s']-old['s']).min()),
        x_predictor_error_scaled_L2=float(np.linalg.norm(new['x']-pred)/p.length_scale),
        total_energy_difference_J=current['metrics']['total_J']-previous['metrics']['total_J'],
        reaction_difference_N=current['reaction_N']-previous['reaction_N'],
        active_set_changes=change,switched_bond_IDs=sorted(switched),switched_bonds=len(switched))
    rules=settings()['continuity'];reasons=[]
    if len(switched)/p.ns>=rules['mass_switch_fraction']:reasons.append('MASS_ACTIVE_SET_SWITCH')
    if event['x_predictor_error_scaled_L2']>rules['maximum_x_predictor_error_scaled_L2']:reasons.append('LARGE_PREDICTOR_DISTANCE')
    if event['maximum_phase_increment']>500*dx+.01:reasons.append('LARGE_PHASE_JUMP')
    event['status']='BRANCH_CONTINUITY_EVENT' if reasons else 'CONTINUOUS_SAMPLED_STEP';event['reasons']=reasons
    return event


def solve_next(label,strain,parent):
    previous=read_npz(OUT/'branch_states'/parent['state']/'STATE.npz')
    p=AuditProblem(physics_config(),strain,previous['s'])
    x=displacement_predictor(p,previous['x'],parent['strain'])
    result=solve_stationary_branch_state(p,x,previous['s'])
    stat=save_solution(label,p,result,parent['state'],'continuation_predictor')
    if not result['accepted']:return None,dict(label=label,strain=strain,stationarity=stat)
    row=audit_stationary_state(label,p,result);row['continuity']=continuity(parent,row)
    save_json(OUT/f'active_sets/{label}_EVENT.json',row['continuity'])
    return row,None


def tracking_schedule(name):
    c=settings();folder=OUT/('coarse_tracking' if name=='A' else 'fine_tracking')
    values=c['schedules'][name];rows=[];interruption=None
    p,d=load_audit_state('A');r=solve_stationary_branch_state(p,d['x'],d['s'])
    label=f'{name}_grid_00';save_solution(label,p,r,d['path'],'saved_seed')
    if not r['accepted']:interruption=dict(label=label,strain=p.strain,reason=r['reason'])
    else:rows.append(audit_stationary_state(label,p,r))
    for i,strain in enumerate(values[1:],start=1):
        if interruption:break
        row,interruption=solve_next(f'{name}_grid_{i:02d}',strain,rows[-1])
        if row:rows.append(row);print('GRID',name,strain,row['classification'],row['stability_indicator'],flush=True)
    summary=dict(schedule=name,requested_strains=values,completed_states=len(rows),
        status='BRANCH_TRACK_INTERRUPTED' if interruption else 'TRACKED',interruption=interruption,states=rows)
    save_json(folder/'TRACKING.json',summary)
    bracket=None
    for lo,hi in zip(rows,rows[1:]):
        if lo['classification']=='STABLE_SAMPLED_STATE' and hi['classification']=='UNSTABLE_SAMPLED_STATE':bracket=(lo,hi);break
    refined=[]
    if bracket:
        lo,hi=bracket;initial=[lo['strain'],hi['strain']];refine_status='BRACKETED'
        for level in range(c['refinement_max_levels']):
            if hi['strain']-lo['strain']<=c['refinement_width']:break
            midpoint=(lo['strain']+hi['strain'])/2
            row,failed=solve_next(f'{name}_refine_{level:02d}',midpoint,lo)
            entry=dict(level=level,before=[lo['strain'],hi['strain']],midpoint=midpoint,failed=failed,state=row)
            if failed:refine_status='BRANCH_TRACK_INTERRUPTED';refined.append(entry);break
            if row['classification']=='STABLE_SAMPLED_STATE':lo=row
            elif row['classification']=='UNSTABLE_SAMPLED_STATE':hi=row
            else:refine_status='REFINEMENT_UNRESOLVED';refined.append(entry);break
            entry['after']=[lo['strain'],hi['strain']];refined.append(entry)
            print('REFINE',name,level,entry['after'],row['classification'],flush=True)
        changed=lo['active']['bond_IDs']!=hi['active']['bond_IDs']
        output=dict(schedule=name,status=refine_status,initial_interval=initial,final_interval=[lo['strain'],hi['strain']],
            width=hi['strain']-lo['strain'],midpoint_estimate=(lo['strain']+hi['strain'])/2,
            method='NONSMOOTH_STABILITY_BRACKET' if changed else 'INTERVAL_BISECTION',
            active_set_changed=changed,lower=lo,upper=hi,refinement=refined)
    else:output=dict(schedule=name,status='NO_STABILITY_SIGN_CHANGE_IN_COARSE_RANGE',initial_interval=None,final_interval=None,refinement=[])
    save_json(OUT/f'critical_bracket/SCHEDULE_{name}.json',output)
    return dict(schedule=name,tracking_status=summary['status'],completed_states=len(rows),bracket_status=output['status'],interval=output['final_interval'])


def main():
    args=argparse.ArgumentParser();args.add_argument('--stage',choices=['reproduce','track','all'],default='all');opt=args.parse_args()
    budget=worker_budget();save_json(OUT/f'provenance/WORKER_BUDGET_{opt.stage.upper()}.json',budget)
    if opt.stage in ['reproduce','all']:
        with ProcessPoolExecutor(max_workers=budget['max_workers']) as pool:
            rows=list(pool.map(reproduction,[(k,start) for k in ['A','B'] for start in ['saved','predictor']]))
        gate=dict(status='VERIFIED' if all(r['passed'] for r in rows) else 'BRANCH_TRACKER_REPRODUCTION_FAILED',cases=rows)
        save_json(OUT/'comparisons/REPRODUCTION.json',gate);print('REPRODUCTION',gate,flush=True)
    gate=read_json(OUT/'comparisons/REPRODUCTION.json')
    if gate['status']!='VERIFIED':raise SystemExit('STOP: BRANCH_TRACKER_REPRODUCTION_FAILED')
    if opt.stage in ['track','all']:
        with ProcessPoolExecutor(max_workers=min(2,budget['max_workers'])) as pool:rows=list(pool.map(tracking_schedule,['A','B']))
        save_json(OUT/'comparisons/SCHEDULE_RUNS.json',rows);print(rows,flush=True)


if __name__=='__main__':main()
