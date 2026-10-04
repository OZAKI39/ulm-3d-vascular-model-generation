"""Saved-state index, exact-map reproduction and full finite-difference spectra."""
from concurrent.futures import ProcessPoolExecutor, as_completed
import time
import os
import numpy as np
import psutil
from scipy.linalg import eig
from pd_clot.pfpd_candidate.coupling_stability_audit import (
    ROOT,OLD,OUT,AuditProblem,read_json,save_json,read_npz,audit_config,physics_config,
    load_saved_iteration,load_audit_state,green_minimum)


def worker_budget():
    c=audit_config()['resource_limits'];m=psutil.virtual_memory()
    available=m.available-c['reserve_total_RAM_fraction']*m.total-c['parent_bytes']
    if available<c['worker_bytes']: raise RuntimeError('Cannot reserve 20% RAM with one worker')
    return dict(max_workers=min(c['worker_cap'],int(available/c['worker_bytes']),psutil.cpu_count(logical=False)),
                RAM_total=m.total,RAM_available=m.available,reserved_total_fraction=c['reserve_total_RAM_fraction'],
                per_worker_bytes=c['worker_bytes'],parent_bytes=c['parent_bytes'],
                BLAS_threads={k:os.environ.get(k) for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']})


def index_saved_states():
    settings=audit_config();index=[];paths={}
    for run in ['tensile_cube','loading_unloading','loading_unloading_extra_iterations']:
        accepted=read_json(OLD/run/'LOAD_PATH.json');hist=read_json(OLD/run/'STAGGERED_HISTORY.json')
        rows=[]
        for load in accepted:
            step=load['load_step'];folder=OLD/run/f'load_{step:03d}'
            a=read_npz(folder/'accepted_state.npz')
            lower=read_npz(OLD/run/f'load_{step-1:03d}/accepted_state.npz')['s_pf_bond'] if step else np.zeros(len(a['s_pf_bond']))
            p=AuditProblem(read_json(OLD/run/'CONFIG.json'),load['strain'],lower)
            index.append(dict(state_name=f'{run}_accepted_{step:03d}',kind='accepted',source=str(folder/'accepted_state.npz'),
                load_step=step,displacement_m=load['displacement_m'],strain=load['strain'],iteration=None,
                previous_accepted_phase_file=str(OLD/run/f'load_{step-1:03d}/accepted_state.npz') if step else None,
                phase_increment=load['maximum_delta_s'],**p.metrics(a['x'],a['s_pf_bond'])))
        for record in hist:
            step=record['load_step'];iteration=record['iteration'];folder=OLD/run/f'load_{step:03d}'
            phase=read_npz(folder/f'iteration_{iteration:03d}_phase.npz')
            lower=read_npz(OLD/run/f'load_{step-1:03d}/accepted_state.npz')['s_pf_bond'] if step else np.zeros(len(phase['s_old']))
            p=AuditProblem(read_json(OLD/run/'CONFIG.json'),record['strain'],lower)
            if iteration: xstart=read_npz(folder/f'iteration_{iteration-1:03d}_mechanical.npz')['x']
            else:xstart=read_npz(OLD/run/f'load_{step-1:03d}/accepted_state.npz')['x'] if step else p.family.X
            start=p.metrics(p.prescribe(xstart),phase['s_old'])
            before=p.metrics(phase['x'],phase['s_old']);after=p.metrics(phase['x'],phase['s_new'])
            increment=phase['s_new']-phase['s_old']
            row=dict(state_name=f'{run}_load_{step:03d}_iteration_{iteration:03d}',kind='phase_after',
                source=str(folder/f'iteration_{iteration:03d}_phase.npz'),load_step=step,iteration=iteration,
                strain=record['strain'],displacement_m=record['strain']*p.config['geometry']['grip_center_separation_m'],
                previous_accepted_phase_file=str(OLD/run/f'load_{step-1:03d}/accepted_state.npz') if step else None,
                phase_increment=float(abs(increment).max()),delta_s_L2=float(np.linalg.norm(increment)),
                mechanical_start_energy_J=start['total_J'],mechanical_end_energy_J=before['total_J'],
                mechanical_energy_change_J=before['total_J']-start['total_J'],
                phase_energy_change_J=after['total_J']-before['total_J'],
                mechanical_force_residual=before['force_residual'],**after)
            rows.append(row);index.append(row)
        for i,row in enumerate(rows):
            row.update(L2_ratio=None,Linf_ratio=None,increment_cosine=None)
            if i and rows[i-1]['load_step']==row['load_step']:
                previous=rows[i-1]
                a=read_npz(previous['source']);b=read_npz(row['source']);da=a['s_new']-a['s_old'];db=b['s_new']-b['s_old']
                if np.linalg.norm(da)>0 and np.linalg.norm(db)>0:
                    row.update(L2_ratio=float(np.linalg.norm(db)/np.linalg.norm(da)),
                        Linf_ratio=float(np.max(abs(db))/np.max(abs(da))),
                        increment_cosine=float(da@db/(np.linalg.norm(da)*np.linalg.norm(db))))
        failure=OLD/run/'FAILED.json'
        if failure.exists():
            event=read_json(failure);folder=OLD/run/f"load_{event['load_step']:03d}"
            last=sorted(folder.glob('iteration_*_mechanical.npz'))[-1];a=read_npz(last)
            lower=read_npz(OLD/run/f"load_{event['load_step']-1:03d}/accepted_state.npz")['s_pf_bond']
            p=AuditProblem(read_json(OLD/run/'CONFIG.json'),event['strain'],lower)
            r=dict(state_name=f'{run}_domain_exit',kind='mechanical_gate_stop',source=str(last),
                load_step=event['load_step'],iteration=int(last.stem.split('_')[1]),phase_increment=None,
                strain=event['strain'],displacement_m=event['strain']*p.config['geometry']['grip_center_separation_m'],
                previous_accepted_phase_file=str(OLD/run/f"load_{event['load_step']-1:03d}/accepted_state.npz"),
                **p.metrics(a['x'],a['s_pf_bond']))
            index.append(r);rows.append(r)
        paths[run]=rows
    save_json(OUT/'state_index/FAILURE_STATE_INDEX.json',dict(selected=settings['states'],states=index))
    save_json(OUT/'energy_path/SAVED_PATHS.json',paths)
    for key in settings['states']:
        p,data=load_audit_state(key);classification=p.classify(data['x'],data['s'],settings)
        np.savez_compressed(OUT/f'state_index/{key}.npz',x=data['x'],s=data['s'],lower=p.lower,
                            x_start=data['x_start'],**classification)
        save_json(OUT/f'state_index/{key}.json',dict(selection=settings['states'][key],source=data['path'],
            metrics=p.metrics(data['x'],data['s']),active_counts={k:int(v.sum()) for k,v in classification.items() if k!='q_scaled'}))
    return paths


def reproduce(task):
    p,data=load_saved_iteration(*task);start=time.perf_counter();r=p.phi(data['s'],data['x_start'])
    row=dict(run=task[0],load_step=task[1],iteration=task[2],status=r['status'],wall_seconds=time.perf_counter()-start,
             ADR_iterations=r['ADR_iterations'],force_residual=r['before']['force_residual'])
    if r['status']=='OK':row.update(max_s_error=float(abs(r['s_output']-data['s_next']).max()),
                                   max_x_error_m=float(abs(r['x_equilibrium']-data['x']).max()))
    return row


def jacobian_column(task):
    key,relative,column=task;p,data=load_audit_state(key);s=data['s'];a=audit_config()
    base=p.phi(s,data['x_start'])
    if base['status']!='OK':return dict(key=key,epsilon=relative,column=column,status=base['status'])
    h=relative*max(.1,abs(s[column]));lower_room=s[column]-p.lower[column];upper_room=p.sc-s[column]
    if min(lower_room,upper_room)>h: coefficients=[(-1.,-1/(2*h)),(1.,1/(2*h))];scheme='central'
    elif upper_room>2*h:coefficients=[(0.,-3/(2*h)),(1.,4/(2*h)),(2.,-1/(2*h))];scheme='forward_second_order'
    elif lower_room>2*h:coefficients=[(0.,3/(2*h)),(-1.,-4/(2*h)),(-2.,1/(2*h))];scheme='backward_second_order'
    else:return dict(key=key,epsilon=relative,column=column,status='NO_FEASIBLE_PERTURBATION')
    derivative=np.zeros(p.ns);trial_derivative=np.zeros(p.ns);switches=np.zeros(p.ns,bool);raw=[]
    for multiplier,weight in coefficients:
        if multiplier==0:r=base
        else:
            pert=s.copy();pert[column]+=multiplier*h;r=p.phi(pert,data['x_start'])
        if r['status']!='OK':return dict(key=key,epsilon=relative,column=column,status=r['status'],scheme=scheme,h=h)
        derivative+=weight*r['s_output'];trial_derivative+=weight*r['s_trial']
        switches|=r['active_set']!=base['active_set']
        raw.append(dict(multiplier=multiplier,ADR_iterations=r['ADR_iterations'],
                        force_residual=r['before']['force_residual'],min_green=r['before']['minimum_Green_strain']))
    return dict(key=key,epsilon=relative,column=column,status='OK',h=h,scheme=scheme,
                derivative=derivative,trial_derivative=trial_derivative,switches=switches,raw=raw)


def main():
    settings=audit_config();index_saved_states();budget=worker_budget()
    save_json(OUT/'provenance/FIXED_POINT_WORKER_BUDGET.json',budget)
    with ProcessPoolExecutor(max_workers=budget['max_workers']) as pool:
        rows=list(pool.map(reproduce,settings['map_reproduction_cases']))
    passed=all(r['status']=='OK' and r['max_s_error']<=settings['map_reproduction_s_atol'] and
               r['max_x_error_m']<=settings['map_reproduction_x_atol_m'] for r in rows)
    save_json(OUT/'fixed_point/REPRODUCTION.json',dict(status='VERIFIED_REPRODUCTION' if passed else 'FIXED_POINT_MAP_REPRODUCTION_FAILED',cases=rows))
    print('REPRODUCTION',passed,'max error',max(r.get('max_s_error',1.) for r in rows),flush=True)
    if not passed:raise SystemExit('STOP: fixed-point reproduction failed')
    tasks=[];bases={}
    for key in settings['states']:
        p,data=load_audit_state(key);r=p.phi(data['s'],data['x_start']);bases[key]=r
        np.savez_compressed(OUT/f'fixed_point/{key}_base.npz',**{k:r[k] for k in ['x_equilibrium','Y','s_input','s_output','s_trial','selector_gap','active_set']})
        save_json(OUT/f'fixed_point/{key}_base.json',dict(status=r['status'],ADR_iterations=r['ADR_iterations'],before=r['before'],after=r.get('after')))
        tasks.extend((key,e,col) for e in settings['jacobian_relative_epsilon'] for col in range(p.ns))
    grouped={};completed=0;began=time.perf_counter()
    with ProcessPoolExecutor(max_workers=budget['max_workers']) as pool:
        futures={pool.submit(jacobian_column,t):t for t in tasks}
        for future in as_completed(futures):
            r=future.result();key=(r['key'],r['epsilon']);grouped.setdefault(key,[]).append(r);completed+=1
            if completed%40==0:print('Jacobian columns',completed,'/',len(tasks),'elapsed',round(time.perf_counter()-began,1),flush=True)
            if len(grouped[key])==185:
                assemble_group(key,grouped[key],bases[key[0]],settings)
    save_json(OUT/'spectrum/COMPLETED.json',dict(columns=completed,wall_seconds=time.perf_counter()-began))


def assemble_group(key,columns,base,settings):
    name,e=key;columns=sorted(columns,key=lambda r:r['column']);p,data=load_audit_state(name)
    folder=OUT/'spectrum'/name/f'eps_{e:g}';folder.mkdir(parents=True,exist_ok=True)
    failures=[{k:v for k,v in r.items() if k not in ['derivative','trial_derivative','switches']} for r in columns if r['status']!='OK']
    save_json(folder/'COLUMN_METADATA.json',[{k:v for k,v in r.items() if k not in ['derivative','trial_derivative','switches']} for r in columns])
    if failures:
        save_json(folder/'SUMMARY.json',dict(status='SPECTRAL_ESTIMATE_UNRESOLVED',failure_columns=failures));return
    raw=np.column_stack([r['derivative'] for r in columns]);trial=np.column_stack([r['trial_derivative'] for r in columns])
    switches=np.column_stack([r['switches'] for r in columns]);update=base['active_set']=='INTERIOR'
    frozen=trial.copy();frozen[~update]=np.eye(p.ns)[~update]
    classification=p.classify(base['x_equilibrium'],data['s'],settings)
    marginal=(abs(base['selector_gap'])<=settings['map_branch_tolerance_s']) & ~classification['LOWER_BOUND_ACTIVE']
    nonsmooth=bool(switches.any() or marginal.any())
    # Do NOT eig the raw switched secant as an ordinary Jacobian. Store the
    # selected-branch matrix separately and mark its spectral interpretation.
    values,vectors=eig(frozen);order=np.argsort(abs(values))[::-1];values=values[order];vectors=vectors[:,order]
    np.savez_compressed(folder/'MATRICES.npz',raw_secant_NOT_ordinary_J=raw,trial_derivative=trial,
        frozen_selected_branch_J=frozen,switch_matrix=switches,marginal_selector=marginal,
        eigenvalues=values,eigenvectors=vectors)
    summary=dict(status='NONSMOOTH_ACTIVE_SET_SWITCH' if nonsmooth else 'ORDINARY_JACOBIAN_RESOLVED',
        ordinary_Phi_spectrum_valid=not nonsmooth,conditional_frozen_branch_only=nonsmooth,
        selected_branch_rho=float(abs(values[0])),dominant=dict(real=float(values[0].real),imag=float(values[0].imag),magnitude=float(abs(values[0]))),
        switched_columns=int(np.any(switches,axis=0).sum()),switched_entries=int(switches.sum()),
        marginal_selector_bonds=int(marginal.sum()),active_phase_relative_to_previous_load={k:int(v.sum()) for k,v in classification.items() if k!='q_scaled'},
        central_columns=sum(r['scheme']=='central' for r in columns),one_sided_columns=sum(r['scheme']!='central' for r in columns))
    save_json(folder/'SUMMARY.json',summary)
    print('SPECTRUM',name,e,summary['status'],'conditional rho',summary['selected_branch_rho'],flush=True)


if __name__=='__main__':main()
