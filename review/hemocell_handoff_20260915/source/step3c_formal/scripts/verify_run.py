#!/usr/bin/python3
"""Offline replay of gates and independent flux/CV reconstruction from raw evidence."""
import csv
import json
import math
from pathlib import Path
import sys
sys.dont_write_bytecode=True
import numpy as np
from remote_common import read,write,contract,verify_bundle,verify_run_inputs,history,fields,field_path,PORTS,OFFSETS
from convergence import COLUMNS,initial_state,evaluate_point,field_residual,check_observed_safety

def equivalent(a,b):
    """Strict roundoff agreement; booleans, gates, keys and iteration IDs exact."""
    if isinstance(a,dict):return isinstance(b,dict) and a.keys()==b.keys() and all(equivalent(a[k],b[k]) for k in a)
    if isinstance(a,list):return isinstance(b,list) and len(a)==len(b) and all(equivalent(x,y) for x,y in zip(a,b))
    if isinstance(a,(float,np.floating)) and isinstance(b,(int,float)) and not isinstance(b,bool):return math.isclose(a,b,rel_tol=2e-11,abs_tol=1e-30)
    return a==b

def reconstruct_flux(quad,sampled,shape,dx,dt,rho_phys):
    """Independent vectorized trilinear reconstruction, no Palabos or solver call."""
    nx,ny,nz=shape
    if np.any(np.diff(sampled['index'])<=0):raise ValueError('Sample coordinates not unique/sorted')
    q=np.zeros(24);m=np.zeros(24)
    for start in range(0,len(quad),16384):
        rows=quad[start:start+16384];base=np.floor(rows[:,1:4]).astype(np.int64);frac=rows[:,1:4]-base
        velocity=np.zeros((len(rows),3));density=np.zeros(len(rows))
        for x in [0,1]:
            for y in [0,1]:
                for z in [0,1]:
                    ijk=base+[x,y,z]
                    if np.any(ijk<0) or np.any(ijk>=shape):raise ValueError('Quadrature leaves frozen lattice')
                    ids=(ijk[:,2]*ny+ijk[:,1])*nx+ijk[:,0]
                    where=np.searchsorted(sampled['index'],ids)
                    if np.any(where>=len(sampled)) or not np.array_equal(sampled['index'][where],ids):raise ValueError('Missing sampled atomic cells')
                    w=(frac[:,0] if x else 1-frac[:,0])*(frac[:,1] if y else 1-frac[:,1])*(frac[:,2] if z else 1-frac[:,2])
                    velocity+=w[:,None]*sampled['u'][where];density+=w*sampled['rho'][where]
        values=rows[:,4]*np.sum(velocity*rows[:,5:8],axis=1)*dx**3/dt
        group=rows[:,0].astype(int)
        q+=np.bincount(group,weights=values,minlength=24)
        m+=np.bincount(group,weights=rho_phys*density*values,minlength=24)
    return q,m

def verify(root,run):
    verify_bundle(root);verify_run_inputs(root,run);c=contract(root)
    status=read(run/'diagnostics/solver_status.json');actual=int(status['timesteps']);maxsteps=c['execution']['max_steps']
    if not (0<actual<=maxsteps and actual%5000==0):raise ValueError('Not a complete long validation result')
    if status['runtime_safety']!='PASS' or status['safety_checks']!=actual+1 or status['mpi_ranks']!=1 or status['cell_count']!=0:raise ValueError('Native safety/MPI/pure-fluid evidence failed')
    if status['all_step_fluid_rho_min']<.99 or status['all_step_fluid_rho_max']>1.01 or status['all_step_ghost_rho_min']<=0 or status['all_step_max_mach_including_ghosts']>.05:raise ValueError('Native all-step safety failed')
    data=history(run);check_observed_safety(data,c)
    if len(data)!=actual//100+1 or int(data[-1]['iteration'])!=actual:raise ValueError('History incomplete')
    with (run/'diagnostics/convergence_history.csv').open() as f:stored_csv=list(csv.DictReader(f))
    if len(stored_csv)!=actual//5000:raise ValueError('Evaluation history incomplete')
    state=initial_state();last_detail=None;last_row=None
    for row_index,step in enumerate(range(5000,actual+1,5000)):
        part=data[data['iteration']<=step]
        fv=field_residual(run,step,c) if step>=c['execution']['first_possible_evaluation_step'] else None
        row,state,detail=evaluate_point(part,step,c,state,fv)
        stored=read(run/f'diagnostics/evaluations/evaluation_{step}.json')
        if not equivalent({'row':row,'state':state,'detail':detail},stored):raise ValueError(f'Offline gate replay differs at {step}')
        for name in COLUMNS:
            text='' if row[name] is None else str(row[name])
            observed=stored_csv[row_index][name]
            if isinstance(row[name],(float,np.floating)):
                if not equivalent(row[name],float(observed)):raise ValueError(f'Convergence CSV differs at {step}/{name}')
            elif observed!=text:raise ValueError(f'Convergence CSV differs at {step}/{name}')
        if detail['stop'] and step!=actual:raise ValueError('Ignored early-stop confirmation')
        last_detail=detail;last_row=row
    confirmed=state['confirmation_iteration'] is not None
    if bool(status['auto_converged'])!=confirmed:raise ValueError('Driver/evaluator convergence disagreement')
    if not confirmed and actual!=maxsteps:raise ValueError('Run stopped before confirmation/hard ceiling')
    expected_status='AUTO_PASS_HUMAN_PENDING' if confirmed else 'NOT_CONVERGED_WITHIN_VALIDATION_HORIZON'
    if status['status']!=expected_status:raise ValueError('Driver final status incorrect')

    mf=np.column_stack([data[f'mass_outward_g{g}'] for g in [2,8,14,20]])
    averaged=.5*(mf[:-1]+mf[1:]);delta_mass=np.diff(data['control_volume_mass'])/(100*c['dt_s'])
    denom=np.maximum.reduce([np.abs(averaged[:,0]),np.abs(averaged[:,1:]).sum(axis=1),np.full(len(averaged),c['rho_kg_m3']*c['Qtarget_m3_s']*1e-12)])
    transient=np.abs(delta_mass+averaged.sum(axis=1))/denom
    if not np.allclose(transient,data['control_volume_mass_balance'][1:],rtol=1e-8,atol=1e-8):raise ValueError('Transient CV residual reconstruction differs')
    if not np.allclose(data['relative_mass_drift'],(data['total_mass']-data['total_mass'][0])/data['total_mass'][0],rtol=1e-11,atol=1e-15):raise ValueError('Mass drift column inconsistent')
    # Verify all reached sparse events, plus final/candidate/confirmation, against
    # independent atomic-sample reconstruction. No requirement to save per-step PDFs.
    with (run/'diagnostics/snapshot_events.csv').open() as f:events=list(csv.DictReader(f))
    event_steps=sorted({int(e['iteration']) for e in events});actual_set=set(event_steps)
    mandatory={x for x in c['execution']['field_snapshots'] if x<=actual}|{actual}
    if not mandatory<=actual_set:raise ValueError('Missing sparse snapshot milestones')
    if not any(int(e['iteration'])==actual and e['reason']=='final' for e in events):raise ValueError('Missing final snapshot marker')
    quad=np.loadtxt(run/'contracts/multiplane_quadrature.tsv')
    geom=read(run/'contracts/geometry_reuse_contract.json');dx=c['dx_effective_m'];dt=c['dt_s'];rho=c['rho_kg_m3']
    with np.load(run/'contracts/control_volume_nodes.npz') as mask:
        physical_ids=mask['flat_index'];cv_ids=physical_ids[mask['control_volume_mask']!=0]
    checks=[]
    for step in event_steps:
        stored=data[data['iteration']==step]
        if len(stored)!=1:raise ValueError('Sparse snapshot lacks scalar sample')
        scalar=stored[0];f=fields(field_path(run,step));sampled=fields(field_path(run,step,True))
        order=np.argsort(f['index']);f=f[order]
        if not np.array_equal(f['index'],physical_ids):raise ValueError('Snapshot physical mask differs')
        where=np.searchsorted(f['index'],cv_ids)
        mass=float(f['rho'][where].sum()*rho*dx**3)
        total=float(f['rho'].sum()*rho*dx**3)
        if not np.isclose(mass,scalar['control_volume_mass'],rtol=2e-11,atol=1e-25):raise ValueError('CV mass snapshot/monitor mismatch')
        if not np.isclose(total,scalar['total_mass'],rtol=2e-11,atol=1e-25):raise ValueError('Total mass snapshot/monitor mismatch')
        q,m=reconstruct_flux(quad,sampled,np.array(geom['lattice_shape']),dx,dt,rho)
        expected_q=np.array([scalar[f'{p}_{k}dx'+('_dx2' if refinement else '')] for p in PORTS for k in OFFSETS for refinement in [0,1]])
        expected_q[:6]*=-1
        expected_m=np.array([scalar[f'mass_outward_g{g}'] for g in range(24)])
        q_error=float(np.max(np.abs(q-expected_q))/c['Qtarget_m3_s'])
        m_error=float(np.max(np.abs(m-expected_m))/(rho*c['Qtarget_m3_s']))
        if q_error>1e-9 or m_error>1e-9:raise ValueError(f'Independent flux reconstruction failed at {step}')
        checks.append({'iteration':step,'volume_flux_max_abs_error_over_Qtarget':q_error,
          'mass_flux_max_abs_error_over_rho_Qtarget':m_error,'control_mass_kg':mass,'status':'PASS'})
    sampler=read(run/'diagnostics/quadrature_ownership_check.json')
    if sampler['status']!='PASS' or sampler['duplicate_owned_points'] or sampler['unowned_points']:raise ValueError('GPU bulk quadrature ownership failed')
    owner=read(run/'diagnostics/control_volume_ownership_check.json')
    if owner['status']!='PASS' or owner['CV_CELL_GLOBAL_COUNT']!=180543 or owner['physical_fluid_global_count']!=182694:raise ValueError('GPU physical/CV ownership failed')
    summary={'status':expected_status,'scope':'REMOTE_NUMERICAL_EVIDENCE_ONLY; local Git/source integrity audit pending',
      'auto_converged':confirmed,'actual_steps':actual,'final_physical_time_s':actual*dt,
      'first_converged_iteration':state['first_converged_iteration'],'candidate_iteration':state['candidate_iteration'],
      'confirmation_iteration':state['confirmation_iteration'],'stop_iteration':actual,'final_row':last_row,'final_evaluation':last_detail,
      'runtime_safety':'PASS','independent_flux_and_mass_checks':checks,'offline_confirmation_replay':'PASS',
      'failed_convergence_gates':last_row['failed_gates'].split(';') if last_row['failed_gates'] else [],
      'persistent_mean_backflow':any(v<0 for v in last_detail['windows']['long']['Qout']),
      'mpi1':'PASS','mpi2':'NOT_VALIDATED','restart':'NOT_IMPLEMENTED','human_review':'PENDING'}
    write(run/'diagnostics/numerical_summary.json',summary)
    return summary

if __name__=='__main__':
    root=Path(sys.argv[1]).resolve();run=Path(sys.argv[2]).resolve()
    try:
        result=verify(root,run);print(json.dumps({'status':result['status'],'steps':result['actual_steps'],'verification':'PASS'},indent=2))
    except Exception as error:
        write(run/'diagnostics/verification_failure.json',{'status':'FAIL_NUMERICAL_EVIDENCE','error':str(error)})
        raise
