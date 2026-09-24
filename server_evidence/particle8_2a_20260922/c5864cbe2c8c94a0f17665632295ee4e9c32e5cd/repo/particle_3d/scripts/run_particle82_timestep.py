#!/usr/bin/env python3
"""Stage E: basin-stratified finite-size dt/dt2/dt4, including stopped IDs."""
from pathlib import Path
import argparse,json,sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
from particle_3d.particle82_batch import run_batch
from particle_3d.particle82_provenance import atomic_json,require_remote


def main(a):
    read=lambda p:json.loads(Path(p).read_text());host=read(a.host_provenance);require_remote(host,host['hostname'])
    previous=Path(a.previous);out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    admission=read(a.admission_audit);scaling=read(a.scaling_benchmark);extended=read(a.extended_audit)
    if not all(x['all_original_prefixes_exact'] for x in extended['factors']):raise ValueError('D prefix gate failed')
    meta={int(p.stem.split('_')[1]):read(p) for p in sorted((previous/'trajectories').glob('*.json'))}
    selected={};groups={}
    for basin in ['OUTLET_01','OUTLET_02','OUTLET_03']:
        ids=sorted(int(pid) for pid,r in admission['per_id'].items() if r['admitted_basin']==basin)
        chosen=ids[:30];groups[basin]=dict(available=len(ids),selected=chosen,target=30,
            target_status='MET' if len(chosen)==30 else 'INSUFFICIENT_ADMITTED_CASES_NO_MANUFACTURED_REPLACEMENTS')
        for pid in chosen:selected[pid]=basin
    stops=[pid for pid,m in meta.items() if m['end_reason']=='INTEGRATION_SAFETY_STOP' and pid not in selected][:30]
    groups['additional_safety_stops']=dict(selected=stops,target=30)
    for pid in stops:selected[pid]=admission['per_id'][str(pid)]['admitted_basin']
    ledger=read(previous/'data/birth_ledger.json');ledger['events']=[meta[pid]['birth_metadata'] for pid in sorted(selected)]
    atomic_json(out/'TIMESTEP_PRESELECTED_IDS.json',dict(groups=groups,criterion='Lowest stable IDs within admitted-position basin; no successful-outcome filtering',selected=selected))
    rows=[]
    for factor in [1,2,4]:
        config=dict(dt_s=.00025/factor,guard_factor=16,dataset_role='BASIN_STRATIFIED_TIMESTEP_DIAGNOSTIC_NOT_PRODUCTION_SELECTION')
        catalog,_=run_batch(ledger,out/f'dt_div_{factor}',scaling['chosen_workers'],1,2200,True,1,config,host)
        for record in catalog['records']:
            pid=record['particle_id'];ref=read(out/f'dt_div_1/trajectories/mb_{pid:06d}.json')
            x=np.load(out/f'dt_div_1/trajectories/mb_{pid:06d}.npz')['samples'];y=np.load(out/f'dt_div_{factor}/trajectories/mb_{pid:06d}.npz')['samples']
            end=min(x[-1,0],y[-1,0]);times=np.linspace(0,end,1001)
            xp=np.column_stack([np.interp(times,x[:,0],x[:,j]) for j in [1,2,3]])
            yp=np.column_stack([np.interp(times,y[:,0],y[:,j]) for j in [1,2,3]])
            rows.append(dict(particle_id=pid,point_basin=selected[pid],dt_divisor=factor,dt_s=config['dt_s'],
                baseline_outlet=ref['exit_outlet'],outlet=record['exit_outlet'],end_reason=record['end_reason'],
                residence_time_s=record['residence_time_s'],path_length_m=record['path_length_m'],
                last_elapsed_s=float(y[-1,0]),trajectory_deviation_max_m=float(np.linalg.norm(xp-yp,axis=1).max()),
                common_comparison_horizon_s=float(end),comparison_interpolation_only_no_extrapolation=True))
    atomic_json(out/'TIMESTEP_AUDIT.json',dict(groups=groups,rows=rows,production_timestep='NOT_FROZEN',guard_factor=16))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for n in ['previous','output-dir','host-provenance','admission-audit','scaling-benchmark','extended-audit']:p.add_argument('--'+n,required=True)
    main(p.parse_args())
