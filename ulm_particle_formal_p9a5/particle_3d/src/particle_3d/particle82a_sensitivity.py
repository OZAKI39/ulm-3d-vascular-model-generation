"""Predeclared numerical checks, without changing the production cohort."""
from pathlib import Path
import argparse,json,multiprocessing as mp,socket,time
import numpy as np
from .particle82a_admission import context,method_c
from .particle82a_trajectories import iter_rows
from .particle82a_pipeline import BASINS,initialize,run_map
from .particle82_provenance import atomic_json,require_remote


def check_event(task):
    row,config=task;c=context();event=row['event']
    path=c.native.trace(event['anchor_m'],step_m=config['point_step_m'])['path']
    result=method_c(event,path,c,tolerance=config['search_tolerance_m'],horizon_factor=config['horizon_factor'])
    base=row['methods']['C']
    return dict(event_id=event['event_id'],anchor_basin=event['point_tracer_basin'],
        base_status=base['status'],new_status=result['status'],status_equal=base['status']==result['status'],
        base_s_birth_m=base.get('s_birth_m'),new_s_birth_m=result.get('s_birth_m'),
        distance_difference_m=(result['s_birth_m']-base['s_birth_m']) if base['accepted'] and result['accepted'] else None,
        radius_m=event['first_radius_m'],anchor_unchanged=result['anchor_m']==event['anchor_m'],
        radius_unchanged=result['radius_m']==event['first_radius_m'])


def main():
    p=argparse.ArgumentParser();p.add_argument('--admission',required=True);p.add_argument('--output',required=True)
    p.add_argument('--provenance',required=True);p.add_argument('--workers',type=int,default=8);a=p.parse_args()
    output=Path(a.output);output.mkdir(parents=True,exist_ok=True)
    prov=json.loads(Path(a.provenance).read_text());require_remote(prov,prov['hostname']);context()
    groups={b:[] for b in BASINS}
    for row in iter_rows(a.admission):
        b=row['event']['point_tracer_basin']
        if len(groups[b])<64:groups[b].append(row)
        if all(len(v)==64 for v in groups.values()):break
    rows=[row for b in BASINS for row in groups[b]]
    selection=dict(rule='FIRST_64_COMMON_EVENT_IDS_PER_ANCHOR_BASIN_INCLUDING_UNRESOLVED_NO_C_OUTCOME_SELECTION',
        ids={b:[r['event']['event_id'] for r in groups[b]] for b in BASINS})
    atomic_json(output/'SELECTION.json',selection)
    variants={'point_step_half':dict(point_step_m=.1e-6,search_tolerance_m=1e-9,horizon_factor=1.),
        'search_tolerance_half':dict(point_step_m=.2e-6,search_tolerance_m=.5e-9,horizon_factor=1.),
        'search_horizon_double':dict(point_step_m=.2e-6,search_tolerance_m=1e-9,horizon_factor=2.)}
    summaries={}
    for name,config in variants.items():
        start=time.time()
        with mp.get_context('fork').Pool(a.workers) as pool:results=pool.map(check_event,[(row,config) for row in rows])
        summaries[name]=dict(config=config,status_changes=sum(not r['status_equal'] for r in results),
            maximum_distance_difference_m=max([abs(r['distance_difference_m']) for r in results if r['distance_difference_m'] is not None],default=None),
            rows=results,seconds=time.time()-start)
        atomic_json(output/(name+'.json'),summaries[name])
    config=dict(point_step_m=.1e-6,search_tolerance_m=1e-9,horizon_factor=1.)
    initialize(output/'refined_basin_map',prov,config);run_map(24,a.workers)
    base=np.load(Path(a.admission)/'map_n24.npz')['rows'];refined=np.load(output/'refined_basin_map/map_n24.npz')['rows']
    if not np.array_equal(base[:,:5],refined[:,:5]):raise ValueError('Basin refinement must use identical quadrature points/weights')
    changed=base[:,5]!=refined[:,5]
    basin_audit=dict(same_points_and_weights=True,changed_labels=int(changed.sum()),
        changed_flux_fraction=float(base[changed,4].sum()/base[:,4].sum()),
        original_unresolved_flux_fraction=float(base[base[:,5]==3,4].sum()/base[:,4].sum()),
        refined_unresolved_flux_fraction=float(refined[refined[:,5]==3,4].sum()/refined[:,4].sum()))
    atomic_json(output/'SENSITIVITY_SUMMARY.json',dict(selection=selection,
        variants={k:{key:value for key,value in v.items() if key!='rows'} for k,v in summaries.items()},
        point_basin_refinement=basin_audit,REMOTE_SERVER_COMPUTE=True,hostname=socket.gethostname(),
        source_commit=prov['source_git_commit'],production_ledger_unchanged=True,
        interpretation='Stratified diagnostic only, not a replacement natural-cohort estimator'))


if __name__=='__main__':main()
