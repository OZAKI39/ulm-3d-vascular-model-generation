#!/usr/bin/env python3
"""Limited P6.5 timestep sensitivity evidence, not a production convergence claim."""
from pathlib import Path
import sys,multiprocessing
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from particle_3d.particle81_simulation import OUTPUT,DT,environment,integrate_one,dump
from particle_3d.particle8_replay import read,digest


def work(task):
    event,factor=task;root=OUTPUT/'data/sensitivity'/f'dt_div_{factor}'
    result=integrate_one(event,dt=DT/factor,output=root)
    print(event['particle_id'],factor,result['end_reason'],flush=True)
    return result


def main():
    births=read(OUTPUT/'data/birth_ledger.json')['events'];ids=[1,4,7,13]
    tasks=[(births[i-1],factor) for i in ids for factor in [2,4]]
    environment()
    with multiprocessing.get_context('fork').Pool(2) as pool:list(pool.imap_unordered(work,tasks))
    rows=[]
    for pid in ids:
        base=read(OUTPUT/f'trajectories/mb_{pid:06d}.json');a=np.load(OUTPUT/base['samples_path'])['samples']
        for factor in [2,4]:
            root=OUTPUT/'data/sensitivity'/f'dt_div_{factor}';meta=read(root/f'trajectories/mb_{pid:06d}.json');b=np.load(root/meta['samples_path'])['samples']
            common=np.linspace(0,min(a[-1,0],b[-1,0]),301)
            pa=np.column_stack([np.interp(common,a[:,0],a[:,j]) for j in [1,2,3]])
            pb=np.column_stack([np.interp(common,b[:,0],b[:,j]) for j in [1,2,3]])
            rows.append(dict(particle_id=pid,dt_factor=factor,dt_s=DT/factor,base_end=base['end_reason'],refined_end=meta['end_reason'],
                base_completed=base['completed'],refined_completed=meta['completed'],same_birth_metadata=base['birth_metadata']==meta['birth_metadata'],
                same_inlet=np.array_equal(a[0,1:4],b[0,1:4]),all_finite=bool(np.isfinite(b).all()),
                base_elapsed_s=float(a[-1,0]),refined_elapsed_s=float(b[-1,0]),
                common_time_position_difference_max_m=float(np.linalg.norm(pa-pb,axis=1).max()),
                residence_difference_s=meta['residence_time_s']-base['residence_time_s'] if meta['completed'] and base['completed'] else None,
                base_samples_sha256=base['samples_sha256'],refined_samples_sha256=meta['samples_sha256'],
                refined_metadata_path=str((root/f'trajectories/mb_{pid:06d}.json').relative_to(OUTPUT)),failure_detail=meta['failure_detail']))
    dump(OUTPUT/'data/timestep_sensitivity.json',dict(scope='FOUR_PRESELECTED_IDS; TWO_COMPLETED_AND_TWO_NUMERICALLY_STOPPED_AT_BASE_DT',
        production_convergence_established=False,rows=rows,physical_cohort_count_excludes_these_reintegrations=True,
        comparison_interpolation='LINEAR_COMMON_ELAPSED_TIME_ONLY_FOR_ERROR_METRIC; NO_EXTRAPOLATION; NOT_NEW_PHYSICAL_DATA'))


if __name__=='__main__':main()
