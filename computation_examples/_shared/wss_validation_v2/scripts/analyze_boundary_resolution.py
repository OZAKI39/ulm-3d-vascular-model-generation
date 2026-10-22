"""Compare actual +/- pressure responses with measured terminal/time-step changes.

These observed differences are not certified error bounds or mesh independence.
Requires all pressure CFD cases and production WSS analyses to be complete.
"""
from pathlib import Path
import csv, json
import numpy as np
import pyvista as pv
from case_common import V, wss
from flow_solver_support.wss_case import coordinate_identity, material
from vessel_regions import masks
from analyze_vessel import stats, csvout

group=json.loads((V/'stage4/boundary_group.json').read_text())
base=Path(group['baseline_path'])
cases=[base,V/'stage4/O2_minus1pct',V/'stage4/O2_plus1pct']
metrics=['mean_Pa','p05_Pa','p50_Pa','p95_Pa','min_Pa','max_Pa','below1_area_pct','below2_area_pct','above30_area_pct']
changes=[]
for case in cases:
    ex=json.loads((case/'reports/execution.json').read_text())
    assert ex['status']=='PASS'
    assert json.loads((case/'reports/flow_quality.json').read_text())['accepted_final_and_log_checks']
    policy=json.loads((case/'policy.json').read_text());mu=material(case)['mu_Pa_s']
    mesh=np.load(case/'SV_MESH/mesh_arrays.npz')
    x,t,b,tags=[mesh[k] for k in ['points_m','tetra','boundary_triangles','facet_tags']]
    wall=b[tags==1];own=wss.boundary_owners(t,wall);c,area,n=wss.wall_geometry(x,t,wall,own);regions=masks(c)
    paths=sorted((case/'run'/f"{policy['MPI_ranks']}-procs").glob('result_*.vtu'),key=lambda p:int(p.stem.rsplit('_',1)[1]))
    full=[p for p in paths if int(p.stem.rsplit('_',1)[1])%5==0][-4:]
    selected=full+([paths[-1]] if paths[-1] not in full else [])
    assert len(full)==4 and int(paths[-1].stem.rsplit('_',1)[1])==ex['final_step']
    previous=None
    for path in selected:
        step=int(path.stem.rsplit('_',1)[1]);g=pv.read(path);coordinate_identity(g.points,x)
        u=np.asarray(g['Velocity']);tau=wss.tangential_traction(wss.p1_gradients(x,t[own],u),n,mu);w=np.linalg.norm(tau,axis=1)
        values={region:stats(w[sel],area[sel],c[sel]) for region,sel in regions.items()}
        if previous:
            old_step,old_values=previous
            for region in regions:
                for metric in metrics:
                    a,z=old_values[region][metric],values[region][metric]
                    changes.append(dict(case=case.name,mesh=group['mesh'],region=region,metric=metric,
                        from_step=old_step,to_step=step,time_interval_s=(step-old_step)*policy['dt_s'],
                        unit='percentage_points' if 'area_pct' in metric else 'Pa',
                        previous_value=a,current_value=z,signed_change=z-a,absolute_change=abs(z-a),
                        statistical_weight='same_wall_triangle_area',reference='actual_last_three_full_save_intervals_plus_final_partial_interval'))
        previous=step,values
csvout(V/'data/boundary_terminal_metric_changes.csv',changes)
responses=list(csv.DictReader((V/'data/boundary_region_responses.csv').open()))
controls=list(csv.DictReader((V/'data/vessel_control_metric_changes.csv').open()))
time_changes={(r['region'],r['metric']):abs(float(r['signed_absolute_change'])) for r in controls if r['comparison_type']=='time_step'}
rows=[]
for r in responses:
    key=r['region'],r['metric']
    terminal=max(z['absolute_change'] for z in changes if (z['region'],z['metric'])==key)
    dt=time_changes[key];minus=abs(float(r['minus_absolute_change']));plus=abs(float(r['plus_absolute_change']));small=min(minus,plus)
    rows.append(dict(mesh=group['mesh'],baseline_case=base.name,region=key[0],metric=key[1],unit=r['unit'],
        absolute_minus_response=minus,absolute_plus_response=plus,
        largest_observed_terminal_change_across_group=terminal,
        observed_coarse_grid_dt_to_half_dt_change_context_only=dt,
        smaller_response_over_terminal_change=small/terminal if terminal else '',
        smaller_medium_response_over_coarse_dt_change_context_only=small/dt if dt else '',
        statistical_weight='wall_triangle_area',
        interpretation='Terminal changes are same medium mesh; dt changes measured only on coarse mesh are context, not a medium-mesh error bound. Zero denominator undefined. No rigorous total-error bound'))
csvout(V/'data/boundary_numerical_scale_comparison.csv',rows)
print('Actual terminal differences:',len(changes),'response/observed numerical-scale comparisons:',len(rows))
