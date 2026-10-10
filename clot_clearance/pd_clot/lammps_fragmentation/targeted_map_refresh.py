"""Merge completed postprocessing and list only real old/new map points."""
import csv,json
from .config import dump
from .failure_map import csv_out
from .targeted_contact_transfer import O,OLD

def collect():
    summaries={p.parent.name:json.loads(p.read_text()) for p in sorted((O/'analysis').glob('*/summary.json'))}
    summary_rows=[]
    for name,q in summaries.items():
        q.update(first_break_times=dict(bulk_s=q['first_bulk_break_time_s'],interface_s=q['first_interface_break_time_s']),broken_bulk_fraction=q['final']['bulk_broken_fraction'],broken_interface_fraction=q['final']['interface_broken_fraction'],debris_fraction=q['final']['debris_fraction'],largest_attached_fraction=q['final']['largest_attached_component_fraction'])
        for d in ['200','300']:
            q['detached_ge'+d+'_count_peak']=q['thresholds'][d]['MAX_SIMULTANEOUS_COUNT'];q['detached_ge'+d+'_count_final']=q['thresholds'][d]['AT_FINAL_FRAME'];q['detached_ge'+d+'_mass_fraction_peak']=q['thresholds'][d]['PEAK_DETACHED_MASS_FRACTION']
        dump(O/'analysis'/name/'summary.json',q)
        if (O/'CASES'/name).exists():dump(O/'CASES'/name/'FRAGMENT_SIZE_RESULT.json',q)
        for d,t in q['thresholds'].items():summary_rows.append(dict(case_id=name,contact_model_id=q['contact_model_id'],pattern_id=q['pattern_id'],tau_budget_Pa=q['tau_budget_Pa'],threshold_um=int(d),**t,largest_ever_detached_d_eq_um=q['largest_ever_detached_d_eq_um'],largest_final_detached_d_eq_um=q['largest_final_detached_d_eq_um'],thin_detached_ge300_episodes=q['thin_detached_ge300_episodes']))
    csv_out(O/'FRAGMENT_SIZE_SUMMARY.csv',summary_rows)
    for source,target,detached_only in [('components.csv','PER_FRAME_DETACHED_COMPONENTS.csv',True),('lifetimes.csv','FRAGMENT_LIFETIMES.csv',False),('timeseries.csv','PER_FRAME_SIZE_AND_CONTACT_METRICS.csv',False),('split_relations.csv','FRAGMENT_SPLIT_RELATIONS.csv',False)]:
        writer=None
        with (O/target).open('w',newline='') as dst:
            for name,q in summaries.items():
                with (O/'analysis'/name/source).open() as src:
                    reader=csv.DictReader(src)
                    for row in reader:
                        if detached_only and row.get('attachment')!='DETACHED':continue
                        row['contact_model_id']=q['contact_model_id']
                        if 'node_count' in row:row['resolution_label']='NUMERICALLY_UNRESOLVED_DEBRIS' if int(row['node_count'])<=2 else 'DISCRETE_COMPONENT_NOT_PHYSICALLY_SIZE_VALIDATED'
                        if writer is None:writer=csv.DictWriter(dst,fieldnames=list(row));writer.writeheader()
                        writer.writerow(row)
    points=[]
    for path in sorted((OLD/'cases').glob('*/result.json')):
        q=json.loads(path.read_text())
        if q['pattern_id'] not in ['LOCAL_GAUSSIAN_SHEAR','OPPOSING_SHEAR'] or q['nominal_duration_factor']!=2:continue
        points.append(dict(case_id=q['case_id'],source_path=str(path.parent),contact_model_id='NATIVE_1P35DX',pattern_id=q['pattern_id'],tau_budget_Pa=q['tau_budget'],duration_over_Tcross=q['duration_over_Tcross'],terminal_mode=q['failure_mode'],bulk_broken_fraction=q['bulk_broken_fraction'],interface_broken_fraction=q['interface_broken_fraction'],debris_fraction=q['debris_fraction'],execution='HISTORICAL_REUSED'))
    for name,q in summaries.items():
        if q['contact_model_id']!='MODIFIED_1P0DX':continue
        points.append(dict(case_id=name,source_path=q['source_path'],contact_model_id=q['contact_model_id'],pattern_id=q['pattern_id'],tau_budget_Pa=q['tau_budget_Pa'],duration_over_Tcross=q['pulse_duration_s']/0.0005801744510205713,terminal_mode=q['final']['failure_mode'],bulk_broken_fraction=q['broken_bulk_fraction'],interface_broken_fraction=q['broken_interface_fraction'],debris_fraction=q['debris_fraction'],execution='NEW_THIS_ROUND' if (O/'CASES'/name).exists() else 'HISTORICAL_REUSED'))
    csv_out(O/'TARGETED_OLD_NEW_MAP_COMPARISON.csv',points)
    dump(O/'LOCAL_CONTACT_MODEL_TRANSITION_COMPARISON.json',dict(points=points,unmeasured_regions='UNKNOWN',interpolation=False,old_map_overwritten=False,scope='Only the two selected patterns at 2Tcross; no new duration boundary inferred'))
    dump(O/'ALL_CASE_SUMMARIES.json',summaries)
    perf={}
    for name,q in summaries.items():
        perf[name]=dict(execution='NEW_THIS_ROUND' if (O/'CASES'/name).exists() else 'HISTORICAL_REUSED',wall_seconds=q['wall_seconds'],original_postprocess_wall_s=q['original_postprocess_wall_s'],new_size_postprocess_wall_s=q['size_postprocess_wall_s'])
        source=__import__('pathlib').Path(q['source_path']);perf[name]['solver_performance']=json.loads((source/'PERFORMANCE.json').read_text())
        perf[name]['raw_file_bytes']={p.name:p.stat().st_size for p in source.iterdir() if p.name in ['trajectory.lammpstrj','broken_bonds.local','surviving_bonds.local']}
    dump(O/'PERFORMANCE_REAL_CASES.json',dict(cases=perf,new_formal_solves=sum(v['execution']=='NEW_THIS_ROUND' for v in perf.values()),new_solver_wall_sum_s=sum(v['wall_seconds'] for v in perf.values() if v['execution']=='NEW_THIS_ROUND'),size_analysis_wall_sum_s=sum(v['new_size_postprocess_wall_s'] for v in perf.values()),baseline_reruns=0,calibration_runs=0,new_builds=0,performance_benchmarks=0,timing_scope='Historical runtime is not new cost. Some read-only analyses overlap with new case solves; sums are not elapsed round time.'))
    print('COLLECT_COMPLETE',len(summaries),len(points),flush=True)
