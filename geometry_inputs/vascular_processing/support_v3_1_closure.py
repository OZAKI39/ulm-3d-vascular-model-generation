"""Actual Bambu residual -> local preflight -> reslice, at most five repairs."""
from __future__ import annotations
import json,pickle,shutil,copy
from pathlib import Path
from .support_inputs import save,table,digest
from .support_v3_1_pipeline import copy_baseline,evaluate
from .support_v3_1_planning import write_geometry
from .support_v3_1_residual import cluster,iteration_stop,area_tolerance
from .support_v3_1_validation import semantic_fields,gates
from .support_v3_1_repair import historical_repairs,repair_grains


def read_rows(path):
    import csv
    rows=[]
    with Path(path).open() as f:
        for row in csv.DictReader(f):
            for key in ('z_mm','uncovered_area_mm2','area_mm2'):
                if key in row:row[key]=float(row[key])
            if 'covered' in row:row['covered']=row['covered']=='True'
            rows.append(row)
    return rows


def compact(record):
    s=record['slice'];r=record['removal'];combined=s['actual_vascular_coverage']
    return dict(candidate_id=record['candidate_id'],folder=record['folder'],
        coverage_percent=combined['percent'],uncovered_area_mm2=combined['uncovered_area_sum_mm2'],
        residual_regions=combined['uncovered_region_count'],support_self=s['actual_self_buildability']['unsupported_region_count'],
        trapped=r['trapped_pieces'],inaccessible_cuts=r['tool_access_failed'],
        semantic_coverage=record['semantic_coverage'],failed_gates=record['failed_gates'],
        actual_slice_success=s['slicer_success'],auto_support_off=s['auto_support_off'])


def run_closure(base,inputs,profile,discovery,baseline,report_only=False):
    out=base/'v3_1';out.mkdir(exist_ok=True);summary=out/'closure_summary.json'
    if summary.exists():return json.loads(summary.read_text())
    if report_only:raise ValueError('PHASE_B_HAS_NOT_FINISHED_RUN_ALL_TO_RESUME')
    records=[];areas=[];stop=None
    for iteration in range(6):
        folder=out/f'iteration_{iteration:02}';folder.mkdir(exist_ok=True);(folder/'QC').mkdir(exist_ok=True)
        completed=folder/'iteration_complete.json'
        if completed.exists():
            record=json.loads((folder/'candidate.json').read_text());records.append(record);areas.append(record['slice']['actual_vascular_coverage']['uncovered_area_sum_mm2'])
            stop=json.loads(completed.read_text()).get('stop_reason')
            if stop:break
            continue
        if iteration==0:
            record=copy_baseline(baseline['folder'],folder,'V3_1_O000_BASELINE')
            table(folder/'repair_actions.csv',[],fields=['grain_id','repair','accepted','reason'])
        else:
            previous=records[-1];old=Path(previous['folder'])
            record=copy_baseline(old,folder,f'V3_1_O000_ITER{iteration:02}')
            # Preserve evidence needed for local preflight; evaluate replaces
            # it with an independently audited final result after repair.
            record['removal']=previous['removal'];record['slice']=previous['slice']
            for name in ('tool_access.json','support_removal_plan.json'):
                if (old/name).exists():shutil.copy2(old/name,folder/name)
            with (old/'rib_geometry_cache.pkl').open('rb') as f:ribs=pickle.load(f)
            baseline_ribs=copy.deepcopy(ribs)
            with (old/'bambu/actual_support_footprints.pkl').open('rb') as f:actual=pickle.load(f)
            _,_,_,branches=semantic_fields(old,inputs)
            rows=read_rows(old/'vascular_residuals.csv')
            grains=cluster(rows,branches,profile)
            self_rows=list(previous['slice']['actual_self_buildability']['rows'])
            for row in previous.get('geometry_support_self_actual_coverage',{}).get('rows',[]):
                if not row['covered']:
                    self_rows.append(dict(z_mm=profile['first_layer_height']+row['layer_index']*profile['layer_height'],
                        polygon_wkt=row['uncovered_wkt'],area_mm2=row['uncovered_area_mm2'],source='GEOMETRY_RS_MINUS_PREVIOUS_ACTUAL_PATHS'))
            grains+=cluster(self_rows,branches,profile,True)
            table(folder/'residual_grains.csv',grains);save(folder/'residual_grains.json',grains)
            historical=[]
            if iteration==1:ribs,historical=historical_repairs(ribs,record,inputs,profile,folder)
            ribs,actions=repair_grains(ribs,grains,record,inputs,profile,actual,folder)
            table(folder/'repair_actions.csv',historical+actions)
            record['repair_iteration']=dict(index=iteration,grain_count=len(grains),attempts=len(historical)+len(actions),accepted=sum(a['accepted'] for a in historical+actions))
            write_geometry(record,ribs,profile,baseline_ribs)
            reverted={x['rib_id'] for x in record.get('post_export_reverted_ribs',[])}
            for action in historical+actions:
                if action.get('accepted') and action.get('rib_id') in reverted:
                    action.update(accepted=False,accepted_preflight=True,reason='NEW_EDIT_REVERTED_FOR_EXPORT_PRECISION')
            record['repair_iteration']['accepted_preflight']=record['repair_iteration']['accepted']
            record['repair_iteration']['accepted']=sum(a['accepted'] for a in historical+actions)
            table(folder/'repair_actions.csv',historical+actions)
        record=evaluate(record,inputs,profile,discovery,phase_a=False)
        records.append(record);areas.append(record['slice']['actual_vascular_coverage']['uncovered_area_sum_mm2'])
        all_coverage=all(c['all_required_covered'] for c in record['semantic_coverage'].values()) and record['slice']['actual_self_buildability']['unsupported_region_count']==0
        stop=iteration_stop(areas,all_coverage,max_repairs=5)
        save(completed,dict(iteration=iteration,stop_reason=stop,metrics=compact(record),
            support_sha256=digest(folder/'custom_support.stl'),archive_sha256=digest(folder/'X1C_project.3mf')))
        print('V3.1 CLOSURE',iteration,compact(record),'stop',stop,flush=True)
        if stop:break
    first,last=records[0],records[-1];before,after=compact(first),compact(last)
    winner=last if not gates(last,False) else None
    result=dict(executed=True,fallback_orientation='O000',fallback_reason='No inverted candidate passed all hard gates; preserve V3 segmented rib architecture',
        iteration_count=len(records)-1,iterations=[compact(r) for r in records],stop_reason=stop,
        coverage_before=before['coverage_percent'],coverage_after=after['coverage_percent'],
        uncovered_area_before=before['uncovered_area_mm2'],uncovered_area_after=after['uncovered_area_mm2'],
        residual_regions_before=before['residual_regions'],residual_regions_after=after['residual_regions'],
        support_self_before=before['support_self'],support_self_after=after['support_self'],
        trapped_before=before['trapped'],trapped_after=after['trapped'],
        inaccessible_cut_before=before['inaccessible_cuts'],inaccessible_cut_after=after['inaccessible_cuts'],
        slicer_result=last['slice']['slicer_success'],status='READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_PARAMETERS_UNCALIBRATED' if winner else 'NEEDS_ADJUSTMENT',
        final_folder=last['folder'],winner=winner,physical_status=['PHYSICAL_GAP_UNCALIBRATED','BRIDGE_LIMIT_UNCALIBRATED','BREAK_FORCE_UNCALIBRATED','ABS_RIB_STRENGTH_UNCALIBRATED'])
    residuals=read_rows(Path(last['folder'])/'vascular_residuals.csv')
    table(out/'residual_failure_map.csv',residuals)
    save(out/'residual_failure_map.json',dict(required_regions=residuals,support_self=last['slice']['actual_self_buildability'],
        trapped=[r for r in last['removal']['fragments'] if r['status']=='TRAPPED_FRAGMENT'],failed_gates=last['failed_gates']))
    save(summary,result);return result
