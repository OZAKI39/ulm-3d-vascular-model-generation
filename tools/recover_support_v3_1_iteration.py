"""Rebuild the reviewed first edit batch after the support-export guard fix.

The failed real slice is archived, not overwritten or silently discarded.
All formerly accepted grains and all self-support grains are re-screened.
Previously rejected vascular/port edits remain in the archived evidence;
the next normal iteration recomputes ALL actual residuals again.
"""
from pathlib import Path
import sys,json,pickle,shutil,copy
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from vascular_processing.support_v3_1_pipeline import preparation,copy_baseline,evaluate
from vascular_processing.support_v3_1_repair import historical_repairs,repair_grains
from vascular_processing.support_v3_1_planning import write_geometry
from vascular_processing.support_v3_1_closure import compact
from vascular_processing.support_inputs import save,table,digest


def main():
    c=preparation();base,out,inputs,profile,discovery=c[:5];b=base/'v3_1';folder=b/'iteration_01'
    archive=b/'history/export_guard_failed_iteration_01';archive.parent.mkdir(exist_ok=True)
    if not archive.exists():shutil.move(folder,archive)
    partial=b/'iteration_02'
    if partial.exists() and not (partial/'iteration_complete.json').exists():
        target=b/'history/interrupted_before_second_iteration'
        if not target.exists():shutil.move(partial,target)
    baseline=b/'iteration_00';r=copy_baseline(baseline,folder,'V3_1_O000_ITER01')
    old=json.loads((baseline/'candidate.json').read_text());r['removal']=old['removal'];r['slice']=old['slice']
    for name in ('tool_access.json','support_removal_plan.json'):shutil.copy2(baseline/name,folder/name)
    with (baseline/'rib_geometry_cache.pkl').open('rb') as f:ribs=pickle.load(f)
    baseline_ribs=copy.deepcopy(ribs)
    with (baseline/'bambu/actual_support_footprints.pkl').open('rb') as f:actual=pickle.load(f)
    previous=json.loads((archive/'repair_preflight.json').read_text());accepted={x['grain_id'] for x in previous if x['accepted']}
    allgrains=json.loads((archive/'residual_grains.json').read_text())
    grains=[g for g in allgrains if g['grain_id'] in accepted or g['component']=='SUPPORT_SELF']
    save(folder/'residual_grains.json',allgrains);table(folder/'residual_grains.csv',allgrains)
    ribs,historical=historical_repairs(ribs,r,inputs,profile,folder)
    ribs,actions=repair_grains(ribs,grains,r,inputs,profile,actual,folder)
    write_geometry(r,ribs,profile,baseline_ribs)
    reverted={x['rib_id'] for x in r.get('post_export_reverted_ribs',[])}
    for a in historical+actions:
        if a['accepted'] and a['rib_id'] in reverted:a.update(accepted=False,accepted_preflight=True,reason='NEW_EDIT_REVERTED_FOR_EXPORT_PRECISION')
    table(folder/'repair_actions.csv',historical+actions)
    r['repair_iteration']=dict(index=1,grain_count=len(allgrains),rescreened_grains=len(grains),attempts=len(historical)+len(actions),accepted=sum(a['accepted'] for a in historical+actions),
        failed_export_attempt_preserved=str(archive),scope='Rebuild reviewed batch; rejected vascular/port alternatives preserved in archived preflight. All actual residuals recomputed in next iteration.')
    r=evaluate(r,inputs,profile,discovery,False)
    save(folder/'iteration_complete.json',dict(iteration=1,stop_reason=None,metrics=compact(r),support_sha256=digest(folder/'custom_support.stl'),archive_sha256=digest(folder/'X1C_project.3mf')))
    print('RECOVERED_FIRST_ITERATION',compact(r),flush=True)


if __name__=='__main__':main()
