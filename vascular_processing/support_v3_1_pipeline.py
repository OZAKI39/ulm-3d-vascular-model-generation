"""Bounded two-phase V3.1 experiment. Historical directories are read-only."""
from __future__ import annotations
import json, pickle, shutil, traceback, time
from pathlib import Path
import numpy as np
import yaml
from .support_inputs import ROOT, discover, save, table, digest, verify
from .support_v2_inputs import load
from .support_v3_1_orientation import orientations, surrogate, surrogate_key
from .support_v3_1_planning import plan_inverted
from .support_v3_1_validation import semantic_audit, removal_audit, gates


def preparation(source=None,output=None):
    inputs=discover(source=source);base=inputs['root'].parent/'vascular_support_planner'
    if output and Path(output).resolve()!=base:
        raise ValueError('V3_1_OUTPUT_MUST_BE_SUPPORT_PLANNER_PARENT; only inverted_study and v3_1 are writable')
    out=base/'inverted_study';out.mkdir(parents=True,exist_ok=True)
    protected=out/'validation/protected_before.json'
    if not protected.exists():
        hashes=json.loads((base/'v3/validation/protected_history_before.json').read_text())
        for p in base.rglob('*'):
            if p.is_file() and not any(s in ('inverted_study','v3_1') for s in p.relative_to(base).parts):hashes[str(p)]=digest(p)
        save(protected,hashes)
    verify(json.loads(protected.read_text()))
    inputs,profile,discovery,mold,v,box,branches,top=load(out/'validation/inputs',source)
    profile.update(cut_neck_lines=2,gap_status='PHYSICAL_GAP_UNCALIBRATED');save(out/'profile_runtime.json',profile)
    summary=json.loads((base/'v3/vascular_support_v3_summary.json').read_text())
    baseline=summary['diagnostic_selected']
    # Read the independent records and transforms, not just embedded summary.
    bf=Path(baseline['folder'])
    actual=json.loads((bf/'bambu_validation.json').read_text());candidate=json.loads((bf/'candidate.json').read_text())
    transforms={n:json.loads((bf/n).read_text()) for n in ('print_transform_v3.json','casting_restore_transform_v3.json')}
    report=(base/'v3/vascular_support_v3_report.md').read_text()
    save(out/'validation/baseline_read.json',dict(summary_sha=digest(base/'v3/vascular_support_v3_summary.json'),
        report_sha=digest(base/'v3/vascular_support_v3_report.md'),candidate_sha=digest(bf/'candidate.json'),
        candidate_id=baseline['candidate_id'],actual_bambu=actual,transforms=transforms,production_selected=summary['production_selected']))
    rows=orientations(mold,profile,baseline['orientation']);save(out/'orientation_candidates.json',rows);table(out/'orientation_matrix.csv',rows)
    metrics=[]
    for row in rows:
        path=out/'surrogates'/f"{row['orientation_id']}.json"
        if path.exists():r=json.loads(path.read_text())
        else:r=surrogate(row,mold,v,box,branches,inputs,profile);save(path,r)
        metrics.append(r)
    save(out/'surrogate_metrics.json',metrics)
    public=[{k:v for k,v in r.items() if k not in ('visibility','direct_targets_print','direct_clear_mask')} for r in metrics]
    table(out/'surrogate_metrics.csv',public);table(out/'box_ceiling_burden.csv',[{k:v for k,v in r.items() if k.startswith('box_') or k=='orientation_id'} for r in public])
    ids=[r['orientation_id'] for r in sorted(metrics[1:],key=surrogate_key)[:6]]
    ids=list(dict.fromkeys(ids+['INV_0','UPRIGHT_O000']));save(out/'full_candidate_ids.json',ids)
    return base,out,inputs,profile,discovery,mold,v,box,branches,top,rows,baseline,metrics


def planning_worker(args):
    folder,row,*rest=args
    try:
        path=folder/'candidate.json'
        if path.exists():
            old=json.loads(path.read_text())
            if old.get('v3_1_planning_complete'):return old
        r=plan_inverted(folder,row,*rest);r['v3_1_planning_complete']=True;save(path,r)
        print('V3.1 PLAN COMPLETE',r['candidate_id'],flush=True);return r
    except Exception as exc:
        save(folder/'v3_1_planning_failure.json',dict(error=str(exc),traceback=traceback.format_exc()))
        raise


def copy_baseline(source,folder,identity=None):
    """Copy, never hardlink: slicer/auditor writes cannot mutate old artifacts."""
    source=Path(source);folder.mkdir(parents=True,exist_ok=True)
    # Reuse immutable geometry/detector evidence only. Actual slice is run in
    # the new folder so old absolute cache paths cannot cause writes to V3.
    for name in ['mold_without_support.stl','vascular_with_ports.stl','box.stl','custom_support_v3.stl',
                 'mold_with_custom_support_v3.stl','support_ribs.stl','support_roofs.stl','support_cut_necks.stl',
                 'support_exterior.stl','support_pillars.stl','rib_geometry_cache.pkl','frozen_detection_cache.pkl',
                 'cut_points.csv','print_transform_v3.json','casting_restore_transform_v3.json']:
        if (source/name).exists():shutil.copy2(source/name,folder/name)
    r=json.loads((source/'candidate.json').read_text());r['folder']=str(folder)
    if identity:r['candidate_id']=identity
    r.update(slice={},removal={},source_unchanged=True,v3_1_planning_complete=True)
    shutil.copy2(folder/'custom_support_v3.stl',folder/'custom_support.stl')
    save(folder/'candidate.json',r);return r


def evaluate(r,inputs,profile,discovery,phase_a=True):
    from .support_v3_validation import slice_candidate
    from .support_v3_verify import geometry_audit
    from .support_v3_removal import run as normal_removal
    folder=Path(r['folder']);r['source_unchanged']=True
    geometry_audit(r,inputs);save(folder/'candidate.json',r)
    if phase_a:removal_audit(r,inputs)
    else:r['removal']=normal_removal(r,inputs)
    save(folder/'candidate.json',r)
    r['slice']=slice_candidate(r,profile,discovery)
    if r['slice'].get('slicer_success'):semantic_audit(r,inputs,profile)
    r['failed_gates']=gates(r,phase_a);r['status']='NEEDS_ADJUSTMENT' if r['failed_gates'] else 'READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_PARAMETERS_UNCALIBRATED'
    save(folder/'candidate.json',r);return r


def phase_a(context,cfg,stage):
    from concurrent.futures import ProcessPoolExecutor
    import multiprocessing
    base,out,inputs,profile,discovery,mold,v,box,branches,top,rows,baseline,metrics=context
    ids=json.loads((out/'full_candidate_ids.json').read_text());records=[];futures={}
    with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
        for name in ids:
            folder=out/'candidates'/name;row=next(r for r in rows if r['orientation_id']==name)
            if name=='UPRIGHT_O000':
                if (folder/'candidate.json').exists():r=json.loads((folder/'candidate.json').read_text())
                else:r=copy_baseline(baseline['folder'],folder,name);r['orientation']=row;save(folder/'candidate.json',r)
                futures[name]=r
            else:futures[name]=pool.submit(planning_worker,(folder,row,inputs,profile,cfg,mold,v,box,branches,top))
        for name in ids:
            r=futures[name] if isinstance(futures[name],dict) else futures[name].result()
            if stage not in ('plan','geometry'):
                r=evaluate(r,inputs,profile,discovery,True)
                print('PHASE A ACTUAL',name,r['slice'].get('actual_vascular_coverage',{}).get('percent'),r['failed_gates'],flush=True)
            records.append(r);save(out/'full_candidates.json',records)
    if stage in ('plan','geometry'):return records,None
    good=[r for r in records if r['candidate_id']!='UPRIGHT_O000' and not r['failed_gates']]
    winner=min(good,key=lambda r:(-r['direct_buildplate_support_fraction'],r['plate_release']['level'],r['slice'].get('contact_area_mm2',1e9),r['support_volume_mm3'],r.get('box_ceiling_rib_count',0),r['slice'].get('estimated_print_time_seconds',1e12))) if good else None
    save(out/'phase_a_result.json',dict(status='INVERTED_STUDY_PASS' if winner else 'INVERTED_STUDY_FAILED',winner=winner['candidate_id'] if winner else None,records=records))
    return records,winner


def production(record,dest):
    dest.mkdir(parents=True,exist_ok=True);folder=Path(record['folder'])
    names={'mold_without_support.stl':'mold_without_support.stl','custom_support_v3.stl':'custom_support_v3_1.stl',
           'mold_with_custom_support_v3.stl':'mold_with_custom_support_v3_1.stl','X1C_project.3mf':'X1C_project_v3_1.3mf','support_removal_plan.json':'support_removal_plan.json'}
    for old,new in names.items():shutil.copy2(folder/old,dest/new)
    for key,name in [('print_transform','print_transform_v3_1.json'),('casting_restore_transform','casting_restore_transform_v3_1.json')]:
        save(dest/name,dict(transform_4x4=record['orientation'][key],scale_factor=1.))


def run(source=None,output=None,stage='all',config=None):
    began=time.monotonic();context=preparation(source,output);base,out,inputs,profile,discovery,*_=context
    if stage=='orientation':return dict(status='ORIENTATION_STUDY_PREPARED',output=str(out))
    cfg=yaml.safe_load(Path(config or ROOT/'config/vascular_support_planner.yaml').read_text())
    if stage in ('report','resume','refine') and (out/'phase_a_result.json').exists():
        a=json.loads((out/'phase_a_result.json').read_text());records=a['records'];winner=next((r for r in records if r['candidate_id']==a['winner']),None)
    else:records,winner=phase_a(context,cfg,stage)
    if stage in ('plan','geometry'):return dict(status='PHASE_A_GEOMETRY_COMPLETE',output=str(out))
    if winner:
        production(winner,out/'production_selected');b=dict(executed=False,reason='PHASE_A_PASSED')
    else:
        from .support_v3_1_closure import run_closure
        b=run_closure(base,inputs,profile,discovery,context[11],report_only=stage=='report')
        if b.get('winner'):production(b['winner'],base/'v3_1/production_selected')
    protection=verify(json.loads((out/'validation/protected_before.json').read_text()));save(out/'validation/source_protection.json',protection)
    from .support_v3_1_review import generate
    return generate(context,records,winner,b,protection,time.monotonic()-began)
