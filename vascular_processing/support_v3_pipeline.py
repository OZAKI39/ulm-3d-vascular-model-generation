"""Bounded V3 planning, actual local slicing and evidence-gated selection."""
from __future__ import annotations
import json, traceback, time, shutil, hashlib
from pathlib import Path
import numpy as np
import yaml
from scipy.spatial.transform import Rotation
from .support_inputs import ROOT, discover, save, table, verify, digest
from .support_v2_inputs import load
from .support_v2_orientation import bed_transform, fits
from .support_v3_planner import POLICY, plan, fingerprint


def gates(r):
    fail=[];s=r.get('slice',{});rem=r.get('removal',{})
    checks={'BUILD_VOLUME_EXCEEDED':r.get('build_volume_fit'), 'NON_RIGID_TRANSFORM':r.get('scale_factor')==1.,
        'ACTUAL_SUPPORT_COVERAGE_INCOMPLETE':s.get('actual_vascular_coverage',{}).get('all_required_covered'),
        'BOX_OVERHANG_UNRESOLVED':s.get('actual_box_coverage',{}).get('all_required_covered'),
        'SUPPORT_REQUIRES_SUPPORT':r.get('self_buildability',{}).get('unsupported_region_count')==0,
        'ACTUAL_SUPPORT_REQUIRES_SUPPORT':s.get('actual_self_buildability',{}).get('unsupported_region_count')==0,
        'ILLEGAL_VESSEL_PORT_INTERSECTION':r.get('vascular_intersection_mm3',1)>-1e-9 and r.get('vascular_intersection_mm3',1)<1e-5,
        'CUT_SECTION_CROSSES_UNWEAKENED_SUPPORT':r.get('cuts_valid'),
        'CUT_POINT_TOOL_ACCESS_FAILED':rem.get('tool_access_passed'),
        'TRAPPED_FRAGMENT':rem.get('trapped_pieces')==0 and rem.get('all_internal_exit_top'),
        'SLICER_VALIDATION_PENDING':s.get('slicer_success') and s.get('support_only_slicer_success'),
        'AUTO_SUPPORT_NOT_OFF':s.get('auto_support_off') and s.get('native_support_toolpaths')==0,
        'CRITICAL_FEATURE_DROPPED':s.get('all_features_preserved'),
        'UNRESOLVED_CRITICAL_WARNING':s.get('unresolved_warning_count')==0,
        'ASSEMBLY_STL_PRECISION_FAILURE':not r.get('assembly_export_error'),
        'SOURCE_PROTECTION_NOT_VERIFIED':r.get('source_unchanged')}
    checks['STALE_CACHED_GEOMETRY_EVIDENCE']=not r.get('evidence_stale',False)
    if 'exported_geometry_audit' in r:checks['EXPORTED_GEOMETRY_OR_BOX_ROOT_FAILURE']=r['exported_geometry_audit'].get('passed')
    return [k for k,v in checks.items() if not v]


def rank(r):
    s=r.get('slice',{});cov=s.get('actual_vascular_coverage',r.get('planned_coverage',{}))
    return (bool(gates(r)),not r.get('build_volume_fit',False),r.get('vascular_intersection_mm3',1)>1e-5,
        -cov.get('percent',0),r.get('self_buildability',{}).get('unsupported_region_count',9999),
        s.get('actual_near_vessel_projection_mm2',1e9),r.get('removal',{}).get('internal_fragment_count',9999),
        r.get('removal',{}).get('tool_access_failed',9999),r.get('support_volume_mm3',1e9),s.get('estimated_print_time_seconds') or 1e12)


def evidence_signature(inputs,row,profile,geometry_parameters):
    # Exclude descriptive labels such as UNVALIDATED versus UNCALIBRATED.
    # Numerical settings, physical preset contents, policy and code all enter
    # the cache. An old candidate is never trusted merely because it exists.
    keys=['printer','process_name','material','nozzle_diameter','layer_height','first_layer_height','line_width','roof_gap_mm','build_bounds_mm','bridge_settings','speeds','temperatures']
    value=dict(input_sha=inputs['hashes'],orientation=row,profile={k:profile.get(k) for k in keys},policy=POLICY,
        geometry_parameters=geometry_parameters,code_sha={n:digest(Path(__file__).with_name(n)) for n in ['support_v3_geometry.py','support_v3_planner.py']})
    return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()


def orientations(v2,mold,profile):
    prior=json.loads((v2/'orientation_candidates.json').read_text());names=['O000','O008','O096','O104','O092','O072']
    rows=[next(r for r in prior if r['orientation_id']==n) for n in names]
    for n,angle in [('N008',8),('N352',-8)]:
        rotation=Rotation.from_euler('x',angle,degrees=True).as_matrix();m=bed_transform(mold.vertices,rotation,profile['build_bounds_mm'])
        rows.append(dict(orientation_id=n,rotation=rotation.tolist(),print_transform=m.tolist(),casting_restore_transform=np.linalg.inv(m).tolist(),
            gravity_in_casting=rotation[2].tolist(),build_volume_fit=fits(mold.vertices,m,profile['build_bounds_mm']),origin='NEW_SMALL_UPRIGHT_TILT'))
    return rows


def worker(args):
    try:return plan(*args)
    except Exception as exc:
        folder,row=args[:2];r=dict(candidate_id=row['orientation_id'],folder=str(folder),orientation=row,error=str(exc),status='NEEDS_ADJUSTMENT')
        save(folder/'planning_failure.json',dict(**r,traceback=traceback.format_exc()));print('V3 PLANNING FAILURE',r['candidate_id'],str(exc),flush=True);return r


def run(source=None,output=None,stage='all',config=None):
    start=time.monotonic();inputs=discover(source=source);out=Path(output).resolve() if output else inputs['root'].parent/'vascular_support_planner/v3'
    if 'v2' in out.parts or 'legacy_explicit_support_v1' in out.parts:raise ValueError('V1_V2_ARE_READ_ONLY_HISTORY')
    out.mkdir(parents=True,exist_ok=True);v2=inputs['root'].parent/'vascular_support_planner/v2'
    inputs,profile,discovery,mold,v,box,branches,topology=load(out/'validation/inputs',source)
    profile['gap_status']='PHYSICAL_GAP_UNCALIBRATED';save(out/'profile_runtime.json',profile);save(out/'v3_policy_runtime.json',POLICY)
    cfg=yaml.safe_load(Path(config or ROOT/'config/vascular_support_planner.yaml').read_text());save(out/'resolved_config.json',cfg)
    protected=out/'validation/protected_history_before.json'
    if not protected.exists():
        hashes=json.loads((v2/'protected_before.json').read_text())
        for directory in [v2,v2.parent/'legacy_explicit_support_v1']:
            hashes.update({str(p):digest(p) for p in directory.rglob('*') if p.is_file()})
        save(protected,hashes)
    hashes=json.loads(protected.read_text());verify(hashes)
    rows=orientations(v2,mold,profile);save(out/'orientation_candidates.json',rows);table(out/'tables/v3_orientation_candidates.csv',rows)
    if stage=='orientation':return dict(status='NEEDS_ADJUSTMENT',orientation_count=len(rows))
    from .support_v3_validation import thin_wall_trial,slice_candidate
    trialpath=out/'sensitivity/line_widths/thin_wall_result.json'
    if stage in ['all','coupon'] or not trialpath.exists():trial=thin_wall_trial(out/'sensitivity/line_widths',profile,discovery)
    else:trial=json.loads(trialpath.read_text())
    # Actual one-line loss forces an explicit upgrade. A two-line body is the
    # stability search starting point only if it survives the local profile;
    # it never overrides a failed preservation check.
    available=[int(k) for k,val in trial.get('line_preservation',{}).items() if val]
    lines=min((n for n in available if n>=2),default=max(available,default=2))
    save(out/'sensitivity/line_selection.json',dict(available_actual_preserved_lines=available,selected_geometry_lines=lines,
        criterion='Among preserved lines, avoid one-line-only tall curtain where a wider printable alternative exists; strength uncalibrated',
        automatically_upgraded_from_one_line=1 not in available,physical_strength_validated=False))
    results=[];pending=[]
    for row in rows:
        folder=out/'candidates'/row['orientation_id'];path=folder/'candidate.json'
        if stage!='plan' and path.exists():
            r=json.loads(path.read_text());r['source_unchanged']=True
            if 'geometry_parameters' not in r:
                r['geometry_parameters']=dict(body_lines=r.get('line_count',lines),cut_neck_lines=r.get('actual_thin_wall_upgrade',{}).get('current_neck_lines',1))
            signature=evidence_signature(inputs,row,profile,r.get('geometry_parameters',{}))
            if r.get('evidence_signature') and r['evidence_signature']!=signature and stage not in ['report']:
                if stage in ['all','resume','geometry','research']:
                    archive=out/'history'/(row['orientation_id']+'_cache_'+r['evidence_signature'][:12])
                    if not archive.exists():shutil.copytree(folder,archive)
                    pending.append((folder,row,inputs,profile,cfg,mold,v,box,branches,topology,lines));continue
                raise ValueError('V3_EVIDENCE_FINGERPRINT_CHANGED_RUN_RESUME: '+row['orientation_id'])
            r['evidence_stale']=bool(r.get('evidence_signature') and r['evidence_signature']!=signature)
            if r.get('artifact_sha256') and any(digest(folder/n)!=sha for n,sha in r['artifact_sha256'].items()):
                raise ValueError('V3_CACHED_GEOMETRY_ARTIFACT_CHANGED: '+row['orientation_id'])
            results.append(r)
        elif stage in ['all','plan','resume','geometry','research']:
            pending.append((folder,row,inputs,profile,cfg,mold,v,box,branches,topology,lines))
    if pending:
        from concurrent.futures import ProcessPoolExecutor,as_completed
        import multiprocessing
        with ProcessPoolExecutor(max_workers=2,mp_context=multiprocessing.get_context('spawn')) as pool:
            futures=[pool.submit(worker,a) for a in pending]
            for f in as_completed(futures):
                r=f.result();r['source_unchanged']=True;results.append(r);save(out/'planning_results.json',results)
    if not results and stage not in ['all','plan','resume','geometry','research'] and (out/'planning_results.json').exists():results=json.loads((out/'planning_results.json').read_text())
    valid=[r for r in results if (Path(r['folder'])/'mold_with_custom_support_v3.stl').exists()]
    top4path=out/'top4_candidate_ids.json'
    if stage!='plan' and top4path.exists() and not pending:
        ids=json.loads(top4path.read_text());finalists=[next(r for r in valid if r['candidate_id']==identity) for identity in ids]
    else:
        # Compare like-for-like geometry at the selection stage. Mixing an
        # unsliced nominal score with a sliced actual score changes the cohort
        # on every resume and wastes additional full-machine evaluations.
        finalists=sorted(valid,key=lambda r:rank({**r,'slice':{},'removal':{}}))[:4]
        save(top4path,[r['candidate_id'] for r in finalists])
    if stage in ['all','removal','refine']:
        from .support_v3_removal import run as removal
        for r in finalists:
            try:r['removal']=removal(r,inputs)
            except Exception as exc:
                r['removal']=dict(error=str(exc));save(Path(r['folder'])/'removal_failure.json',dict(error=str(exc),traceback=traceback.format_exc()))
            save(Path(r['folder'])/'candidate.json',r);save(out/'planning_results.json',results)
    if stage in ['all','slice']:
        # Sequential: the Windows adapter uses a shared staging namespace.
        for r in finalists:
            try:
                actual_profile={**profile,'cut_neck_lines':r.get('geometry_parameters',{}).get('cut_neck_lines',1)}
                r['slice']=slice_candidate(r,actual_profile,discovery)
                lost=[f for f in r['slice'].get('feature_preservation',{}).get('rows',[]) if not f['preserved']]
                if any(f['kind']=='CUT_NECK' for f in lost) and actual_profile['cut_neck_lines']==1:
                    folder=Path(r['folder']);archive=out/'history'/(r['candidate_id']+'_before_actual_one_line_upgrade')
                    if not archive.exists():shutil.copytree(folder,archive)
                    upgraded_profile={**profile,'cut_neck_lines':2}
                    new=plan(folder,r['orientation'],inputs,upgraded_profile,cfg,mold,v,box,branches,topology,3)
                    new['actual_thin_wall_upgrade']=dict(reason='Actual curved one-line neck preservation failed',failed_features=lost,previous_neck_lines=1,current_neck_lines=2,body_lines=3)
                    from .support_v3_removal import run as removal
                    new['removal']=removal(new,inputs);new['slice']=slice_candidate(new,upgraded_profile,discovery)
                    r.clear();r.update(new)
            except Exception as exc:
                r['slice']=dict(slicer_success=False,error=str(exc));save(Path(r['folder'])/'slice_failure.json',dict(error=str(exc),traceback=traceback.format_exc()))
            save(Path(r['folder'])/'candidate.json',r);save(out/'planning_results.json',results)
    for r in results:
        if (Path(r['folder'])/'custom_support_v3.stl').exists():
            from .support_v3_verify import geometry_audit
            geometry_audit(r,inputs)
        path=Path(r['folder'])/'support_removal_plan.json'
        if path.exists():r['removal']=json.loads(path.read_text())
        path=Path(r['folder'])/'bambu_validation.json'
        if path.exists():r['slice']=json.loads(path.read_text())
        r['source_unchanged']=True;r['failed_gates']=gates(r)
        if not r.get('evidence_stale'):r['evidence_signature']=evidence_signature(inputs,r['orientation'],profile,r.get('geometry_parameters',{}))
        r['artifact_sha256']={n:digest(Path(r['folder'])/n) for n in ['mold_without_support.stl','custom_support_v3.stl','mold_with_custom_support_v3.stl'] if (Path(r['folder'])/n).exists()}
        save(Path(r['folder'])/'candidate.json',r)
    protection=verify(hashes);save(out/'validation/source_protection.json',protection)
    save(out/'planning_results.json',results)
    from .support_v3_review import generate
    return generate(out,results,profile,inputs,protection,time.monotonic()-start)
