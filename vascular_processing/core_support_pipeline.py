"""s1-11 orchestration; all s1-10 files remain read-only."""
from __future__ import annotations
import json,pickle,shutil,traceback
from pathlib import Path
import numpy as np
import trimesh
from .core_support_inputs import default_output,discover_source,freeze,current_profiles,branches,search
from .core_support_geometry import generate,bed_roots
from .core_support_slicing import evaluate
from .support_inputs import save,table,digest,verify,ROOT


def hard_gates(record,evaluation,removal,protection):
    c=evaluation.get('coverage',{});s=evaluation.get('support_self',{});l=evaluation.get('layer_coverage',{})
    gates=dict(SOURCE_UNCHANGED=protection['all_unchanged'],RIGID_ONLY=np.allclose(np.asarray(record['orientation']['rotation']).T@np.asarray(record['orientation']['rotation']),np.eye(3),atol=1e-9),CENTERLINE_100=c.get('passed',False),PORTS_100=c.get('port_underbody_percent')==100.,ZERO_UNCOVERED_ARC=c.get('max_uncovered_arc_mm')==0.,LAYER_DOMAIN_100=l.get('passed',False),SUPPORT_SELF_ZERO=s.get('self_unsupported_regions')==0,NO_VASCULAR_INTERSECTION=record['intersection_mm3']<1e-7,NO_WRAP=removal.get('wrap_traps',1)==0,ALL_NOTCHES_ACCESSIBLE=removal.get('accessible_notch_fraction')==1.,ALL_SEGMENTS_REMOVABLE=removal.get('removable_segment_fraction')==1.,SAFE_VIRTUAL_CUTS=all(r['passed'] for r in removal.get('cut_checks',[{'passed':False}])),AUTO_SUPPORT_OFF=evaluation.get('auto_support_off',False),REAL_BAMBU_SLICE=evaluation.get('slice',{}).get('slicer_success',False),VASCULAR_FEATURES_PRESERVED=evaluation.get('feature_preservation',{}).get('passed',False),NATIVE_GEOMETRY_UNCHANGED=evaluation.get('archive_geometry',{}).get('passed',False),WARNINGS_RESOLVED=not evaluation.get('unresolved_warnings',['MISSING_SLICE']))
    gates['RIGID_ONLY']=bool(gates['RIGID_ONLY'] and np.isclose(np.linalg.det(np.asarray(record['orientation']['rotation'])),1.,atol=1e-9))
    gates['SEPARATE_VASCULAR_AND_SUPPORT_PARTS']=evaluation.get('independent_part_count')==2
    gates['SUPPORT_ROOTS_ON_BED']=bool(record.get('bed_roots')) and all(x['bed_rooted'] for x in record['bed_roots'])
    return gates,[k for k,v in gates.items() if not v]


def run(output=None,stage='all'):
    out=Path(output or default_output()).resolve()
    for name in ['source','orientations','candidates','selected','tables','QC','bambu','reports','validation']:(out/name).mkdir(parents=True,exist_ok=True)
    inp=discover_source();hashes=freeze(inp,out)
    if stage in ['all','prepare'] or not (out/'source/context.pkl').exists():
        profile,actual,engine=current_profiles(inp,out/'bambu/profiles');bs,top=branches(inp,min(profile['line_width'],profile['support_line_width'])/2)
        c=dict(inputs=inp,profile=profile,actual=actual,engine=engine,branches=bs,topology=top)
        with (out/'source/context.pkl').open('wb') as f:pickle.dump(c,f)
        ranking=search(trimesh.load_mesh(inp['core']),bs,profile,out)
    else:c=pickle.load((out/'source/context.pkl').open('rb'));ranking=json.loads((out/'orientations/ranking.json').read_text())
    if stage=='prepare':return dict(status='PREPARED')
    records=[];errors=[];mesh=trimesh.load_mesh(inp['core'])
    for row in ranking[:6]:
        for suffix,gap,lines,roof in [('R20_C2',.2,2,2),('R20_C3',.2,3,2),('R25_C2',.25,2,2)]:
            folder=out/'candidates'/(row['orientation_id']+'_'+suffix)
            try:
                if (folder/'geometry.json').exists():r=json.loads((folder/'geometry.json').read_text())
                elif stage in ['all','geometry']:r=generate(mesh,c['branches'],row,c['profile'],folder,gap,lines,roof)
                else:continue
                records.append(r)
            except Exception as e:
                errors.append(dict(candidate_id=folder.name,error=str(e),traceback=traceback.format_exc()));print('CANDIDATE ERROR',folder.name,str(e),flush=True)
    error_path=out/'reports/generation_errors.json'
    prior_errors=json.loads(error_path.read_text()) if error_path.exists() else []
    error_map={x['candidate_id']:x for x in prior_errors}
    error_map.update({x['candidate_id']:x for x in errors})
    for r in records:error_map.pop(r['candidate_id'],None)
    save(error_path,list(error_map.values()));save(out/'tables/geometry_candidates.json',records)
    if stage=='geometry':return dict(status='GEOMETRY_COMPLETE',candidate_count=len(records),errors=errors)
    # Top three coarse orientations, all required line/gap variants, receive
    # real native slices; additionally evaluate any adaptive repair candidate.
    ids={r['orientation_id'] for r in ranking[:3]};todo=[r for r in records if r['orientation']['orientation_id'] in ids]
    for p in sorted((out/'candidates').glob('*/geometry.json')):
        r=json.loads(p.read_text())
        if r['candidate_id'] not in {x['candidate_id'] for x in records}:todo.append(r)
    evaluations=[]
    for r in todo:
        folder=Path(r['folder']);path=folder/'evaluation.json'
        if stage=='review' and path.exists():e=json.loads(path.read_text())
        elif stage in ['all','slice']:e=evaluate(r,c['actual'],c['engine'])
        else:continue
        if not e.get('coverage'):continue
        evaluations.append((r,e))
    if not evaluations:raise RuntimeError('NO_SUCCESSFUL_NATIVE_SLICE')
    for r,e in evaluations:
        if 'bed_roots' not in r or any('signed_volume_mm3' not in x for x in r['bed_roots']):
            support=trimesh.load_mesh(Path(r['folder'])/'custom_support.stl')
            r['bed_roots']=bed_roots(support)
            save(Path(r['folder'])/'geometry.json',r)
    # Bounded repair studies are part of the reproducible entry point, not
    # one-off manual output files. Both opposite build directions are tried.
    if stage in ['all','slice'] and not any(e.get('passed') for _,e in evaluations):
        directions=[]
        for row in ranking[:6]:
            normal=np.asarray(row['rotation'])[2]
            if not any(np.allclose(normal,n,atol=1e-8) for n in directions):
                directions.append(normal)
                index=len(directions)-1
                suffix='R20_R3_C3_LAYER' if index==0 else 'R25_C2_LAYER_BASE'
                gap,lines=(.2,3) if index==0 else (.25,2)
                folder=out/'candidates'/(row['orientation_id']+'_'+suffix)
                if folder.name not in {r['candidate_id'] for r,e in evaluations}:
                    rr=json.loads((folder/'geometry.json').read_text()) if (folder/'geometry.json').exists() else generate(mesh,c['branches'],row,c['profile'],folder,gap,lines,lines,reinforce=True,wide_foot=index>0)
                    ee=evaluate(rr,c['actual'],c['engine']);evaluations.append((rr,ee))
                if len(directions)==2:break
        adaptive=[(r,e) for r,e in evaluations if r.get('layer_reinforcement_count',0)>0 and e.get('coverage')]
        if adaptive:
            parent,pe=min(adaptive,key=lambda x:(x[1]['layer_coverage']['uncovered_region_count'],x[1]['support_self']['self_unsupported_regions'],x[1]['coverage']['unsupported_samples']))
            repairfolder=out/'candidates'/(parent['orientation']['orientation_id']+f"_R{int(parent['gap_mm']*100)}_C{parent['curtain_lines']}_REPAIR2")
            if repairfolder.name not in {r['candidate_id'] for r,e in evaluations}:
                from .core_support_repair import repair
                rr=json.loads((repairfolder/'geometry.json').read_text()) if (repairfolder/'geometry.json').exists() else repair(parent['folder'],repairfolder)
                ee=evaluate(rr,c['actual'],c['engine']);evaluations.append((rr,ee))
    # Pre-sort diagnostics by uncovered targets/regions. Final readiness is
    # solely the hard gates below, never this diagnostic ranking.
    evaluations.sort(key=lambda x:(x[1]['coverage']['unsupported_samples'],x[1]['layer_coverage']['uncovered_region_count'],x[1]['support_self']['self_unsupported_regions'],x[0]['intersection_mm3'],x[0]['support_volume_mm3']))
    from .core_support_removal import audit as removal_audit
    from .core_support_travel import audit as travel_audit
    detailed=[]
    review_jobs=evaluations[:3]
    adaptive=[x for x in evaluations if x[0].get('layer_reinforcement_count',0)>0]
    if adaptive:
        best=min(adaptive,key=lambda x:(x[1]['coverage']['unsupported_samples'],x[1]['layer_coverage']['uncovered_region_count'],x[1]['support_self']['self_unsupported_regions']))
        if best[0]['candidate_id'] not in {r['candidate_id'] for r,e in review_jobs}:review_jobs.append(best)
    for r,e in review_jobs:
        folder=Path(r['folder']);rp=folder/'removal.json';tp=folder/'travel_near_miss.json'
        removal=json.loads(rp.read_text()) if rp.exists() else removal_audit(r)
        travel=json.loads(tp.read_text()) if tp.exists() else travel_audit(e['slice']['archive'],folder,c['profile']['line_width'])
        gates,failures=hard_gates(r,e,removal,verify(hashes));detailed.append((r,e,removal,travel,gates,failures))
    detailed.sort(key=lambda x:(bool(x[5]),len(x[5]),x[1]['coverage']['unsupported_samples'],x[3]['near_miss_segments'],x[0]['support_volume_mm3'],x[1]['slice'].get('estimated_print_time_seconds',1e9)))
    r,e,removal,travel,gates,failures=detailed[0];folder=Path(r['folder']);status='NEEDS_ADJUSTMENT' if failures else 'READY_FOR_HUMAN_REVIEW_WITH_REMOVAL_UNCALIBRATED'
    for src,dest in [('vascular_core.stl','frozen_vascular_core.stl'),('custom_support.stl','custom_vascular_support.stl'),('vascular_core_plus_support_review.stl','vascular_core_plus_support_review.stl'),('vascular_print_transform.json','print_transform.json'),('inverse_transform.json','inverse_transform.json')]:shutil.copy2(folder/src,out/'selected'/dest)
    shutil.copy2(folder/'vascular_core.stl',out/'selected/vascular_core.stl');shutil.copy2(folder/'custom_support.stl',out/'selected/custom_support.stl');shutil.copy2(folder/'vascular_print_transform.json',out/'selected/vascular_print_transform.json')
    filename='BG001_RMCA_vascular_core_supported_X1C.3mf';shutil.copy2(e['slice']['archive'],out/'selected'/filename)
    protection=verify(hashes);protection.update(mismatch_count=len(protection['changed']));save(out/'validation/source_protection.json',protection)
    s=dict(status=status,production_selected=not failures,diagnostic_selected=bool(failures),selected_candidate=r['candidate_id'],input=dict(vascular_path=str(inp['core']),vascular_hash=digest(inp['core']),SWC_hash=digest(inp['smoothed_swc']),compensated_SWC_hash=digest(inp['compensated_swc']),port_manifest_hash=digest(inp['geometry'])),orientation=dict(candidate_count=len(ranking),full_candidate_count=len(records),sliced_candidate_count=len(evaluations),selected=r['orientation']['orientation_id'],transform=r['orientation']['transform'],scale_factor=1.),coverage=e['coverage'],layer_coverage=e['layer_coverage'],support=dict(roof_gap=r['gap_mm'],roof_width_lines=r['roof_width_lines'],curtain_lines=r['curtain_lines'],segment_count=r['segment_count'],support_volume=r['support_volume_mm3'],self_unsupported_regions=e['support_self']['self_unsupported_regions'],intersection_mm3=r['intersection_mm3'],gap_status='PHYSICAL_REMOVAL_GAP_UNCALIBRATED'),removal={k:v for k,v in removal.items() if k not in ['fragments','access','cut_checks']},bambu=dict(slicer_success=e['slice']['slicer_success'],auto_support_off=e['auto_support_off'],unresolved_warnings=e['unresolved_warnings'],feature_preservation=e['feature_preservation'],estimated_time=e['slice'].get('estimated_print_time_seconds'),filament=e['slice'].get('filament_used_g'),archive=str(out/'selected'/filename)),travel={k:v for k,v in travel.items() if k!='rows'},gates=gates,failures=failures,source_protection=protection,no_printer_upload_or_start=True,physical_validation=False)
    save(out/'vascular_core_support_summary.json',s);save(out/'tables/candidate_comparison.json',[dict(candidate=r['candidate_id'],orientation=r['orientation']['orientation_id'],gap_mm=r['gap_mm'],roof_lines=r['roof_width_lines'],curtain_lines=r['curtain_lines'],volume_mm3=r['support_volume_mm3'],intersection_mm3=r['intersection_mm3'],**{k:e[k] for k in ['coverage','layer_coverage','support_self','feature_preservation']}) for r,e in evaluations])
    save(out/'tables/removal_scorecard.json',[dict(candidate=rr['candidate_id'],gap_mm=rr['gap_mm'],roof_width_mm=rr['roof_width_mm'],roof_width_lines=rr['roof_width_lines'],curtain_lines=rr['curtain_lines'],support_volume_mm3=rr['support_volume_mm3'],wrap_angle_p95=rm['wrap_angle_p95'],segment_count=rm['segment_count'],median_segment_length_mm=rm['median_segment_length_mm'],max_segment_length_mm=rm['max_segment_length_mm'],accessible_notch_fraction=rm['accessible_notch_fraction'],removable_segment_fraction=rm['removable_segment_fraction'],physical_validation=False,hard_gate_failures=ff) for rr,ee,rm,tt,gg,ff in detailed])
    from .core_support_review import render,report
    render(out,r,e,removal,travel,ranking);report(out,s,r,e,removal)
    (out/'selected/README.md').write_text(f"# 当前候选\n\n状态：**{status}**。\n\n候选：{r['candidate_id']}。\n\n"+('全部软件门槛通过，仍需实物验证。' if not failures else '这是一套诊断文件，尚未通过生产门槛，不建议直接进行完整实物打印。')+'\n\n详见上一级 vascular_core_support_report.md 和 QC/。3MF 中血管及自定义支撑保留为独立 part，Auto Support 已关闭。\n')
    save(out/'reports/code_manifest.json',{str(p.relative_to(ROOT)):digest(p) for p in [ROOT/'s1-11_vascular_core_support.py',*sorted((ROOT/'vascular_processing').glob('core_support_*.py'))]})
    print('s1-11 RESULT',status,failures,flush=True);return s
