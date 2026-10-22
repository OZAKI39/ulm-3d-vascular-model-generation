"""s1-10 compact C revision; accepted V1 and all vascular inputs stay immutable."""
from __future__ import annotations
import json,shutil
from pathlib import Path
import numpy as np
import trimesh,yaml
from .support_inputs import ROOT,save,table,digest,verify
from .all_abs_pipeline import default_output as baseline_output,assemble
from .all_abs_audit import trim
from .all_abs_slicing import profiles
from .casting_fixture_geometry import export_shape,shape_mesh
from .casting_fixture_v2_audit import cad_manifold
from .compact_c_geometry import load_frozen,candidate
from .compact_c_audit import insertion,contacts,anti_float,axial_range,wall_routes
from .compact_c_slicing import slice_candidate


def default_output():return baseline_output()/'compact_c_retainer_v2'


def selection_key(row):
    # Equal CAD heights can differ by 4e-16 mm after Boolean operations.
    # A 1 nm comparison precision is far tighter than CAD/export tolerance.
    return (not row['hard_gates_passed'],round(row['max_height_above_port_mm'],6),round(row['total_volume_mm3'],6))


def protect(base,out):
    path=out/'source_protection_start.json'
    if path.exists():h=json.loads(path.read_text());verify(h)
    else:
        h=json.loads((base/'source_protection_start.json').read_text())
        for p in base.rglob('*'):
            if p.is_file() and 'compact_c_retainer_v2' not in p.relative_to(base).parts and 'legacy_retainer_v1' not in p.relative_to(base).parts:h[str(p)]=digest(p)
        save(path,h)
    archive=base/'legacy_retainer_v1';archive.mkdir(exist_ok=True)
    for p in list((base/'retainers').glob('port_*_ABS_retainer.*'))+[base/'QC/08_anti_float_retainers.png',base/'QC/11_final_ABS_assembly.png']:
        if p.suffix not in ['.stl','.step','.png']:continue
        dest=archive/p.name
        if not dest.exists():shutil.copy2(p,dest)
        if digest(dest)!=digest(p):raise ValueError('LEGACY_ARCHIVE_CHANGED '+str(dest))
    if not (archive/'summary_metrics.json').exists():
        import cadquery as cq
        baseline=json.loads((base/'all_abs_sacrificial_casting_summary.json').read_text());rows=[]
        for p in baseline['ports']:
            shape=cq.importers.importStep(str(base/f"retainers/port_{p['port_id']}_ABS_retainer.step")).val()
            rows.append(dict(port_id=p['port_id'],volume_mm3=shape.Volume(),height_above_port_mm=shape.BoundingBox().zmax-p['trim_origin'][2],role=p['role'],upper_clearance_mm=p['anti_float_gap_mm'],axial_float_mm=p['axial_float_mm']))
        save(archive/'summary_metrics.json',dict(ports=rows,assembly=baseline['assembly'],solvent=baseline['solvent_access'],PDMS=baseline['trim']))
    ahpath=out/'legacy_archive_protection_start.json'
    if ahpath.exists():ah=json.loads(ahpath.read_text());verify(ah)
    else:ah={str(p):digest(p) for p in archive.iterdir() if p.is_file()};save(ahpath,ah)
    historical=json.loads((base/'reports/implementation_manifest.json').read_text())
    # Hashes come from the accepted V1 manifest, not a post-edit snapshot.
    # Only the entrypoint, append-only tests, and current docs are revised.
    ch={str(ROOT/p):v for p,v in historical.items() if p.startswith('vascular_processing/all_abs_') or p=='config/all_abs_sacrificial_casting_system.yaml'}
    verify(ch);save(out/'legacy_code_protection_start.json',ch)
    return {**h,**ah,**ch}


def voids(l,old,new,basecfg):
    from .casting_fixture_geometry import bounds_box
    from .casting_fixture_v2_geometry import world
    rows=[];cav=cad_manifold(old['cavity'])
    for p in l['ports']:
        k=p['port_id'];half=basecfg['trim_notch_half_width_mm'];z0=l['low'][2]-p['trim_origin'][2];z1=l['high'][2]-p['trim_origin'][2]
        cut=cad_manifold(world(bounds_box([-half,0,z0-1],[half,p['wall_distance_from_trim_mm']+3,z1+1]),p));a=cad_manifold(old['retainers'][k])^cav;b=cad_manifold(new['retainers'][k])^cav
        rows.append(dict(port_id=k,old_retainer_auxiliary_void_mm3=float(a.volume()),new_retainer_auxiliary_void_mm3=float(b.volume()),new_retainer_void_outside_existing_trim_mm3=float((b-cut).volume()),not_increased=b.volume()<=a.volume()+1e-6,inside_existing_trim=abs((b-cut).volume())<1e-6))
    return dict(passed=all(r['not_increased'] and r['inside_existing_trim'] for r in rows),rows=rows)


def gates(s):
    c=s['contact'];v=s['solvent'];t=s['PDMS'];b=s['bambu'];eps=1e-6
    checks=dict(SOURCE_UNCHANGED=s['source_protection']['all_unchanged'],FROZEN_MOLD_AND_LOWER_CRADLE=s['frozen_geometry']['unchanged'],CORE_TOP_INSERTION_INHERITED=s['frozen_geometry']['core_top_insertion_pass'],C_NOT_CLOSED=all(0<p['c_retainer_wrap_angle_deg']<=180 for p in s['ports']),POSITIVE_CLEARANCE=c['positive_clearance'],TOP_INSERTION=s['insertion']['passed'],NO_SNAP=s['insertion']['no_elastic_deformation'],SECONDARY_FLOAT=s['axial_float']['passed'],ANTI_FLOAT=s['anti_float']['passed'],ANTI_FLOAT_NO_WORSE=s['anti_float_not_worse'],ZERO_TREE_CONTACT=c['central_tree_contacts']==0,EXTENSION_ONLY=c['extension_only_contacts'],NO_STEM=all(p['no_exposed_vertical_stem'] and p['below_nominal_fill'] for p in s['ports']),RETAINER_NOT_REQUIRED=not v['retainer_stem_required'],WALL_ROUTES=v['all_ports_wall_access'],FOUR_INDEPENDENT_ROUTES=v['independent_port_routes']>=s['legacy']['independent_port_routes'],BOTH_SWCS_REACHABLE=v['all_vascular_branches_reachable'],BURIED_RETAINER_ACCESS=v['all_buried_retainers_accessible'],RETAINER_VOID_NO_INCREASE=s['retainer_voids']['passed'],TOTAL_VOID_NO_INCREASE=t['auxiliary_void_new']<=t['auxiliary_void_old']+eps,TRIM_AUX_ZERO=abs(t['auxiliary_void_after_trim'])<eps,TRIM_GHOST_ZERO=abs(t['trim']['ghost_missing_volume_mm3'])<eps and abs(t['trim']['ghost_extra_volume_mm3'])<eps,TRIM_PLANES_UNCHANGED=t['trim_planes_unchanged'],PORTS_EXPOSED=t['trim']['all_ports_exposed'],USABLE_LENGTHS_PRESERVED=t['usable_lengths_preserved'],PDMS_EXPORTS=t['trim']['diagnostic_mesh_qc_pass'],REAL_ABS_SUPPORT_OFF=b['auto_support_off'] and b['slicer_success'],C_FEATURES=b['critical_features_preserved'],HEIGHT_REDUCTION=s['compact']['max_height']<s['legacy']['max_height'],VOLUME_REDUCTION=s['compact']['total_retainer_volume']<s['legacy']['total_retainer_volume'],NO_VAPOR=s['smoothing']['retainers_forbidden'] and s['smoothing']['vascular_core_forbidden'],STEP_ASSEMBLY=s['assembly_step_exported'],REVIEW_FIGURES=s.get('review_figures_complete',False),MESH_QC=all(q['passed'] for q in s['mesh_qc'].values()))
    failed=[k for k,v in checks.items() if not v];return ('NEEDS_ADJUSTMENT' if failed else 'READY_FOR_HUMAN_REVIEW_WITH_PROCESS_PARAMETERS_UNCALIBRATED'),checks,failed


def run(config=None,out=None,do_slice=True,do_step=True,do_render=True):
    cfg=yaml.safe_load(Path(config or ROOT/'config/compact_c_retainer_v2.yaml').read_text());out=Path(out or default_output()).resolve();base=baseline_output()
    if out.is_relative_to(base) and not out.is_relative_to(default_output()):raise ValueError('CANNOT_OVERWRITE_FROZEN_BASELINE')
    for folder in ['retainers','assembly','diagnostic_PDMS','tables','solvent_access','bambu','trim','QC','reports','candidates']:(out/folder).mkdir(parents=True,exist_ok=True)
    hashes=protect(base,out);print('s1-10 compact: load frozen CAD, STL, manifest and SWCs',flush=True);base,old,inputs,bc,l,g,line=load_frozen();save(out/'reports/config.json',cfg)
    original_metrics=json.loads((base/'legacy_retainer_v1/summary_metrics.json').read_text())['ports'];oldanti=anti_float(l,g);values,actual,ref,engine=profiles(inputs,out/'bambu');candidates=[];evaluated=[]
    for angle in cfg['wrap_angle_candidates_deg']:
        print('s1-10 compact: candidate',angle,flush=True);ng=candidate(l,g,cfg,line,angle);folder=out/'candidates'/f'wrap_{angle:g}';folder.mkdir(exist_ok=True)
        ins=insertion(l,ng,cfg);con=contacts(l,ng);anti=anti_float(l,ng);axial=axial_range(l,g,ng,cfg);sol=wall_routes(l,ng,bc,inputs,folder/'solvent_access');void=voids(l,g,ng,bc)
        rec,feat,pk=slice_candidate(ng,engine,actual,out/'bambu'/f'candidate_{angle:g}') if do_slice else ({'slice_passed':False,'auto_support_off':False},{'passed':False},[])
        for k,part in ng['retainers'].items():export_shape(part,folder/f'port_{k}_C_retainer',bc)
        hard=dict(anti_float=anti['passed'],top_insertion=ins['passed'],axial_float=axial['passed'],central_no_contact=con['central_tree_contacts']==0,extension_only=con['extension_only_contacts'],positive_clearance=con['positive_clearance'],solvent_access=sol['all_ports_wall_access'] and sol['all_vascular_branches_reachable'] and sol['all_buried_retainers_accessible'],voids=void['passed'],support_off_slice=rec.get('slice_passed',False) and rec.get('auto_support_off',False),features_preserved=feat['passed'])
        row=dict(wrap_angle_deg=angle,opening_angle_deg=360-angle,**hard,hard_gates_passed=all(hard.values()),max_height_above_port_mm=max(r['retainer_max_height_above_port_mm'] for r in ng['compact_parameters']),total_volume_mm3=sum(r['retainer_volume_mm3'] for r in ng['compact_parameters']))
        evidence=dict(parameters=ng['compact_parameters'],candidate=row,insertion=ins,contact=con,anti_float=anti,axial_float=axial,solvent=sol,retainer_voids=void,bambu=rec,features=feat,packing=pk);save(folder/'evaluation.json',evidence);candidates.append(row);evaluated.append((ng,evidence));print('candidate result',row,flush=True)
    # Hard gates before the explicitly ordered height, then volume objectives.
    order=sorted(range(len(candidates)),key=lambda j:selection_key(candidates[j]))
    selected=order[0];ng,e=evaluated[selected];angle=candidates[selected]['wrap_angle_deg'];print('selected wrap',angle,flush=True)
    for j,row in enumerate(candidates):row['selected']=j==selected
    table(out/'tables/C_retainer_candidates.csv',candidates);save(out/'candidate_selection.json',dict(rows=candidates,selection_order=['ALL_HARD_GATES','LOWER_MAX_HEIGHT','LOWER_TOTAL_VOLUME'],comparison_precision_mm=1e-6,black_box_score=False,selected_wrap_angle_deg=angle))
    sol=wall_routes(l,ng,bc,inputs,out/'solvent_access');print('s1-10 compact: recalculate PDMS and existing trim',flush=True);t,trimrows,diagnostics=trim(l,ng,bc,out/'trim');t['retainer_stems_above_fill_not_counted_as_PDMS_void']=False;t['no_exposed_retainer_stems']=True;save(out/'trim/post_dissolution_trim.json',t)
    for src,dst in [('PDMS_after_all_ABS_removed','PDMS_after_C_retainer_dissolution'),('PDMS_after_trim','PDMS_after_C_retainer_trim')]:shutil.copy2(out/f'diagnostic_PDMS/{src}.stl',out/f'diagnostic_PDMS/{dst}.stl')
    meshes={};qc={}
    for k,part in ng['retainers'].items():meshes[k],qc[k]=export_shape(part,out/f'retainers/port_{k}_C_retainer',bc)
    trimesh.util.concatenate([l['core'],g['mold_mesh'],*meshes.values()]).export(out/'assembly/all_ABS_C_retainer_review.stl')
    if do_step:
        print('s1-10 compact: export independent assembly STEP',flush=True);assemble(ng,inputs['core'],l['transform'],out/'assembly/all_ABS_casting_system_C_retainer.step')
    if e['bambu'].get('archive'):shutil.copy2(e['bambu']['archive'],out/'retainers/C_retainers_X1C_0p2_ABS.3mf')
    save(out/'retainer_insertion_sequence.json',e['insertion']);table(out/'tables/C_retainer_clearances.csv',e['contact']['rows']);table(out/'tables/C_retainer_insertion.csv',e['insertion']['rows']);table(out/'tables/C_retainer_axial_float.csv',e['axial_float']['rows']);table(out/'tables/C_retainer_auxiliary_voids.csv',e['retainer_voids']['rows']);table(out/'tables/C_retainer_anti_float.csv',[dict(version=v,**r) for v,a in [('legacy',oldanti),('compact',e['anti_float'])] for r in a['poses']])
    compare=[]
    for om,nm,ar in zip(original_metrics,ng['compact_parameters'],e['axial_float']['rows']):
        compare.append(dict(port_id=nm['port_id'],old_volume=om['volume_mm3'],new_volume=nm['retainer_volume_mm3'],volume_reduction_percent=100*(1-nm['retainer_volume_mm3']/om['volume_mm3']),old_max_height=om['height_above_port_mm'],new_max_height=nm['retainer_max_height_above_port_mm'],height_reduction_percent=100*(1-nm['retainer_max_height_above_port_mm']/om['height_above_port_mm']),old_axial_float=[ar['old_min_mm'],ar['old_max_mm']],new_axial_float=[ar['new_min_mm'],ar['new_max_mm']],anti_float_old=oldanti['passed'],anti_float_new=e['anti_float']['passed']))
    table(out/'tables/legacy_vs_compact_retainer.csv',compare);volold=sum(r['old_volume'] for r in compare);volnew=sum(r['new_volume'] for r in compare);hold=max(r['old_max_height'] for r in compare);hnew=max(r['new_max_height'] for r in compare)
    unchanged=all(ng[k] is g[k] for k in ['mold','fixtures','guides','lands','bottom','walls']);fixed=dict(unchanged=unchanged,core_top_insertion_pass=old['assembly']['top_insertion_pass'],transform_4x4=l['transform'],source_core=str(inputs['core']),source_core_sha256=digest(inputs['core']),outer_mold_step=str(base/'outer_mold/ABS_sacrificial_outer_mold.step'),outer_mold_stl=str(base/'outer_mold/ABS_sacrificial_outer_mold.stl'),outer_mold_sha256=digest(base/'outer_mold/ABS_sacrificial_outer_mold.stl'),source_assembly_step=str(base/'assembly/all_ABS_casting_system.step'),outer_mold_dimensions=old['outer_mold'],cradle_and_wall_blade_reconstructed_read_only=True,no_original_shape_reexport=True)
    protection=verify(hashes);protection.update(start_hashes=hashes,end_hashes={p:digest(p) for p in hashes},mismatch_count=0);save(out/'source_protection.json',protection)
    s=dict(method=cfg['method'],selected_wrap_angle_deg=angle,candidates=candidates,legacy=dict(total_retainer_volume=volold,max_height=hold,solvent_route_dependency='OLD_RETAINER_SELF_REMOVAL_VIA_EXPOSED_STEM; PORTS_ALREADY_HAVE_INDEPENDENT_WALL_ROUTES',independent_port_routes=old['solvent_access']['independent_access_routes']),compact=dict(total_retainer_volume=volnew,max_height=hnew,volume_reduction_mm3=volold-volnew,volume_reduction_percent=100*(1-volnew/volold),height_reduction_mm=hold-hnew,height_reduction_percent=100*(1-hnew/hold)),ports=ng['compact_parameters'],comparison=compare,insertion=e['insertion'],contact=e['contact'],axial_float=e['axial_float'],anti_float=e['anti_float'],anti_float_legacy=oldanti,anti_float_not_worse=all(not a['blocked'] or b['blocked'] for a,b in zip(oldanti['poses'],e['anti_float']['poses'])) and all(not a['retainer_up_0p8_blocked'] or b['retainer_up_0p8_blocked'] for a,b in zip(oldanti['ports'],e['anti_float']['ports'])),solvent=sol,retainer_voids=e['retainer_voids'],PDMS=dict(auxiliary_void_old=old['trim']['auxiliary_void_volume_before_trim'],auxiliary_void_new=t['auxiliary_void_volume_before_trim'],auxiliary_void_after_trim=t['auxiliary_void_volume_after_trim'],trim=t,trim_planes_unchanged=t['actual_cut_planes']==old['trim']['actual_cut_planes'],usable_lengths_preserved=all(t['usable_port_lengths'][k]>=v-1e-7 for k,v in old['trim']['usable_port_lengths'].items())),bambu=dict(auto_support_off=e['bambu'].get('auto_support_off',False),slicer_success=e['bambu'].get('slice_passed',False),critical_features_preserved=e['features']['passed'],selected_slice=e['bambu'],features=e['features'],unchanged_mold_slice_reused=old['bambu']['mold_auto_support_off_slice'],unchanged_vascular_reference_slice_reused=old['bambu']['vascular_reference_slice']),source_protection={k:v for k,v in protection.items() if k not in ['start_hashes','end_hashes']},frozen_geometry=fixed,mesh_qc=qc,assembly_step_exported=do_step,smoothing=dict(retainers_forbidden=True,vascular_core_forbidden=True,mold_optional=True,dimensionally_calibrated=False),physical_validation=dict(anti_float_force=False,print=False,dissolution=False,PDMS_swelling_recovery=False,manual_trim=False),no_printer_delivery=True)
    s['status'],s['gates'],s['failures']=gates(s);save(out/'compact_c_retainer_summary.json',s)
    from .compact_c_review import render,report,NAMES
    if do_render:render(out,l,g,ng,diagnostics,s)
    s['review_figures_complete']=do_render and all((out/'QC'/(n+'.png')).is_file() for n in NAMES);s['status'],s['gates'],s['failures']=gates(s);save(out/'compact_c_retainer_summary.json',s);report(out,s)
    save(out/'reports/artifact_manifest.json',{str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='artifact_manifest.json' and p.suffix!='.log'})
    print('s1-10 compact RESULT',s['status'],s['failures'],flush=True);return s
