"""s1-10 dimensional revision: 2 mm walls/base and moderate top/bottom cover.

Vascular geometry and C-cap geometry stay rigid. The box and vertical fixture
stems adapt to the new floor/rim. Previous accepted output is immutable history.
"""
from __future__ import annotations
import copy,json,shutil,zipfile
from pathlib import Path
import numpy as np
import trimesh,yaml
from scipy.spatial import cKDTree
from .support_inputs import ROOT,save,table,digest,verify,numeric
from .all_abs_pipeline import default_output as base_output,assemble
from .compact_c_pipeline import default_output as compact_output
from .compact_c_geometry import load_frozen,candidate
from .all_abs_geometry import build
from .all_abs_audit import assembly,trim
from .compact_c_audit import insertion,contacts,anti_float,axial_range,wall_routes
from .all_abs_slicing import profiles,slice_abs
from .compact_c_slicing import slice_candidate
from .casting_fixture_geometry import export_shape,shape_mesh
from .casting_fixture_v2_geometry import port_transform
from .bambu_casting_mold import parse_paths,sample_segments


def default_output():return compact_output()/'mold_2mm_clearance_5mm'


def snapshot(out):
    path=out/'source_protection_start.json'
    if path.exists():h=json.loads(path.read_text());verify(h);return h
    previous=json.loads((compact_output()/'source_protection.json').read_text());h=dict(previous['start_hashes']);verify(h)
    for p in compact_output().rglob('*'):
        if p.is_file() and not p.is_relative_to(default_output()) and not p.is_relative_to(out):h[str(p)]=digest(p)
    save(path,h);return h


def revised_layout(old,cfg,reference):
    for k in ['wall_thickness_mm','bottom_thickness_mm','bottom_clearance_mm','top_clearance_mm','minimum_fixture_to_rim_mm']:
        if not np.isfinite(cfg[k]) or cfg[k]<=0:raise ValueError('NONPOSITIVE_MOLD_DIMENSION_'+k)
    l=copy.deepcopy(old);delta=cfg['bottom_thickness_mm']+cfg['bottom_clearance_mm']-old['core'].bounds[0,2];shift=np.array([0.,0.,delta]);m=np.eye(4);m[:3,3]=shift;l['transform']=m@old['transform']
    for mesh in [l['core'],l['central'],*l['source_ports'].values()]:mesh.apply_translation(shift)
    l['bifurcations']+=shift
    for p in l['ports']:
        for key in ['straight_start','end','attachment','trim_origin','approved_contact_start','approved_contact_end']:p[key]=np.asarray(p[key])+shift
    l['low'][2]=cfg['bottom_thickness_mm'];requested_top=l['core'].bounds[1,2]+cfg['top_clearance_mm']
    # Guide height comes from the accepted production kernel. Cap height comes
    # from its accepted summary; neither geometry is shortened to fit the box.
    fixed_tops=[]
    for p in l['ports']:
        k=p['port_id'];d=next(r for r in reference['caps'] if r['port_id']==k)
        fixed_tops.extend([p['trim_origin'][2]+d['retainer_max_height_above_port_mm'],p['trim_origin'][2]+reference['guides'][k]])
    l['high'][2]=max(requested_top,max(fixed_tops)+cfg['minimum_fixture_to_rim_mm'])
    l['vertical_revision']=dict(translation_mm=shift,requested_top_clearance_mm=cfg['top_clearance_mm'],actual_top_clearance_mm=l['high'][2]-l['core'].bounds[1,2],actual_bottom_clearance_mm=l['core'].bounds[0,2]-l['low'][2],minimum_fixed_fixture_to_rim_mm=l['high'][2]-max(fixed_tops),top_definition='OPEN_RIM_AND_NOMINAL_FILL_PLANE',xy_placement_and_rotation_unchanged=True)
    # These records describe the inherited XY decision, not a new height search.
    for key in ['placements','selected_placement']:l.pop(key,None)
    return l,shift


def mold_toolpaths(l,g,rec,out):
    if not rec.get('slice_passed'):return dict(passed=False,reason='MOLD_SLICE_MISSING')
    paths=parse_paths(rec['archive']);seg=paths['segments'][paths['kind']==0];tree=cKDTree(sample_segments(seg,.08));M=np.asarray(json.loads((out/'print_transform.json').read_text())['transform_4x4']);settings=rec['actual_profile'];layer=numeric(settings['layer_height']);first=numeric(settings['initial_layer_print_height']);t=g['blade_thickness_mm'];rows=[]
    def check(k,name,points):
        qs=np.asarray(points)@M[:3,:3].T+M[:3,3];dist=tree.query(qs)[0];tol=g['line_width_mm']*.8
        rows.append(dict(port_id=k,feature=name,query_count=len(qs),missing_query_count=int((dist>tol).sum()),max_distance_mm=float(dist.max()),tolerance_mm=tol,passed=bool(len(qs) and (dist<=tol).all())))
    for p in l['ports']:
        k=p['port_id'];R=p['radius_mm']+p['cradle_clearance_mm'];inner=R+.55;w=p['wall_distance_from_trim_mm'];top=g['local_parts'][k]['guide_top_local_z'];levels=np.arange(first,l['high'][2],layer);zc=levels-layer/2-p['trim_origin'][2];zf=l['low'][2]-p['trim_origin'][2];P=port_transform(p)
        def vertical(name,xy,zmax,zmin=zf+.03):
            qs=[];mask=(zc>zmin)&(zc<zmax-.04)
            for x,y in xy:qs.extend(np.c_[np.full(mask.sum(),x),np.full(mask.sum(),y),levels[mask]-p['trim_origin'][2]])
            if qs:check(k,name,np.asarray(qs)@P[:3,:3].T+P[:3,3])
        ys=np.linspace(min(g['local_parts'][k]['cradle_ys'])-.2,w-.35,20)
        for sign in [-1,1]:vertical('wall_fin_'+str(sign),[(sign*(inner+t/2),y) for y in ys],top)
        vertical('solvent_web',[(0,y) for y in np.linspace(2.3,w-.4,15)],-R-.08)
        access=next(r for r in g['access_records'] if r['port_id']==k);cz=min(np.asarray(access['contact_witnesses'])[:,2])-p['trim_origin'][2]
        vertical('conformal_access_land',[(0,y) for y in [2.,2.2,2.4]],cz,zmin=-R-.12+.03)
        vertical('transferred_trim_guide',[(x,t/2) for x in np.linspace(-3.8,-2.8,5)],zf+.2,zmin=zf+.02)
    result=dict(passed=bool(rows) and all(r['passed'] for r in rows),rows=rows,scope='Actual positive extrusion on every occupied thin-feature layer; physical contact uncalibrated')
    save(out/'mold_feature_preservation.json',result);table(out/'mold_feature_preservation.csv',rows);return result


def checks(s):
    a=s['assembly'];c=s['contact'];t=s['trim'];v=s['solvent'];b=s['bambu'];d=s['dimensions'];eps=1e-6
    gate=dict(SOURCE_HISTORY_IMMUTABLE=s['source_protection']['all_unchanged'],RIGID_CORE_ONLY=s['core_rigid_only'],C_CAP_SHAPE_UNCHANGED=s['cap_shape_unchanged'],WALLS_AND_BOTTOM_2MM=abs(d['wall_thickness_mm']-2)<eps and abs(d['bottom_thickness_mm']-2)<eps,REQUESTED_BOTTOM_GAP=abs(d['bottom_clearance_mm']-s['config']['bottom_clearance_mm'])<eps,SUFFICIENT_TOP_GAP=d['top_clearance_mm']>=s['config']['top_clearance_mm']-eps,FIXTURE_BELOW_RIM=d['fixture_to_rim_mm']>=s['config']['minimum_fixture_to_rim_mm']-eps,XY_UNCHANGED=s['xy_unchanged'],CLOSED_MOLD=a['walls_closed'] and a['bottom_closed'] and a['top_open'],CORE_INSERTION=a['top_insertion_pass'],C_INSERTION=s['insertion']['passed'],NO_SNAP=s['insertion']['no_elastic_deformation'],POSITIVE_CAP_GAP=c['positive_clearance'],NO_CENTRAL_CONTACT=a['central_tree_contacts']==0 and c['central_tree_contacts']==0,EXTENSION_ONLY=a['all_contacts_on_approved_extensions'] and c['extension_only_contacts'],AXIAL_FLOAT=s['axial_float']['passed'] and s['previous_C_fit_preserved'],ANTI_FLOAT=s['anti_float']['passed'],SOLVENT_WALL_ROUTES=v['all_ports_wall_access'] and v['independent_port_routes']==4,SOLVENT_C_CAP_ACCESS=v['all_buried_retainers_accessible'],BOTH_SWCS_REACHABLE=v['all_vascular_branches_reachable'],AUX_VOID_ZERO=abs(t['auxiliary_void_volume_after_trim'])<eps,PDMS_GHOST_ZERO=abs(t['ghost_missing_volume_mm3'])<eps and abs(t['ghost_extra_volume_mm3'])<eps,CENTRAL_TRIM_UNTOUCHED=abs(t['central_vascular_removed_mm3'])<eps,PORTS_EXPOSED=t['all_ports_exposed'],TRIM_PLANE_GEOMETRY_UNCHANGED=s['trim_planes_same_geometric_planes'],PORT_LENGTHS_UNCHANGED=s['port_lengths_unchanged'],PDMS_EXPORTS=t['diagnostic_mesh_qc_pass'],MESH_EXPORTS=all(x['passed'] for x in s['mesh_qc'].values()),MOLD_REAL_SUPPORT_OFF=b['mold'].get('slice_passed',False) and b['mold'].get('auto_support_off',False),MOLD_THIN_FEATURES=b['mold_features']['passed'],C_REAL_SUPPORT_OFF=b['retainers'].get('slice_passed',False) and b['retainers'].get('auto_support_off',False),C_FEATURES=b['retainer_features']['passed'],ASSEMBLY_STEP=s['assembly_step_exported'],QC_IMAGES=s.get('review_figures_complete',False))
    fail=[k for k,v in gate.items() if not v];return ('NEEDS_ADJUSTMENT' if fail else 'READY_FOR_HUMAN_REVIEW_WITH_PROCESS_PARAMETERS_UNCALIBRATED'),gate,fail


def run(config=None,out=None,do_slice=True,do_step=True,do_render=True):
    cfg=yaml.safe_load(Path(config or ROOT/'config/compact_abs_mold_2mm.yaml').read_text());out=Path(out or default_output()).resolve()
    if out.is_relative_to(base_output()) and not out.is_relative_to(default_output()):raise ValueError('CANNOT_OVERWRITE_ACCEPTED_HISTORY')
    for f in ['reports','outer_mold','retainers','assembly','vascular_core','tables','solvent_access','trim','diagnostic_PDMS','bambu','QC']:(out/f).mkdir(parents=True,exist_ok=True)
    h=snapshot(out);print('s1-10: load accepted C caps and frozen vascular source',flush=True);base,bs,inputs,bc,old,oldg,line=load_frozen();cs=json.loads((compact_output()/'compact_c_retainer_summary.json').read_text());capcfg=yaml.safe_load((ROOT/cfg['retainer_config']).read_text());reference=dict(caps=cs['ports'],guides={k:d['guide_top_local_z'] for k,d in oldg['local_parts'].items()});l,shift=revised_layout(old,cfg,reference)
    nc=dict(bc);nc.update({k:cfg[k] for k in ['wall_thickness_mm','bottom_thickness_mm','bottom_clearance_mm','top_clearance_mm']});save(out/'reports/config.json',cfg);save(out/'reports/mold_geometry_config.json',nc)
    print('s1-10: new surface clearances',l['vertical_revision'],flush=True);raw=build(l,nc,line);g=candidate(l,raw,capcfg,line,cs['selected_wrap_angle_deg']);g['mold_mesh'],moldqc=export_shape(g['mold'],out/'outer_mold/ABS_sacrificial_outer_mold_2mm',nc)
    # Reuse the exact accepted C-cap CAD; only its common assembly translation
    # changes. The print-local shape and resulting printing orientation remain.
    import cadquery as cq
    capshape=True;meshqc={'mold':moldqc};capmeshes={}
    for p in l['ports']:
        k=p['port_id'];accepted=cq.importers.importStep(str(compact_output()/f'retainers/port_{k}_C_retainer.step')).val().translate(tuple(shift));difference=abs(accepted.cut(g['retainers'][k]).Volume())+abs(g['retainers'][k].cut(accepted).Volume());capshape&=difference<1e-6;g['retainers'][k]=accepted;capmeshes[k],meshqc[k]=export_shape(accepted,out/f'retainers/port_{k}_C_retainer',nc)
    l['core'].export(out/'vascular_core/assembly_placed_ABS_core.stl');shutil.copy2(inputs['core'],out/'vascular_core/BG001_RMCA_frozen_ABS_core.stl')
    save(out/'assembly/layout.json',{k:v for k,v in l.items() if k not in ['core','central','source_ports']});save(out/'assembly/casting_core_transform.json',dict(transform_4x4=l['transform'],inverse_4x4=np.linalg.inv(l['transform']),previous_transform_4x4=old['transform'],common_translation_mm=shift,scale_factor=1.))
    a=assembly(l,g,nc,out/'tables');ins=insertion(l,g,capcfg);a.update(retainer_insertions=ins['rows'],retainer_insertion_pass=ins['passed'],retainer_insertion_algorithm='C_CAP_CONTINUOUS_WITH_PLANAR_WALL_SEATING');save(out/'tables/assembly_audit.json',a)
    con=contacts(l,g);anti=anti_float(l,g);axial=axial_range(l,raw,g,capcfg);fit=all(abs(n['new_min_mm']-o['new_min_mm'])<.002 and abs(n['new_max_mm']-o['new_max_mm'])<.002 for n,o in zip(axial['rows'],cs['axial_float']['rows']));sol=wall_routes(l,g,nc,inputs,out/'solvent_access')
    print('s1-10: recalculate PDMS domain and trimming for reduced height',flush=True);t,tr,d=trim(l,g,nc,out/'trim');t['retainer_stems_above_fill_not_counted_as_PDMS_void']=False;t['no_exposed_retainer_stems']=True;save(out/'trim/post_dissolution_trim.json',t)
    planes=t['actual_cut_planes'];oldplanes=cs['PDMS']['trim']['actual_cut_planes'];sameplanes=len(planes)==len(oldplanes) and all(a['port']==b['port'] and a['kind']==b['kind'] and np.allclose(a['normal'],b['normal']) and abs((np.asarray(a['origin'])-b['origin'])@a['normal'])<1e-6 for a,b in zip(planes,oldplanes))
    print('s1-10: actual X1 Carbon 0.2 mm ABS mold and cap slices',flush=True);values,actual,ref,engine=profiles(inputs,out/'bambu');mr=slice_abs(g['mold_mesh'],engine,actual,out/'bambu/outer_mold',True) if do_slice else dict(slice_passed=False);mf=mold_toolpaths(l,g,mr,out/'bambu/outer_mold')
    rr,rf,packing=slice_candidate(g,engine,actual,out/'bambu/retainers') if do_slice else (dict(slice_passed=False),dict(passed=False),[])
    for rec,dest in [(mr,out/'outer_mold/ABS_sacrificial_outer_mold_2mm_X1C_ABS.3mf'),(rr,out/'retainers/C_retainers_X1C_0p2_ABS.3mf')]:
        if rec.get('archive'):shutil.copy2(rec['archive'],dest)
    trimesh.util.concatenate([l['core'],g['mold_mesh'],*capmeshes.values()]).export(out/'assembly/all_ABS_C_retainer_mold_2mm_review.stl')
    if do_step:assemble(g,inputs['core'],l['transform'],out/'assembly/all_ABS_casting_system_C_retainer_2mm.step')
    proof=verify(h);proof.update(start_hashes=h,end_hashes={p:digest(p) for p in h},mismatch_count=0);save(out/'source_protection.json',proof)
    dimensions=dict(wall_thickness_mm=nc['wall_thickness_mm'],bottom_thickness_mm=nc['bottom_thickness_mm'],bottom_clearance_mm=l['core'].bounds[0,2]-l['low'][2],top_clearance_mm=l['high'][2]-l['core'].bounds[1,2],inner_dimensions_mm=l['high']-l['low'],outer_dimensions_mm=g['outer_bounds'][1]-g['outer_bounds'][0],fixture_to_rim_mm=l['vertical_revision']['minimum_fixed_fixture_to_rim_mm'],old_bottom_clearance_mm=old['core'].bounds[0,2]-old['low'][2],old_top_clearance_mm=old['high'][2]-old['core'].bounds[1,2],old_outer_height_mm=oldg['outer_bounds'][1,2]-oldg['outer_bounds'][0,2],old_bottom_thickness_mm=bc['bottom_thickness_mm'])
    s=dict(method=cfg['method'],config=cfg,dimensions=dimensions,common_translation_mm=shift,core_rigid_only=bool(np.array_equal(old['core'].faces,l['core'].faces) and np.allclose(l['core'].vertices,old['core'].vertices+shift,atol=1e-12,rtol=0)),cap_shape_unchanged=bool(capshape),xy_unchanged=bool(np.array_equal(old['low'][:2],l['low'][:2]) and np.array_equal(old['high'][:2],l['high'][:2])),source_protection={k:v for k,v in proof.items() if k not in ['start_hashes','end_hashes']},assembly=a,insertion=ins,contact=con,anti_float=anti,axial_float=axial,previous_C_fit_preserved=fit,solvent=sol,trim=t,trim_planes_same_geometric_planes=bool(sameplanes),port_lengths_unchanged=t['usable_port_lengths']==cs['PDMS']['trim']['usable_port_lengths'],bambu=dict(mold=mr,mold_features=mf,retainers=rr,retainer_features=rf,vascular_reference=bs['bambu']['vascular_reference_slice']),mesh_qc=meshqc,assembly_step_exported=do_step,review_figures_complete=False,physical_process_calibrated=False,no_printer_delivery=True,source_core_sha256=digest(inputs['core']))
    table(out/'tables/dimension_changes.csv',[dict(parameter=k,value=v) for k,v in dimensions.items()]);table(out/'tables/C_retainer_clearances.csv',con['rows']);table(out/'tables/C_retainer_axial_float.csv',axial['rows']);save(out/'retainer_insertion_sequence.json',ins)
    s['status'],s['gates'],s['failures']=checks(s);save(out/'mold_revision_summary.json',s)
    from .compact_mold_review import render,report,NAMES
    if do_render:render(out,old,l,oldg,g,d,s)
    s['review_figures_complete']=do_render and all((out/'QC'/(n+'.png')).is_file() for n in NAMES);s['status'],s['gates'],s['failures']=checks(s);save(out/'mold_revision_summary.json',s);report(out,s)
    save(out/'reports/artifact_manifest.json',{str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='artifact_manifest.json' and p.suffix!='.log'})
    print('s1-10 dimension revision RESULT',s['status'],s['failures'],flush=True);return s
