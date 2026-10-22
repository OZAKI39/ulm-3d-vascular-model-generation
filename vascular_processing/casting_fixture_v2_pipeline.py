"""Current s1-9 production: closed walls, internal wall fins and trim guides."""
from __future__ import annotations
import json,shutil
from pathlib import Path
import numpy as np
import trimesh
import cadquery as cq
import yaml
from scipy.spatial.transform import Rotation
from .support_inputs import ROOT,discover,digest,save,table,numeric,verify
from .casting_fixture_pipeline import default_output as legacy_output,dish_check
from .casting_fixture_geometry import shape_mesh,export_shape,mesh_qc
from .casting_fixture_profiles import profiles,slice_part
from .casting_fixture_v2_geometry import layout,build,port_transform,clearance_study
from .casting_fixture_v2_audit import audit,trim
from .casting_fixture_v2_slicing import validate_blades


def default_output():return legacy_output()/'internal_wall_supported_fixture_v2'


def gates(s):
    a=s['assembly'];t=s['trim'];b=s['bambu'];q=s['mesh_qc']
    v=dict(FULL_CORE_INSIDE=t['core_outside_cavity_mm3']<1e-7,FOUR_WALLS_CLOSED=a['four_walls_closed'],BOTTOM_OPEN=a['bottom_open'],TOP_OPEN=a['top_open'],CENTRAL_NO_CONTACT=a['central_contact_count']==0,EXTENSION_CONTACTS_ONLY=a['all_design_contacts_on_extensions'],CORE_TOP_INSERTION=a['top_insertion']['passed'],RETAINER_INSERTION=a['retainer_insertion_pass'],ANTI_FLOAT=a['anti_float']['geometric_Z_constraint'],TRIM_IMPRINT_REMOVED=t['fixture_imprint_remaining_after_trim_mm3']<1e-7 and t['ghost_channel_missing_volume_mm3']<1e-7,CENTRAL_TRIM_UNTOUCHED=t['central_vascular_removed_mm3']<1e-7,PORTS_EXPOSED=t['all_ports_exposed'],TRIM_NONPORT_CLEARANCE=t['minimum_nonport_vascular_trim_clearance_mm']>0,FRAME_REAL_SLICE=b.get('frame_slice',{}).get('support_free_slice_passed',False),BLADES_PRESERVED=b.get('blade_preservation',{}).get('passed',False),RETAINERS_REAL_SLICE=b.get('retainer_slice',{}).get('support_free_slice_passed',False),MESHES_VALID=q['frame']['passed'] and q['vascular']['passed'] and all(x['passed'] for x in q['retainers'].values()),SOURCE_UNCHANGED=s['source_protection']['all_unchanged'],STEP_ASSEMBLY=s['assembly_step_exported'])
    v['NO_UNRELATED_PORT_TRIM']=t.get('maximum_unrelated_port_cut_volume_mm3',1.)<1e-7
    v['NO_OBVIOUS_RELEASE_INTERLOCK']=t.get('no_obvious_vertical_interlock',False)
    v['PDMS_STL_ROUNDTRIP']=t.get('diagnostic_mesh_qc_passed',False)
    v['DISH_FIT_IF_MEASURED']=not s['dish']['actual_dimensions_known'] or s['dish']['fit_verified']
    failed=[k for k,ok in v.items() if not ok];return ('NEEDS_ADJUSTMENT' if failed else 'READY_FOR_HUMAN_REVIEW_WITH_PHYSICAL_SEAL_PENDING'),v,failed


def assemble(g,source,m,path):
    from OCP.StlAPI import StlAPI_Reader
    from OCP.TopoDS import TopoDS_Shape
    raw=TopoDS_Shape()
    if not StlAPI_Reader().Read(raw,str(source)):raise ValueError('FROZEN_STL_STEP_READ_FAILED')
    core=cq.Shape.cast(raw);angle=np.degrees(np.arctan2(m[1,0],m[0,0]));core=core.rotate((0,0,0),(0,0,1),angle).translate(tuple(m[:3,3]))
    assembly=cq.Assembly(name='INTERNAL_WALL_FIXTURE_V2_SEPARATE_PARTS');assembly.add(core,name='FROZEN_ABS_CORE',color=cq.Color(.85,.24,.16));assembly.add(g['frame'],name='PLA_CLOSED_WALL_FRAME',color=cq.Color(.15,.5,.8))
    for k,v in g['retainers'].items():assembly.add(v,name='PLA_RETAINER_'+k,color=cq.Color(.2,.75,.4))
    assembly.export(str(path))


def run(config=None,out=None,do_slice=True,do_step=True,do_render=True):
    cfg=yaml.safe_load(Path(config or ROOT/'config/internal_wall_supported_fixture_v2.yaml').read_text());out=Path(out or default_output());out.mkdir(parents=True,exist_ok=True)
    inputs=discover(source=cfg.get('source'));save(out/'config/source_discovery.json',inputs);save(out/'config/fixture_config.json',cfg)
    source=out/'vascular_core/frozen_ABS_core_with_extensions.stl';source.parent.mkdir(exist_ok=True);shutil.copy2(inputs['core'],source)
    if digest(source)!=digest(inputs['core']):raise ValueError('FROZEN_COPY_HASH_MISMATCH')
    _,pla,engine=profiles(inputs,cfg,out/'config/profiles');line=numeric(pla['process']['line_width'],numeric(pla['machine']['nozzle_diameter'])) if pla else .22
    l=layout(inputs,cfg,line);save(out/'assembly/layout.json',{k:v for k,v in l.items() if k not in ['core','central','source_ports']});save(out/'assembly/casting_core_transform.json',dict(transform_4x4=l['transform'],inverse_4x4=np.linalg.inv(l['transform']),source_sha256=digest(source),scale_factor=1.,frame='DISH_Z0'))
    table(out/'tables/rigid_placement_comparison.csv',l['placements']);table(out/'tables/extension_fixture_zones.csv',l['ports'])
    clearances=[dict(port_id=p['port_id'],port_diameter_mm=p['diameter_mm'],clearance_mm=line*f,clearance_lines=f,selected=abs(line*f-p['cradle_clearance_mm'])<1e-8,physical_tolerance_calibrated=False,criterion='Positive slip clearance; master one line, secondary 1.5 lines to reduce preload') for p in l['ports'] for f in cfg['radial_clearance_line_multipliers']];table(out/'tables/cradle_clearance_candidates.csv',clearances)
    clearance_study(l,cfg,line,out/'config/clearance_study')
    candidates=[];chosen=None
    for lines in cfg['blade_line_candidates']:
        print('V2 building / actual slicing',lines,'line fins',flush=True);g=build(l,cfg,line,lines);m=shape_mesh(g['frame']);folder=out/'bambu'/f'frame_{lines}lines';folder.mkdir(parents=True,exist_ok=True)
        record=slice_part(m,engine,pla,folder) if do_slice and pla and engine else dict(support_free_slice_passed=False)
        preservation=validate_blades(record,l,g,folder) if do_slice and record.get('slicer_success') else dict(passed=False)
        candidates.append(dict(lines=lines,thickness_mm=line*lines,slice=record,blade_preservation=preservation))
        if chosen is None and record['support_free_slice_passed'] and preservation['passed']:chosen=(g,record,preservation)
    if chosen is None:g=build(l,cfg,line,2);frame_record=candidates[0]['slice'];preservation=candidates[0]['blade_preservation']
    else:g,frame_record,preservation=chosen
    save(out/'bambu/blade_line_comparison.json',candidates)
    meshes={'retainers':{}};q={'vascular':l['source_qc'],'retainers':{}}
    meshes['frame'],q['frame']=export_shape(g['frame'],out/'frame/closed_wall_bottomless_frame_with_support_blades',cfg)
    for k,v in g['retainers'].items():meshes['retainers'][k],q['retainers'][k]=export_shape(v,out/'retainers'/f'port_{k}_anti_float_retainer',cfg)
    l['core'].export(out/'vascular_core/assembly_placed_ABS_core.stl');trimesh.util.concatenate([meshes['frame'],l['core'],*meshes['retainers'].values()]).export(out/'assembly/assembled_review_only.stl')
    a=audit(l,g,meshes,cfg,out/'tables');t,trimrows,diagnostics=trim(l,g,cfg,out/'trim')
    pack=[];packing=[];cursor=0.
    for p in l['ports']:
        k=p['port_id'];part=shape_mesh(g['local'][k]['retainer']);r=Rotation.from_euler('y',90,degrees=True).as_matrix();mat=np.eye(4);mat[:3,:3]=r;part.apply_transform(mat);shift=np.array([cursor,0,0])-part.bounds[0];part.apply_translation(shift);cursor=part.bounds[1,0]+5;pack.append(part);packing.append(dict(port_id=k,local_rotation=r,translation=shift,assembly_transform=port_transform(p)))
    packed=trimesh.util.concatenate(pack);packed.export(out/'retainers/retainers_print_layout.stl');save(out/'retainers/print_layout.json',packing)
    print('V2 separate flat PLA retainer slice',flush=True);retrecord=slice_part(packed,engine,pla,out/'bambu/retainers_v2') if do_slice and pla and engine else dict(support_free_slice_passed=False)
    if frame_record.get('archive'):shutil.copy2(frame_record['archive'],out/'frame/frame_with_wall_supports.3mf')
    if retrecord.get('archive'):shutil.copy2(retrecord['archive'],out/'retainers/anti_float_retainers.3mf')
    if do_step:
        print('V2 independent-part STEP assembly',flush=True);assemble(g,source,l['transform'],out/'assembly/assembled_casting_fixture.step')
    dish=dish_check(ROOT,cfg,float(np.linalg.norm(np.diff(g['flange_bounds'][:,:2],axis=0))))
    save(out/'assembly/mesh_qc.json',q);table(out/'tables/wall_support_design.csv',[dict(**d,support_printable=frame_record.get('support_free_slice_passed',False) and preservation['passed']) for d in g['designs']])
    s=dict(architecture=cfg['architecture'],vascular=dict(source=str(inputs['core']),source_hash=digest(source),port_count=len(l['ports']),scale=1.,independent_frozen_copy=str(source)),frame=dict(four_walls_closed=a['four_walls_closed'],top_open=a['top_open'],bottom_open=a['bottom_open'],dimensions=(l['high']-l['low']).tolist(),wall_thickness_mm=cfg['wall_thickness_mm'],outer_flange_dimensions_mm=(g['flange_bounds'][1]-g['flange_bounds'][0]).tolist()),placement=dict(full_core_inside=t['core_outside_cavity_mm3']<1e-7,minimum_wall_clearance=float(min(np.r_[l['core'].bounds[0,:2]-l['low'][:2],l['high'][:2]-l['core'].bounds[1,:2]])),minimum_dish_clearance=float(l['core'].bounds[0,2]),selected_yaw_degrees=l['selected_placement']['yaw_degrees']),supports=dict(count=len(l['ports']),support_only_on_extensions=a['all_design_contacts_on_extensions'],central_contact_count=a['central_contact_count'],max_wall_span=max(x['wall_to_extension_span_mm'] for x in g['designs']),blade_line_count=g['lines'],blade_thickness_mm=line*g['lines'],real_slice_pass=frame_record.get('support_free_slice_passed',False) and preservation['passed']),ports=l['ports'],assembly=a,anti_float=a['anti_float'],trim=t,dish=dict(**dish,bottom_is_dish=True,physical_bottom_seal_validated=False),mesh_qc=q,bambu=dict(frame_slice=frame_record,retainer_slice=retrecord,blade_preservation=preservation),source_protection=verify(inputs['hashes']),assembly_step_exported=do_step,physical_status=['BOTTOM_SEAL_PHYSICALLY_UNVALIDATED','ANTI_FLOAT_FORCE_NOT_VALIDATED','MANUAL_CUT_ACCURACY_NOT_VALIDATED','PLA_PDMS_RELEASE_NOT_VALIDATED'])
    s['status'],s['gates'],s['failures']=gates(s);save(out/'internal_wall_supported_fixture_summary.json',s)
    if do_render:
        from .casting_fixture_v2_review import render,report
        render(out,l,g,meshes,diagnostics,s,trimrows);report(out,s,trimrows)
    print('V2 RESULT',s['status'],s['failures'],flush=True);return s
