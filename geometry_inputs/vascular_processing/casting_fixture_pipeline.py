"""BCF production orchestration. Anatomy inputs are copied byte-for-byte."""
from __future__ import annotations
import copy,json,shutil,time
from pathlib import Path
import numpy as np
import trimesh
import cadquery as cq
import yaml
from scipy.spatial.transform import Rotation
from .support_inputs import ROOT,discover,digest,save,table,numeric,verify
from .casting_fixture_geometry import source_layout,build_frame,export_shape,mesh_qc,shape_mesh,rigid,port_matrix
from .casting_fixture_profiles import profiles,slice_part
from .casting_fixture_audit import audit,perturbation_fit


def default_output():return ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/vascular_casting_fixture'


def source_copy(inputs,out):
    dest=Path(out)/'vascular_core/BG001_RMCA_core_with_extensions.stl';dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(inputs['core'],dest)
    if digest(dest)!=digest(inputs['core']):raise RuntimeError('FROZEN_COPY_HASH_MISMATCH')
    return dest


def dish_check(root,cfg,min_diameter):
    path=Path(cfg['dish_measurement_file']);path=path if path.is_absolute() else Path(root)/path
    if not path.exists():return dict(actual_dimensions_known=False,fit_verified=False,status='DISH_FIT_NOT_VERIFIED',measurement_file=str(path),minimum_required_flat_dish_diameter_mm=min_diameter)
    values=json.loads(path.read_text());required=['dish_inner_flat_diameter','dish_wall_height','dish_floor_thickness'];known=all(isinstance(values.get(k),(int,float)) and values[k]>0 for k in required)
    fits=bool(known and values['dish_inner_flat_diameter']>=min_diameter)
    return dict(actual_dimensions_known=known,fit_verified=fits,status='DISH_FIT_VERIFIED' if fits else 'DISH_FIT_NOT_VERIFIED',measurements=values,measurement_file=str(path),measurement_sha256=digest(path),minimum_required_flat_dish_diameter_mm=min_diameter)


def summary_gates(summary):
    a=summary['assembly'];f=summary['frame'];z=summary['anti_float'];b=summary['bambu'];q=summary['mesh_qc'];failed=[]
    gates={'DISH_FIT_IF_MEASURED':not summary['dish']['actual_dimensions_known'] or summary['dish']['fit_verified'],'FROZEN_CORE_QC':q['vascular']['passed'],'BOTTOMLESS':f['bottom_open'],'NO_CENTRAL_CONTACT':a['central_vascular_contact_count']==0,'NO_CAVITY_INTRUSION':a['illegal_fixture_intrusion_volume']<1e-8,'TOP_LOADING':a['top_loading_path']['passed'],'CLOSURE_INSERTION':a['closure_insertion_pass'],'ANTI_FLOAT':z['geometric_Z_lock'] and z['pitch_lock'] and z['roll_lock'] and z['closures_positive_Z_locked'],'FRAME_MANIFOLD':q['frame']['passed'],'CLOSURES_MANIFOLD':all(v['passed'] for v in q['closures'].values()),'PINS_MANIFOLD':all(v['passed'] for v in q['pins'].values()),'BOTTOM_CLEARANCE':a['core_to_dish_clearance_mm']>=9.999,'SOURCE_UNCHANGED':summary['source_protection']['all_unchanged'],'FRAME_PLA_SUPPORT_FREE':b.get('frame_slice',{}).get('support_free_slice_passed',False),'CLOSURE_PLA_SUPPORT_FREE':b.get('closure_slice',{}).get('support_free_slice_passed',False)}
    for k,ok in gates.items():
        if not ok:failed.append(k)
    return ('NEEDS_ADJUSTMENT' if failed else 'READY_FOR_HUMAN_REVIEW' if summary['dish']['fit_verified'] else 'READY_FOR_HUMAN_REVIEW_WITH_DISH_FIT_PENDING'),gates,failed


def assemble_step(geometry,core_path,matrix,path):
    # STL has no analytic centerline CAD. Preserve its existing triangle faces
    # as a separate faceted BRep in the STEP assembly; no refit/Boolean union.
    from OCP.StlAPI import StlAPI_Reader
    from OCP.TopoDS import TopoDS_Shape
    shape=TopoDS_Shape();reader=StlAPI_Reader()
    if not reader.Read(shape,str(core_path)):raise ValueError('FROZEN_STL_STEP_READ_FAILED')
    core=cq.Shape.cast(shape).translate(tuple(matrix[:3,3]))
    assembly=cq.Assembly(name='BCF_SEPARATE_PARTS');assembly.add(geometry['frame'],name='PLA_BOTTOMLESS_FRAME',color=cq.Color(.15,.55,.8));assembly.add(core,name='FROZEN_ABS_CORE_FACETED',color=cq.Color(.9,.2,.15))
    for k,s in geometry['closures'].items():assembly.add(s,name='PLA_CLOSURE_'+k,color=cq.Color(.1,.8,.45))
    for k,s in geometry['pins'].items():assembly.add(s,name='PLA_LOCK_PIN_'+k,color=cq.Color(.6,.2,.8))
    assembly.export(str(path));return dict(core_representation='FACETED_BREP_FROM_UNALTERED_STL_TRIANGLES',boolean_union=False)


def run(config=None,out=None,do_slice=True,do_step=True,do_render=True):
    cfg=yaml.safe_load(Path(config or ROOT/'config/vascular_casting_fixture.yaml').read_text());out=Path(out or default_output());out.mkdir(parents=True,exist_ok=True)
    gate=json.loads((ROOT/'research_retrospective/cleanup_validation.json').read_text())
    if gate['status']!='PASS':raise RuntimeError('CLEANUP_RETAINED_TEST_GATE_FAILED')
    inputs=discover(source=cfg.get('source'));save(out/'config/source_discovery.json',inputs);save(out/'config/fixture_config.json',cfg)
    core_path=source_copy(inputs,out);abs_discovery,pla,engine=profiles(inputs,cfg,out/'config/profiles')
    if pla:line=numeric(pla['process']['line_width'],numeric(pla['machine']['nozzle_diameter']))
    else:line=.22
    save(out/'config/clearance_candidates.json',dict(line_width_mm=line,source='ACTUAL_INSTALLED_PLA_PROFILE' if pla else 'GEOMETRY_ONLY_FALLBACK',radial_clearances_mm=[line*x for x in cfg['port_radial_clearance_line_multipliers']],axial_float_candidates_mm=cfg['axial_float_candidates_mm'],physical_fit_validated=False))
    from .casting_fixture_tolerance import evaluate
    tolerance=evaluate(inputs,cfg,line,out/'config/clearance_study')
    print('BCF frozen input:',inputs['core'],flush=True);layout=source_layout(inputs,cfg,line);geometry=build_frame(layout,cfg,line)
    save(out/'assembly/casting_core_transform.json',dict(transform_4x4=layout['transform'],scale_factor=1.,source_sha256=digest(core_path),dish_contact_plane=[0,0,1,0],inverse_4x4=np.linalg.inv(layout['transform'])))
    meshes={'closures':{},'pins':{}};q={'vascular':layout['source_qc'],'closures':{},'pins':{}}
    meshes['frame'],q['frame']=export_shape(geometry['frame'],out/'frame/bottomless_PLA_casting_frame',cfg)
    for k in geometry['closures']:
        meshes['closures'][k],q['closures'][k]=export_shape(geometry['closures'][k],out/'closures'/f'port_{k}_closure',cfg)
        meshes['pins'][k],q['pins'][k]=export_shape(geometry['pins'][k],out/'closures'/f'port_{k}_locking_pin',cfg)
    review=trimesh.util.concatenate([layout['core'],meshes['frame'],*meshes['closures'].values(),*meshes['pins'].values()]);review.export(out/'assembly/full_fixture_review.stl')
    save(out/'assembly/layout.json',{k:v for k,v in layout.items() if k not in ('core','central')});save(out/'assembly/mesh_qc.json',q)
    a=audit(layout,geometry,meshes,cfg,out/'tables');save(out/'assembly/assembly_sequence.json',dict(sequence=['Fix/seal frame to dish','Lower independent cleaned ABS core in -Z',*[f'Insert {p} closure then external locking pin' for p in [layout['master_port']]+[x['port_id'] for x in layout['ports'] if x['port_id']!=layout['master_port']]],'Inspect external port seal regions and bottom perimeter; physical leak test before PDMS'],no_central_contact=True))
    if do_step:
        print('BCF exporting separate-part faceted STEP assembly',flush=True);save(out/'assembly/step_representation.json',assemble_step(geometry,core_path,layout['transform'],out/'assembly/full_fixture_assembled.step'))
    # Default footprint uses the smallest searched flange that contains the
    # continuous channel with > 2 mm edge material on both sides.
    dish=dish_check(ROOT,cfg,geometry['minimum_dish_diameter'])
    slices={}
    if do_slice and pla and engine:
        print('BCF actual PLA frame slice',flush=True);slices['frame_slice']=slice_part(meshes['frame'],engine,pla,out/'bambu/frame')
        # Normal closures lie on their external key face (-90 deg about local X).
        # The datum closure uses the opposite, inner face (+90 deg about X),
        # selected by actual per-part slices; its outer end stop then grows up.
        pack=[];packing=[];cursor=0.
        for i,p in enumerate(layout['ports']):
            k=p['port_id']
            for kind,xoff in [('closure',0),('pin',15)]:
                local=shape_mesh(geometry['local_parts'][k][kind]);rot=Rotation.from_euler('x',90,degrees=True).as_matrix() if kind=='closure' and p['role']=='MASTER_DATUM' else Rotation.from_euler('x',-90,degrees=True).as_matrix()
                local.apply_transform(np.block([[rot,np.zeros((3,1))],[np.zeros((1,3)),np.ones((1,1))]]));shift=np.array([cursor,0,0])-local.bounds[0];local.apply_translation(shift);cursor=float(local.bounds[1,0]+5.);pack.append(local);packing.append(dict(part=kind,port_id=k,rotation=rot,translation=shift))
        packed=trimesh.util.concatenate(pack);save(out/'closures/print_layout.json',packing);packed.export(out/'closures/closures_and_pins_print_layout.stl')
        print('BCF actual PLA closure/pin slice',flush=True);slices['closure_slice']=slice_part(packed,engine,pla,out/'bambu/closures')
        from .support_slicer_adapter import engine_profiles
        abs_profiles=engine_profiles(abs_discovery,'Bambu');print('BCF actual independent ABS core reference slice',flush=True);slices['core_reference_slice']=slice_part(layout['core'],engine,abs_profiles,out/'bambu/ABS_core_reference');slices['core_reference_slice']['role']='REFERENCE_ONLY_EXISTING_ABS_AUTO_SUPPORT'
    summary=dict(tolerance_study=tolerance,architecture=cfg['architecture'],cleanup=json.loads((ROOT/'cleanup_result.json').read_text()),vascular=dict(source=str(inputs['core']),source_hash=digest(core_path),independent_core_copy=str(core_path),port_count=len(layout['ports']),scale_factor=1.,material='ABS',read_only=True),frame=dict(inner_dimensions=layout['high']-layout['low'],outer_dimensions=np.array([*(geometry['flange_bounds'][1,:2]-geometry['flange_bounds'][0,:2]),layout['high'][2]]),height=layout['high'][2],wall_thickness=cfg['wall_thickness_mm'],flange_width=geometry['flange_width'],flange_width_candidates_mm=cfg['flange_width_candidates_mm'],bottom_open=a['cavity']['bottom_open'],minimum_required_dish_flat_diameter=geometry['minimum_dish_diameter'],PDMS_bottom='PETRI_DISH_SMOOTH_FLOOR',floor_or_grid=False),ports=dict(master_port=layout['master_port'],secondary_ports=[p['port_id'] for p in layout['ports'] if p['role']!='MASTER_DATUM'],clearances={p['port_id']:p['radial_clearance_mm'] for p in layout['ports']},axial_float={p['port_id']:p['axial_float_mm'] for p in layout['ports']},closure_status=a['closure_insertions'],designs=layout['ports']),assembly=dict(top_loading_path=a['top_loading_path'],closure_insertion_pass=a['closure_insertion_pass'],central_vascular_contact_count=a['central_vascular_contact_count'],illegal_fixture_intrusion_volume=a['cavity']['illegal_fixture_intrusion_volume'],core_to_dish_clearance_mm=a['core_to_dish_clearance_mm'],central_min_clearance_mm=a['central_min_clearance_mm'],nominal_capture_sites=len(layout['ports']),nominal_hard_contact_count=a['nominal_hard_contact_count'],Boolean_union_of_core_and_fixture=False),anti_float=a['anti_float'],dish=dish,seal=a['seal'],mesh_qc=q,bambu=slices,source_protection=verify(inputs['hashes']),physical_status=['PHYSICAL_SEAL_NOT_VALIDATED','PHYSICAL_LEAK_TEST_NOT_DONE','ANTI_FLOAT_FORCE_NOT_PHYSICALLY_VALIDATED','PLA_ABS_FIT_NOT_PHYSICALLY_VALIDATED'])
    summary['status'],summary['gates'],summary['failures']=summary_gates(summary);save(out/'vascular_casting_fixture_summary.json',summary)
    if do_render:
        from .casting_fixture_review import render,report
        render(out,layout,geometry,meshes,summary);report(out,summary)
    print('BCF RESULT',summary['status'],summary['failures'],flush=True);return summary
