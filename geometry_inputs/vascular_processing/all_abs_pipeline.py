"""s1-10 all-ABS sacrificial production stage, with immutable s1-9 history."""
from __future__ import annotations
import json,shutil
from pathlib import Path
import numpy as np
import trimesh,cadquery as cq,yaml
from scipy.spatial.transform import Rotation
from .support_inputs import ROOT,discover,digest,save,table,verify
from .casting_fixture_geometry import shape_mesh,export_shape
from .casting_fixture_v2_geometry import port_transform
from .all_abs_geometry import layout,build
from .all_abs_audit import assembly,solvent,trim
from .all_abs_slicing import profiles,slice_abs,feature_preservation


def default_output():return ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/all_abs_sacrificial_casting_system'


def protected_snapshot(inputs,out):
    path=out/'source_protection_start.json'
    if path.exists():return json.loads(path.read_text())
    hashes=dict(inputs['hashes']);old=default_output().parent/'vascular_casting_fixture'
    paths=list(old.rglob('*'))+list(ROOT.glob('s1-9*.py'))+list((ROOT/'vascular_processing').glob('casting_fixture*.py'))+list((ROOT/'tests').glob('test_s1_9*.py'))+[ROOT/'config/internal_wall_supported_fixture_v2.yaml',ROOT/'config/vascular_casting_fixture.yaml']
    for p in paths:
        if p.is_file():hashes[str(p)]=digest(p)
    save(path,hashes);return hashes


def gates(s):
    a=s['assembly'];t=s['trim'];v=s['solvent_access'];b=s['bambu'];q=s['mesh_qc'];eps=1e-6
    checks=dict(SOURCE_IMMUTABLE=s['source_protection']['all_unchanged'],ALL_ABS=s['all_components_ABS'],RIGID_ONLY=s['vascular']['scale_factor']==1.,WALLS_CLOSED=a['walls_closed'],BOTTOM_CLOSED=a['bottom_closed'],TOP_OPEN=a['top_open'],FULL_CORE_INSIDE=t['core_outside_cavity_mm3']<eps,CORE_INSERTION=a['top_insertion_pass'],RETAINER_INSERTION=a['retainer_insertion_pass'],EXTENSION_CONTACTS_ONLY=a['all_contacts_on_approved_extensions'],NO_CENTRAL_CONTACT=a['central_tree_contacts']==0,NO_PRESS_FIT=a['core_interference_mm3']<eps,ANTI_FLOAT=a['anti_float_geometric_pass'],ALL_PORTS_ACCESS=v['all_ports_accessible'],ALL_RETAINERS_ACCESS=v['all_retainers_accessible'],SWC_ALL_BRANCHES_ACCESS=v['all_vascular_branches_reachable'],AUX_VOID_REMOVED=abs(t['auxiliary_void_volume_after_trim'])<eps and abs(t['ghost_missing_volume_mm3'])<eps and abs(t['ghost_extra_volume_mm3'])<eps,CENTRAL_TRIM_UNTOUCHED=abs(t['central_vascular_removed_mm3'])<eps,UNRELATED_PORT_TRIM_UNTOUCHED=abs(t['maximum_unrelated_port_cut_mm3'])<eps,PORTS_EXPOSED=t['all_ports_exposed'],NONPORT_CLEARANCE=t['minimum_nonport_clearance']>0,PDMS_EXPORT_VALID=t['diagnostic_mesh_qc_pass'],PRINT_MESHES_VALID=q['mold']['passed'] and q['vascular']['passed'] and all(x['passed'] for x in q['retainers'].values()),MOLD_REAL_ABS_SUPPORT_OFF=b['mold_auto_support_off_slice'].get('slice_passed',False),RETAINER_REAL_ABS_SUPPORT_OFF=b['retainer_auto_support_off_slice'].get('slice_passed',False),CORE_REFERENCE_REAL_ABS=b['vascular_reference_slice'].get('slice_passed',False),THIN_FEATURE_TOOLPATHS=b['feature_preservation'].get('passed',False),ASSEMBLY_STEP=s['assembly_step_exported'],NO_CORE_VAPOR=s['smoothing']['vascular_core_forbidden'],NO_UNCALIBRATED_COMPENSATION=not s['smoothing']['software_compensation_applied'])
    checks['REVIEW_FIGURES_COMPLETE']=s.get('review_figures_complete',False)
    checks['EXTERNAL_MOLD_SOLVENT_ACCESS']=v['external_mold_access']
    failed=[k for k,v in checks.items() if not v];return ('NEEDS_ADJUSTMENT' if failed else 'READY_FOR_HUMAN_REVIEW_WITH_PROCESS_PARAMETERS_UNCALIBRATED'),checks,failed


def assemble(g,source,m,path):
    from OCP.StlAPI import StlAPI_Reader
    from OCP.TopoDS import TopoDS_Shape
    raw=TopoDS_Shape()
    if not StlAPI_Reader().Read(raw,str(source)):raise ValueError('FROZEN_STL_STEP_READ_FAILED')
    core=cq.Shape.cast(raw);angle=float(np.degrees(np.arctan2(m[1,0],m[0,0])));core=core.rotate((0,0,0),(0,0,1),angle).translate(tuple(m[:3,3]));a=cq.Assembly(name='ALL_ABS_SEPARATELY_PRINTED_PARTS');a.add(core,name='FROZEN_ABS_CORE',color=cq.Color(.85,.24,.16));a.add(g['mold'],name='ABS_SACRIFICIAL_MOLD',color=cq.Color(.2,.55,.75))
    for k,v in g['retainers'].items():a.add(v,name='ABS_RETAINER_'+k,color=cq.Color(.25,.75,.4))
    a.export(str(path))


def smoothing_dimensions(l,g,cfg,out):
    rows=[]
    for p,d in zip(l['ports'],g['design']):
        for feature,value in [('cradle_radial_clearance',p['cradle_clearance_mm']),('retainer_guide_clearance',d['retainer_guide_clearance_mm']),('retainer_top_clearance',p['anti_float_gap_mm']),('fixture_blade_thickness',d['blade_thickness_mm']),('access_channel_minimum_width',d['access_channel_minimum_width_mm']),('conformal_land_nominal_gap',0.)]:
            rows.append(dict(port_id=p['port_id'],feature=feature,nominal_mm=value,dimensionally_calibrated=False,software_compensation_applied=False,check_required='Dry-fit and confirm solvent contact; a PDMS membrane at a printed contact invalidates that access edge' if feature=='conformal_land_nominal_gap' else 'Measure after complete drying; no negative fit or blocked access allowed',status='VAPOR_SMOOTHING_DIMENSIONAL_EFFECT_UNCALIBRATED'))
    for feature,value in [('wall_thickness',cfg['wall_thickness_mm']),('bottom_thickness',cfg['bottom_thickness_mm']),('shallow_guide_depth',cfg['trim_guide_depth_mm'])]:rows.append(dict(port_id='ALL',feature=feature,nominal_mm=value,dimensionally_calibrated=False,software_compensation_applied=False,check_required='Measure printed and optionally treated mold; uncalibrated dimensions',status='VAPOR_SMOOTHING_DIMENSIONAL_EFFECT_UNCALIBRATED'))
    table(out/'vapor_smoothing_critical_dimensions.csv',rows)


def run(config=None,out=None,do_slice=True,do_step=True,do_render=True):
    cfg=yaml.safe_load(Path(config or ROOT/'config/all_abs_sacrificial_casting_system.yaml').read_text());out=Path(out or default_output()).resolve()
    for name in ['vascular_core','outer_mold','retainers','solvent_access','trim','diagnostic_PDMS','assembly','QC','bambu','tables','reports']: (out/name).mkdir(parents=True,exist_ok=True)
    inputs=discover(source=cfg.get('source'));hashes=protected_snapshot(inputs,out);verify(hashes);save(out/'reports/source_discovery.json',inputs);save(out/'reports/config.json',cfg)
    source=out/'vascular_core/BG001_RMCA_frozen_ABS_core.stl';shutil.copy2(inputs['core'],source)
    if digest(source)!=digest(inputs['core']):raise ValueError('FROZEN_COPY_HASH_MISMATCH')
    values,actual,ref,engine=profiles(inputs,out/'bambu');line=values['line_width'];l=layout(inputs,cfg,line);print('s1-10: build new all-ABS solids',flush=True);g=build(l,cfg,line)
    save(out/'assembly/layout.json',{k:v for k,v in l.items() if k not in ['core','central','source_ports']});save(out/'assembly/casting_core_transform.json',dict(transform_4x4=l['transform'],inverse_4x4=np.linalg.inv(l['transform']),scale_factor=1.,source_hash=digest(source),frame='OUTER_BOTTOM_Z0'))
    table(out/'tables/rigid_placement_comparison.csv',l['placements']);table(out/'tables/extension_fixture_zones.csv',l['ports']);table(out/'tables/fixture_design.csv',g['design']);smoothing_dimensions(l,g,cfg,out/'tables')
    meshes={'retainers':{}};q={'retainers':{},'vascular':l['source_qc']};meshes['mold'],q['mold']=export_shape(g['mold'],out/'outer_mold/ABS_sacrificial_outer_mold',cfg)
    for k,part in g['retainers'].items():meshes['retainers'][k],q['retainers'][k]=export_shape(part,out/'retainers'/f'port_{k}_ABS_retainer',cfg)
    l['core'].export(out/'vascular_core/assembly_placed_ABS_core.stl');trimesh.util.concatenate([meshes['mold'],l['core'],*meshes['retainers'].values()]).export(out/'assembly/all_ABS_review.stl')
    a=assembly(l,g,cfg,out/'tables');v=solvent(l,g,cfg,inputs,out/'solvent_access');t,trimrows,diagnostics=trim(l,g,cfg,out/'trim')
    pack=[];packing=[];cursor=0.
    for p in l['ports']:
        k=p['port_id'];part=shape_mesh(g['local_parts'][k]['retainer']);r=Rotation.from_euler('y',90,degrees=True).as_matrix();m=np.eye(4);m[:3,:3]=r;part.apply_transform(m);shift=np.array([cursor,0,0])-part.bounds[0];part.apply_translation(shift);cursor=part.bounds[1,0]+5;pack.append(part);packing.append(dict(port_id=k,local_rotation=r,translation=shift,assembly_transform=port_transform(p)))
    packed=trimesh.util.concatenate(pack);packed.export(out/'retainers/retainers_print_layout.stl');save(out/'retainers/print_layout.json',packing)
    records=[]
    for name,part,profile,off in [('outer_mold',meshes['mold'],actual,True),('retainers',packed,actual,True),('vascular_reference',l['core'],ref,False)]:
        print('s1-10: actual ABS slice',name,flush=True);records.append(slice_abs(part,engine,profile,out/'bambu'/name,off) if do_slice else dict(slice_passed=False,reason='GEOMETRY_ONLY'))
    mr,rr,vr=records;pres=feature_preservation(l,g,mr,rr,packing,out/'bambu') if do_slice else dict(passed=False)
    for rec,dest in [(mr,out/'outer_mold/ABS_sacrificial_outer_mold.3mf'),(rr,out/'retainers/ABS_retainers.3mf'),(vr,out/'vascular_core/ABS_core_reference.3mf')]:
        if rec.get('archive'):shutil.copy2(rec['archive'],dest)
    if do_step:
        print('s1-10: export independent-component faceted STEP assembly',flush=True);assemble(g,source,l['transform'],out/'assembly/all_ABS_casting_system.step')
    protection=verify(hashes);protection.update(start_hashes=hashes,end_hashes={p:digest(p) for p in hashes},mismatch_count=0);save(out/'source_protection.json',protection)
    s=dict(architecture=cfg['architecture'],all_components_ABS=True,vascular=dict(source=str(inputs['core']),source_hash=digest(source),scale_factor=1.,port_count=len(l['ports']),source_copy=str(source),never_vapor_smoothed=True),outer_mold=dict(walls_closed=a['walls_closed'],bottom_closed=a['bottom_closed'],top_open=a['top_open'],wall_thickness=cfg['wall_thickness_mm'],bottom_thickness=cfg['bottom_thickness_mm'],inner_dimensions=l['high']-l['low'],outer_dimensions=np.diff(g['outer_bounds'],axis=0)[0]),assembly=a,ports=l['ports'],smoothing=dict(outer_mold_targeted=True,optional_post_print_only=True,target_surfaces=['INNER_BOTTOM','INNER_SIDE_WALLS'],vascular_core_forbidden=True,dimensionally_calibrated=False,software_compensation_applied=False,status='VAPOR_SMOOTHING_NOT_DIMENSIONALLY_CALIBRATED'),solvent_access=v,voids={k:t[k] for k in ['auxiliary_void_volume_before_trim','auxiliary_void_volume_after_trim']},trim=t,bambu=dict(mold_auto_support_off_slice=mr,retainer_auto_support_off_slice=rr,vascular_reference_slice=vr,feature_preservation=pres),mesh_qc=q,source_protection={k:val for k,val in protection.items() if k not in ['start_hashes','end_hashes']},assembly_step_exported=do_step,physical_validation=dict(print=False,dissolution=False,vapor_dimensions=False,anti_float_force=False,manual_trim=False,PDMS_swelling_recovery=False),no_printer_delivery=True)
    s['status'],s['gates'],s['failures']=gates(s);save(out/'all_abs_sacrificial_casting_summary.json',s)
    from .all_abs_review import report,render,NAMES
    if do_render:render(out,l,g,meshes,diagnostics,s)
    s['review_figures_complete']=do_render and all((out/'QC'/(n+'.png')).is_file() for n in NAMES)
    s['status'],s['gates'],s['failures']=gates(s);save(out/'all_abs_sacrificial_casting_summary.json',s)
    report(out,s,trimrows)
    save(out/'reports/artifact_manifest.json',{str(p.relative_to(out)):digest(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='artifact_manifest.json' and p.suffix!='.log'})
    print('s1-10 RESULT',s['status'],s['failures'],flush=True);return s
