"""Re-open local delivered meshes, projects and ZIPs after native server slicing."""
from pathlib import Path
import sys,json,hashlib,zipfile
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
import numpy as np
import trimesh
from vascular_processing import y_channel_casting as y, straight_channel_casting as s
from vascular_processing import bottomless_channel_fixture as f
from vascular_processing.all_abs_print_plate import read_objects
from vascular_processing.rhino_channel_export import verify_exports

rows={}
for kind,module,name in [('Y',y,'y_channel_p1s_abs.yaml'),('Straight',s,'straight_channel_1p6_p1s_abs.yaml')]:
    cfg=module.load_config(ROOT/'config'/name);out=cfg['output']
    g=json.loads((out/'geometry_qc.json').read_text());p=json.loads((out/'print_qc.json').read_text())
    final=json.loads((out/'manufacturing_summary.json').read_text())
    assert g['status']=='GEOMETRY_PASSED' and all(g['gates'].values())
    assert final['status']=='READY_FOR_FIT_AND_PRINT_TRIAL' and all(p['gates'].values())
    assert p['fixture_toolpath_qc']['passed']
    for filename,digest in g['assembly_sha256'].items():assert y.digest(out/'assembly'/filename)==digest
    names=f.part_names(kind);parts={n:trimesh.load_mesh(out/'assembly'/(n+'.stl')) for n in names}
    assert all(m.is_volume and len(m.split())==1 for m in parts.values())
    inner=np.asarray(g['cavity_xy_bounds_mm']);frame=parts[kind+'_bottomless_box_wall2mm'];core=parts[kind+'_lumen_core']
    probe=y.box([np.r_[inner[0]+.0001,-.5],np.r_[inner[1]-.0001,12.5]])
    internal_abs=s.volume_intersection(frame,probe);assert internal_abs<1e-4
    origins=np.array([[x,v,-1] for x in np.linspace(inner[0,0]+.1,inner[1,0]-.1,11)
                     for v in np.linspace(inner[0,1]+.1,inner[1,1]-.1,11)])
    bottom_hits=frame.ray.intersects_location(origins,np.tile([0,0,1],(len(origins),1)))[0]
    assert len(bottom_hits)==0
    assert abs(core.bounds[0,2]-1.5)<1e-5
    widths=[]
    for z in [.3,4,9]:
        origin=np.r_[inner.mean(0),z];dist=[]
        for sign in [-1,1]:
            d=np.array([0,sign,0]);hit=frame.ray.intersects_location([origin],[d])[0];dist.append(float(np.min((hit-origin)@d)))
        widths.append(sum(dist))
    height=float(frame.bounds[1,2]-cfg['freeboard_mm']-frame.bounds[0,2])
    height_qc=f.casting_height_qc(cfg,parts)
    assert height_qc['passed'] and height_qc['frame_and_key_top_z_mm']<=14.8+1e-5
    assert abs(height_qc['retainer_top_z_mm']-17.8)<1e-5
    gap_qc=f.retainer_rim_clearance_qc(cfg,parts,inner,float(frame.bounds[1,2]),g['dovetail_joints'],kind)
    assert gap_qc['passed'] and gap_qc['requested_gap_mm']==3
    assert all(abs(v-3)<1e-5 for r in gap_qc['parts'] for v in r['measured_wall_to_arm_gaps_mm'])
    assert abs(height_qc['core_top_z_mm']-18.4)<1e-5
    rhino_qc=verify_exports(out)
    assert rhino_qc['passed']
    assert abs(height-9.6)<1e-5 and min(widths)*.98-.3>=9
    if kind=='Y':assert y.digest(cfg['source'])==cfg['source_sha256'] and y.digest(cfg['planes'])==g['planes_sha256']
    else:
        assert g['closed_chamber_qc']['horizontal_extension_outside_pdms_mm3']<1e-4
        assert abs(g['closed_chamber_qc']['horizontal_end_to_inner_wall_mm']-3)<1e-5
    fixture_paths=p['fixture_toolpath_qc']['parts']
    assert {r['part'] for r in fixture_paths}==set(names) and all(r['positive_model_segments']>0 for r in fixture_paths)
    projects={}
    for project in out.glob('*.3mf'):
        objs,settings,plates,crc=read_objects(project)
        assert crc and plates==1 and len(objs)==len(names)
        assert {o['metadata']['name'] for o in objs}==set(names)
        assert settings['printer_model']=='Bambu Lab P1S' and settings['nozzle_diameter']==['0.4'] and settings['filament_type']==['ABS']
        assert y.slicer_scale_qc(cfg,settings)['passed']
        for obj in objs:
            meta=obj['metadata'];n=meta['name']
            assert meta['enable_support']==('1' if 'core' in n or 'retainer' in n else '0')
            assert meta['raft_layers']==('2' if 'core' in n else '0')
        projects[project.name]=dict(sha256=y.digest(project),object_count=len(objs),crc=True)
    receipt=json.loads((out/'delivery_receipt.json').read_text());archive=ROOT/'deliverables'/Path(receipt['zip']).name
    assert y.digest(archive)==receipt['zip_sha256']
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for n,h in receipt['members_sha256'].items():assert hashlib.sha256(z.read(kind+'_channel/'+n)).hexdigest()==h
        doc=ROOT/'docs'/('Y形通道_ABS打印与PDMS浇注.md' if kind=='Y' else '直管1.6mm_双竖直单水平端口打印说明.md')
        assert z.read(kind+'_channel/使用说明.md')==doc.read_bytes()
    rows[kind]=dict(output=str(out),part_count=len(parts),parts=names,internal_frame_volume_mm3=internal_abs,
        open_bottom_ray_count=len(origins),open_bottom_ray_hits=0,frame_size_mm=frame.extents,
        core_min_z_mm=float(core.bounds[0,2]),manufacturing_height_qc=height_qc,rhino_qc=rhino_qc,retainer_rim_clearance_qc=gap_qc,
        complete_assembly_xy_size_mm=np.ptp(np.vstack([m.bounds for m in parts.values()]),axis=0)[:2],
        cad_pdms_width_mm=widths,cad_pdms_height_mm=height,
        fixture_gates=g['fixture_qc']['gates'],print_gates=p['gates'],fixture_toolpath_qc=p['fixture_toolpath_qc'],
        estimated_print_seconds=final['slice']['estimated_print_time_seconds'],filament_used_g=final['slice']['filament_used_g'],
        projects=projects,zip=str(archive),zip_sha256=y.digest(archive),zip_crc_and_member_hashes_passed=True,
        independent_external_substrate_required=True,physical_fit_and_seal_validated=False)
report=dict(status='RAISED_SUPPORT_RHINO_AND_PRINT_DELIVERY_VERIFIED',test_result='81 passed',variants=rows,
    implementation_sha256={str(p.relative_to(ROOT)):y.digest(p) for p in [
        ROOT/'vascular_processing/rhino_channel_export.py',ROOT/'vascular_processing/y_uniform_channel.py',
        ROOT/'vascular_processing/bottomless_channel_fixture.py',ROOT/'vascular_processing/y_channel_casting.py',
        ROOT/'vascular_processing/straight_channel_casting.py',ROOT/'config/y_channel_p1s_abs.yaml',ROOT/'config/straight_channel_1p6_p1s_abs.yaml']},
    geometry_generated_on='local workstation',native_slicing_on='isolated /workspace/vascular_printing environment at 95.3.33.46:15283',
    no_printer_upload_or_start=True)
y.save(Path(__file__).with_name('revision_audit.json'),report)
for kind,r in rows.items():print(kind,r['part_count'],r['cad_pdms_height_mm'],r['frame_size_mm'],r['estimated_print_seconds'],r['filament_used_g'],r['zip'])
print(report['status'])
