"""Independent delivery check; run with the project Python after sync."""
from pathlib import Path
import sys, json, hashlib, zipfile
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT))
import numpy as np
import trimesh, yaml
from vascular_processing import y_channel_casting as y
from vascular_processing import straight_channel_casting as s
from vascular_processing.y_uniform_channel import local_section_diameter
from vascular_processing.all_abs_print_plate import read_objects

rows={}
for kind,name,prefix in [('Y','y_channel_p1s_abs.yaml','Y'),('Straight','straight_channel_1p6_p1s_abs.yaml','Straight')]:
    cfg=yaml.safe_load((ROOT/'config'/name).read_text());out=ROOT/cfg['output']
    g=json.loads((out/'geometry_qc.json').read_text())
    p=json.loads((out/'print_qc.json').read_text());final=json.loads((out/'manufacturing_summary.json').read_text())
    assert g['status']=='GEOMETRY_PASSED' and all(g['gates'].values())
    assert final['status']=='READY_FOR_FIT_AND_PRINT_TRIAL' and all(p['gates'].values())
    for file,expected in g['assembly_sha256'].items():assert y.digest(out/'assembly'/file)==expected
    mold=trimesh.load_mesh(out/'assembly'/(prefix+'_casting_box_base_2mm.stl'))
    inner=np.asarray(g['cavity_xy_bounds_mm'] if kind=='Y' else g['design']['inner'])
    widths=[]
    for height in [3.,6.,11.1]:
        origin=np.r_[inner.mean(0),height];dist=[]
        for sign in [-1,1]:
            direction=np.array([0,sign,0]);hits=mold.ray.intersects_location([origin],[direction])[0]
            dist.append(float(np.min((hits-origin)@direction)))
        widths.append(sum(dist))
    floor_hits=mold.ray.intersects_location([np.r_[inner.mean(0),1]],[[0,0,1]])[0]
    floor=float(floor_hits[:,2].min());h=float(mold.bounds[1,2]-cfg['freeboard_mm']-floor)
    assert abs(floor-2)<1e-5 and abs(h-9.6)<1e-5
    assert min(widths)*.98-.3 >= 9 and h*.98-.3 >= 9
    for dimension,qckey in [('width','pdms_width_qc'),('height','pdms_height_qc')]:
        assert g[qckey]['allowance']['finished_minimum_mm']==9
        assert g[qckey]['allowance']['required_cad_minimum_mm']==9.6
        assert not g[qckey]['allowance']['physical_finished_dimensions_validated']
    if kind=='Y':
        channel=trimesh.load_mesh(out/'assembly/Y_channel_without_positioning_collars_reference.stl')
        diameters=[local_section_diameter(channel,r['center_mm'],r['normal']) for r in g['section_checks']]
        assert y.digest(ROOT/cfg['source'])==cfg['source_sha256']
        pdms_ml=g['pdms_fill_volume_ml']
    else:
        channel=trimesh.load_mesh(out/'assembly/Straight_channel_without_collars_reference.stl')
        diameters=[s.section_diameter(channel,[r['x_mm'],0,4.3],[1,0,0]) for r in g['straight_section_checks']]
        diameters += [local_section_diameter(channel,r['center_mm'],r['normal']) for r in g['elbow_riser_section_checks']]
        assert g['chamber_mode']=='closed_wall_internal_end'
        assert abs(g['closed_chamber_qc']['horizontal_end_to_inner_wall_mm']-3)<1e-5
        assert g['closed_chamber_qc']['horizontal_extension_outside_pdms_mm3']<1e-4
        pdms_ml=g['pdms_net_fill_ml']
    assert max(abs(np.array(diameters)-1.6))<.003
    projects={}
    for project in out.glob('*.3mf'):
        objects,settings,plates,crc=read_objects(project)
        assert crc and plates==1 and len(objects)==(5 if kind=='Y' else 4)
        assert settings['printer_model']=='Bambu Lab P1S'
        assert settings['nozzle_diameter']==['0.4'] and settings['filament_type']==['ABS']
        assert y.slicer_scale_qc(cfg,settings)['passed']
        projects[project.name]=dict(sha256=y.digest(project),object_count=len(objects),no_duplicate_compensation=True)
    receipt=json.loads((out/'delivery_receipt.json').read_text());archive=ROOT/'deliverables'/Path(receipt['zip']).name
    assert y.digest(archive)==receipt['zip_sha256']
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for file,expected in receipt['members_sha256'].items():
            assert hashlib.sha256(z.read(prefix+'_channel/'+file)).hexdigest()==expected
        doc=ROOT/'docs'/('Y形通道_ABS打印与PDMS浇注.md' if kind=='Y' else '直管1.6mm_双竖直单水平端口打印说明.md')
        assert z.read(prefix+'_channel/使用说明.md')==doc.read_bytes()
    rows[kind]=dict(output=str(out),mold_outer_size_mm=mold.extents,
        inner_width_ray_measurements_mm=widths,cad_fill_height_mm=h,
        budgeted_finished_minimum_width_mm=min(widths)*.98-.3,
        budgeted_finished_minimum_height_mm=h*.98-.3,
        cad_diameter_measurement_range_mm=[min(diameters),max(diameters)],diameter_sections_count=len(diameters),
        pdms_net_volume_ml=pdms_ml,projects=projects,zip=str(archive),zip_sha256=y.digest(archive),zip_bytes=archive.stat().st_size,
        all_geometry_gates=g['gates'],all_print_gates=p['gates'],zip_crc_and_hashes_passed=True,
        slicer_estimated_seconds=final['slice']['estimated_print_time_seconds'],
        slicer_filament_g=final['slice']['filament_used_g'],physical_dimensions_validated=False)
report=dict(status='DIGITAL_DELIVERY_VERIFIED',tests='57 passed',variants=rows,
    sources=y.dimension_budget(cfg,'height')['sources'],
    implementation_sha256={str(p.relative_to(ROOT)):y.digest(p) for p in [
        ROOT/'vascular_processing/abs_dimensional_allowance.py',ROOT/'vascular_processing/y_channel_casting.py',
        ROOT/'vascular_processing/y_uniform_channel.py',ROOT/'vascular_processing/straight_channel_casting.py',
        ROOT/'config/y_channel_p1s_abs.yaml',ROOT/'config/straight_channel_1p6_p1s_abs.yaml']},
    method='2% cooling shrinkage design budget + 0.30 mm total dimension reserve, round CAD upward to 0.20 mm',
    limitations='Uncalibrated ABS design budget; physical printed/cured dimensions and PDMS cure shrinkage have not been measured.')
y.save(Path(__file__).with_name('revision_audit.json'),report)
print(json.dumps(report,ensure_ascii=False,indent=2,default=lambda v:v.tolist() if isinstance(v,np.ndarray) else v))
