"""Reopen delivered files and compare with the protected pre-magnet variants."""
from pathlib import Path
import sys,json,zipfile,hashlib
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT))
import numpy as np
import trimesh
from vascular_processing import y_channel_casting as y,straight_channel_casting as s
from vascular_processing import channel_magnet_slots as m
from vascular_processing.rhino_channel_export import verify_exports
from vascular_processing.all_abs_print_plate import read_objects

baseline=json.loads(Path(__file__).with_name('baseline.json').read_text())
rows={};old_count=0
for config,old in baseline.items():
    kind='Y' if config.startswith('y_') else 'Straight';module=y if kind=='Y' else s
    cfg=module.load_config(ROOT/'config'/config);out=cfg['output'];prior=ROOT/old['original_output']
    for rel,h in old['files'].items():
        assert y.digest(ROOT/rel)==h,rel+' original changed';old_count+=1
    assert y.digest(ROOT/'config/history/before_magnet_slots_20261012'/config)==old['config_sha256']
    geom=json.loads((out/'geometry_qc.json').read_text());printed=json.loads((out/'print_qc.json').read_text())
    final=json.loads((out/'manufacturing_summary.json').read_text());slots=json.loads((out/'magnet_slots_qc.json').read_text())
    assert geom['status']=='GEOMETRY_PASSED' and all(geom['gates'].values())
    assert final['status']=='READY_FOR_FIT_AND_PRINT_TRIAL' and all(printed['gates'].values())
    assert printed['magnet_toolpath_qc']['passed'] and printed['fixture_toolpath_qc']['passed']
    assert slots==geom['magnet_slots'] and slots['count']==2
    assert y.digest(out/'reference/magnet_slot.3dm')==cfg['magnet_slots']['source_sha256']
    assert y.digest(ROOT/cfg['magnet_slots']['source'])==cfg['magnet_slots']['source_sha256']
    parts={name:trimesh.load_mesh(out/'assembly'/(name+'.stl')) for name in geom['production_part_names']}
    assert len(parts)==(8 if kind=='Y' else 6)
    protected=[]
    for name,mesh in parts.items():
        assert mesh.is_volume and len(mesh.split())==1
        if 'box' not in name:
            assert y.digest(out/'assembly'/(name+'.stl'))==y.digest(prior/'assembly'/(name+'.stl'))
            protected.append(name)
    frame_name=kind+'_bottomless_box_wall2mm';frame=parts[frame_name]
    old_frame=trimesh.load_mesh(prior/'assembly'/(frame_name+'.stl'))
    actual_qc=m.check_slots(frame,old_frame,{n:o for n,o in parts.items() if n!=frame_name},
        np.array(geom['cavity_xy_bounds_mm']),cfg,slots['source'],slots['placements'])
    assert actual_qc['passed']
    for filename,h in geom['assembly_sha256'].items():assert y.digest(out/'assembly'/filename)==h
    rhino=verify_exports(out)
    objects_audit={}
    layouts=json.loads((out/'print_layout.json').read_text())['records']
    T=np.array(next(r for r in layouts if r['name']==frame_name)['casting_to_plate_transform'])
    for project in out.glob('*.3mf'):
        objects,settings,plates,crc=read_objects(project)
        assert crc and plates==1 and len(objects)==len(parts)
        assert {o['metadata']['name'] for o in objects}==set(parts)
        assert settings['printer_model']=='Bambu Lab P1S' and settings['nozzle_diameter']==['0.4'] and settings['filament_type']==['ABS']
        assert y.slicer_scale_qc(cfg,settings)['passed']
        obj=next(o for o in objects if o['metadata']['name']==frame_name)
        assert obj['metadata']['enable_support']=='0'
        for slot in slots['placements']:
            c=(T@np.r_[slot['pocket_floor_center_mm'],1])[:3]
            assert not obj['mesh'].contains([c+[0,0,1.]])[0]
            assert obj['mesh'].contains([c-[0,0,.5]])[0]
        objects_audit[project.name]=dict(sha256=y.digest(project),count=len(objects),crc=True,
            two_pockets_open=True,actual_printer=settings['printer_model'],actual_nozzle=settings['nozzle_diameter'],material=settings['filament_type'])
    assert len(objects_audit)==2
    receipt=json.loads((out/'delivery_receipt.json').read_text());archive=ROOT/'deliverables'/Path(receipt['zip']).name
    assert y.digest(archive)==receipt['zip_sha256']
    with zipfile.ZipFile(archive) as z:
        assert z.testzip() is None
        for member,h in receipt['members_sha256'].items():
            assert hashlib.sha256(z.read(kind+'_channel/'+member)).hexdigest()==h
        assert z.read(kind+'_channel/磁铁槽使用说明.md')==(ROOT/'docs/三版模具_外置磁铁槽使用说明.md').read_bytes()
        assert z.read(kind+'_channel/reference/magnet_slot.3dm')==(ROOT/cfg['magnet_slots']['source']).read_bytes()
    rows[config]=dict(output=str(out),protected_parts_byte_identical=protected,actual_magnet_qc=actual_qc,
        cavity_xy_size_mm=geom['cavity_xy_size_mm'],frame_size_mm=frame.extents,
        heights=geom['manufacturing_height_qc'],pocket_positions=slots['placements'],
        original_frame_missing_volume_mm3=actual_qc['original_frame_missing_volume_mm3'],
        geometry_gates=geom['gates'],native_print_gates=printed['gates'],
        native_magnet_toolpath_qc=printed['magnet_toolpath_qc'],rhino_files=rhino['exported_files'],
        projects=objects_audit,zip_sha256=y.digest(archive),zip_bytes=archive.stat().st_size,zip=str(archive),
        estimated_print_seconds=final['slice']['estimated_print_time_seconds'],filament_g=final['slice']['filament_used_g'])
    print(config,'DIGITAL_PASS',frame.extents,rows[config]['estimated_print_seconds'],rows[config]['filament_g'],flush=True)
code=['vascular_processing/channel_magnet_slots.py','vascular_processing/bottomless_channel_fixture.py',
      'vascular_processing/y_channel_casting.py','vascular_processing/straight_channel_casting.py',
      'vascular_processing/rhino_channel_export.py',*['config/'+n for n in baseline],
      'data/fixtures/magnet_slot.3dm']
report=dict(status='THREE_MAGNET_FIXTURES_DIGITALLY_VERIFIED',variants=rows,
    protected_original_files_verified=old_count,test_result='107 passed',
    implementation_sha256={n:y.digest(ROOT/n) for n in code},
    physical_fit_tested=False,magnetic_force_tested=False,no_printer_upload_or_start=True)
y.save(Path(__file__).with_name('delivery_audit.json'),report)
print(report['status'],old_count,flush=True)
