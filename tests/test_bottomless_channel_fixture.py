"""Actual joint geometry, open interior, core preservation and native layout."""
import copy
import json
import numpy as np
import pytest
import trimesh
import yaml
from vascular_processing import y_channel_casting as y, straight_channel_casting as s
from vascular_processing import bottomless_channel_fixture as f
from vascular_processing.all_abs_print_plate import read_objects


@pytest.fixture(scope='module',params=['Y','Straight'])
def real(request,tmp_path_factory):
    kind=request.param;module=y if kind=='Y' else s
    name='y_channel_p1s_abs.yaml' if kind=='Y' else 'straight_channel_1p6_p1s_abs.yaml'
    cfg=module.load_config(y.ROOT/'config/history/before_magnet_slots_20261012'/name);out=tmp_path_factory.mktemp('bottomless_'+kind)
    parts,g=module.geometry(cfg,out)
    return kind,module,name,cfg,out,parts,g


def test_full_inner_prism_has_no_floor_and_four_walls_are_2mm(real):
    kind,_,_,cfg,_,parts,g=real;inner=np.array(g['cavity_xy_bounds_mm'])
    frame=parts[kind+'_bottomless_box_wall2mm']
    probe=y.box([np.r_[inner[0]+.001,-1],np.r_[inner[1]-.001,11.9]])
    assert s.volume_intersection(frame,probe)<1e-4
    origins=np.array([[x,v,-1] for x in np.linspace(*inner[:,0]+[.1,-.1],11) for v in np.linspace(*inner[:,1]+[.1,-.1],11)])
    assert len(frame.ray.intersects_location(origins,np.tile([0,0,1],(len(origins),1)))[0])==0
    for z in [.1,5,float(frame.bounds[1,2])-.2]:
        hits=frame.ray.intersects_location([[*inner.mean(0),z]],[[1,0,0]])[0]
        assert sorted(hits[:,0])==pytest.approx([inner[1,0],inner[1,0]+2],abs=2e-5)
    assert cfg['base_mm']==0 and g['base_mm']==0
    assert all(g['fixture_qc']['gates'].values())


def test_casting_channel_is_preserved_below_extended_upper_stems(real):
    kind,_,name,cfg,_,parts,g=real
    old_cfg=yaml.safe_load((y.ROOT/'config/history/before_raised_support_20261012'/name).read_text())
    old=trimesh.load_mesh(y.ROOT/old_cfg['output']/'assembly'/(kind+'_lumen_core.stl'))
    new=parts[kind+'_lumen_core']
    assert np.allclose(old.bounds[:,:2],new.bounds[:,:2],atol=1e-5)
    assert old.bounds[1,2]==pytest.approx(15.4)
    assert new.bounds[1,2]==pytest.approx(18.4)
    old_lower,new_lower=[trimesh.intersections.slice_mesh_plane(m,[0,0,-1],[0,0,10],cap=True) for m in (old,new)]
    assert old_lower.volume==pytest.approx(new_lower.volume,abs=.003)
    for a,b in [(old_lower,new_lower),(new_lower,old_lower)]:
        idx=np.linspace(0,len(a.faces)-1,300,dtype=int)
        assert trimesh.proximity.closest_point(b,a.triangles_center[idx])[1].max()<.003
    assert new.bounds[0,2]==pytest.approx(1.5,abs=1e-5)
    assert g['pdms_fill_height_mm']==pytest.approx(9.6)
    assert g['pdms_width_qc']['allowance']['budgeted_finished_lower_bound_mm']>=9
    assert g['pdms_height_qc']['allowance']['budgeted_finished_lower_bound_mm']>=9
    assert len(g['section_checks'])==(42 if kind=='Y' else 18)
    assert max(abs(r['equivalent_diameter_mm']-1.6) for r in g['section_checks'])<.003


def test_assembled_height_and_print_height_use_correct_datums(real):
    kind,module,_,cfg,out,parts,g=real
    h=g['manufacturing_height_qc']
    assert h['passed'] and h['enforced']
    assert h['fixture_top_z_mm']==pytest.approx(17.8)
    assert h['retainer_top_z_mm']==pytest.approx(17.8)
    assert h['frame_and_key_top_z_mm']==pytest.approx(10.4)
    assert h['core_top_z_mm']==pytest.approx(18.4)
    assert parts[kind+'_bottomless_box_wall2mm'].bounds[1,2]==pytest.approx(10.4)
    assert min(r['dovetail_engagement_mm'] for r in g['dovetail_joints'])==pytest.approx(8.4)
    meshes,records=module.arrange(parts,out)
    for mesh,row in zip(meshes,records):
        assert mesh.bounds[0,2]==pytest.approx(0)
        assert mesh.bounds[1,2] <= f.part_height_limit(cfg,row['name'])
    # A displaced installed retainer keeps its own height but must FAIL the
    # dish envelope check. This protects against checking extents alone.
    displaced=parts.copy();name=kind+'_dovetail_retainer_INLET'
    displaced[name]=f.shifted(parts[name],[0,0,3.])
    assert not f.casting_height_qc(cfg,displaced)['passed']
    displaced=parts.copy();name=kind+'_bottomless_box_wall2mm'
    displaced[name]=f.shifted(parts[name],[0,0,5.])
    assert not f.casting_height_qc(cfg,displaced)['passed']


def test_every_retainer_clears_entire_box_rim_and_catches_old_touching_arms(real):
    kind,_,_,cfg,_,parts,g=real
    qc=g['fixture_qc']['retainer_rim_clearance_qc']
    assert qc['passed'] and qc['requested_gap_mm']==3.
    for row in qc['parts']:
        assert row['measured_wall_to_arm_gaps_mm']==pytest.approx([3.,3.,3.],abs=1e-5)
        assert row['retainer_in_required_air_gap_mm3']<1e-4
    displaced={n:f.shifted(m,[0,0,-3]) if 'retainer' in n else m for n,m in parts.items()}
    failed=f.retainer_rim_clearance_qc(cfg,displaced,np.asarray(g['cavity_xy_bounds_mm']),
        g['heights']['rim_z_mm'],g['dovetail_joints'],kind)
    assert not failed['passed']
    assert all(r['retainer_in_required_air_gap_mm3']>1 for r in failed['parts'])


def test_raising_support_keeps_existing_frame_and_keys(real):
    kind,_,name,_,_,parts,_=real
    old_cfg=yaml.safe_load((y.ROOT/'config/history/before_raised_support_20261012'/name).read_text())
    for n,m in parts.items():
        if 'box' not in n and 'key' not in n:continue
        old=trimesh.load_mesh(y.ROOT/old_cfg['output']/'assembly'/(n+'.stl'))
        assert m.bounds==pytest.approx(old.bounds,abs=1e-5)
        assert m.volume==pytest.approx(old.volume,abs=.002)
        assert trimesh.proximity.closest_point(old,m.vertices)[1].max()<1e-5


def test_real_rhino_files_preserve_parts_units_and_height(real):
    from vascular_processing import rhino_channel_export as r
    kind,module,_,cfg,out,parts,g=real
    module.arrange(parts,out)
    qc=r.export_models(cfg,parts,g,out)
    assert r.verify_exports(out)['passed']
    for row in qc['exported_files']:
        assert row['archive_version']==70 and row['units']=='mm'
        assert row['part_count']==(8 if kind=='Y' else 6)
        assert row['height_qc']['passed']
        assert all(p['closed'] and p['vertex_max_error_mm']<2e-5 for p in row['parts'])


def test_rhino_rejects_wrong_units_missing_parts_and_oversize_geometry(real,tmp_path):
    import rhino3dm as rhino
    from vascular_processing import rhino_channel_export as r
    _,_,_,cfg,_,parts,_=real
    file=tmp_path/'invalid.3dm';r.write_model(file,parts,'casting_assembly',cfg)
    doc=rhino.File3dm.Read(str(file));doc.Settings.ModelUnitSystem=rhino.UnitSystem.Meters;doc.Write(str(file),7)
    with pytest.raises(ValueError,match='RHINO_MM_UNITS_REQUIRED'):
        r.read_and_check(file,parts,cfg,'casting_assembly')
    r.write_model(file,dict(list(parts.items())[1:]),'casting_assembly',cfg)
    with pytest.raises(ValueError,match='RHINO_PART_SET_MISMATCH'):
        r.read_and_check(file,parts,cfg,'casting_assembly')
    tall={n:f.shifted(m,[0,0,3]) for n,m in parts.items()}
    r.write_model(file,tall,'casting_assembly',cfg)
    with pytest.raises(ValueError,match='RHINO_HEIGHT_LIMIT_FAILED'):
        r.read_and_check(file,tall,cfg,'casting_assembly')


def test_lateral_dovetail_and_cross_key_physically_block_motion(real):
    kind,_,_,_,_,parts,g=real;frame=parts[kind+'_bottomless_box_wall2mm']
    for r in g['dovetail_joints']:
        post=parts[kind+'_dovetail_retainer_'+r['port']];key=parts[kind+'_dovetail_key_'+r['port']]
        for direction in [r['outward_axis'],r['key_insertion_axis']]:
            for sign in [-1,1]:
                assert s.volume_intersection(f.shifted(post,np.array(direction)*sign*.8),frame)>.01
        lifted=f.shifted(post,[0,0,.7])
        assert s.volume_intersection(lifted,key)>.01
        assert s.volume_intersection(lifted,frame)<1e-4
        assert s.volume_intersection(f.shifted(post,[0,0,-.5]),frame)>.01
        assert s.volume_intersection(f.shifted(parts[kind+'_lumen_core'],[0,0,.5]),post)>.01


def test_each_interface_is_assembleable_and_outside_pdms(real):
    kind,_,_,_,_,parts,g=real
    for r in g['fixture_qc']['joints']:
        for k in ['collar_lateral_insertion_intersection_mm3','vertical_rail_insertion_intersection_mm3','key_insertion_intersection_mm3']:
            assert np.isfinite(r[k]).all() and max(r[k])<1e-4
    inner=np.asarray(g['cavity_xy_bounds_mm']);fill=y.box([np.r_[inner[0],0],np.r_[inner[1],9.6]])
    for name,m in parts.items():
        if 'core' not in name:assert s.volume_intersection(m,fill)<1e-4
        assert m.is_volume and len(m.split())==1
    assert not g['physical_fit_validated']


def test_straight_horizontal_end_still_embedded(real):
    kind,_,_,_,_,parts,g=real
    if kind=='Y':
        assert all(p['minimum_cap_normal_dot']>.99999 for p in g['ports']);return
    assert g['design']['horizontal_end']==pytest.approx(23.15)
    assert g['closed_chamber_qc']['horizontal_end_to_inner_wall_mm']==pytest.approx(3.,abs=1e-5)
    assert g['closed_chamber_qc']['horizontal_extension_outside_pdms_mm3']<1e-4
    assert g['ports'][-1]['outward_axis']==[1,0,0]


def test_all_parts_in_native_project_with_correct_supports_and_raft(real):
    kind,module,_,cfg,out,parts,_=real
    meshes,records=module.arrange(parts,out);dest=out/'test.3mf'
    module.write_project(dest,meshes,records,{k:{'name':k} for k in ['machine','process','filament']},cfg)
    objects,_,plates,crc=read_objects(dest)
    assert crc and plates==1 and len(objects)==(8 if kind=='Y' else 6)
    assert {o['metadata']['name'] for o in objects}==set(f.part_names(kind))
    for obj in objects:
        meta=obj['metadata'];name=meta['name']
        assert meta['enable_support']==('1' if 'core' in name or 'retainer' in name else '0')
        assert meta['raft_layers']==('2' if 'core' in name else '0')


def test_internal_base_cannot_reappear_in_bottomless_mode(real,tmp_path):
    _,module,_,cfg,_,_,_=real
    cfg=copy.deepcopy(cfg);cfg['base_mm']=2
    file=tmp_path/'invalid.yaml'
    file.write_text(yaml.safe_dump({k:str(v) if hasattr(v,'read_text') else v for k,v in cfg.items()}))
    with pytest.raises(ValueError,match='BOTTOMLESS_REQUIRES_ZERO_INTERNAL_BASE'):module.load_config(file)
