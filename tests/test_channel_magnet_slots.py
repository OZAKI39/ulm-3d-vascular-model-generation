"""Real supplied Rhino geometry, integrated frames and usable magnet pockets."""
import copy
import json
from pathlib import Path

import numpy as np
import pytest
import trimesh

from vascular_processing import y_channel_casting as y, straight_channel_casting as s
from vascular_processing import bottomless_channel_fixture as f, channel_magnet_slots as m
from vascular_processing.rhino_channel_export import export_models, verify_exports
from vascular_processing.all_abs_print_plate import read_objects

NAMES=['y_channel_p1s_abs.yaml','straight_channel_1p6_p1s_abs.yaml','straight_channel_two_top_ports_1p6_p1s_abs.yaml']


@pytest.fixture(scope='module',params=NAMES)
def real(request,tmp_path_factory):
    name=request.param;kind='Y' if name.startswith('y_') else 'Straight';mod=y if kind=='Y' else s
    cfg=mod.load_config(y.ROOT/'config'/name);out=tmp_path_factory.mktemp('magnets_'+name)
    parts,g=mod.geometry(cfg,out)
    old_cfg=mod.load_config(y.ROOT/'config/history/before_magnet_slots_20261012'/name)
    original={n:trimesh.load_mesh(old_cfg['output']/'assembly'/(n+'.stl')) for n in parts}
    return kind,mod,cfg,out,parts,g,original


def test_source_is_read_from_real_trimmed_brep_not_underlying_surface_bbox(real):
    _,_,cfg,out,_,g,_=real
    slot,meta=m.read_source(cfg)
    assert meta['pocket_diameter_mm']==pytest.approx(6.4,abs=2e-6)
    assert meta['pocket_depth_mm']==pytest.approx(2.,abs=2e-6)
    assert meta['pad_size_mm']==pytest.approx([8.,8.60641266,3.],abs=2e-6)
    assert meta['pocket_floor_mm']==1 and slot.euler_number==2
    assert (out/'reference/magnet_slot.3dm').read_bytes()==(y.ROOT/cfg['magnet_slots']['source']).read_bytes()
    assert meta['source_solid_count']==4 and meta['ignored_construction_object_count']==10
    assert not g['magnet_slots']['source']['source_dimensions_rescaled']


def test_all_original_core_support_key_geometry_and_inner_chamber_preserved(real):
    kind,_,_,_,parts,g,original=real
    for n,a in parts.items():
        b=original[n]
        if 'box' in n:
            assert max(0,b.volume-s.volume_intersection(a,b))<1e-3
            assert a.volume>b.volume+200
            continue
        assert a.bounds==pytest.approx(b.bounds,abs=1e-5)
        assert abs(a.volume-b.volume)<.001
        assert trimesh.proximity.closest_point(b,a.vertices)[1].max()<1e-5
    assert all(g['fixture_qc']['gates'].values()) and all(g['gates'].values())
    assert g['pdms_minimum_floor_mm']==pytest.approx(1.5,abs=1e-5)
    assert g['pdms_fill_height_mm']==pytest.approx(9.6)
    assert g['manufacturing_height_qc']['retainer_top_z_mm']==pytest.approx(17.8)
    assert g['manufacturing_height_qc']['core_top_z_mm']==pytest.approx(18.4)


def test_actual_pockets_are_opposed_open_from_above_and_have_one_mm_floor(real):
    kind,_,_,_,parts,g,_=real;r=g['magnet_slots'];frame=parts[kind+'_bottomless_box_wall2mm']
    assert r['count']==2 and r['qc']['passed']
    centers=np.array([p['pocket_floor_center_mm'] for p in r['placements']])
    assert abs(centers[0,0]-centers[1,0])<3e-6
    assert centers[0,1]<g['cavity_xy_bounds_mm'][0][1]
    assert centers[1,1]>g['cavity_xy_bounds_mm'][1][1]
    for p in r['placements']:
        c=np.array(p['pocket_floor_center_mm'])
        hit=frame.ray.intersects_location([c+[0,0,21]],[[0,0,-1]])[0]
        assert sorted(hit[:,2])==pytest.approx([0,1],abs=2e-5)
        assert not frame.contains([c+[0,0,.2],c+[0,0,1.8]]).any()
        assert frame.contains([c-[0,0,.5]])[0]
        assert np.linalg.det(np.array(p['slot_to_casting_transform'])[:3,:3])==pytest.approx(1)
    assert frame.bounds[1,2]==pytest.approx(10.4,abs=2e-5)


def test_pocket_check_rejects_filled_hole_removed_floor_and_obstructed_access(real):
    kind,_,cfg,_,parts,g,old=real;name=kind+'_bottomless_box_wall2mm'
    report=g['magnet_slots'];c=np.array(report['placements'][0]['pocket_floor_center_mm'])
    others={n:o for n,o in parts.items() if n!=name}
    check=lambda mesh,others: m.check_slots(mesh,old[name],others,np.array(g['cavity_xy_bounds_mm']),cfg,
                                          report['source'],report['placements'])
    assert check(parts[name],others)['passed']
    filled=y.union([parts[name],y.cylinder(c[:2],3.,.8,3.)])
    assert not check(filled,others)['passed']
    no_floor=trimesh.boolean.difference([parts[name],y.cylinder(c[:2],3.,-.1,1.1)],engine='manifold')
    assert not check(no_floor,others)['passed']
    assert not check(parts[name],dict(others,obstacle=y.cylinder(c[:2],2.,5.,6.)))['passed']


def test_magnet_slots_are_integrated_in_native_and_rhino_frame_not_extra_parts(real):
    kind,mod,cfg,out,parts,g,_=real
    meshes,rows=mod.arrange(parts,out);export_models(cfg,parts,g,out)
    assert verify_exports(out)['passed']
    project=out/'test.3mf';mod.write_project(project,meshes,rows,{k:{'name':k} for k in ['machine','process','filament']},cfg)
    objects,_,plates,crc=read_objects(project)
    assert len(objects)==len(parts)==(8 if kind=='Y' else 6) and crc and plates==1
    name=kind+'_bottomless_box_wall2mm';obj=next(o for o in objects if o['metadata']['name']==name)
    assert obj['metadata']['enable_support']=='0'
    transform=np.array(next(r for r in rows if r['name']==name)['casting_to_plate_transform'])
    for p in g['magnet_slots']['placements']:
        c=(transform@np.r_[p['pocket_floor_center_mm'],1])[:3]
        assert not obj['mesh'].contains([c+[0,0,1]])[0]
        assert obj['mesh'].contains([c-[0,0,.5]])[0]


def test_source_hash_and_impossible_clearance_fail_instead_of_silently_changing_slot(real):
    _,_,cfg,out,parts,g,old=real
    bad=copy.deepcopy(cfg);bad['magnet_slots']['source_sha256']='0'*64
    with pytest.raises(ValueError,match='MAGNET_SOURCE_SHA256_MISMATCH'):m.read_source(bad)
    name=next(n for n in parts if 'box' in n)
    obstacle=y.box([[-60,-60,0],[60,60,22]])
    with pytest.raises(ValueError,match='NO_CLEAR_OPPOSING_MAGNET_SLOT_POSITION'):
        m.add_slots(cfg,old[name],{'obstacle':obstacle},np.array(g['cavity_xy_bounds_mm']),out)
