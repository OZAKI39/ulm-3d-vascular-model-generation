"""Real geometry checks: continuous chamber, internal end and later PDMS cut."""
import copy
from pathlib import Path
import numpy as np
import pytest
import trimesh
import yaml
from scipy.spatial import cKDTree
from vascular_processing import straight_channel_casting as s


@pytest.fixture(scope='module')
def real(tmp_path_factory):
    cfg = s.load_config(s.ROOT/'config/history/before_bottomless_dovetail_20261012/straight_channel_1p6_p1s_abs.yaml')
    out = tmp_path_factory.mktemp('closed_straight')
    parts, report = s.geometry(cfg, out)
    return cfg, parts, report, out


def test_four_parts_and_no_gate_or_hole(real):
    _, parts, report, out = real
    assert len(parts) == 4
    assert not any('gate' in k for k in parts)
    assert not list((out/'assembly').glob('*gate*'))
    assert not report['side_gate_external_seal_required']
    assert all(report['gates'].values())
    for m in parts.values():
        assert m.is_volume and len(m.split()) == 1
        assert np.all(np.bincount(m.edges_unique_inverse) == 2)


def test_actual_end_wall_is_continuous_2mm_through_old_hole_region(real):
    _, parts, report, _ = real
    d = report['design'];mold=parts['Straight_casting_box_base_2mm']
    # Probe the whole former hole/slot height, rather than accepting STL closure.
    for y in np.linspace(-3.5, 3.5, 7):
        for z in np.linspace(2.2, 13.2, 12):
            hits = mold.ray.intersects_location([[d['horizontal_end'], y, z]], [[1,0,0]])[0]
            assert len(hits) == 2
            assert np.sort(hits[:,0]) == pytest.approx([26.15,28.15], abs=2e-5)
    assert max(report['closed_chamber_qc']['wall_and_bottom_missing_volume_mm3'].values()) < 1e-3


def test_horizontal_end_is_inside_pdms_with_3mm_end_cover(real):
    cfg, parts, report, _ = real
    d = report['design'];core=parts['Straight_lumen_core']
    assert core.bounds[1,0] == pytest.approx(23.15,abs=1e-5)
    assert d['inner'][1,0]-core.bounds[1,0] == pytest.approx(3.,abs=1e-5)
    rows=report['closed_chamber_qc']
    assert rows['horizontal_extension_outside_pdms_mm3'] < 1e-4
    assert min(*rows['extension_clearance_to_pdms_bounds_mm']['lower_xyz'],
               *rows['extension_clearance_to_pdms_bounds_mm']['upper_xyz']) >= 1.5-1e-5


def test_low_channel_unchanged_only_upper_stems_are_lengthened(real):
    cfg, parts, _, _ = real
    legacy = s.load_config(s.ROOT/'config/history/before_closed_chamber_20261012/straight_channel_1p6_p1s_abs.yaml')
    old_cad, _ = s.build_core(legacy, s.design(legacy))
    old=s.cad_mesh(old_cad,legacy);new=parts['Straight_lumen_core']
    # Raising fill/rim necessarily lengthens the upper stems and moves collars.
    # Verify the entire low channel using independent volumetric intersections.
    clip=s.box([[-30,-5,0],[30,5,12]])
    a=trimesh.boolean.intersection([old,clip],engine='manifold')
    b=trimesh.boolean.intersection([new,clip],engine='manifold')
    common=s.volume_intersection(a,b)
    assert a.volume+b.volume-2*common < 1e-3
    assert new.bounds[1,2]-old.bounds[1,2] == pytest.approx(.6,abs=1e-5)


def test_manual_cut_opens_1p6_channel_away_from_tee(real):
    _, _, report, out=real
    trim=report['closed_chamber_qc']['manual_cut_reference']
    assert trim['plane_x_mm'] == pytest.approx(22.65)
    assert trim['pdms_to_remove_from_end_face_mm'] == pytest.approx(3.5)
    m=trimesh.load_mesh(out/'assembly/Straight_channel_without_collars_reference.stl')
    assert s.section_diameter(m,[trim['plane_x_mm'],0,4.3],[1,0,0]) == pytest.approx(1.6,abs=.003)


def test_dimensions_floor_fill_wall_and_vertical_ports(real):
    cfg, parts, report, _=real
    core,mold=parts['Straight_lumen_core'],parts['Straight_casting_box_base_2mm']
    assert np.allclose(mold.extents,[50.95,13.6,14.0],atol=2e-5)
    floor=mold.ray.intersects_location([[0,0,1]],[[0,0,1]])[0][:,2].min()
    assert floor == pytest.approx(2)
    assert core.bounds[0,2]-floor == pytest.approx(1.5,abs=.003)
    assert mold.bounds[1,2]-cfg['freeboard_mm']-floor == pytest.approx(9.6,abs=2e-5)
    assert [p['outward_axis'] for p in report['ports']] == [[0,0,1],[0,0,1],[1,0,0]]
    for p in report['ports']: assert p['measured_cap_diameter_mm'] == pytest.approx(1.6,abs=.003)


def test_actual_pdms_width_is_9p6_excluding_abs_walls(real):
    _,parts,report,_=real
    mold=parts['Straight_casting_box_base_2mm']
    for z in [3.,6.,10.5]:
        for sign in [-1,1]:
            hits=mold.ray.intersects_location([[0,0,z]],[[0,sign,0]])[0]
            assert sorted(abs(hits[:,1])) == pytest.approx([4.8,6.8],abs=1e-5)
    assert report['pdms_width_qc']['measured_inner_width_mm'] == pytest.approx(9.6,abs=1e-5)
    assert report['gates']['pdms_width_minimum']


def test_widened_outlet_retainer_straddles_actual_wall(real):
    _,parts,_,_=real
    m=parts['Straight_retainer_OUTLET1']
    # Cross its two locating feet below the horizontal support plate.
    hits=m.ray.intersects_location([[12,0,13]],[[0,1,0]])[0]
    assert sorted(hits[:,1]) == pytest.approx([3.35,4.55,7.05,8.25],abs=2e-5)


def test_core_assembly_and_fixtures_clear(real):
    _, parts, report, _=real
    for row in report['closed_chamber_qc']['core_insertion_path']:
        assert row['core_mold_intersection_mm3'] < 1e-4
    assert max(report['assembly_intersections_mm3'].values()) < 1e-4
    for name,m in parts.items():
        if 'retainer' in name: assert m.bounds[0,2]-report['design']['fill'] > .39


def test_adaptive_end_wall_does_not_change_core_endpoint(real):
    cfg, _, report, _=real
    alt=copy.deepcopy(cfg);alt['horizontal_end_clearance_mm']=4.5
    d=s.design(alt)
    assert d['horizontal_end'] == pytest.approx(report['design']['horizontal_end'])
    assert d['inner'][1,0] == pytest.approx(27.65)


def test_four_independent_native_objects_and_core_only_support(real):
    from vascular_processing.all_abs_print_plate import read_objects
    cfg,parts,_,out=real
    meshes,records=s.arrange(parts,out)
    p=out/'test.3mf'
    s.write_project(p,meshes,records,{k:{'name':k} for k in ['machine','process','filament']},cfg)
    objects,_,plates,crc=read_objects(p)
    assert crc and plates == 1 and len(objects)==4
    assert {o['metadata']['name'] for o in objects} == set(parts)
    enabled=[o for o in objects if o['metadata']['enable_support']=='1']
    assert len(enabled)==1 and enabled[0]['metadata']['name']=='Straight_lumen_core'


def test_sidegate_geometry_cannot_be_loaded_as_closed_chamber(real,tmp_path):
    cfg,_,_,out=real
    cfg=copy.deepcopy(cfg);cfg['output']=out;cfg['chamber_mode']='side_gate'
    old=s.load_config(s.ROOT/'config/history/before_closed_chamber_20261012/straight_channel_1p6_p1s_abs.yaml')
    for k,v in old.items():
        if k.startswith('side_gate_') or k=='horizontal_projection_beyond_gate_mm':cfg[k]=v
    p=tmp_path/'wrong.yaml'
    p.write_text(yaml.safe_dump({k:str(v) if isinstance(v,Path) else v for k,v in cfg.items()}))
    with pytest.raises(ValueError,match='CONFIG_CHANGED_REGENERATE_GEOMETRY'):s.run(p,'slice')
