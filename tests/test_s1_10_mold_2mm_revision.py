"""Actual 2 mm shell, surface-cover and rigid-placement regression checks."""
import copy,json,zipfile
import numpy as np
import pytest,trimesh
from vascular_processing.compact_mold_revision import default_output,checks,revised_layout
from vascular_processing.compact_c_pipeline import default_output as previous_output
from vascular_processing.support_inputs import verify,digest

@pytest.fixture(scope='module')
def out():return default_output()

@pytest.fixture(scope='module')
def summary(out):return json.loads((out/'mold_revision_summary.json').read_text())

@pytest.fixture(scope='module')
def mold(out):return trimesh.load_mesh(out/'outer_mold/ABS_sacrificial_outer_mold_2mm.stl')


def test_actual_mold_base_and_four_walls_are_2mm(mold):
    loc,_,_=mold.ray.intersects_location([[0,0,-1]],[ [0,0,1] ],multiple_hits=True)
    heights=np.unique(np.round(loc[:,2],5));assert heights==pytest.approx([0,2],abs=1e-5)
    for axis in [0,1]:
        start=np.zeros(3);start[axis]=-100;start[2]=10;direction=np.eye(3)[axis]
        pts,_,_=mold.ray.intersects_location([start],[direction],multiple_hits=True);v=np.unique(np.round(pts[:,axis],5));assert len(v)==4;assert [v[1]-v[0],v[3]-v[2]]==pytest.approx([2,2],abs=2e-5)
    assert mold.is_watertight and mold.is_winding_consistent


def test_clearance_measured_from_actual_exported_vascular_surface(out,mold,summary):
    core=trimesh.load_mesh(out/'vascular_core/assembly_placed_ABS_core.stl');assert core.bounds[0,2]-2==pytest.approx(5,abs=1e-5);assert mold.bounds[1,2]-core.bounds[1,2]==pytest.approx(5,abs=1e-5)
    assert summary['dimensions']['outer_dimensions_mm'][2]==pytest.approx(31.01095581)


def test_complete_core_is_only_translated_not_refitted(out):
    transform=json.loads((out/'assembly/casting_core_transform.json').read_text());a=np.asarray(transform['previous_transform_4x4']);b=np.asarray(transform['transform_4x4']);assert b[:3,:3]==pytest.approx(a[:3,:3]);assert b[:3,3]-a[:3,3]==pytest.approx([0,0,-6]);assert transform['scale_factor']==1
    original=trimesh.load_mesh(out/'vascular_core/BG001_RMCA_frozen_ABS_core.stl');original.apply_transform(b);placed=trimesh.load_mesh(out/'vascular_core/assembly_placed_ABS_core.stl');assert original.bounds==pytest.approx(placed.bounds,abs=1e-5);assert original.volume==pytest.approx(placed.volume,abs=.0001)


def test_previous_outputs_preserved(out):
    p=json.loads((out/'source_protection.json').read_text());assert p['mismatch_count']==0;assert p['start_hashes']==p['end_hashes'];assert verify(p['start_hashes'])['all_unchanged'];assert str(previous_output()/'compact_c_retainer_summary.json') in p['start_hashes']


def test_retainer_and_fit_logic_preserved(summary):
    assert summary['cap_shape_unchanged'] and summary['previous_C_fit_preserved'];assert summary['axial_float']['passed'];assert summary['insertion']['passed'];assert summary['anti_float']['passed'];assert summary['contact']['positive_clearance'];assert summary['assembly']['top_insertion_pass'];assert summary['xy_unchanged']


def test_reduced_height_does_not_block_solvent_or_damage_PDMS(summary):
    s=summary['solvent'];t=summary['trim'];assert s['independent_port_routes']==4 and s['all_vascular_branches_reachable'] and s['all_buried_retainers_accessible'];assert abs(t['auxiliary_void_volume_after_trim'])<1e-6;assert abs(t['ghost_missing_volume_mm3'])<1e-6;assert abs(t['ghost_extra_volume_mm3'])<1e-6;assert summary['trim_planes_same_geometric_planes'] and summary['port_lengths_unchanged']


@pytest.mark.parametrize('part',['mold','retainers'])
def test_new_actual_abs_slices_support_off(summary,part):
    r=summary['bambu'][part];assert r['slice_passed'] and r['auto_support_off'] and r['native_support_extrusion_segments']==0;assert r['actual_profile']['filament_type']==['ABS'];assert r['actual_profile']['nozzle_diameter']==['0.2'];assert r['archive_geometry']['passed']
    with zipfile.ZipFile(r['archive']) as z:assert z.testzip() is None
    assert summary['bambu']['mold_features' if part=='mold' else 'retainer_features']['passed']


def test_current_model_and_review_exports(out):
    from vascular_processing.compact_mold_review import NAMES
    assert (out/'assembly/all_ABS_casting_system_C_retainer_2mm.step').stat().st_size>1000
    for name in NAMES:assert (out/'QC'/(name+'.png')).stat().st_size>10000
    for name in ['PDMS_after_all_ABS_removed','PDMS_after_trim']:
        mesh=trimesh.load_mesh(out/f'diagnostic_PDMS/{name}.stl');assert mesh.is_watertight and mesh.is_winding_consistent


def test_missing_slice_or_wrong_base_cannot_be_ready(summary):
    assert checks(summary)[0]=='READY_FOR_HUMAN_REVIEW_WITH_PROCESS_PARAMETERS_UNCALIBRATED'
    bad=copy.deepcopy(summary);bad['dimensions']['bottom_thickness_mm']=3;assert checks(bad)[0]=='NEEDS_ADJUSTMENT'
    bad=copy.deepcopy(summary);bad['bambu']['mold']['slice_passed']=False;assert checks(bad)[0]=='NEEDS_ADJUSTMENT'


def test_invalid_clearance_rejected_before_geometry():
    cfg=dict(wall_thickness_mm=2,bottom_thickness_mm=2,bottom_clearance_mm=-1,top_clearance_mm=5,minimum_fixture_to_rim_mm=1)
    with pytest.raises(ValueError,match='NONPOSITIVE_MOLD_DIMENSION_bottom_clearance_mm'):revised_layout({},cfg,{})
