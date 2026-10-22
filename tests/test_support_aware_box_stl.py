import json
import numpy as np
import pytest
import trimesh

from vascular_processing.support_aware_box_stl import resized_box,measured_thickness
from vascular_processing.support_aware_port_design import configuration
from vascular_processing.support_aware_roi import DEFAULT_OUTPUT,sha256


@pytest.mark.parametrize('thickness',[2.,3.])
def test_exported_box_has_requested_thickness_and_preserved_cavity(tmp_path,thickness):
    cfg=configuration();inner=np.array([[89.,98.,10.],[167.,169.,39.]])
    original=inner.copy();box=resized_box(inner,thickness,cfg)
    file=tmp_path/'box.stl';box['mesh'].export(file);mesh=trimesh.load_mesh(file)
    assert np.array_equal(inner,original) and np.array_equal(box['inner'],original)
    assert all(v==pytest.approx(thickness,abs=2e-5) for v in measured_thickness(mesh).values())
    assert mesh.is_watertight and mesh.is_volume
    expected=np.prod([78+2*thickness,71+2*thickness,29+thickness])-78*71*29
    assert mesh.volume==pytest.approx(expected,rel=1e-6)


@pytest.mark.parametrize('thickness',[0.,-1.,float('nan'),float('inf')])
def test_invalid_wall_thickness_rejected(thickness):
    with pytest.raises(ValueError,match='THICKNESS'):
        resized_box([[0,0,0],[20,20,20]],thickness,configuration())


def test_real_2mm_stl_preserves_selected_core_and_connections():
    out=DEFAULT_OUTPUT/'casting_mold_2mm';path=out/'geometry_summary.json'
    if not path.exists():pytest.skip('Run s1-7 --stage box-stl --box-thickness-mm 2 first')
    result=json.loads(path.read_text());source=DEFAULT_OUTPUT/'candidates'/result['selected_candidate']/'manufacturing'
    assert result['status']=='GEOMETRY_STL_CHECKS_PASSED' and not result['failures']
    assert result['source_core_sha256']==sha256(source/'core_with_ports.stl')
    assert result['union_qc']['core_material_loss_mm3']<=.001
    assert result['source_protection']['all_unchanged']
    assert all(r['status']=='PASS' and r['intersection_axial_span_mm']>=1.8 for r in result['port_overlaps'])
    assert all(r['passed'] for r in result['interfaces'])
    mesh=trimesh.load_mesh(out/'intact_open_top_box_2mm.stl')
    assert all(v==pytest.approx(2.,abs=2e-5) for v in measured_thickness(mesh).values())
    native=trimesh.load_mesh(out/'final_abs_casting_mold_2mm.stl')
    placed=trimesh.load_mesh(out/'final_abs_casting_mold_2mm_oriented.stl')
    assert native.is_volume and placed.is_volume and len(native.split())==1
    assert placed.bounds[0,2]==pytest.approx(0.,abs=1e-6)
    assert placed.volume==pytest.approx(native.volume,rel=1e-5)
    assert result['slicing_status']=='NOT_RESLICED_FOR_NEW_THICKNESS'
