import copy
import csv
import json
from pathlib import Path
import zipfile

import networkx as nx
import numpy as np
import pyvista as pv
import pytest
from scipy.spatial.transform import Rotation

from vascular_processing.boundary_review_export import HUMAN_SHA, ROOT
from vascular_processing.manufacturing_roi import (
    profile,compensate_radius,export_manufacturing,minimal_subtree,extend_proximal,
    prune_twigs,geometry_stats,roi_size_gate,load_semantic)
from vascular_processing.print_orientation import (
    cap_native_ports,transform_points,bed_transform,fits_volume,score_orientation,
    orientation_rotations,triangle_geometry,branch_descriptors)
from vascular_processing.bambu_manufacturing import parse_sliced_archive
from vascular_processing.swc_export import Source,read_source
from vascular_processing.topbrain_qc import sha256


@pytest.fixture
def config():
    return profile(ROOT/'config/bambu_manufacturing_profile.yaml')


def line(radii,spacing=1.):
    graph=nx.DiGraph()
    for i,r in enumerate(radii):
        graph.add_node(i,coords=np.array([i*spacing,0.,0.,r]),swc_type=3)
        if i:graph.add_edge(i-1,i)
    return graph


def branching():
    g=line([.7,.3,.3],spacing=5.)
    for n,point,parent in [(3,[15,0,0,.3],2),(4,[20,0,0,.3],3),
                           (5,[10,4,0,.3],2),(6,[5,15,0,.3],1),
                           (7,[5,-15,0,.7],1)]:
        g.add_node(n,coords=np.array(point,float),swc_type=3);g.add_edge(parent,n)
    return g


def test_nozzle_thresholds_and_detail_profile(config):
    assert config['effective']['hard_min_diameter_mm']==pytest.approx(1.)
    assert config['effective']['preferred_min_diameter_mm']==pytest.approx(1.2)
    detail=profile(ROOT/'config/bambu_abs_0p2_detail.yaml')
    assert detail['effective']['hard_min_diameter_mm']==pytest.approx(.6)
    assert detail['effective']['preferred_min_diameter_mm']==pytest.approx(.8)
    assert detail['effective']['radius_floor_mm']==pytest.approx(.4)


def test_floor_changes_only_radius():
    raw=line([.31,.31,.62,.93,.31,.31])
    before=copy.deepcopy(raw)
    out,_=compensate_radius(raw,.6,5)
    assert set(out.edges)==set(raw.edges)
    for n in raw:
        np.testing.assert_array_equal(out.nodes[n]['coords'][:3],raw.nodes[n]['coords'][:3])
        np.testing.assert_array_equal(before.nodes[n]['coords'],raw.nodes[n]['coords'])
        assert out.nodes[n]['coords'][3]>=max(.6,raw.nodes[n]['coords'][3])
        assert out.nodes[n]['swc_type']==raw.nodes[n]['swc_type']


def test_mapping_preserves_original_radius_and_print_header(tmp_path):
    raw=line([.31,.31,.62,.93,.31,.31])
    source_path=tmp_path/'source.swc'
    source_path.write_text(''.join(f'{n} 3 {n} 0 0 {d["coords"][3]:.17g} {n-1}\n' for n,d in raw.nodes(data=True)))
    source=read_source(source_path);digest=sha256(source_path)
    out,_=compensate_radius(raw,.6,5);path=tmp_path/'print.swc'
    info=export_manufacturing(path,out,source,{2,3,4,5},{e:'b' for e in raw.edges},True)
    rows=list(csv.DictReader(path.with_name('print_mapping.csv').open()))
    assert info['round_trip'] and sha256(source_path)==digest
    assert 'MANUFACTURING RADIUS COMPENSATION APPLIED' in path.read_text()
    for r in rows:
        n=int(r['original_swc_id'])
        assert float(r['original_radius_mm'])==raw.nodes[n]['coords'][3]
        assert float(r['diameter_manufacturing_mm'])==2*float(r['manufacturing_radius_mm'])
        assert float(r['radius_added_mm'])==pytest.approx(float(r['manufacturing_radius_mm'])-float(r['original_radius_mm']))


def test_local_blend_reduces_floor_activation_step():
    raw=line([.4]*8+[.8]*8)
    out,windows=compensate_radius(raw,.6,5.)
    radius=np.array([out.nodes[n]['coords'][3] for n in out])
    assert np.all(np.diff(radius)>=-1e-12)
    assert np.max(np.diff(radius))<.2
    assert radius[0]==.6 and radius[-1]==.8 and windows
    assert set(out)==set(raw)


def test_short_terminal_twig_removed_but_long_thin_branch_retained(config):
    graph=branching();kept,removed=prune_twigs(graph,config)
    assert 5 not in kept and 6 in kept
    assert len(removed)==2  # 2->3->4 is 10 mm; 1->6 is 15 mm and stays.
    assert any(r['original_node_ids']==[5] for r in removed)


def test_bifurcation_connector_not_recursively_removed(config):
    graph=branching();kept,_=prune_twigs(graph,config)
    assert 1 in kept and 2 in kept and kept.has_edge(1,2)
    assert nx.is_arborescence(kept)
    assert kept.has_edge(1,6) and kept.has_edge(1,7)


def test_minimal_connecting_subtree_uses_native_edges():
    graph=branching();core={4,5,6}
    connected=minimal_subtree(graph,core)
    assert set(connected)=={1,2,3,4,5,6}
    assert set(connected.edges)<=set(graph.edges) and nx.is_arborescence(connected)
    assert 0 not in connected and 7 not in connected


def test_proximal_context_native_only_and_length_bounded(config):
    graph=line([.7]*20,spacing=2.)
    roi=copy.deepcopy(config['roi']);roi['proximal_context_max_mm']=7
    extended,report=extend_proximal(graph.subgraph(range(10,20)).copy(),graph,roi)
    assert report['extended_length_mm']==6 and report['final_root']==7
    assert set(extended.edges)<=set(graph.edges)
    for n in extended:np.testing.assert_array_equal(extended.nodes[n]['coords'],graph.nodes[n]['coords'])


def test_roi_size_gate_checks_geometry_and_ports(config):
    stats=geometry_stats(branching())
    assert not roi_size_gate(stats,config)['passed']
    stats['bbox'].update(longest_dimension_mm=110.,diagonal_mm=140.)
    assert roi_size_gate(stats,config)['passed']
    stats['inlet_count']=2
    assert not roi_size_gate(stats,config)['passed']


def test_build_volume_rejection_does_not_scale(config):
    assert fits_volume([226,226,240],[226,226,240])
    assert not fits_volume([226.1,226,240],[226,226,240])
    mesh=pv.Cube(x_length=400,y_length=30,z_length=40).triangulate()
    areas,normals,centers=triangle_geometry(mesh)
    result=score_orientation(mesh.points,areas,normals,centers,np.eye(3),config,[])
    assert not result['fits_build_volume'] and result['status']=='ROI_EXCEEDS_BUILD_VOLUME'
    assert result['scale']==1.0 and result['bbox_extents_mm'][0]==400.


def test_rotation_transform_is_rigid_and_invertible():
    points=np.array([[0,0,0],[1,0,0],[0,2,0],[0,0,3]],float)
    rotation=Rotation.from_euler('xyz',[15,-30,45],degrees=True).as_matrix()
    matrix,extent,low,high=bed_transform(points,rotation,[256,256,256],15)
    new=transform_points(points,matrix)
    np.testing.assert_allclose(transform_points(new,np.linalg.inv(matrix)),points,atol=1e-12)
    np.testing.assert_allclose(matrix[:3,:3].T@matrix[:3,:3],np.eye(3),atol=1e-12)
    assert low[2]==pytest.approx(0) and np.linalg.det(matrix[:3,:3])==pytest.approx(1.)


def test_orientation_score_and_search_are_deterministic(config):
    mesh=pv.Sphere(theta_resolution=12,phi_resolution=12).triangulate()
    mesh.points[:,0]*=3;mesh.points[:,1]*=2
    rotations=orientation_rotations(mesh.points,config)
    assert 100<=len(rotations)<=300
    areas,normals,centers=triangle_geometry(mesh)
    args=(mesh.points,areas,normals,centers,rotations[10][1],config,branch_descriptors(branching()))
    assert score_orientation(*args)==score_orientation(*args)
    for first,second in zip(rotations,orientation_rotations(mesh.points,config)):
        np.testing.assert_array_equal(first[1],second[1])


def test_all_frozen_sources_and_semantic_layers_unchanged():
    _,_,_,_,snapshot=load_semantic(ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/refined_roi')
    assert len(snapshot)>1000
    assert all(sha256(path)==digest for path,digest in snapshot.items())


def test_human_s1_2_unchanged():
    assert sha256(ROOT/'s1-2_swc_roi_generate_human.py')==HUMAN_SHA


def test_print_port_caps_close_mesh_without_moving_sidewalls():
    native=pv.Cylinder(capping=False,resolution=16).triangulate()
    original=native.points.copy();closed,qc=cap_native_ports(native)
    assert qc['cap_count']==2 and closed.n_open_edges==0 and qc['nonmanifold_edges']==0
    np.testing.assert_array_equal(native.points,original)


def test_bambu_export_without_gcode_is_not_slice_success(tmp_path):
    path=tmp_path/'unsliced.3mf'
    with zipfile.ZipFile(path,'w') as archive:
        archive.writestr('Metadata/slice_info.config','<config><header/></config>')
    assert not parse_sliced_archive(path)['slice_success']
    with zipfile.ZipFile(path,'w') as archive:
        archive.writestr('Metadata/slice_info.config','<config><plate><metadata key="prediction" value="120"/><filament used_g="1.25"/></plate></config>')
        archive.writestr('Metadata/plate_1.gcode','; FEATURE: Support\nG1 X1 Y2 E0.03\n')
    report=parse_sliced_archive(path)
    assert report['slice_success'] and report['estimated_print_time_seconds']==120
    assert report['filament_used_g']==1.25 and report['support_information']['filament_usage_g'] is None
