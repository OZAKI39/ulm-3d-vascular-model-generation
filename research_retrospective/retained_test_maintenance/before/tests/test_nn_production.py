"""Core support rejection and exact BraVa subset invariants; no matching benchmark."""
from pathlib import Path
from types import SimpleNamespace
import numpy as np
import pytest
from vascular_processing.semantic_ensemble import nearest_support
from vascular_processing.nn_support_gate import apply_gate
from vascular_processing.brava_branch_labels import (load_exact_graph,topology_branches,aggregate_branches,
    color_source_mapping)
from vascular_processing.brava_mevo_roi import strict_components,export_component,extract_and_model
from vascular_processing.swc_export import read_source
from vascular_processing.topbrain_qc import sha256

ROOT=Path(__file__).resolve().parents[1]


def donor(label,case,scale=1.):
    sample=SimpleNamespace(geometry=SimpleNamespace(points=np.array([[0.,0,0]]),case_id=case),labels=np.array([label]))
    transform=np.diag([scale,scale,scale,1.])
    return sample,dict(status='VALID',transform=transform,inlier_rmse=scale)


def test_binary_votes_merge_native_m2_m3_before_voting():
    f=nearest_support([donor(1,'a'),donor(2,'b'),donor(3,'c')],np.array([[0.,0,0]]))
    assert f['per_donor_native_vote'].ravel().tolist()==[1,2,3]
    assert f['winner'].tolist()==[2]
    assert f['vote_MeVO'].tolist()==[2]
    assert f['p_MeVO'][0]==pytest.approx(2/3)
    assert f['valid_donor_count'].tolist()==[3]
    assert f['vote_margin'][0]==pytest.approx(1/3)


def test_agreement_and_valid_donor_minimum_reject_without_relabeling():
    features=dict(winner=np.array([2,2,1]),agreement=np.array([.59,.60,1.]),
        median_normalized_distance=np.zeros(3),valid_donor_count=np.array([5,5,2]))
    gate=dict(agreement_threshold=.6,distance_threshold=3.,minimum_donors=3)
    assert apply_gate(features,gate).tolist()==[0,2,0]
    assert features['winner'].tolist()==[2,2,1]


def test_distance_gate_is_normalized_by_registration_rmse():
    points=np.array([[.2,0,0],[10.,0,0]])
    a=nearest_support([donor(2,str(i)) for i in range(3)],points)
    b=nearest_support([donor(2,str(i),10.) for i in range(3)],points*10)
    np.testing.assert_allclose(a['median_normalized_distance'],b['median_normalized_distance'])
    gate=dict(agreement_threshold=.6,distance_threshold=2.,minimum_donors=3)
    assert apply_gate(a,gate).tolist()==[2,0]
    assert apply_gate(b,gate).tolist()==[2,0]


@pytest.fixture
def tree(tmp_path):
    path=tmp_path/'tree.swc'
    path.write_text('1 8 0 0 0 1 -1\n2 8 5 0 0 .7 1\n3 8 10 0 0 .6 2\n'
        '5 9 20 0 0 .5 3\n6 9 25 0 0 .4 5\n7 9 30 0 0 .3 6\n'
        '8 9 20 5 0 .3 5\n9 8 10 5 0 .4 3\n')
    exact=load_exact_graph(path);branches=topology_branches(exact)
    return exact,branches


def test_branch_aggregation_uses_arc_length_and_known_length(tree):
    exact,branches=tree;branch=next(b for b in branches.values() if b['node_ids']==[1,2,3])
    ids=np.array([1,2,3]);features=dict(p_M1=np.array([0.,0.,1.]),p_MeVO=np.array([1.,1.,0.]),
        median_normalized_distance=np.ones(3),valid_donor_count=np.full(3,5))
    rows=aggregate_branches(exact,{branch['branch_id']:branch},ids,features,np.array([2,2,1]))
    assert rows[0]['mean_p_MeVO']==.75 and rows[0]['label']=='MeVO'
    rows=aggregate_branches(exact,{branch['branch_id']:branch},ids,features,np.array([2,0,0]))
    assert rows[0]['known_fraction']==.25 and rows[0]['label']=='UNKNOWN'


def choose_labels(branches,by_start_end):
    return [dict(branch_id=k,label=by_start_end.get((b['node_ids'][0],b['node_ids'][-1]),'UNKNOWN')) for k,b in branches.items()]


def test_unknown_gap_remains_disconnected_even_when_endpoints_are_in_roi(tree):
    exact,branches=tree;rows=choose_labels(branches,{(1,3):'MeVO',(5,7):'MeVO'})
    parts=strict_components(exact,branches,rows)
    assert len(parts)==2
    assert any(3 in g for g in parts) and any(5 in g for g in parts)
    assert all(not g.has_edge(3,5) for g in parts)
    assert exact.graph.has_edge(3,5)


def test_roi_geometry_and_type_are_exact_after_export(tree,tmp_path):
    exact,branches=tree;before=sha256(exact.source.path)
    rows=choose_labels(branches,{(1,3):'MeVO'})
    graph,=strict_components(exact,branches,rows)
    result=export_component(tmp_path/'roi.swc',graph,exact,branches,set(graph))
    saved=read_source(Path(result['path']))
    for mapping in result['node_mapping']:
        a=exact.graph.nodes[mapping['original_swc_id']];b=saved.graph.nodes[mapping['roi_node_id']]
        np.testing.assert_array_equal(a['coords'],b['coords'])
        assert a['swc_type']==b['swc_type']
    assert result['directed_edges_preserved'] and sha256(exact.source.path)==before


def test_modelable_roundtrip_adds_only_one_real_m1_context_branch(tree,tmp_path):
    exact,branches=tree;rows=choose_labels(branches,{(1,3):'M1',(3,5):'MeVO',(5,7):'MeVO'})
    out=tmp_path/'export';out.mkdir()
    results=extract_and_model(exact,branches,rows,out,model=False)
    assert len(results)==1
    item=results[0]
    assert item['strict']['node_count']==4 and item['modelable']['node_count']==6
    assert len(item['context_branch_ids'])==1
    assert item['strict']['root_original_id']==3 and item['modelable']['root_original_id']==1
    assert item['modelable']['geometry_exact_preserved']
    assert 'True,False' in (out/'node_mapping.csv').read_text()


def test_real_version_mapping_and_human_hash_protected():
    raw=ROOT/'vessel_model/T - Brava/swc_files/BG001.CNG.swc'
    color=raw.with_name('BG001_ColorCoded.CNG.swc')
    expected=(sha256(raw),sha256(color));exact=load_exact_graph(raw)
    labels,mapping,audit=color_source_mapping(exact,color)
    assert len(mapping)==len(set(mapping.values()))==2810
    assert audit['undirected_topology_identical']
    assert set(labels.values())=={0,2,3,4,5,6,7}
    assert (sha256(raw),sha256(color))==expected
    assert sha256(ROOT/'s1-2_swc_roi_generate_human.py')=='3361cea44a97e8bfa907ff1504b0f977e329d35eb59cb9f57cc1b3a2a55cd729'
