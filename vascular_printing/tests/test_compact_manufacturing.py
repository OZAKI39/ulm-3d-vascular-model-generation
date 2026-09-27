import copy
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pyvista as pv
import pytest

from vascular_processing.boundary_review_export import HUMAN_SHA, ROOT, read_csv
from vascular_processing.compact_roi import (MODES, add_context, compact_subtree, distances,
    effective_limits, generate_candidate, metrics, radius_metrics, branch_audit)
from vascular_processing.compact_orientation import optimize_compact
from vascular_processing.manufacturing_roi import compensate_radius, profile, export_manufacturing, topology_paths
from vascular_processing.swc_export import read_source, write_swc, validate_tree
from vascular_processing.topbrain_qc import sha256


@pytest.fixture
def cfg():
    return profile(ROOT/'config/compact_manufacturing_profile.yaml')


def tree():
    graph = nx.DiGraph()
    values = [(0,-5,0,0,-1), (1,0,0,0,0), (2,10,0,5,1), (3,20,5,10,2),
        (4,20,-5,10,2), (5,35,10,20,3), (6,30,0,20,3), (7,35,-10,20,4),
        (8,30,0,15,4), (9,48,15,22,5), (10,46,7,27,5)]
    for n,x,y,z,parent in values:
        graph.add_node(n, coords=np.array([x,y,z,.65]), swc_type=3)
        if parent >= 0: graph.add_edge(parent,n)
    return graph


def support(graph):
    return {edge:.95 for edge in graph.edges}


def test_single_component_constraint(cfg):
    raw=tree(); core=raw.subgraph({3,5,6,9,10}).copy()
    out,context=add_context(core,raw,{2,4,7,8},cfg)
    assert not set(out)&{2,4,7,8}
    assert set(out)==set(core)  # Native ancestry may not bridge into another semantic core.
    assert context['length_mm']==0


def test_branch_budget(cfg):
    for mode in MODES:
        graph,_=compact_subtree(tree(),cfg,mode,support(tree()))
        assert len(topology_paths(graph))<=effective_limits(cfg,mode)['max_branches']
    mini,_=compact_subtree(tree(),cfg,'MINI',support(tree()))
    assert metrics(mini)['branch_count']==3


def test_anatomical_bifurcation_depth_budget(cfg):
    cfg['compact_phantom']['rich'].update(max_branches=20,max_bifurcations=10,max_outlets=10,max_path_mm=200)
    cfg['compact_roi'].update(max_topology_branches=20,hard_max_bifurcations=10,hard_max_outlets=10,hard_max_path_length_mm=200)
    graph,audit=compact_subtree(tree(),cfg,'RICH',support(tree()))
    assert graph.out_degree(5)==0 and 9 not in graph and 10 not in graph
    assert audit['stops'][5,9]=='EXCEEDS_BIFURCATION_DEPTH'


def test_outlet_budget_and_daughter_ranking(cfg):
    graph=tree(); graph.add_node(11,coords=np.array([22,0,-15,.9]),swc_type=3);graph.add_edge(2,11)
    weights=support(graph);weights[2,11]=.01
    out,audit=compact_subtree(graph,cfg,'MINI',weights)
    assert 11 not in out and sum(out.out_degree(n)==0 for n in out)==2
    assert audit['stops'][2,11]=='LOW_MANUFACTURING_SCORE'


def test_native_path_clipping_and_adaptive_lengths(cfg):
    g=tree(); cfg['compact_phantom']['balanced']['max_path_mm']=25
    out,_=compact_subtree(g,cfg,'BALANCED',support(g))
    d=distances(g,0)
    assert max(d[n] for n in out)<=25
    assert all(n in g for n in out)
    # Adaptive search may exceed the former fixed path threshold, never the
    # total centerline or branch/size caps, and still uses original samples.
    source=dict(graph=g,component_support=.95)
    candidate=generate_candidate(source,set(g),type('Raw',(),{'graph':g})(),cfg,'BALANCED',support(g),set())
    assert candidate['adaptive_search']['selected_cutoff_mm']>25
    assert candidate['stats']['centerline_length_mm']<=200


def test_true_connector_preservation(cfg):
    g=tree(); core,_=compact_subtree(g,cfg,'BALANCED',support(g))
    root=validate_tree(core)
    for tip in core:
        if core.out_degree(tip)==0:
            assert nx.shortest_path(core,root,tip)==nx.shortest_path(g,root,tip)
    assert set(core.edges)<=set(g.edges)


def test_three_candidates_deterministic(cfg):
    g=tree();source=dict(graph=g,component_support=.95);raw=type('Raw',(),{'graph':g})()
    for mode in MODES:
        first=generate_candidate(source,set(g),raw,cfg,mode,support(g),set())
        second=generate_candidate(source,set(g),raw,cfg,mode,support(g),set())
        assert set(first['graph'].edges)==set(second['graph'].edges)
        assert first['adaptive_search']==second['adaptive_search']
        assert first['gate']['structure_passed']


def test_balanced_default_and_real_single_source(cfg):
    assert cfg['compact_phantom']['default_mode']=='BALANCED'
    base=ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA/compact_manufacturing_roi'
    m=json.loads((base/'compact_manifest.json').read_text())
    assert m['default_print_candidate']=='BALANCED'
    assert len({x['source_component'] for x in m['candidates'].values()})==1
    assert len({x['refined_root_original_id'] for x in m['candidates'].values()})==1
    for item in m['candidates'].values():
        assert item['gate']['structure_passed']
        assert item['stats']['branch_count']<=6 and item['stats']['outlet_count']<=4
        assert all(row['discard_reasons'] for row in item['discarded_eligible_branches'])


def test_radius_compensation_provenance(tmp_path,cfg):
    g=tree();g.nodes[3]['coords'][3]=.31
    src=tmp_path/'source.swc';write_swc(src,g,tmp_path/'raw-source.swc','test','original')
    raw=read_source(src);comp,_=compensate_radius(raw.graph,.6,5)
    export_manufacturing(tmp_path/'comp.swc',comp,raw,set(raw.graph),{e:'b' for e in raw.graph.edges},True)
    rows=read_csv(tmp_path/'comp_mapping.csv')
    for row in rows:
        n=int(row['original_swc_id'])
        assert float(row['original_radius_mm'])==raw.graph.nodes[n]['coords'][3]
        assert float(row['manufacturing_radius_mm'])==comp.nodes[n]['coords'][3]
    info=radius_metrics(raw.graph,comp,cfg)
    assert info['warning']=='LARGE_MANUFACTURING_COMPENSATION'
    assert 0<info['centerline_compensated_fraction']<1


def test_source_geometry_preserved(cfg):
    g=tree(); before=copy.deepcopy(g);source=dict(graph=g,component_support=.95)
    result=generate_candidate(source,set(g),type('Raw',(),{'graph':g})(),cfg,'BALANCED',support(g),set())
    for n in g: np.testing.assert_array_equal(g.nodes[n]['coords'],before.nodes[n]['coords'])
    for n in result['graph']:
        np.testing.assert_array_equal(result['graph'].nodes[n]['coords'],g.nodes[n]['coords'])
        np.testing.assert_array_equal(result['compensated'].nodes[n]['coords'][:3],g.nodes[n]['coords'][:3])


def test_compact_orientation_deterministic(tmp_path,cfg):
    mesh=pv.Sphere(theta_resolution=12,phi_resolution=12).scale([15,20,30],inplace=False).rotate_vector((1,2,3),17,inplace=False)
    graph,_=compact_subtree(tree(),cfg,'BALANCED',support(tree()))
    first=optimize_compact(mesh,graph,cfg,tmp_path/'first')
    second=optimize_compact(mesh,graph,cfg,tmp_path/'second')
    assert 50<=first['tested']<=150 and len(first['top_candidates'])==3
    assert [x['score'] for x in first['all_candidates']]==[x['score'] for x in second['all_candidates']]
    np.testing.assert_array_equal(first['top_candidates'][0]['transform_4x4'],second['top_candidates'][0]['transform_4x4'])


def test_human_visualization_hash_unchanged():
    assert sha256(ROOT/'s1-2_swc_roi_generate_human.py')==HUMAN_SHA
