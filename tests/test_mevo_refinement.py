from pathlib import Path

import networkx as nx
import numpy as np
import pytest

from vascular_processing.boundary_review_export import HUMAN_SHA, ROOT, load_cache
from vascular_processing.mevo_refinement import (
    MIN_DIAMETER_MM, MIN_RADIUS_MM, arc, branch_edges, distal_candidates, export_subset,
    frozen_upstream, proximal_start, prune_diameter, refine_component)
from vascular_processing.swc_export import Source, read_source, validate_tree
from vascular_processing.topbrain_qc import sha256


def line(radii, spacing=1.):
    graph = nx.DiGraph()
    for i,r in enumerate(radii):
        graph.add_node(i,coords=np.array([i*spacing,0.,0.,r]),swc_type=3)
        if i:
            graph.add_edge(i-1,i)
    coordinates = {n: v['coords'] for n,v in graph.nodes(data=True)}
    return graph, coordinates, {e:'branch' for e in graph.edges}


def prediction(probabilities):
    return {i: dict(p_M1=1-p,p_MeVO=p,valid_donor_count=5,label='MeVO' if p>=.6 else 'M1',
                   agreement=max(p,1-p), normalized_support_distance=.8)
            for i,p in enumerate(probabilities)}


def test_diameter_conversion_and_threshold():
    assert MIN_DIAMETER_MM == .75 and MIN_RADIUS_MM == .375
    assert MIN_DIAMETER_MM == 2 * MIN_RADIUS_MM
    g,c,e = line([.375]*9)
    kept,rows,_ = prune_diameter(g,c,e,1,'test')
    assert set(kept)==set(g) and not rows
    g,c,e = line([.3749]+[.5]*8)
    kept,rows,_ = prune_diameter(g,c,e,1,'test')
    assert not kept and rows[0]['reason']=='COMPONENT_BELOW_DIAMETER_THRESHOLD'
    assert rows[0]['diameter_mm']==.7498


def test_sustained_proximal_uses_physical_length_and_reversion():
    g,c,_ = line([.5]*20)
    p = prediction([.2]*3+[.8]*17)
    cut,_ = proximal_start(list(g),c,p)
    assert cut['node']==3 and cut['window_length_mm']==5
    p = prediction([.2]*3+[.8]*7+[.2]*5+[.8]*5)
    cut,audit = proximal_start(list(g),c,p)
    assert cut['node']==15 and cut['window_length_mm']==4
    assert any(a['reversion_detected'] for a in audit)
    # A low-support UNKNOWN label must not hide finite M1 reversion evidence.
    for i in range(10,15):
        p[i]['label']='UNKNOWN'
        p[i]['valid_donor_count']=2
    cut,_ = proximal_start(list(g),c,p)
    assert cut['node']==15


def test_short_spike_and_insufficient_samples_rejected():
    g,c,_ = line([.5]*10)
    cut,_ = proximal_start(list(g),c,prediction([.2]*3+[.8]+[.2]*6))
    assert cut['node'] is None
    cut,_ = proximal_start([8,9],c,prediction([.8]*10))
    assert cut['status']=='PROXIMAL_REFINEMENT_INSUFFICIENT_SAMPLES'
    p = prediction([.8]*10)
    for i in range(8):
        p[i]['label']='UNKNOWN'
    cut,_ = proximal_start(list(g),c,p)
    assert cut['node'] is None


def test_branch_interior_cut_preserves_only_suffix(tmp_path):
    g,c,_ = line([.5]*12,spacing=1.7)
    cut,_ = proximal_start(list(g),c,prediction([.2]*4+[.8]*8))
    assert cut['node']==4
    suffix=g.subgraph({4}|nx.descendants(g,4)).copy()
    assert validate_tree(suffix)==4 and set(suffix)==set(range(4,12))
    assert arc(list(suffix),c)[-1]==pytest.approx(11.9)


def test_sustained_distal_cut_keeps_endpoint():
    g,c,e = line([.5]*4+[.3]*9)
    kept,rows,_ = prune_diameter(g,c,e,1,'test')
    assert set(kept)==set(range(5)) and rows[0]['cut_original_swc_id']==4
    assert rows[0]['removed_downstream_nodes']==8
    assert rows[0]['window_span_mm']==3.


def test_single_dip_not_pruned_even_with_sparse_sampling():
    g,c,e = line([.5]*4+[.3]+[.5]*8,spacing=4.)
    kept,rows,_ = prune_diameter(g,c,e,1,'test')
    assert set(kept)==set(g) and not rows


def test_rebound_veto_in_physical_five_mm_window():
    # Dense low samples, then a long high interval: count alone would prune.
    positions=np.r_[0.,np.linspace(1.,2.9,20),3.,4.,5.,6.,7.]
    g,c,e=line([.5]+[.3]*20+[.5]*5)
    for n,x in enumerate(positions):
        c[n][0]=x
    candidates=distal_candidates(list(g),c)
    assert candidates[0]['node']==1 and candidates[0]['status']=='DIAMETER_REBOUND'
    kept,rows,_=prune_diameter(g,c,e,1,'test')
    assert set(kept)==set(g)
    assert rows and all(r['reason']=='DIAMETER_REBOUND' for r in rows)


def test_cut_removes_entire_descendant_subtree():
    g,c,e=line([.5]*4+[.3]*6)
    g.add_node(10,coords=np.array([8.,2.,0.,.8]),swc_type=3)
    g.add_edge(7,10)
    c[10]=g.nodes[10]['coords'];e[7,10]='daughter'
    kept,rows,_=prune_diameter(g,c,e,1,'test')
    assert set(kept)==set(range(5))
    assert 10 not in kept and nx.is_arborescence(kept)
    assert rows[0]['removed_downstream_branches']==2


def test_exact_geometry_type_and_swc_roundtrip(tmp_path):
    graph,coords,edges=line([.52345678912345,.62,.31,.44,.93,.31])
    graph.nodes[3]['swc_type']=7
    source=tmp_path/'raw.swc'
    source.write_text(''.join(f'{i} {graph.nodes[i]["swc_type"]} {i} 0 0 {r[3]:.17g} {i-1}\n' for i,r in coords.items()))
    raw=read_source(source)
    subset=raw.graph.subgraph([2,3,4,5]).copy()
    output=tmp_path/'derived.swc'
    info=export_subset(output,subset,raw,edges,tmp_path/'maps')
    readback=read_source(output)
    assert info['geometry_exact'] and info['type_exact'] and info['round_trip']
    assert validate_tree(readback.graph)==1
    for new,old in enumerate([2,3,4,5],1):
        np.testing.assert_array_equal(readback.graph.nodes[new]['coords'],coords[old])
        assert readback.graph.nodes[new]['swc_type']==graph.nodes[old]['swc_type']
    assert len(list((tmp_path/'maps').glob('*.csv')))==2


def test_real_frozen_refinement_part03_and_split_entrances():
    cache=load_cache(ROOT/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA')
    source=read_source(Path(cache.manifest['source']['raw_source']))
    with frozen_upstream():
        _,graphs,rows,_=refine_component(cache,3,source)
        assert len(graphs)==1 and validate_tree(graphs[0])==1697
        assert graphs[0].out_degree(1697)==1
        expected_stem = sum(np.linalg.norm(source.graph.nodes[n+1]['coords'][:3] -
                                           source.graph.nodes[n]['coords'][:3]) for n in range(1697,1702))
        assert rows[0]['distance_from_original_boundary_mm']==pytest.approx(-expected_stem)
        _,graphs,_,_=refine_component(cache,4,source)
        assert {validate_tree(g) for g in graphs}=={2045,2133}
        edges=branch_edges(cache)
        f,rows,_=prune_diameter(graphs[1],cache.coordinates,edges,4,'part04_02')
        assert not f and rows[0]['reason']=='COMPONENT_BELOW_DIAMETER_THRESHOLD'
    cache.verify_unchanged()


def test_upstream_guard_rejects_queries_and_does_not_allow_silent_failure():
    from vascular_processing import semantic_ensemble
    with pytest.raises(AssertionError,match='FORBIDDEN_UPSTREAM'):
        with frozen_upstream():
            semantic_ensemble.nearest_support()


def test_human_visualization_unchanged():
    assert sha256(ROOT/'s1-2_swc_roi_generate_human.py')==HUMAN_SHA
