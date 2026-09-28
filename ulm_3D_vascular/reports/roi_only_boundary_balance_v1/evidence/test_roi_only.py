from pathlib import Path
from dataclasses import replace
import ast
import inspect
import json
import numpy as np
import pytest
from network_1d0d import roi_only_hydraulics as hydraulic,roi_boundary_design as design
from network_1d0d.hydraulic_resistance import ROI_TARGET_Q_M3_S as Q
from network_1d0d.fem_h0_case import serialize_fixed_pressure_xml,audit_xml_change

ROOT=Path(__file__).resolve().parents[2]
H0=ROOT/'reports/a_network_1d0d_boundary_v2_idealized/data'
GRAPH=H0/'analysis_A_H0_graph_si.npz'
PORTS=ROOT/'reports/a_network_1d0d_boundary_v1/data/roi_ports_in_a.json'

@pytest.fixture(scope='module')
def cache():return hydraulic.load_roi_cache(GRAPH,PORTS)

@pytest.fixture(scope='module')
def regression(cache):return design.regress_h0(cache,json.loads((H0/'roi_operating_point_H0.json').read_text()))

@pytest.fixture(scope='module')
def solution(cache,regression):
    return design.design(cache,regression,json.loads((ROOT/'reports/balanced_three_outlet_forward_v1/feasibility_certificate.json').read_text()))[0]

def test_roi_only_geometry_contains_only_internal_edges(cache):
    with np.load(GRAPH) as g:
        np.testing.assert_array_equal(cache.roi_edge_indices,np.flatnonzero(g['roi_internal_edge_mask']))
        np.testing.assert_array_equal(cache.source_node_indices[cache.edge_nodes],g['edges'][g['roi_internal_edge_mask']])
    assert len(cache.node_ids)==109 and len(cache.edge_nodes)==108

def test_roi_only_ports_are_exactly_four(cache):
    saved=json.loads((H0/'roi_port_flow_sign_convention.json').read_text())['ports']
    for row,p in zip(saved,cache.summary()['ports']):
        assert row['port']==p['port'] and row['node_id']==p['node_id']
        assert row['roi_segment']==p['source_edge_index'] and row['coefficient']==p['coefficient']
    assert len(cache.port_node_indices)==4

def test_roi_mixed_bc_mass_balance(cache):
    r=hydraulic.solve_roi_mixed_bc(cache,[-800,3700,0])
    assert r['mass_audit']['status']=='PASS'
    assert abs(sum(r['outlet_fraction'])-1)<1e-10

def test_roi_mixed_bc_regresses_fullA_port_state(regression):
    assert regression['maximum_port_flow_error_relative_Qin']<1e-10

def test_roi_pressure_gauge_invariance(cache):
    a=hydraulic.solve_roi_mixed_bc(cache,[-800,3700,0]);b=hydraulic.solve_roi_mixed_bc(cache,[-800+2000,5700,2000])
    np.testing.assert_allclose(a['port_flow_m3_s'],b['port_flow_m3_s'],rtol=1e-10,atol=1e-26)
    np.testing.assert_allclose(b['pressure_pa']-a['pressure_pa'],2000,rtol=1e-11)

def test_roi_basis_map_is_linear(cache):
    b=design.pressure_basis(cache)
    for dp in [[-1200,1000],[340,3065],[300,-1700]]:
        actual=hydraulic.solve_roi_mixed_bc(cache,[*dp,0])['port_flow_m3_s'][1:]
        np.testing.assert_allclose(np.array(b['q0_m3_s'])+np.array(b['B_m3_s_per_pa'])@dp,actual,rtol=1e-9,atol=1e-26)

def test_roi_basis_solution_matches_iterative_solution(solution):
    assert max(abs(np.array(solution['direct_minus_iterative_pa'])))<1e-6
    assert solution['max_abs_deviation_from_one_third']<1e-8

def test_roi_equal_flow_requires_o1_lower_pressure_than_o3(solution):
    assert solution['relative_realcut_pressure_pa']['O1']<solution['relative_realcut_pressure_pa']['O3']

def test_roi_equal_flow_path_identity_matches_certificate(solution):
    assert solution['relative_realcut_pressure_pa']['O1']==pytest.approx(-795.60984355,abs=1e-7)
    assert solution['relative_realcut_pressure_pa']['O1']==pytest.approx(solution['fixed_path_identity']['equal_split_dP_O1_O3_pa'],abs=1e-7)

def test_roi_design_does_not_import_fullA_downstream_masks():
    for module in (hydraulic,design):
        source=inspect.getsource(module);tree=ast.parse(source)
        imports=[n.module or '' for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]
        assert not any('parameterized' in n or 'idealized_h0' in n or 'balance_design' in n for n in imports)
        assert 'downstream_node_mask' not in source and 'downstream_edge_mask' not in source

def test_roi_design_does_not_modify_radius(cache,solution):
    with np.load(GRAPH) as g:np.testing.assert_array_equal(cache.radius_m,g['radius_m'][cache.source_node_indices])
    with pytest.raises(ValueError):cache.radius_m.setflags(write=True)

def test_roi_design_has_zero_cfd_dependency():
    for module in (hydraulic,design):
        tree=ast.parse(inspect.getsource(module))
        names=[n.id for n in ast.walk(tree) if isinstance(n,ast.Name)]
        assert not set(names)&{'subprocess','pyvista','SolutionMeasurements','svmultiphysics'}

def test_roi_frozen_design_roundtrip(solution,tmp_path):
    p=tmp_path/'frozen.yaml';design.freeze_design(p,solution)
    assert design.load_frozen_design(p)['outlet_fraction']==solution['outlet_fraction']
    with pytest.raises(FileExistsError):design.freeze_design(p,solution)
    p.write_text(p.read_text().replace('particle_RBC_calls: 0','particle_RBC_calls: 1'))
    with pytest.raises(ValueError,match='modified'):design.load_frozen_design(p)

def test_only_three_fem_outlet_values_change(solution):
    path=ROOT.parent/'ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1/run/solver.xml'
    old=path.read_bytes();values=[0.,4500.,850.]
    new=serialize_fixed_pressure_xml(old,values)
    assert len(audit_xml_change(old,new,values))==3
    with pytest.raises(ValueError):audit_xml_change(old,new.replace(b'1056',b'1060'),values)

def test_reversed_edge_storage_preserves_flow(cache):
    reverse=replace(cache,edge_nodes=cache.edge_nodes[:,::-1],port_coefficients=-cache.port_coefficients)
    a=hydraulic.solve_roi_mixed_bc(cache,[-800,3700,0]);b=hydraulic.solve_roi_mixed_bc(reverse,[-800,3700,0])
    np.testing.assert_allclose(a['port_flow_m3_s'],b['port_flow_m3_s'],rtol=1e-11,atol=1e-26)

def test_reverse_flow_is_signed_not_absoluted(cache):
    r=hydraulic.solve_roi_mixed_bc(cache,[20000,0,0])
    assert r['port_flow_m3_s'][1]<0 and r['outlet_fraction'][0]<0

def test_invalid_inlet_sign_is_rejected(cache):
    with pytest.raises(ValueError):hydraulic.solve_roi_mixed_bc(cache,[0,0,0],-Q)
