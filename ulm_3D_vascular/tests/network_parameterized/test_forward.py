from dataclasses import replace
from pathlib import Path
import json
import numpy as np
import pytest
from network_1d0d.idealized_h0 import solve_operating_point, external_radius_variant, port_signs
from network_1d0d.hydraulic_resistance import linear_radius_resistance, ROI_TARGET_Q_M3_S
from network_1d0d.parameterized_hydraulics import (
    GeometryHydraulicCache, ParameterizedHydraulicSpec as Spec,
    solve_parameterized_operating_point as solve, downstream_node_mask, downstream_edge_mask)


def test_legacy_h0_parameterized_regression(geometry_cache,geometry_domain):
    saved=json.loads((Path(__file__).resolve().parents[2]/'reports/a_network_1d0d_boundary_v2_idealized/data/roi_operating_point_H0.json').read_text())
    actual=solve(geometry_cache,Spec())
    legacy=solve_operating_point(geometry_domain)
    np.testing.assert_allclose(actual.signed_outlet_fractions,[r['signed_fraction'] for r in saved['ports'][1:]],rtol=2e-11,atol=2e-13)
    np.testing.assert_allclose(actual.port_pressure_pa,[r['pressure_realcut_Pa'] for r in saved['ports']],rtol=2e-11,atol=1e-9)
    assert actual.port_flow_m3_s[0]==pytest.approx(ROI_TARGET_Q_M3_S,rel=2e-15,abs=0)
    assert actual.mass_audit['max_relative_residual']<1e-10
    np.testing.assert_allclose(actual.pressure_pa,legacy['pressure'],rtol=2e-12,atol=1e-9)
    np.testing.assert_allclose(actual.edge_flow_m3_s,legacy['flow'],rtol=2e-10,atol=1e-26)


@pytest.mark.parametrize('port',['O1','O2'])
def test_downstream_mask_matches_existing_behavior(geometry_domain,geometry_cache,port):
    r,oldmask=external_radius_variant(geometry_domain,port,.95)
    mask=downstream_node_mask(geometry_domain,port)
    np.testing.assert_array_equal(mask,oldmask)
    actual=solve(geometry_cache,Spec(**{'s_'+port:.95}))
    np.testing.assert_array_equal(actual.effective_radius_m,r)
    np.testing.assert_array_equal(r[geometry_cache.port_indices],geometry_domain.radius_m[geometry_cache.port_indices])
    edge_mask=downstream_edge_mask(geometry_domain,port)
    transition=np.flatnonzero(edge_mask & (mask[geometry_domain.edges].sum(axis=1)==1))
    assert len(transition)==1
    i=transition[0]
    assert not np.isclose(actual.edge_resistance_pa_s_m3[i]/geometry_cache.baseline_resistance_pa_s_m3[i],.95**-4,rtol=1e-3)
    np.testing.assert_allclose(actual.edge_resistance_pa_s_m3,geometry_domain.resistances(r),rtol=1e-15)


def test_o3_has_no_downstream_geometry(geometry_domain):
    with pytest.raises(ValueError,match='NOT_AVAILABLE'): downstream_node_mask(geometry_domain,'O3')


def test_geometry_snapshot_is_not_mutated(geometry_domain,geometry_cache):
    before=geometry_domain.radius_m.copy()
    solve(geometry_cache,Spec(s_O1=.92,s_O2=1.07))
    np.testing.assert_array_equal(before,geometry_domain.radius_m)
    np.testing.assert_array_equal(before,geometry_cache.baseline_radius_m)
    with pytest.raises(ValueError): geometry_cache.baseline_radius_m[0]=1
    with pytest.raises(ValueError): geometry_cache.baseline_radius_m.setflags(write=True)


def test_port_orientation_independent_of_edge_storage(geometry_domain,geometry_cache):
    d=geometry_domain
    edges=d.edges[:,::-1].copy()
    signs=port_signs(d.ids,d.xyz_m,edges,d.internal_roi,d.ports)
    reversed_cache=GeometryHydraulicCache.from_domain(replace(d,edges=edges,signs=signs))
    original=solve(geometry_cache,Spec());rev=solve(reversed_cache,Spec())
    np.testing.assert_allclose(original.port_flow_m3_s,rev.port_flow_m3_s,rtol=2e-10,atol=1e-26)
    np.testing.assert_allclose(original.edge_flow_m3_s,-rev.edge_flow_m3_s,rtol=2e-10,atol=1e-26)


def test_common_pressure_reference_is_added_not_scaled(geometry_cache):
    a=solve(geometry_cache,Spec());b=solve(geometry_cache,Spec(distal_reference_pa=123.))
    np.testing.assert_array_equal(a.edge_flow_m3_s,b.edge_flow_m3_s)
    np.testing.assert_allclose(b.pressure_pa-a.pressure_pa,123.,rtol=1e-13,atol=1e-12)


@pytest.mark.parametrize('kw',[{'s_O1':0},{'s_O2':-1},{'mu_pa_s':float('nan')},{'mu_pa_s':0}])
def test_invalid_physics_is_rejected(geometry_cache,kw):
    with pytest.raises(ValueError):solve(geometry_cache,Spec(**kw))


@pytest.mark.parametrize('q',[0,-1,float('nan')])
def test_invalid_inlet_target_is_not_absoluted(geometry_cache,q):
    with pytest.raises(ValueError,match='never take abs'):solve(geometry_cache,Spec(),q)


def test_geometry_resistance_formula_regression():
    from utils.cfd_preprocess.one_d_flow import edge_resistance
    for ratio in [1.,1+1e-12,1+1e-8,.3,2.,10.]:
        r,L,mu=1.4e-6,37e-6,.00345312
        a=float(linear_radius_resistance(L,r,r*ratio,mu))
        assert a==pytest.approx(edge_resistance(L,r,r*ratio,mu),rel=2e-13)
        if ratio==1:assert a==pytest.approx(8*mu*L/(np.pi*r**4),rel=2e-15)


def test_o2_downstream_radius_monotonicity(geometry_cache):
    results=[solve(geometry_cache,Spec(s_O2=s)) for s in [.95,1.,1.05]]
    f=[r.signed_outlet_fractions[1] for r in results]
    assert f[0]<f[1]<f[2]
    for r in results:
        assert r.mass_audit['status']=='PASS'
        assert np.all(np.isfinite(r.effective_radius_m)) and np.all(r.effective_radius_m>0)
        assert np.all(np.isfinite(r.edge_resistance_pa_s_m3)) and np.all(r.edge_resistance_pa_s_m3>0)


@pytest.mark.parametrize('s1,s2',[(.8,1.2),(1.15,.85),(.92,1.07)])
def test_mass_balance_parameterized_network(geometry_cache,s1,s2):
    r=solve(geometry_cache,Spec(s_O1=s1,s_O2=s2))
    assert r.mass_audit['max_relative_residual']<1e-10
    assert r.signed_outlet_fractions.sum()==pytest.approx(1.,abs=1e-11)


@pytest.mark.parametrize('resistance',[1e16,1e17,5e17])
def test_o3_terminal_resistance_closure(geometry_cache,resistance):
    pd=87.
    r=solve(geometry_cache,Spec(o3_mode='terminal_resistance',terminal_resistance_O3_pa_s_m3=resistance,distal_reference_pa=pd))
    assert r.port_pressure_pa[3]-pd==pytest.approx(resistance*r.port_flow_m3_s[3],rel=1e-10)
    assert r.pressure_pa[-1]==pd and r.port_pressure_pa[3]>pd
    assert int(geometry_cache.port_indices[3]) not in r.solver_audit['fixed_pressure_node_indices']
    assert len(geometry_cache.node_ids) in r.solver_audit['fixed_pressure_node_indices']
    assert len(r.pressure_pa)==len(geometry_cache.node_ids)+1
    assert r.closure_provenance['source']=='MODEL_CLOSURE_PARAMETER'
    assert r.closure_provenance['virtual_nodes'][0]['original_swc_node_id'] is None
    assert r.mass_audit['max_relative_residual']<1e-10


@pytest.mark.parametrize('resistance',[None,0,-1,float('inf')])
def test_invalid_closure_is_rejected(geometry_cache,resistance):
    with pytest.raises(ValueError):solve(geometry_cache,Spec(o3_mode='terminal_resistance',terminal_resistance_O3_pa_s_m3=resistance))


def test_source_component_and_geometry_guards(geometry_domain):
    d=geometry_domain
    with pytest.raises(ValueError,match='connected'):
        GeometryHydraulicCache.from_domain(replace(d,edges=d.edges[1:],internal_roi=d.internal_roi[1:]))
    r=d.radius_m.copy();r[0]=0.
    with pytest.raises(ValueError,match='radius'):GeometryHydraulicCache.from_domain(replace(d,radius_m=r))


def test_wrong_inlet_orientation_is_not_hidden(geometry_cache):
    coefficients=geometry_cache.port_coefficients.copy();coefficients[0]*=-1
    bad=replace(geometry_cache,port_coefficients=coefficients)
    with pytest.raises(ValueError,match='never take abs'):solve(bad,Spec())


def test_nonfinite_port_normal_is_rejected(geometry_domain):
    import copy
    ports=copy.deepcopy(geometry_domain.ports);ports[0]['outward_normal'][0]=float('nan')
    with pytest.raises(ValueError,match='normal'):
        GeometryHydraulicCache.from_domain(replace(geometry_domain,ports=ports))
