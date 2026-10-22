"""Refactor-specific scientific identity, boundary and output protection checks."""
from pathlib import Path
import json
import numpy as np
import pytest
from human_vessel_flow.paths import ROOT, read_json, fresh_directory

def test_default_sonovue_distribution_contract_unchanged():
    from particle_3d.injection_method_c import TruncatedSonoVue
    d=TruncatedSonoVue(ROOT/'inputs/distributions/sonovue')
    c=read_json(ROOT/'configs/population/FIXED_WINDOW_POPULATION_EXPERIMENT_V1.json')
    assert d.contract_sha256==c['source_distribution_sha256']

def test_centreline_pressure_recipe_matches_original_numbers():
    from human_vessel_flow.boundary import estimate
    actual,edges=estimate();expected=read_json(ROOT/'inputs/geometry/boundary_design.json')
    assert actual['outlet_pressures_Pa']==expected['outlet_pressures_Pa']
    assert len(edges)>0
    assert actual['outlet_pressures_Pa']['O1']!=556.9706555488448

def test_new_case_uses_centreline_pressures_and_no_links(tmp_path):
    from human_vessel_flow.fem import prepare_case
    import xml.etree.ElementTree as ET
    case=prepare_case(tmp_path/'new_case')
    xml=ET.parse(case/'run/solver.xml')
    values={e.get('name'):float(e.findtext('Value')) for e in xml.findall('.//Add_BC')}
    assert values['INLET']==-3e-7
    assert values['OUTLET_01']==508.233436510425
    assert values['OUTLET_02']==0 and values['OUTLET_03']==24.309360181533066
    assert all(not p.is_symlink() for p in case.rglob('*'))
    assert not (case/'run/solver.log').exists()
    with pytest.raises(FileExistsError):prepare_case(case)

def test_frozen_output_protection():
    for relative in ['cases/balanced_pressure_final/new','inputs/new','diagnostics/new','results/population/new']:
        with pytest.raises(ValueError):fresh_directory(ROOT/relative)

def test_rebuild_accepted_surface_does_not_change_arrays(tmp_path):
    from human_vessel_flow.geometry import prepare_surface
    output=tmp_path/'geometry';prepare_surface(output)
    with np.load(output/'source_surface.npz') as actual,np.load(ROOT/'inputs/geometry/source_surface.npz') as expected:
        for key in expected.files:np.testing.assert_array_equal(actual[key],expected[key])

def test_future_population_plan_keeps_fixed_protocol(tmp_path):
    from human_vessel_flow.particles import prepare_run
    job=tmp_path/'future';command=prepare_run(job,'cpu')
    c=read_json(job/'contract.json')
    assert c['T_observation_s']==10 and c['T_end_s']==11.240815660
    assert not c['outcome_dependent_stop'] and not (job/'results').exists()
    assert 'particle_3d.population_v1.run_fixed_window_population' in command
    assert '--minimum-exits' not in command

def test_p2_memory_gate_keeps_resource_stop():
    from human_vessel_flow.p2_gpu_budget import gate
    assert not gate(None,8192)['PASS']
    budget=dict(status='BOUNDED',MPI_ranks=8,source='native_workspace_query_without_numeric_allocation',
                sum_rank_simultaneous_workspace_bound_MiB=4096,resident_increment_bound_MiB=4096)
    assert not gate(budget,8192)['PASS']

def test_p2_curved_reference_remains_separate_from_native_solver():
    from human_vessel_flow.p2_geometry import geometry, RN10
    result=geometry(np.array([.2,.15,.1]),RN10.copy())
    np.testing.assert_allclose(result['J'],np.eye(3),atol=1e-15)
    np.testing.assert_allclose(result['N'].sum(),1.,atol=1e-15)
    np.testing.assert_allclose(result['Dp'].sum(axis=0),0.,atol=1e-15)
