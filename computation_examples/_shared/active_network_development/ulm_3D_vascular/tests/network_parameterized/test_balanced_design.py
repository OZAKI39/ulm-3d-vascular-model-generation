"""The requested equal-split design must fail closed when physics forbids it."""
from copy import deepcopy
from pathlib import Path
import ast
import json
import numpy as np
import pytest
from network_1d0d.parameter_config import read_configuration,configuration_targets
from network_1d0d.parameterized_hydraulics import ParameterizedHydraulicSpec,solve_parameterized_operating_point
from network_1d0d.parameter_freeze import freeze_parameters
from network_1d0d.balanced_feasibility import passive_o1_o3_certificate,require_balanced_fit,balance_metrics

ROOT=Path(__file__).resolve().parents[2]
CONFIG=ROOT/'configs/parameterized_hydraulics/balanced_three_outlet_fit.yaml'
FIT=ROOT/'reports/parameterized_0d_v1/balanced_three_outlet'


def test_balanced_target_has_two_independent_parameters(geometry_cache):
    loaded=read_configuration(CONFIG)
    target,truth=configuration_targets(loaded,geometry_cache)
    assert loaded['free']==('s_O1','s_O2') and not loaded['priors']
    assert truth is None and target.kind=='SYNTHETIC_TARGET'
    assert all(o.value==1/3 and o.sigma==.005 for o in target.observations)


def test_equal_split_requires_negative_pressure_relative_to_passive_distal(geometry_cache):
    c=passive_o1_o3_certificate(geometry_cache,mu_pa_s=.00345312,target_roi_flow_m3_s=1.551359160440232e-14)
    assert c['status']=='EQUAL_SPLIT_IMPOSSIBLE_UNDER_PRESCRIBED_MODEL'
    assert c['junction_original_id']==3274
    assert c['required_P_O1_minus_Pd_for_equal_pa']<0
    k=c['ratio_bound_Q_O1_over_Q_O3']
    assert k/(1+k)<1/3
    assert c['necessary_max_abs_deviation_lower_bound']>1e-5


@pytest.mark.parametrize('s1,s2',[(.7,.8),(1.,1.),(2.,1.1),(10.,.893),(10000.,.893)])
def test_fixed_path_certificate_agrees_with_full_network(geometry_cache,s1,s2):
    c=passive_o1_o3_certificate(geometry_cache,mu_pa_s=.00345312,target_roi_flow_m3_s=1.551359160440232e-14)
    r=solve_parameterized_operating_point(geometry_cache,ParameterizedHydraulicSpec(s_O1=s1,s_O2=s2))
    expected=c['paths']['O3']['resistance_pa_s_m3']*r.port_flow_m3_s[3]-c['paths']['O1']['resistance_pa_s_m3']*r.port_flow_m3_s[1]
    assert abs((r.port_pressure_pa[1]-r.port_pressure_pa[3])-expected)<1e-7
    assert r.signed_outlet_fractions[0]<=c['ratio_bound_Q_O1_over_Q_O3']*r.signed_outlet_fractions[2]+1e-11
    assert r.mass_audit['status']=='PASS'


def test_certificate_explicit_viscosity_and_gauge(geometry_cache):
    args=dict(cache=geometry_cache,target_roi_flow_m3_s=1.551359160440232e-14)
    a=passive_o1_o3_certificate(mu_pa_s=.00345312,**args)
    b=passive_o1_o3_certificate(mu_pa_s=.00690624,**args)
    assert b['required_P_O1_minus_Pd_for_equal_pa']==2*a['required_P_O1_minus_Pd_for_equal_pa']
    assert b['ratio_bound_Q_O1_over_Q_O3']==a['ratio_bound_Q_O1_over_Q_O3']
    r=solve_parameterized_operating_point(geometry_cache,ParameterizedHydraulicSpec(distal_reference_pa=500.))
    assert r.port_pressure_pa[3]==500. and r.port_pressure_pa[1]>500.


def test_recorded_balanced_fit_cannot_be_frozen(tmp_path):
    summary=json.loads((FIT/'fit_summary.json').read_text())
    with pytest.raises(ValueError,match='BALANCED_0D_TARGET_NOT_REACHED'):require_balanced_fit(summary)
    output=tmp_path/'must_not_exist.yaml'
    with pytest.raises(ValueError,match='validated forward or converged fit'):
        freeze_parameters(CONFIG,FIT/'fit_summary.json',output)
    assert not output.exists()


def test_balance_tolerance_is_independent_of_success_flag():
    summary=json.loads((FIT/'fit_summary.json').read_text())
    summary=deepcopy(summary);summary['status']='FIT_CONVERGED'
    summary['identifiability_audit']['final_data_jacobian']['rank']=2
    with pytest.raises(ValueError,match='balance tolerance'):require_balanced_fit(summary)


def test_final_failed_fit_matches_analytic_least_squares_limit(geometry_cache):
    s=json.loads((FIT/'fit_summary.json').read_text())
    c=passive_o1_o3_certificate(geometry_cache,mu_pa_s=.00345312,target_roi_flow_m3_s=1.551359160440232e-14)
    f=[s['prediction']['outlet_flow_fraction'][p] for p in ('O1','O2','O3')]
    np.testing.assert_allclose(f,c['equal_sigma_least_squares_boundary_projection'],rtol=0,atol=1e-7)
    assert balance_metrics(f)['max_abs_deviation']>.17


def test_no_cfd_feedback_or_automatic_launch_path():
    for path in (ROOT/'network_1d0d/balanced_feasibility.py',ROOT/'scripts/diagnose_balanced_three_outlet.py'):
        tree=ast.parse(path.read_text())
        imports=[a.name for n in ast.walk(tree) if isinstance(n,(ast.Import,ast.ImportFrom)) for a in n.names]
        assert not any(x.startswith(('subprocess','particle','svMultiPhysics')) for x in imports)
    status=json.loads((ROOT/'reports/balanced_three_outlet_forward_v1/workflow_status.json').read_text())
    assert status['scientific_CFD_cases_created']==status['new_scientific_CFD_run_count']==0
    assert status['automatic_CFD_launch_registered'] is False
    assert not (FIT/'frozen_parameterized_boundary.yaml').exists()
