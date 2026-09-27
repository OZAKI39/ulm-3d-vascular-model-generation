from dataclasses import replace
import csv
import numpy as np
import pytest
from network_1d0d.parameterized_hydraulics import ParameterizedHydraulicSpec as Spec, solve_parameterized_operating_point as solve
from network_1d0d.parameter_fit import (
    Observation, Targets, LogPrior, UnidentifiableParameterizationError,fit_parameters)


def synthetic(cache,spec=None):
    true=spec or Spec(s_O1=.92,s_O2=1.07)
    r=solve(cache,true)
    targets=Targets('SYNTHETIC_TARGET',tuple(Observation('flow_fraction',p,float(f),.005)
        for p,f in zip(('O1','O2','O3'),r.signed_outlet_fractions)), 'pure_0D_generated_test')
    return true,targets


def test_synthetic_two_parameter_recovery(geometry_cache,tmp_path):
    true,target=synthetic(geometry_cache)
    fit=fit_parameters(geometry_cache,Spec(),target,trace_path=tmp_path/'trace.csv')
    assert fit['status']=='FIT_CONVERGED'
    assert fit['final_parameters']['s_O1']==pytest.approx(true.s_O1,rel=1e-7)
    assert fit['final_parameters']['s_O2']==pytest.approx(true.s_O2,rel=1e-7)
    assert fit['final_objective_norm']<1e-8
    assert fit['initial_objective_norm']>1
    assert fit['identifiability_audit']['final_data_jacobian']['rank']==2
    assert fit['identifiability_audit']['fixed_anchor_parameters']==[]
    assert fit['identifiability_audit']['inactive_parameters']==['terminal_resistance_O3_pa_s_m3']
    trace=list(csv.DictReader((tmp_path/'trace.csv').open()))
    assert len(trace)==fit['total_forward_evaluations'] and len(trace)>fit['nfev']
    assert all(float(row['mass_residual_dimensionless'])<1e-9 for row in trace)


def test_underdetermined_fit_is_rejected(geometry_cache):
    spec=Spec(o3_mode='terminal_resistance',terminal_resistance_O3_pa_s_m3=1e17)
    _,target=synthetic(geometry_cache,spec)
    with pytest.raises(UnidentifiableParameterizationError) as caught:
        fit_parameters(geometry_cache,spec,target,free_parameters=('s_O1','s_O2','terminal_resistance_O3_pa_s_m3'))
    a=caught.value.audit
    assert a['number_of_free_parameters']==3
    assert a['number_of_independent_flow_constraints']==2
    assert a['pressure_constraints']==a['prior_constraints']==0


def test_flow_and_Q_are_not_double_counted(geometry_cache):
    spec=Spec(o3_mode='terminal_resistance',terminal_resistance_O3_pa_s_m3=1e17)
    _,target=synthetic(geometry_cache,spec);r=solve(geometry_cache,spec)
    obs=target.observations+tuple(Observation('outlet_flow_m3_s',p,float(q),1e-16) for p,q in zip(('O1','O2','O3'),r.port_flow_m3_s[1:]))
    with pytest.raises(UnidentifiableParameterizationError):
        fit_parameters(geometry_cache,spec,replace(target,observations=obs),free_parameters=('s_O1','s_O2','terminal_resistance_O3_pa_s_m3'))


def test_fixed_reference_pressure_does_not_identify_extra_parameter(geometry_cache):
    target=Targets('SYNTHETIC_TARGET',(Observation('real_cut_pressure_pa','O3',0.,1.),),'fixed_reference_is_not_data')
    with pytest.raises(UnidentifiableParameterizationError):fit_parameters(geometry_cache,Spec(),target,free_parameters=('s_O1',))


def test_pressure_observation_can_anchor_third_parameter(geometry_cache):
    true=Spec(s_O1=.92,s_O2=1.07,o3_mode='terminal_resistance',terminal_resistance_O3_pa_s_m3=1.2e17)
    _,target=synthetic(geometry_cache,true);r=solve(geometry_cache,true)
    target=replace(target,observations=target.observations+(Observation('real_cut_pressure_pa','O3',float(r.port_pressure_pa[3]),10.),))
    fit=fit_parameters(geometry_cache,replace(true,s_O1=1,s_O2=1,terminal_resistance_O3_pa_s_m3=1e17),target,
                       free_parameters=('s_O1','s_O2','terminal_resistance_O3_pa_s_m3'))
    assert fit['status']=='FIT_CONVERGED'
    for name in ['s_O1','s_O2','terminal_resistance_O3_pa_s_m3']:
        assert fit['final_parameters'][name]==pytest.approx(getattr(true,name),rel=1e-7)
    assert fit['identifiability_audit']['final_data_jacobian']['rank']==3


def test_prior_regularization_is_not_called_data_identification(geometry_cache):
    spec=Spec(o3_mode='terminal_resistance',terminal_resistance_O3_pa_s_m3=1e17)
    _,target=synthetic(geometry_cache,spec)
    fit=fit_parameters(geometry_cache,spec,target,free_parameters=('s_O1','s_O2','terminal_resistance_O3_pa_s_m3'),
                       priors=(LogPrior('terminal_resistance_O3_pa_s_m3',0.,.1),))
    assert fit['identifiability_audit']['final_data_jacobian']['rank']==2
    assert fit['identifiability_audit']['final_augmented_jacobian']['rank']==3
    assert fit['identifiability_audit']['interpretation']=='PRIOR_REGULARIZED_NOT_DATA_IDENTIFIED'


def test_parameter_trace_is_deterministic(geometry_cache,tmp_path):
    _,target=synthetic(geometry_cache)
    fits=[]
    for i in range(2):fits.append(fit_parameters(geometry_cache,Spec(),target,trace_path=tmp_path/f'trace{i}.csv'))
    assert (tmp_path/'trace0.csv').read_bytes()==(tmp_path/'trace1.csv').read_bytes()
    assert fits[0]['final_parameters']==fits[1]['final_parameters']


@pytest.mark.parametrize('kind',['CFD_TARGET','PARTICLE_TARGET','CFD_SURROGATE',''])
def test_prohibited_target_kind_is_rejected(geometry_cache,kind):
    _,target=synthetic(geometry_cache)
    with pytest.raises(ValueError,match='forbidden'):fit_parameters(geometry_cache,Spec(),replace(target,kind=kind))


def test_experimental_observation_residual_units(geometry_cache):
    r=solve(geometry_cache,Spec())
    target=Targets('EXPERIMENTAL_TARGET',(
        Observation('outlet_flow_m3_s','O1',float(r.port_flow_m3_s[1]),1e-16),
        Observation('pressure_difference_pa','INLET',float(r.port_pressure_pa[0]-r.port_pressure_pa[2]),10.,'O2')),
        'unit_test_only_not_real_experimental_data')
    fit=fit_parameters(geometry_cache,Spec(),target)
    assert fit['final_objective_norm']==0


def test_legacy_solver_is_not_used_as_new_forward_return(geometry_cache,monkeypatch):
    import network_1d0d.idealized_h0 as legacy
    def forbidden(*a,**k):raise AssertionError('Legacy forward return must not be used')
    monkeypatch.setattr(legacy,'solve_operating_point',forbidden)
    monkeypatch.setattr(legacy,'load_h0',forbidden)
    _,target=synthetic(geometry_cache)
    assert fit_parameters(geometry_cache,Spec(),target)['status']=='FIT_CONVERGED'


def test_numerical_rank_deficiency_is_rejected_even_when_counts_match(geometry_cache,monkeypatch):
    import network_1d0d.parameter_fit as pf
    _,target=synthetic(geometry_cache)
    original=pf.solve_parameterized_operating_point
    # Deliberately remove one parameter's effect to exercise the rank guard.
    monkeypatch.setattr(pf,'solve_parameterized_operating_point',lambda c,s,q:original(c,replace(s,s_O2=1.),q))
    with pytest.raises(UnidentifiableParameterizationError,match='rank deficient'):
        fit_parameters(geometry_cache,Spec(),target)
