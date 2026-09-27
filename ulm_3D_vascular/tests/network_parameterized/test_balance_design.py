from pathlib import Path
from dataclasses import replace
import ast
import json
import numpy as np
import pytest
import yaml
from network_1d0d.audit import sha256
from network_1d0d.balance_design import (DesignEnvelope,DesignEvaluator,metrics,run_design,
    acceptance,freeze_design,load_frozen_design)

ROOT=Path(__file__).resolve().parents[2]


@pytest.fixture(scope='module')
def design_run(tmp_path_factory):
    output=tmp_path_factory.mktemp('bounded-design')
    original=ROOT/'configs/parameterized_hydraulics/best_feasible_balance_design.yaml'
    c=yaml.safe_load(original.read_text())
    c['geometry']={k:str((original.parent/v).resolve()) for k,v in c['geometry'].items()}
    c['infeasibility_certificate']=str((original.parent/c['infeasibility_certificate']).resolve())
    c['grid_size_per_axis']=5
    config=output/'config.yaml';config.write_text(yaml.safe_dump(c))
    result,cache=run_design(config,output/'results')
    return result,cache,output


def test_equal_split_remains_certified_infeasible(design_run):
    s,_,_=design_run
    assert s['equal_split_feasible'] is False and s['equal_split_required'] is False
    assert np.isclose(s['ratio_bound_Q_O1_over_Q_O3'],.35320878510001563)
    assert s['optimized']['metrics']['max_abs_deviation']>1e-5


def test_calibration_identifiability_guard_unchanged():
    assert sha256(ROOT/'network_1d0d/parameter_fit.py')=='4bad7b5760c6976acfa5bb9550c477101dc4a9c59e0d38e5e55011cbe5b56610'


def test_design_optimizer_respects_bounds(design_run):
    import csv
    s,_,out=design_run
    rows=list(csv.DictReader((out/'results/design_trace.csv').open()))
    assert all(.5<=float(r[p])<=2. for r in rows for p in ('s_O1','s_O2'))
    assert all(acceptance(s).values())


@pytest.mark.parametrize('upper',[float('inf'),float('nan'),-.1,.5])
def test_design_does_not_allow_infinite_radius_escape(upper):
    with pytest.raises(ValueError):DesignEnvelope((.5,.5),(upper,2.))


def test_out_of_envelope_is_rejected_not_clipped(geometry_cache):
    evaluator=DesignEvaluator(geometry_cache,DesignEnvelope())
    with pytest.raises(ValueError,match='outside'):evaluator.evaluate((1e8,1.))
    assert evaluator.trace[-1]['status']=='FAILED'


def test_design_improves_balance_over_baseline(design_run):
    s,_,_=design_run
    assert s['optimized']['metrics']['J_balance']<s['baseline']['metrics']['J_balance']
    assert s['optimized']['metrics']['range']<s['baseline']['metrics']['range']


def test_design_objective_matches_fraction_variance(design_run):
    s,_,_=design_run
    for result in (s['baseline'],s['optimized']):
        f=list(result['prediction']['outlet_flow_fraction'].values())
        np.testing.assert_allclose(metrics(f)['J_balance'],3*np.var(f),rtol=1e-13,atol=1e-15)


def test_design_trace_is_deterministic(design_run):
    s,_,_=design_run
    assert s['deterministic_reproducibility']['trace_payload_exact_equal'] is True


def test_grid_and_local_design_agree(design_run):
    s,_,_=design_run
    assert s['crosscheck']['status']=='CONSISTENT' and s['crosscheck']['grid_objective_not_better']


def test_bound_active_is_reported_not_rejected(design_run):
    s,_,_=design_run
    assert s['status']=='BEST_FEASIBLE_BALANCE_DESIGN_PASS'
    assert s['optimum_classification']=='BOUND_ACTIVE_OPTIMUM'
    assert s['optimized']['active_bounds']['s_O1']=='UPPER'


def test_frozen_design_roundtrip_and_no_data_identification_claim(design_run):
    s,cache,out=design_run
    path=out/'frozen_balance_design.yaml'
    f=freeze_design(s,cache,path)
    assert load_frozen_design(path)==f
    assert f['data_identification_claim']=='NONE' and f['status']=='FROZEN_DESIGN'
    assert f['workflow_kind']=='DESIGN_OPTIMIZATION' and f['parameter_interpretation']=='DESIGN_VARIABLES'
    assert path.stat().st_mode & 0o222==0
    with pytest.raises(FileExistsError):freeze_design(s,cache,path)
    d=yaml.safe_load(path.read_text());d['optimized']['parameters']['s_O1']=1.
    bad=out/'tampered.yaml';bad.write_text(yaml.safe_dump(d))
    with pytest.raises(ValueError,match='integrity'):load_frozen_design(bad)


def test_sensitivity_does_not_replace_primary(design_run):
    s,_,_=design_run
    assert [r['s_max_dimensionless'] for r in s['bound_sensitivity']]==[1.25,1.5,2.,3.]
    assert s['design_bounds']=={'s_O1':[.5,2.],'s_O2':[.5,2.]}


def test_design_has_zero_cfd_dependency():
    for path in (ROOT/'network_1d0d/balance_design.py',ROOT/'scripts/optimize_balanced_0d_design.py'):
        source=path.read_text();tree=ast.parse(source)
        modules=[n.module for n in ast.walk(tree) if isinstance(n,ast.ImportFrom)]+[a.name for n in ast.walk(tree) if isinstance(n,ast.Import) for a in n.names]
        assert not any(m and any(token in m.lower() for token in ('subprocess','particle','fem','hemocell')) for m in modules)
