from pathlib import Path
import importlib.util
import json
import pytest

ROOT=Path(__file__).resolve().parents[2]
REPORT=ROOT/'reports/best_feasible_balance_forward_v1'


def module(name):
    spec=importlib.util.spec_from_file_location(name,ROOT/'scripts'/(name+'.py'))
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def test_exactly_one_scientific_cfd_case_created():
    r=json.loads((REPORT/'scientific_case_registry.json').read_text())
    assert r['scientific_case_count']==1
    p=json.loads((REPORT/'final_preflight.json').read_text())
    assert p['planned_CFD_run_count']==1 and p['design_bounds']=={'s_O1':[.5,2.],'s_O2':[.5,2.]}
    assert len(p['configuration_diff']['fields'])==3
    with pytest.raises(FileExistsError):
        module('prepare_best_balance_fem').prepare(
            ROOT/'reports/parameterized_0d_v1/best_feasible_balance/frozen_balance_design.yaml',
            Path(p['configuration_diff']['original_case']),Path(p['local_case']),REPORT)


def test_dispatch_claim_is_exclusive_and_cannot_be_replaced(tmp_path):
    m=module('launch_best_balance_once');p=tmp_path/'claim.json'
    m.claim_once(p,{'scientific_case_count':1})
    before=p.read_bytes()
    with pytest.raises(FileExistsError):m.claim_once(p,{'scientific_case_count':2})
    assert p.read_bytes()==before
