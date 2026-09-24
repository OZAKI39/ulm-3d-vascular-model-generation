import json
from pathlib import Path
from fem3d.audit import sha256
ROOT=Path(__file__).resolve().parents[1]


def test_frozen_stage1_baseline_and_acceptance_policy():
    b=json.loads((ROOT/'inputs/stage01_5/baseline_contract.json').read_text())
    assert b['fixed_baseline']['total_lt_0_1']==153
    assert b['fixed_baseline']['cap_adjacent_lt_0_1']==129
    for p,digest in b['verified_stage1_files'].items(): assert sha256(ROOT/p)==digest
    for category in ('reports','inputs'):
        assert sha256(ROOT/category/'stage01_5/acceptance_policy.json')==b['acceptance_policy_sha256']
