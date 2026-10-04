from pathlib import Path
import hashlib
import json
import numpy as np
import pytest
from scripts.verify_at1_evolution import OUT, config
from pd_clot.pfpd_candidate.adr_adapter_at1 import fixture, boundary_view, tensile_domain, FrozenPhaseForce
from pd_clot.pfpd_candidate.baqp import evaluate


def test_multibond_energy_and_bitwise_kinematic_identity():
    summary = json.loads((OUT/'fixed_x/SUMMARY.json').read_text())
    assert all(summary['bitwise_identical'].values())
    e = summary['energies_J']
    assert e['elastic_after'] < e['elastic_before']
    assert e['crack_after'] > e['crack_before']
    assert e['total_after'] < e['total_before']
    assert summary['existing_crack_energy_difference_J'] == pytest.approx(0.,abs=1e-26)


def test_frozen_phase_callback_is_gradient_of_degraded_elastic_energy():
    cloud,f,material,_ = fixture(config())
    x = f.X.copy(); x[:,2] *= 1.012
    s = np.linspace(.1,.3,len(f.pairs)); callback = FrozenPhaseForce(f,material,s)
    state = callback(x); h = cloud.spacing*1e-5
    plus,minus = x.copy(),x.copy(); plus[4,2] += h; minus[4,2] -= h
    gradient = (callback(plus).stored_energy_J-callback(minus).stored_energy_J)/(2*h)
    assert gradient == pytest.approx(-state.internal_force_N[4,2],rel=2e-6,abs=1e-15)
    assert not np.any(state.external_force_N)


def test_boundary_view_never_moves_original_baqp_reference():
    cloud,f,_,top = fixture(config()); original = f.X.copy()
    view = boundary_view(cloud,top,.02)
    assert np.array_equal(f.X,original) and np.array_equal(cloud.X,original)
    assert np.all(view.X[top,2] > original[top,2])
    assert np.array_equal(view.X[~cloud.fixed],original[~cloud.fixed])


def test_compression_is_rejected_without_inventing_energy_split():
    _,f,material,_ = fixture(config()); x = f.X.copy(); x[:,2] *= .99
    with pytest.raises(ValueError,match='PURE_TENSION_DOMAIN_VIOLATION'):
        tensile_domain(evaluate(f,x,material),-1e-9)


def test_all_preexisting_python_modules_and_tests_are_unchanged():
    manifest = json.loads((OUT/'provenance/PRE_RUN_MANIFEST.json').read_text())
    root = Path(manifest['repo_root'])
    protected = [(name,record) for name,record in manifest['files'].items()
                 if name.endswith('.py') and ('/pd_clot/' in name or '/tests/' in name)]
    assert len(protected) > 50
    for name,record in protected:
        assert hashlib.sha256((root/name).read_bytes()).hexdigest() == record['sha256'],name
