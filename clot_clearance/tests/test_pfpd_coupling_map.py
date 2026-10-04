"""Read-only regression checks for the separate coupling diagnostic API."""
import hashlib
from pathlib import Path
import numpy as np
from pd_clot.pfpd_candidate.coupling_stability_audit import (
    OUT, audit_config, read_json, load_saved_iteration)


def test_phi_reproduces_saved_next_iterate_without_mutating_inputs_or_file():
    p,d=load_saved_iteration('loading_unloading',10,6)
    source=Path(d['path']); before=hashlib.sha256(source.read_bytes()).hexdigest()
    s=d['s'].copy(); x=d['x_start'].copy(); lower=p.lower.copy()
    result=p.phi(s,x)
    assert result['status']=='OK'
    np.testing.assert_allclose(result['s_output'],d['s_next'],rtol=0,atol=audit_config()['map_reproduction_s_atol'])
    np.testing.assert_allclose(result['x_equilibrium'],d['x'],rtol=0,atol=audit_config()['map_reproduction_x_atol_m'])
    np.testing.assert_array_equal(s,d['s']); np.testing.assert_array_equal(x,d['x_start'])
    np.testing.assert_array_equal(p.lower,lower)
    assert hashlib.sha256(source.read_bytes()).hexdigest()==before


def test_all_predeclared_reproduction_cases_pass():
    result=read_json(OUT/'fixed_point/REPRODUCTION.json'); a=audit_config()
    assert len(result['cases'])==len(a['map_reproduction_cases'])
    assert all(r['status']=='OK' and r['max_s_error']<=a['map_reproduction_s_atol']
               and r['max_x_error_m']<=a['map_reproduction_x_atol_m'] for r in result['cases'])


def test_preexisting_production_candidate_and_tests_have_original_hashes():
    manifest=read_json(OUT/'provenance/PRE_RUN_MANIFEST.json'); repo=Path(manifest['repo_root'])
    checked=0
    for name,record in manifest['files'].items():
        if name.startswith(('clot_clearance/pd_clot/','clot_clearance/tests/')) and name.endswith('.py'):
            assert hashlib.sha256((repo/name).read_bytes()).hexdigest()==record['sha256'],name
            checked+=1
    assert checked>20


def test_predeclared_physics_config_unchanged():
    a=audit_config(); root=OUT.parents[1]
    assert hashlib.sha256((root/a['physics_config_source']).read_bytes()).hexdigest()==a['physics_config_sha256']
    assert hashlib.sha256((OUT/'provenance/PREDECLARED_AUDIT_CONFIG.json').read_bytes()).hexdigest()==(OUT/'provenance/PREDECLARED_AUDIT_CONFIG.sha256').read_text().strip()
