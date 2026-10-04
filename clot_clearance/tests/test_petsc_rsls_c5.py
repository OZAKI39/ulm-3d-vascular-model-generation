import hashlib
import inspect
from scripts import petsc_vi_audit as audit
from pd_clot.pfpd_candidate import petsc_vi_candidate as vi


def all_roots():return audit.read(audit.OUT/'continuation/PETSC_RSLS/ALL.json')['states']


def test_first_activation_discovered_automatically():
    r=all_roots()[0];a=audit.prev.ca.read_npz(audit.OUT/'continuation/PETSC_RSLS/first_A1/STATE.npz')
    assert (a['z_initial'][225:]==0).all()
    assert r['interior_IDs']==[4488,5059] and r['scientifically_accepted']
    source=inspect.getsource(vi)
    assert '4488' not in source and '5059' not in source


def test_original_failure_load_recovered():
    r=next(r for r in all_roots() if r['label']=='same_load_as_original')
    assert r['scientifically_accepted'] and r['agrees'] and r['reason']>0


def test_all_seven_roots_and_inactive_inequalities():
    rows=all_roots();assert len(rows)==7
    for r in rows:
        assert r['scientifically_accepted'] and r['agrees']
        assert r['metrics']['force_residual']<=1e-9
        assert r['metrics']['projected_phase_KKT']<=1e-8
        assert r['gates']['inactive_inequalities'] and r['gates']['inactive_exact']


def test_legacy_projected_and_production_sources_unchanged():
    m=audit.read(audit.OUT/'provenance/PRE_RUN_MANIFEST.json')
    # Inventory schema is inspected by the audit; compare every old Python source.
    for name,digest in m['files'].items():
        if name.endswith('.py') and ('pd_clot/' in name or 'scripts/' in name or 'tests/' in name):
            path=audit.ROOT.parent/name
            if not path.exists():path=audit.ROOT/name
            assert hashlib.sha256(path.read_bytes()).hexdigest()==digest
