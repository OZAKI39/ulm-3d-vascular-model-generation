"""Candidate-only domain admission and preserved historical rejection."""
import hashlib
import numpy as np
import pytest
from scripts.spectral_split_audit import *

def test_mixed_state_probe_is_admitted_but_old_candidate_rejects():
    from pd_clot.pfpd_candidate.adr_adapter_at1 import tensile_domain
    p=SpectralProblem(3,.005);x=v2.affine(p);x[p.free,0]*=.999;z=p.pack(x,p.lower)
    state=p.state(x,p.lower);assert ca.green_minimum(state)<-1e-9
    new=mixed_probe_point(p,z,np.zeros_like(z),0.)
    assert new['feasible'] and new['green_strain_policy']=='DIAGNOSTIC_ONLY'
    with pytest.raises(ValueError,match='PURE_TENSION_DOMAIN_VIOLATION'):tensile_domain(state,-1e-9)
    with pytest.raises(tracker.BranchDomainExit,match='PURE_TENSION_DOMAIN_VIOLATION'):tracker.stationary_residual(p,z)

def test_split_metrics_keep_negative_strain_and_energy_fractions():
    p=SpectralProblem(3,.005);x=v2.affine(p);x[p.free,0]*=.999;m=p.metrics(x,p.lower)
    assert m['min_principal_green_strain']<0 and m['negative_principal_fraction']>0
    assert m['Psi_minus_fraction']>0
    np.testing.assert_allclose(m['Psi_plus_fraction']+m['Psi_minus_fraction'],1,rtol=1e-15)

def test_candidate_cache_is_invalidated_by_phase_or_displacement():
    p=SpectralProblem(3,.005);x=v2.affine(p);a=p.state(x,p.lower)
    assert p.state(x.copy(),p.lower.copy()) is a
    b=p.state(x,np.full(p.ns,.01));assert b is not a and b.energy<a.energy
    xx=x.copy();xx[p.free,2]+=1e-9;assert p.state(xx,np.full(p.ns,.01)) is not b

def test_all_preexisting_production_and_candidate_sources_unchanged():
    manifest=read(OUT/'provenance/PRE_RUN_MANIFEST.json')['files']
    for name,digest in manifest.items():
        if name.startswith(('pd_clot/','scripts/','tests/')):
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
