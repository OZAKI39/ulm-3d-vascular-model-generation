import numpy as np
from pd_clot.pfpd_candidate.coupling_stability_audit import (
    OUT,audit_config,read_json,load_audit_state,deterministic_starts,projected_kkt)


def test_joint_objective_and_gradient_finite():
    for key in ['A','B','C']:
        p,d=load_audit_state(key); f,g=p.objective_gradient(p.pack(d['x'],d['s']))
        assert np.isfinite(f) and np.isfinite(g).all() and len(g)==p.nx+p.ns
        np.testing.assert_allclose(f*p.energy_scale,p.metrics(d['x'],d['s'])['total_J'],rtol=1e-14)


def test_predeclared_directional_gradient_gate_all_directions_and_epsilons():
    a=audit_config()
    for key in ['A','B','failed_load_initial']:
        r=read_json(OUT/f'monolithic/gradient_{key}.json')
        assert len(r['rows'])==a['gradient_directions']*len(a['gradient_epsilon'])
        assert set(x['epsilon'] for x in r['rows'])==set(a['gradient_epsilon'])
        assert max(x['normalized_error'] for x in r['rows'])<=a['gradient_relative_tolerance']


def test_force_sign_and_phase_gradient_against_live_directional_energy():
    p,d=load_audit_state('C'); z=p.pack(d['x'],d['s']); _,g=p.objective_gradient(z)
    # Displace a free coordinate and a phase bond separately away from stationarity.
    z[0]+=1e-3; z[p.nx+10]+=2e-3; _,g=p.objective_gradient(z)
    for j in [0,p.nx+10]:
        e=np.zeros(len(z)); e[j]=3e-6
        fd=(p.objective_gradient(z+e)[0]-p.objective_gradient(z-e)[0])/(6e-6)
        np.testing.assert_allclose(fd,g[j],rtol=1e-3,atol=3e-9)


def test_projected_kkt_distinguishes_lower_upper_interior_and_wrong_signs():
    s=np.array([.2,.2,.9,.9,.5]); q=np.array([.3,-.1,-.3,.1,.0])
    np.testing.assert_allclose(projected_kkt(s,q,np.full(5,.2),.9),[0,-.1,0,.1,0],atol=1e-15)


def test_joint_lower_bound_is_previous_load_not_inner_iterate():
    p,starts=deterministic_starts()
    np.testing.assert_array_equal(p.lower,starts['previous_accepted_with_affine_boundary_lift'][1])
    assert np.all(starts['failed_iteration_023_phase_after'][1]>p.lower)
    assert not p.lower.flags.writeable
    for x,s in starts.values():
        assert np.all(s>=p.lower) and np.all(s<p.sc)
        xx,ss=p.unpack(p.pack(x,s))
        np.testing.assert_array_equal(xx[p.cloud.fixed],p.boundary[p.cloud.fixed])
