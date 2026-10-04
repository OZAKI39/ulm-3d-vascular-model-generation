"""V2 geometry and audit-evidence regression checks, separate from scientific PASS gates."""
import hashlib
import numpy as np
from scripts.spatial_audit_v2 import (ROOT,OUT,L,DELTA,LEVELS,EPS,fixture,Problem,affine,
    evaluate,read,participation,hessian_operator,aggregate)

def test_coordinates_hit_exact_physical_boundaries():
    for n in LEVELS:
        X=fixture(n)[1].X
        assert np.array_equal(X.min(axis=0),np.full(3,-L/2))
        assert np.array_equal(X.max(axis=0),np.full(3,L/2))
        assert len(np.unique(X[:,0]))==n

def test_conservative_total_volume():
    for n in LEVELS:np.testing.assert_allclose(fixture(n)[1].volume.sum(),L**3,rtol=2e-15)

def test_surface_weights_conserve_both_clamp_areas():
    for n in LEVELS:
        c=fixture(n)[0]
        for side in [-1,1]:np.testing.assert_allclose(c.face_area[c.face_normal[:,2]==side].sum(),L**2,rtol=2e-15)

def test_physical_horizon_is_fixed():
    assert all(fixture(n)[1].horizon==DELTA for n in LEVELS)

def test_horizon_sampling_ratio_increases():
    np.testing.assert_allclose([DELTA/fixture(n)[0].spacing for n in LEVELS],[2.01,4.02,6.03],rtol=1e-15)

def test_Gc_not_recalibrated():
    expected=read(ROOT/'verification/pfpd_at1_evolution_v1/provenance/PREDECLARED_CONFIG.json')['phase']['Gc_demo_J_per_m2']
    assert all(Problem(n,.01).config['phase']['Gc_demo_J_per_m2']==expected for n in LEVELS)

def test_Yc_independent_of_spacing():
    yc=[Problem(n,.01).yc for n in LEVELS]
    assert len(set(yc))==1
    np.testing.assert_allclose(yc,.1010025,rtol=2e-15)

def test_affine_F_reproduction_all_levels():
    F=np.array([[1.01,.007,0],[0,1.002,.004],[0,0,.997]])
    for n in LEVELS:
        _,f,m,_=fixture(n);state=evaluate(f,f.X@F.T,m)
        assert np.max(abs(state.Fbar-F))<2e-12
        assert np.max(abs(state.Ftilde-F))<2e-12

def test_affine_energy_density_fixed_volume():
    for n in LEVELS:
        p=Problem(n,.01);state=p.state(affine(p),p.lower)
        np.testing.assert_allclose(state.energy/L**3,.1010025,rtol=2e-11)

def test_nominal_reaction_uses_physical_area_and_records_sampling_error():
    rows=read(OUT/'mechanics/AFFINE_REACTIONS.json')
    for r in rows:np.testing.assert_allclose(r['nominal_stress_Pa'],r['reaction_N']/L**2,rtol=1e-15)
    fine=[next(r for r in rows if r['level']==name and r['strain']==.01)['reaction_N'] for name in ['C5','C7']]
    assert fine[0]>0 and fine[1]>0
    # Fixture comparison is diagnostic, not an onset-consistency assertion.
    assert np.isfinite(abs(fine[1]-fine[0])/max(fine))

def test_corrected_C3_recomputes_prior_proposal():
    p=Problem(3,.01);st=p.state(affine(p),p.lower)
    np.testing.assert_allclose(st.energy,1.5781640625000852e-12,rtol=1e-10)
    np.testing.assert_allclose(-st.force[p.top,2].sum(),1.2688125000000346e-6,rtol=1e-10)

def test_legacy_files_are_unchanged():
    manifest=read(OUT/'provenance/PRE_RUN_MANIFEST.json')['files']
    prefixes=('results/pfpd_at1_evolution_v1/','verification/pfpd_at1_branch_stability_tracking_v1/','verification/pfpd_at1_spatial_discretization_audit_v1/')
    for name,digest in manifest.items():
        if name.startswith(prefixes):assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name

def test_C3_dense_Hv_validation_includes_required_vectors():
    r=read(OUT/'matrix_free/C3_DENSE_VS_HVP.json');assert r['passed'] and r['dense_reference_valid']
    assert sum(t['direction'].startswith('random_') for t in r['tests'])>=30
    assert {t['direction'] for t in r['tests']} >= {'lowest_coupled_mode','uniform_phase_opening','random_phase_opening'}

def test_matrix_free_three_epsilon_repeatability():
    p=Problem(3,.01,common_scale=True);z=p.pack(affine(p),p.lower);v=np.r_[np.zeros(p.nx),np.ones(p.ns)];v/=np.linalg.norm(v)
    values=[hessian_operator(p,z,e)@v for e in EPS]
    for a in values:np.testing.assert_allclose(a,values[1],rtol=5e-5,atol=2e-9)

def test_critical_cone_witness_phase_and_domain_feasibility():
    r=read(OUT/'critical_brackets/C3_A.json')['first_unstable'];w=r['negative_witness']
    assert w['phase_feasible'] and w['critical_cone']
    assert abs(w['linear_variation'])<1e-8 and any(p['feasible'] for p in w['probes'])

def test_negative_witness_repeats_on_both_schedules_and_all_eps():
    for schedule in ['A','B']:
        w=read(OUT/f'critical_brackets/C3_{schedule}.json')['first_unstable']['negative_witness']
        assert w['verified_negative']
        assert all(v<0 for v in w['curvatures']+w['matrix_free_curvatures'])

def test_accepted_stationary_states_obey_all_original_gates():
    for path in (OUT/'branch_tracking').glob('*/STATIONARY.json'):
        r=read(path)
        if r['accepted']:
            assert all(r['gates'].values())
            assert r['metrics']['force_residual']<=1e-9
            assert r['metrics']['projected_phase_KKT']<=1e-8

def test_active_set_counts_partition_every_bond():
    for path in (OUT/'active_sets').glob('*.json'):
        r=read(path);ids=[i for v in r['bond_IDs'].values() for i in v]
        assert len(ids)==len(set(ids))==sum(r['counts'].values())
        np.testing.assert_allclose(sum(r['fractions'].values()),1,rtol=1e-15)

def test_participation_ratio_known_limits():
    assert participation(np.ones(100))==dict(PR=1.,top10_fraction=.1)
    assert participation(np.r_[1.,np.zeros(99)])==dict(PR=.01,top10_fraction=1.)
    assert participation(np.zeros(3))['PR'] is None

def test_legacy_excluded_from_refinement_data():
    rows=read(OUT/'geometry/LEVELS.json')
    assert [r['name'] for r in rows]==[f'corrected_boundary_node_C{n}' for n in LEVELS]
    assert all(abs(r['volume_m3']/L**3-1)<2e-15 for r in rows)

def test_historical_sources_and_tests_unchanged():
    manifest=read(OUT/'provenance/PRE_RUN_MANIFEST.json')['files']
    for name,digest in manifest.items():
        if name.startswith(('pd_clot/','tests/','scripts/')):
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name

def test_endpoint_face_edge_corner_weights():
    for n in LEVELS:
        c,f,_,_=fixture(n);touch=(abs(f.X)==L/2).sum(axis=1)
        np.testing.assert_allclose(f.volume,c.spacing**3/2.**touch,rtol=2e-15)

def test_symmetric_quadrature_first_moments_and_linear_surface_integrals():
    for n in LEVELS:
        c=fixture(n)[0]
        assert np.max(abs(c.volume@c.X))<1e-29
        for side in [-1,1]:
            sel=c.face_normal[:,2]==side
            np.testing.assert_allclose(c.face_area[sel]@c.face_position[sel,2],side*L**3/2,rtol=2e-15)

def test_rank_gate_passes_without_regularization():
    for r in read(OUT/'geometry/LEVELS.json'):assert r['rank_histogram']=={'3':r['particles']}

def test_C7_dense_memory_allocation_forbidden():
    r=read(OUT/'corrected_C7/MEMORY_ESTIMATE.json')
    assert not r['dense_permitted'] and r['dense_Hessian_bytes']>.4*r['memory']['available']

def test_fixed_version_external_oracle_all_levels():
    r=read(OUT/'external_reference/COMPARISON.json')
    assert r['status']=='VERIFIED'
    assert r['commit']=='b6e86563845a83b5cb06c805af19a817bde6b106'
    assert len(r['cases'])==6 and all(t['passed'] for t in r['cases'])

def test_analytic_threshold_mesh_independence():
    rows=read(OUT/'geometry/LEVELS.json')
    assert np.ptp([r['analytic_threshold'] for r in rows])<1e-13

def test_bound_aware_operator_never_samples_negative_phase():
    p=Problem(3,.01,common_scale=True);z=p.pack(affine(p),p.lower);original=p.objective_gradient;seen=[]
    def checked(q):
        assert np.all(q[p.nx:]>=p.lower) and np.all(q[p.nx:]<p.sc)
        seen.append(q.copy());return original(q)
    p.objective_gradient=checked
    v=np.random.default_rng(951).normal(size=len(z));hessian_operator(p,z,3e-6)@v
    assert len(seen)>=5

def test_weighted_bond_aggregation_preserves_constant_amplitude():
    f=fixture(3)[1];np.testing.assert_allclose(aggregate(f,np.ones(len(f.pairs))),1,rtol=1e-15)

def test_refined_nonstationary_affine_state_not_accepted_as_stability_endpoint():
    for n in [5,7]:
        p=Problem(n,.009);assert p.metrics(affine(p),p.lower)['force_residual']>1e-9
        for name in ['A','B']:
            r=read(OUT/f'critical_brackets/C{n}_{name}.json')
            for key in ['last_stable','first_unstable']:
                if r[key] is not None:assert r[key]['stationary'] and all(r[key]['gates'].values())
