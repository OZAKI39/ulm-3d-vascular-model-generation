import json
import numpy as np
import pytest
from particle_3d.particle_shapes import Sphere,Ellipsoid,Capsule,unit,EPS
from particle_3d.hydrodynamic_resistance import PhysicalNearField,sphere_self_diagonal,viscosity_from_frozen
from particle_3d.resistance_assembly import assemble_resistance_system
from particle_3d.resistance_solver import solve_resistance,ResistanceSystemIllConditioned
from particle_3d.particle5_cases import scalar_pair_solution,sparse_scene,algebraic_pair_block


@pytest.mark.parametrize('radius',[.2e-6,1e-6,5e-6])
@pytest.mark.parametrize('viscosity',[.001,.00345312,.02])
def test_sphere_self_stokes_translation_and_rotation(radius,viscosity):
    d=sphere_self_diagonal(Sphere([0,0,0],radius),viscosity)
    np.testing.assert_array_equal(d[:3],np.full(3,6*np.pi*viscosity*radius))
    np.testing.assert_array_equal(d[3:],np.full(3,8*np.pi*viscosity*radius**3))


@pytest.mark.parametrize('seed',range(5))
def test_isolated_sphere_free_motion(seed,mu):
    rng=np.random.default_rng(seed);free=rng.normal(size=6)
    system=assemble_resistance_system({9:Sphere([0,0,0],(seed+1)*1e-6)},{9:free},mu=mu)
    np.testing.assert_allclose(solve_resistance(system).velocity,free,rtol=64*EPS,atol=0)
    zero=assemble_resistance_system({9:Sphere([0,0,0],1e-6)},{9:np.zeros(6)},mu=mu)
    np.testing.assert_array_equal(solve_resistance(zero).velocity,0)


@pytest.mark.parametrize('ratio',[.1,.01,.001,.0001])
@pytest.mark.parametrize('direction',[[0,0,1],[1,2,-3]])
def test_sphere_wall_lubrication_scaling_velocity_analytic_and_tangent(ratio,direction,mu):
    a=1e-6;h=a*ratio;n=unit(direction);free=np.array([2e-6,-3e-6,4e-6,1.,2.,3.])
    system=assemble_resistance_system({1:Sphere([0,0,0],a)},{1:free},[PhysicalNearField(1,None,h,n,0,'VALIDATION_GAPS_ONLY')],mu=mu)
    z=system.blocks[0]['coefficient_kg_s']
    assert abs(z*h/(6*np.pi*mu*a*a)-1)<8*EPS
    solved=solve_resistance(system);expected=free.copy();expected[:3]-=(free[:3]@n)*n*a/(a+h)
    budget=64*EPS*solved.record['condition_estimate']*np.linalg.norm(free[:3])
    np.testing.assert_allclose(solved.velocity[:3],expected[:3],rtol=0,atol=budget)
    t=np.eye(3)-np.outer(n,n)
    np.testing.assert_allclose(t@solved.velocity[:3],t@free[:3],rtol=0,atol=budget)
    np.testing.assert_allclose(solved.velocity[3:],free[3:],rtol=64*EPS)


@pytest.mark.parametrize('radii',[(1e-6,1e-6),(.5e-6,2e-6)])
@pytest.mark.parametrize('ratio',[.1,.01,.001,.0001])
def test_unequal_sphere_pair_scaling_and_symmetry(radii,ratio,mu):
    ai,aj=radii;reff=ai*aj/(ai+aj);h=ratio*reff
    shapes={3:Sphere([0,0,0],ai),8:Sphere([ai+aj+h,0,0],aj)}
    free={3:np.array([1e-6,2e-6,0,1,2,3]),8:np.array([-2e-6,2e-6,0,4,5,6])}
    p=PhysicalNearField(3,8,h,[-1,0,0],0,'VALIDATION_GAPS_ONLY')
    system=assemble_resistance_system(shapes,free,pair_blocks=[p],mu=mu);solved=solve_resistance(system)
    reverse=assemble_resistance_system(dict(reversed(list(shapes.items()))),free,pair_blocks=[PhysicalNearField(8,3,h,[1,0,0],0,'VALIDATION_GAPS_ONLY')],mu=mu)
    np.testing.assert_array_equal(system.matrix.toarray(),reverse.matrix.toarray())
    np.testing.assert_array_equal(solved.velocity,solve_resistance(reverse).velocity)
    assert abs(system.blocks[0]['coefficient_kg_s']*h/(6*np.pi*mu*reff**2)-1)<8*EPS
    expected=scalar_pair_solution(ai,aj,mu,h,1e-6,-2e-6)
    np.testing.assert_allclose(solved.velocity[[0,6]],expected,rtol=0,atol=64*EPS*solved.record['condition_estimate']*3e-6)
    pair=system.matrix.toarray()-np.diag(system.self_diagonal)
    np.testing.assert_allclose((pair@solved.velocity)[:3],-(pair@solved.velocity)[6:9],rtol=0,atol=64*EPS*np.linalg.norm(pair)*np.linalg.norm(solved.velocity))


def test_pair_common_translation_null_mode(mu):
    a=Sphere([0,0,0],1e-6);b=Sphere([2.0001e-6,0,0],1e-6)
    free={1:np.array([1.,2,3,4,5,6]),2:np.array([1.,2,3,7,8,9])}
    system=assemble_resistance_system({1:a,2:b},free,pair_blocks=[PhysicalNearField(1,2,1e-10,[1,2,3],0)],mu=mu)
    np.testing.assert_allclose(solve_resistance(system).velocity,np.r_[free[1],free[2]],rtol=1e-11)
    assert system.block_rows[0]@system.free == 0.


@pytest.mark.parametrize('count',[5,12,20])
def test_sparse_all_pairs_equivalence_and_matrix_symmetry(count,mu):
    _,a,candidates=sparse_scene(mu,count);_,b,all_pairs=sparse_scene(mu,count,candidates=False)
    assert len(candidates)<len(all_pairs)
    np.testing.assert_array_equal(a.matrix.toarray(),b.matrix.toarray())
    np.testing.assert_array_equal(a.matrix.toarray(),a.matrix.toarray().T)
    np.testing.assert_array_equal(solve_resistance(a).velocity,solve_resistance(b).velocity)


def test_resistance_positive_dissipation_many_vectors(mu):
    _,s,_=sparse_scene(mu);rng=np.random.default_rng(2026092506)
    for _ in range(512):
        v=rng.normal(size=len(s.free))*np.tile([1e-5]*3+[10.]*3,len(s.ids))
        assert min(s.dissipation(v).values())>=0
        assert v@(s.matrix@v)>=0
    assert solve_resistance(s).record['minimum_scaled_eigenvalue']>0


def test_scaled_vs_unscaled_unchanged_solution(mu):
    _,s,_=sparse_scene(mu,5)
    np.testing.assert_allclose(solve_resistance(s).velocity,solve_resistance(s,scaled=False).velocity,rtol=5e-12,atol=1e-20)


def test_extreme_condition_reports_failure_without_floor(mu):
    s=assemble_resistance_system({1:Sphere([0,0,0],1.)},{1:np.ones(6)},[PhysicalNearField(1,None,1e-18,[0,0,1],0,'VALIDATION_GAPS_ONLY')],mu=mu)
    assert s.blocks[0]['gap_m']==1e-18
    with pytest.raises(ResistanceSystemIllConditioned,match='RESISTANCE_SYSTEM_ILL_CONDITIONED'):
        solve_resistance(s)


def test_negative_eigenvalue_is_rejected(mu):
    s=assemble_resistance_system({1:Sphere([0,0,0],1.)},{1:np.ones(6)},mu=mu)
    s.matrix[0,0]*=-1
    with pytest.raises(ResistanceSystemIllConditioned):solve_resistance(s)


@pytest.mark.parametrize('gap,roundoff,reason,active',[(1e-8,1e-17,'VALIDATION_NEAR_FIELD',True),(2e-8,1e-17,'OUTSIDE_VALIDATION_ASYMPTOTIC_RANGE',False),(1e-17,1e-17,'CONTACT_REGIME',False),(0.,1e-17,'CONTACT_REGIME',False),(-1e-18,1e-17,'CONTACT_REGIME',False)])
def test_activation_no_gap_floor(gap,roundoff,reason,active,mu):
    s=assemble_resistance_system({1:Sphere([0,0,0],1e-6)},{1:np.ones(6)},[PhysicalNearField(1,None,gap,[0,0,1],roundoff)],mu=mu)
    assert s.blocks[0]['reason']==reason and s.blocks[0]['active']==active
    assert s.blocks[0]['gap_m']==gap


@pytest.mark.parametrize('shape',[Ellipsoid([0,0,0],[2e-6,2e-6,1e-6],np.eye(3)),Capsule([0,0,0],[0,0,1],1e-6,2e-6)])
def test_nonspherical_lubrication_not_frozen(shape,mu):
    from particle_3d.particle5_motion import assemble_scene
    with pytest.raises(ValueError,match='NONSPHERICAL_LUBRICATION_NOT_FROZEN') as exc:
        assemble_scene({1:shape},{1:np.zeros(6)},mu)
    assert exc.value.record['reason']=='NONSPHERICAL_MODEL_NOT_FROZEN'
    assert not exc.value.record['active'] and exc.value.record['relevant_radius_m'] is None


def test_rbc_mb_algebraic_block_symmetry():
    r=algebraic_pair_block();permutation=np.r_[np.arange(6,12),np.arange(6)]
    np.testing.assert_array_equal(r,r.T)
    np.testing.assert_array_equal(r,r[np.ix_(permutation,permutation)])
    assert np.linalg.eigvalsh(r).min()>=-64*EPS*np.linalg.norm(r,2)
    for x in np.random.default_rng(17).normal(size=(128,12)):
        assert x@r@x>=-64*EPS*np.linalg.norm(x)**2*np.linalg.norm(r,2)


def test_validation_block_cannot_enter_physical_solver(mu):
    with pytest.raises(TypeError,match='VALIDATION_BLOCK_CANNOT_ENTER_PHYSICAL_SOLVER'):
        assemble_resistance_system({1:Sphere([0,0,0],1e-6)},{1:np.zeros(6)},pair_blocks=[algebraic_pair_block()],mu=mu)


def test_viscosity_provenance(repo,mu):
    _,p=viscosity_from_frozen(repo/'formal_3D_flow_solver/FEM_SimVascular')
    assert mu==.00345312 and p['status']=='PASS' and len(p['sources'])==2


def test_viscosity_mismatch_stops(tmp_path):
    root=tmp_path/'frozen_reference';(root/'run').mkdir(parents=True)
    (root/'baseline_summary.json').write_text(json.dumps({'physics':{'density_kg_m3':1056,'kinematic_viscosity_m2_s':3.27e-6,'dynamic_viscosity_pa_s':.004}}))
    (root/'run/solver.xml').write_text('<root><Add_equation><Density>1056</Density><Viscosity><Value>0.004</Value></Viscosity></Add_equation></root>')
    with pytest.raises(ValueError,match='VISCOSITY_PROVENANCE_MISMATCH'):viscosity_from_frozen(tmp_path)
