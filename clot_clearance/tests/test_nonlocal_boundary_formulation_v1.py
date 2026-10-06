"""Independent geometry, continuum work and force-interface regression checks."""
import numpy as np
import pytest

from pd_clot.pfpd_candidate.baqp import evaluate
from pd_clot.pfpd_candidate.energy import Material
from pd_clot.pfpd_candidate.kernel import Family
from pd_clot.pfpd_candidate.nonlocal_boundary import (
    BlockBoundary, AffineSurrogate, block_family, convert_volume_force)
from pd_clot.pfpd_candidate.nonlocal_boundary.diagnostics import load_checks
from pd_clot.pfpd_candidate.nonlocal_boundary.geometry import FACE_NAMES
from pd_clot.pfpd_candidate.nonlocal_boundary.local_surrogate import pressure_shear_traction
from pd_clot.pfpd_candidate.nonlocal_boundary.mechanics_adapter import IntactForce, solve
from pd_clot.pfpd_candidate.nonlocal_boundary.volume_dirichlet_reference import constrained_support, impose


@pytest.fixture(scope='module')
def block():
    return block_family(7, 6*62.5e-6, 2.1*62.5e-6)


def affine(kind):
    F = np.eye(3)
    F[(2, 2) if kind == 'normal' else (0, 2)] += .003
    return AffineSurrogate(F, Material(model='svk', mu=1000., bulk=2e3/3))


def test_boundary_normal_projection_and_partition():
    g = BlockBoundary(np.array([1., 2., 3.]), np.array([3., 6., 9.]), ('x-',))
    X = np.array([[1., 4., 6.], [2., 4., 8.5]])
    projected, normals, face = g.project(X)
    np.testing.assert_allclose(projected, [[1., 4., 6.], [2., 4., 9.]])
    np.testing.assert_array_equal(face, [0, 5])
    np.testing.assert_allclose(normals, [[-1., 0., 0.], [0., 0., 1.]])
    for name, n in zip(FACE_NAMES, g.normals):
        t1, t2 = g.tangent_basis(name)
        np.testing.assert_allclose(np.cross(t1, t2), n)
    assert g.dirichlet_mask(X).tolist() == [True, False]
    empty = BlockBoundary(g.lower, g.upper, ())
    assert empty.dirichlet_mask(X).tolist() == [False, False]
    masks = g.masks(np.array([[.9, 4., 6.], [2., 4., 6.]]), .5)
    assert masks['exterior_interaction'].tolist() == [True, False]


@pytest.mark.parametrize('kind', ['normal', 'shear'])
def test_affine_continuum_face_quadrature(block, kind):
    g, f = block; u = affine(kind); P = u.response()['P']
    load = u.surface_nodal_load(g, f.X)
    np.testing.assert_allclose(load.sum(axis=0), 0., atol=1e-20)
    np.testing.assert_allclose(load.T@(f.X-g.center), g.volume*P, rtol=1e-12, atol=1e-24)
    # Check every face independently, including lateral shear components.
    for q, exact in zip(g.face_quadrature(f.X), u.exact_face_integrals(g)):
        np.testing.assert_allclose(q['area'].sum(), exact['area_m2'], rtol=1e-14)
        np.testing.assert_allclose(u.surface_nodal_load(g, f.X, (q['face'],)).sum(axis=0), exact['force'], rtol=1e-13, atol=1e-23)


@pytest.mark.parametrize('kind', ['normal', 'shear'])
def test_volume_force_resultant_moment_and_work(block, kind):
    g, f = block; u = affine(kind)
    state = evaluate(f, u.position(f.X), u.material)
    converted = convert_volume_force(g, f, u, state)
    checks = load_checks(g, f.X, converted.nodal_force_N, u)
    assert checks['resultant_PASS'] and checks['moment_PASS'] and checks['virtual_work_PASS']
    assert len(checks['work']) == 12
    np.testing.assert_allclose(converted.force_density_N_m3*f.volume[:, None], converted.nodal_force_N)
    np.testing.assert_allclose(state.force+converted.nodal_force_N,
                               converted.omitted_internal_force_N, rtol=0., atol=0.)
    assert not converted.nodal_force_N.flags.writeable
    # Rigid translations / net zero do not suffice: zero loading fails strain work.
    empty = load_checks(g, f.X, np.zeros_like(f.X), u)
    assert empty['resultant_PASS'] and not empty['virtual_work_PASS']


def test_shear_reference_moment_and_lateral_traction(block):
    g, f = block; u = affine('shear'); P = u.response()['P']
    expected = 1000.*.003**3*g.volume
    moment = sum(row['moment'] for row in u.exact_face_integrals(g))
    np.testing.assert_allclose(moment, [0., expected, 0.], rtol=1e-10, atol=1e-25)
    assert abs(P[2, 0]) > 1.
    side = u.surface_nodal_load(g, f.X, ('x-', 'x+'))
    check = load_checks(g, f.X, side, u)
    assert not check['virtual_work_PASS']  # Side faces alone also incomplete.


def test_volume_dirichlet_preserves_constraint_without_new_material(block):
    g, f = block; u = affine('normal'); fixed, prescribed = constrained_support(g, f, u)
    assert fixed.any() and (~fixed).any()
    x = f.X.copy(); x += 1e-7
    result = impose(x, fixed, prescribed)
    np.testing.assert_array_equal(result[fixed], prescribed[fixed])
    np.testing.assert_array_equal(result[~fixed], x[~fixed])
    np.testing.assert_allclose(f.volume.sum(), g.volume)


def test_particle_relabeling_changes_no_physical_load(block):
    g, f = block; u = affine('shear')
    a = convert_volume_force(g, f, u)
    perm = np.random.default_rng(51).permutation(len(f.X)); inv = np.argsort(perm)
    pairs = np.sort(inv[f.pairs], axis=1)
    other = Family.from_arrays(f.X[perm], f.volume[perm], pairs, f.horizon)
    b = convert_volume_force(g, other, u)
    assert np.linalg.norm(b.nodal_force_N[inv]-a.nodal_force_N)/np.linalg.norm(a.nodal_force_N) < 1e-10


def test_matrix_free_tangent_is_force_derivative_and_reciprocal(block):
    _, f = block; u = affine('shear'); x = u.position(f.X)
    callback = IntactForce(f, u.material, np.zeros_like(x)); callback(x)
    state = callback.state; rng = np.random.default_rng(33)
    v, w = rng.normal(size=(2, *x.shape))*1e-5
    jv, jw = callback.tangent(v, state), callback.tangent(w, state)
    original = IntactForce(f, u.material, np.zeros_like(x), assembly='numpy')
    original(x)
    np.testing.assert_allclose(jv, original.tangent(v), rtol=1e-11, atol=1e-20)
    step = 1e-5
    fd = (evaluate(f, x+step*v, u.material).force-evaluate(f, x-step*v, u.material).force)/(2*step)
    assert np.linalg.norm(jv-fd)/np.linalg.norm(fd) < 1e-7
    assert abs(np.sum(w*jv)-np.sum(v*jw))/max(abs(np.sum(w*jv)), 1e-30) < 1e-10


def test_existing_solver_adapter_returns_valid_summary():
    g, f = block_family(3, 125e-6, 130e-6); u = affine('normal')
    x = u.position(f.X); load = convert_volume_force(g, f, u).nodal_force_N
    cfg = dict(solver_max_wall_seconds_per_case=60, solver_max_nfev=10,
               mechanical_residual_gate=1e-8, solver_lsmr_maxiter=30)
    result, info = solve(f, u.material, x, g.dirichlet_mask(f.X), x, load, 62.5e-6,
                         np.linalg.norm(u.response()['P'])*125e-6**2, cfg)
    assert info['converged'] and info['no_new_solver']
    assert info['final_jacobian_condition'] is None
    np.testing.assert_allclose(result, x, rtol=0., atol=1e-18)


def test_pressure_shear_input_rejects_normal_shear():
    np.testing.assert_allclose(pressure_shear_traction(2., [0., 0., 1.], [3., 0., 0.]), [3., 0., -2.])
    with pytest.raises(ValueError):
        pressure_shear_traction(2., [0., 0., 1.], [0., 0., 1.])
