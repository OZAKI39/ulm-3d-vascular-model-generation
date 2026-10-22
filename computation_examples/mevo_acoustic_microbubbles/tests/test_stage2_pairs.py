import numpy as np
import pytest

from mevo_mb.bubbles.forces import ClassicalFarFieldBjerknes
from mevo_mb.bubbles.validity import PhysicsValidityMonitor, PhysicsValidityError
from mevo_mb.validation.zhang2011 import (
    Zhang2011EncapsulatedPairReference,
    Zhang2011Parameters,
)
from mevo_mb.validation.doinikov_bouakaz_pair import (
    DoinikovBouakazPairReferenceModel,
)


def test_far_field_harmonics_and_rejection():
    x = np.array([[0.0, 0.0, 0.0], [30e-6, 0.0, 0.0]])
    r = np.full(2, 1e-6)
    v = np.array([[1e-18, 2e-19j], [1e-18, 2e-19j]])
    f = ClassicalFarFieldBjerknes().forces(x, r, v, 1e6, 1000.0)
    expected = 1000 / (4 * np.pi * (30e-6) ** 2) * 0.5 * (2 * np.pi * 1e6) ** 2 * (1e-36 + 4 * 4e-38)
    np.testing.assert_allclose(f[:, 0], [expected, -expected], rtol=2e-14, atol=0)
    with pytest.raises(ValueError, match="far-field"):
        ClassicalFarFieldBjerknes().forces(x * 0.1, r, v, 1e6, 1000.0)


def test_validity_strict_and_legacy_labels():
    x = np.array([[0.0, 0.0, 0.0], [3e-6, 0.0, 0.0]])
    for policy in ("strict", "legacy_demo"):
        m = PhysicsValidityMonitor("table", policy=policy)
        def call():
            return m.observe(
                    x, np.full(2, 1e-6), 0.0, radial_fraction=np.zeros(2), frequency=1e6, sound_speed=1500.0
                )
        if policy == "strict":
            with pytest.raises(PhysicsValidityError):
                call()
        else:
            call()
        assert m.summary()["near_field_events"] == 1
        assert m.summary()["physics_status"] == ("outside validity" if policy == "strict" else "demo only")


def test_doinikov_mass_limits_relabel_symmetry():
    m = DoinikovBouakazPairReferenceModel(multipole_order=8)
    q = np.array([1.5e-6, 2e-6, 0.0, 1.0])
    mass = m.mass_matrix(q)
    expected = 4 * np.pi * m.p.density_kg_m3 * np.r_[q[:2] ** 3, q[:2] ** 3 / 6]
    np.testing.assert_allclose(np.diag(mass), expected, rtol=1e-12, atol=0)
    q[3] = 7e-6
    mass = m.mass_matrix(q)
    assert np.all(np.linalg.eigvalsh(mass) > 0)
    # Exchange bubble labels AND reflect the common z-axis.
    swapped = m.mass_matrix(np.array([q[1], q[0], 0.0, q[3]]))
    perm = [1, 0, 3, 2]
    sign = np.array([1, 1, -1, -1])
    np.testing.assert_allclose(
        swapped, mass[np.ix_(perm, perm)] * sign[:, None] * sign[None, :], rtol=1e-10, atol=1e-28
    )


@pytest.mark.literature
@pytest.mark.deprecated_reference
def test_zhang_two_initial_translation_signs():
    distances = []
    for r in (3e-6, 4e-6):
        p = Zhang2011Parameters(radius1_m=r)
        z = Zhang2011EncapsulatedPairReference(p).solve(cycles=24)
        distances.append(np.mean(z["distance_m"][-64:]))
    assert distances[0] > 50e-6 > distances[1]


@pytest.mark.gpu
def test_cupy_same_harmonic_physics():
    cp = pytest.importorskip("cupy")
    if cp.cuda.runtime.getDeviceCount() == 0:
        pytest.skip("No CUDA device")
    rng = np.random.default_rng(4)
    x = np.column_stack((np.arange(64) * 25e-6, np.zeros((64, 2))))
    r = np.full(64, 1e-6)
    v = (rng.normal(size=(64, 4)) + 1j * rng.normal(size=(64, 4))) * 1e-18
    a = ClassicalFarFieldBjerknes().forces(x, r, v, 1e6, 1000.0)
    b = ClassicalFarFieldBjerknes("cupy").forces(x, r, v, 1e6, 1000.0).get()
    np.testing.assert_allclose(a, b, rtol=2e-12, atol=1e-25)


def test_time_varying_vtu_and_wall_refusal(tmp_path):
    import pyvista as pv
    from mevo_mb.flow.fields import VTKFlowField
    from mevo_mb.bubbles.wall import DoinikovWallReference

    points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    paths = []
    for i in (0, 1):
        grid = pv.UnstructuredGrid(np.array([4, 0, 1, 2, 3]), np.array([pv.CellType.TETRA]), points.copy())
        grid["flow"] = points * (1 + i)
        path = tmp_path / f"{i}.vtu"
        grid.save(path)
        paths.append(path)
    flow = VTKFlowField(paths, [0.0, 2.0], velocity_array="flow")
    np.testing.assert_allclose(flow.velocity([[0.1, 0.2, 0.3]], 1.0), [[0.15, 0.3, 0.45]])
    with pytest.raises(ValueError, match="time outside"):
        flow.velocity([[0.1, 0.2, 0.3]], 3.0)
    with pytest.raises(NotImplementedError):
        DoinikovWallReference().apply(points, np.ones(4), points, 0.0)
