import numpy as np
import pytest

from mevo_mb.bubbles.contact import constrain_pair_velocities, project_contacts
from mevo_mb.bubbles.forces import NumpyPairBackend, primary_bjerknes
from mevo_mb.bubbles.population import BubblePopulation, BubbleSizeDistribution
from mevo_mb.bubbles.response import LinearResonantBubbleModel, MarmottantBubbleModel, MarmottantParameters
from mevo_mb.coupling.geometry import Cylinder
from mevo_mb.flow.fields import PoiseuilleCylinderFlowField, SimVascularFlowField, ZeroFlowField
from mevo_mb.utils.config import load_config


def test_poiseuille_and_zero():
    x = np.array([[0., 0., 0.], [0., 1., 0.], [0., 0.5, 0.]])
    field = PoiseuilleCylinderFlowField(1., 2.)
    assert np.allclose(field.velocity(x, 0)[:, 0], [2, 0, 1.5])
    assert np.array_equal(ZeroFlowField().velocity(x, 0), x*0)


def test_primary_sign_scaling_phase_and_units():
    model = LinearResonantBubbleModel(2.5e6, 0.15)
    medium = {"density_kg_m3": 1060.}
    p = np.array([1e4+2e3j])
    g = np.array([[1e7+2e6j, 0, 0]])
    volume = model.response(np.array([1.25e-6]), 1e6, p, medium)
    force = primary_bjerknes(volume, g)
    assert force[0, 0] > 0  # below resonance: toward increasing pressure magnitude
    assert np.allclose(primary_bjerknes(volume, g*0), 0, atol=0)
    assert np.allclose(primary_bjerknes(3*volume, 3*g), 9*force, rtol=1e-12, atol=0)
    assert np.allclose(primary_bjerknes(volume*np.exp(0.7j), g*np.exp(0.7j)), force, rtol=1e-12, atol=0)
    assert np.allclose(primary_bjerknes(-volume, g), -force, rtol=1e-12, atol=0)
    # m^3 * (Pa/m) = N; explicit numerical value independent of the response model.
    assert primary_bjerknes(np.array([-2e-18]), np.array([[1e8, 0, 0]]))[0, 0] == 1e-10


def test_secondary_attraction_repulsion_reciprocity():
    backend = NumpyPairBackend(2e-6, 1e-8, 1e-3)
    x = np.array([[0., 0., 0.], [1e-5, 0., 0.]])
    r = np.full(2, 1e-6)
    v = np.full(2, 1e-18+0j)
    force = backend.forces(x, r, v, 1e6, 1000.)
    expected = 1000/(4*np.pi*1e-10)*0.5*(2*np.pi*1e6*1e-18)**2
    assert np.isclose(force[0, 0], expected, rtol=1e-12, atol=0)
    assert np.allclose(force[0], -force[1], rtol=1e-12, atol=0)
    v[1] *= -1
    assert backend.forces(x, r, v, 1e6, 1000.)[0, 0] < 0


def test_contact_and_initialization():
    cyl = Cylinder(1., 4.)
    x = np.array([[0., 0., 0.], [0.1, 0., 0.]])
    projected, _, events, _ = project_contacts(x, np.full(2, 0.1), cyl, 1e-12, 100)
    assert np.linalg.norm(projected[1]-projected[0]) >= 0.2-1e-12
    assert events == [(0, 1)]
    c = load_config("configs/demo_aggregation.yaml")
    a, b = BubblePopulation.initialize(c), BubblePopulation.initialize(c)
    assert np.array_equal(a.positions_m, b.positions_m)
    from scipy.spatial.distance import pdist
    assert pdist(a.positions_m).min() >= 2.5e-6
    d = dict(c["particles"]["size_distribution"], kind="lognormal")
    samples = BubbleSizeDistribution(d).sample(1000, np.random.default_rng(1))
    assert samples.min() >= d["minimum_radius_m"] and samples.max() <= d["maximum_radius_m"]


def test_synthetic_vtu_linear_velocity(tmp_path):
    import pyvista as pv
    points = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    grid = pv.UnstructuredGrid(np.array([4, 0, 1, 2, 3]), np.array([pv.CellType.TETRA]), points)
    grid["Velocity"] = points*np.array([2., 3., 4.])
    path = tmp_path/"synthetic.vtu"
    grid.save(path)
    flow = SimVascularFlowField(path)
    assert np.allclose(flow.velocity([[0.1, 0.2, 0.3]], 0), [[0.2, 0.6, 1.2]])
    with pytest.raises(ValueError, match="outside"):
        flow.velocity([[2., 2., 2.]], 0)


def test_unimplemented_shell_is_explicit():
    params = MarmottantParameters(1, 1e-8, 1e-6, 2e-6, 0.072, 1.4, 101325)
    with pytest.raises(NotImplementedError):
        MarmottantBubbleModel(params).response(1e-6, 1e6, 1e4, {})


def test_simultaneous_contact_chain():
    count = 30
    x = np.zeros((count, 3))
    x[:, 0] = np.arange(count)*2e-6
    radii = np.full(count, 1e-6)
    v = np.zeros_like(x)
    v[:, 0] = np.linspace(1., -1., count)
    projected = constrain_pair_velocities(x, radii, v, 1e-12, 100)
    assert np.max(np.abs(projected)) < 1e-12
    assert np.all(np.diff(projected[:, 0]) >= -1e-12)
    assert np.sum(projected**2) <= np.sum(v**2)
