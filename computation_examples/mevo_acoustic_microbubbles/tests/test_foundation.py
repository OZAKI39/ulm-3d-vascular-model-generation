"""Small deterministic checks run before any production acoustics."""
import copy

import numpy as np
import pytest

from mevo_mb.acoustics.basis import AcousticBasisDatabase, ComplexAcousticField
from mevo_mb.control.pzt import PZTArray
from mevo_mb.coupling.geometry import Cylinder, grid_axes
from mevo_mb.utils.config import load_config, validate_config


def test_config_and_nonlinearity_guard():
    c = load_config("configs/demo_aggregation.yaml")
    assert c["acoustics"]["frequency_Hz"] == 1e6
    assert c["medium"]["viscosity_Pa_s"] == 0.0035
    d = copy.deepcopy(c)
    d["acoustics"]["BonA"] = 6
    with pytest.raises(ValueError, match="Nonlinear"):
        validate_config(d)


def test_cylinder_mask_and_finite_sphere_projection():
    cylinder = Cylinder(1.0, 4.0)
    axes = grid_axes([8, 8, 8], [0.5]*3)
    assert cylinder.mask(axes)[4, 4, 4]
    assert not cylinder.mask(axes)[4, 0, 0]
    x, hit = cylinder.project(np.array([[0., 1.2, 0.]]), np.array([0.1]))
    assert hit[0] and np.isclose(x[0, 1], 0.9)
    assert cylinder.contains(x, 0.1)[0]
    assert cylinder.constrain_velocity(x, np.array([0.1]), np.array([[1., 2., 0.]]))[0, 1] == 0


def test_pzt_controls():
    array = PZTArray.from_config(load_config("configs/demo.yaml"))
    array.set_amplitudes([1., 2., 3., 4.])
    array.set_phases([0., np.pi/2, np.pi, 0.])
    array.set_enabled([True, True, True, False])
    assert np.allclose(array.get_control_vector(), [1, 2j, -3, 0])


def test_complex_affine_interpolation_and_gradient(tmp_path):
    axes = tuple(np.linspace(-1, 1, 5) for _ in range(3))
    x, y, z = np.meshgrid(*axes, indexing="ij")
    p = (2+3j)*x + (4-1j)*y + z
    g = np.stack(np.gradient(p, *axes, edge_order=2), axis=-1)
    db = AcousticBasisDatabase(axes, np.array([p, 2*p]), np.array([g, 2*g]), {"frequency_Hz": 1e6})
    field = db.combine([1, 1j])
    query = np.array([[0.123, -0.231, 0.345]])
    exact = (1+2j) * ((2+3j)*query[:, 0] + (4-1j)*query[:, 1] + query[:, 2])
    assert np.allclose(field.interpolate_complex_pressure(query), exact)
    assert np.allclose(field.interpolate_complex_gradient(query), (1+2j)*np.array([2+3j, 4-1j, 1]))
    db.save_basis(tmp_path/"basis.h5")
    assert np.allclose(AcousticBasisDatabase.load_basis(tmp_path/"basis.h5").combine([1, 1j]).pressure_Pa, field.pressure_Pa)
    field.save(tmp_path/"combined.h5")
    loaded = ComplexAcousticField.load(tmp_path/"combined.h5")
    assert np.array_equal(loaded.pressure_Pa, field.pressure_Pa)
    assert np.array_equal(loaded.gradient_Pa_m, field.gradient_Pa_m)
    assert loaded.frequency_Hz == field.frequency_Hz
    with pytest.raises(ValueError):
        field.interpolate_complex_pressure([[2, 0, 0]])
