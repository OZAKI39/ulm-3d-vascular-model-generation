import numpy as np
import pytest
from particle_3d.field import FrozenFEMField


def test_input_and_output_arrays_cannot_mutate_sampler(affine):
    source = affine["field"]
    points = source.points_m.copy()
    velocity = source.velocity_nodes_m_s.copy()
    pressure = source.pressure_nodes_pa.copy()
    field = FrozenFEMField(points, source.tetra, velocity, pressure)
    position = affine["positions"][0]
    before = field.sample(position)
    points[:] = 999; velocity[:] = 999; pressure[:] = 999
    after = field.sample(position)
    np.testing.assert_array_equal(before.velocity_m_s, after.velocity_m_s)
    assert before.pressure_pa == after.pressure_pa
    for array in [field.points_m, field.tetra, field.velocity_nodes_m_s, field.pressure_nodes_pa,
                  after.velocity_m_s, after.velocity_gradient_s_inv, after.vorticity_s_inv, after.strain_rate_s_inv]:
        with pytest.raises(ValueError):
            array.flat[0] = 0
