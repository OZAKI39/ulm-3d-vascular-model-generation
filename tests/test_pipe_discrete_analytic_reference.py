import numpy as np
import pytest

from utils.cfd_flow.musubi_wall_force_diagnostics import discrete_poiseuille_reference


def test_discrete_reference_averages_actual_cell_centres() -> None:
    points = np.asarray([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]])
    result = discrete_poiseuille_reference(
        points,
        center_m=(0.0, 0.0, 0.0),
        axis=(0.0, 0.0, 1.0),
        radius_m=1.0,
        continuum_mean_m_s=2.0,
    )
    assert result["velocity_m_s"].tolist() == pytest.approx([4.0, 3.0])
    assert result["discrete_analytic_mean_m_s"] == pytest.approx(3.5)
