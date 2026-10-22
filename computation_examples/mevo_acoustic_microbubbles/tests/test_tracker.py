import h5py
import numpy as np
import pytest

from mevo_mb.acoustics.basis import ComplexAcousticField
from mevo_mb.bubbles.population import BubblePopulation
from mevo_mb.bubbles.tracker import ParticleTracker
from mevo_mb.flow.fields import PoiseuilleCylinderFlowField
from mevo_mb.io.artifacts import new_run
from mevo_mb.utils.config import load_config


@pytest.mark.parametrize("integrator", ["RK2", "RK4"])
def test_rk_flow_exact_and_complete_recording(tmp_path, integrator):
    c = load_config("configs/demo.yaml")
    c["particles"].update(count=3, duration_s=0.001, dt_s=0.0001, primary_enabled=False, secondary_enabled=False, integrator=integrator)
    axes = (np.linspace(-0.011, 0.011, 3), np.linspace(-0.002, 0.002, 3), np.linspace(-0.002, 0.002, 3))
    field = ComplexAcousticField(axes, np.zeros((3, 3, 3), complex), np.zeros((3, 3, 3, 3), complex), 2e6)
    flow = PoiseuilleCylinderFlowField(0.0015, 0.02)
    population = BubblePopulation.initialize(c)
    run = new_run(tmp_path, "test", c)
    summary = ParticleTracker(c, field, flow).run(population, run)
    with h5py.File(run/"particles/trajectories.h5") as f:
        exact = population.positions_m + 0.001*flow.velocity(population.positions_m, 0)
        assert np.allclose(f["positions_m"][-1], exact, atol=1e-14, rtol=0)
        assert len(f["time_s"]) == summary["accepted_substeps"]+1
        assert f["time_s"][-1] == 0.001
