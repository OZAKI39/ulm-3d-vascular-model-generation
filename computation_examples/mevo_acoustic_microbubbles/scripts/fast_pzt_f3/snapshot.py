"""Evaluate exact production force/mobility prefix without advancing particle positions."""

import numpy as np
from mevo_mb.fast_sim import kernels as k


def evaluate(sim, event):
    sim.reset()
    sim.control(event)
    sim.fields()
    primary = sim.force.numpy().copy()
    sim.pairs()
    total = sim.force.numpy().copy()
    p = sim.p
    sim.launch(k.free_velocity, sim.total, [sim.force, sim.free, 6 * np.pi * p.viscosity_Pa_s * p.radius_m])
    velocity = sim.free
    for it in range(p.mobility_iterations):
        out = sim.vel1 if it % 2 == 0 else sim.vel2
        sim.launch(
            k.local_mobility,
            sim.total,
            [
                sim.x,
                sim.free,
                velocity,
                out,
                sim.n,
                p.radius_m,
                p.lubrication_activation,
                p.lubrication_floor,
            ],
        )
        velocity = out
    sim.launch(k.component_init, sim.total, [sim.component_parent, sim.component_delta, sim.component_count])
    sim.launch(
        k.component_links,
        (sim.batch, sim.n, sim.n),
        [sim.x, sim.component_parent, sim.n, (2 + p.lubrication_activation) * p.radius_m],
    )
    sim.launch(
        k.component_balance,
        sim.total,
        [sim.component_parent, sim.free, velocity, sim.component_delta, sim.component_count],
    )
    parent = sim.component_parent.numpy()
    root = np.arange(sim.total)
    for _ in range(sim.n):
        root = parent[root]
    u = velocity.numpy()
    delta = sim.component_delta.numpy()
    count = sim.component_count.numpy()
    balanced = u + delta[root] / count[root, None]
    assert int(sim.error.numpy()[0]) == 0
    np.testing.assert_array_equal(sim.state(), sim.initial)
    return dict(
        primary_N=primary,
        secondary_N=total - primary,
        total_N=total,
        free_velocity_m_s=sim.free.numpy(),
        local_velocity_m_s=u,
        mobility_velocity_m_s=balanced,
        precontact_relative_step_m=1e-5 * (balanced - balanced.mean(0)),
        pressure_complex_Pa=sim.pressure.numpy().view(np.complex128).ravel(),
        component_roots=root,
        component_delta_m_s=delta,
        component_counts=count,
        harmonics_m3=sim.harmonics.numpy().view(np.complex128).reshape(sim.total, 8),
    )
