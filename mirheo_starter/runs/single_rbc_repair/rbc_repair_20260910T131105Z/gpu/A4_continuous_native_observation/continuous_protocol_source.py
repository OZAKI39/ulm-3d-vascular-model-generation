"""Observation protocol for the isolated single-RBC diagnostic."""


def evolve_stage(coordinator, count, dt, sample):
    """Run once; native plugins observe at `sample` without reinitialization.

    The caller registers the observers first. The pinned native API rebuilds
    tasks, splits particle vectors, and clears object forces on every run().
    """
    if count <= 0 or sample <= 0 or dt <= 0:
        raise ValueError('POSITIVE_CONTINUOUS_PROTOCOL_REQUIRED')
    coordinator.run(count, dt=dt)
    return count
