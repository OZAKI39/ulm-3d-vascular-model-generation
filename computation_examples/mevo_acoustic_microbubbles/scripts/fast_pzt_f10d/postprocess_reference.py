"""Pure-CPU reconstruction of the unchanged F10C analysis definitions."""

from time import perf_counter
import numpy as np
from mevo_mb.fast_sim.metrics import trajectory_metrics
from mevo_mb.fast_sim.morphodynamics import morphology_series, topology


def process(x, radius, times, edge_ids, edge_times):
    started = perf_counter()
    metrics, lineage = trajectory_metrics(x, radius, times)
    metrics_s = perf_counter() - started
    morph, top = {}, []
    morphology_s = topology_s = 0.0
    for i in range(x.shape[1]):
        t = perf_counter()
        m = morphology_series(x[:, i], np.arange(x.shape[2]), metrics["labels"][:, i])
        for key, value in m.items():
            morph.setdefault(key, []).append(value)
        morphology_s += perf_counter() - t
        t = perf_counter()
        take = (edge_ids[:, 0] == i) & (edge_times <= times[-1] + 1e-14)
        ei = edge_ids[take].copy()
        ei[:, 0] = 0
        tp, _, _ = topology(x[0, i], times, metrics["labels"][:, i], edge_times[take], ei, times[-1])
        top.append(tp)
        topology_s += perf_counter() - t
    morph = {key: np.stack(value, axis=1) for key, value in morph.items()}
    return (
        metrics,
        lineage,
        morph,
        top,
        dict(
            metrics_and_sampled_lineage_s=metrics_s,
            fixed_dynamic_shape_kabsch_s=morphology_s,
            accepted_step_topology_s=topology_s,
            total_s=perf_counter() - started,
            scope="Uninstrumented unchanged reference; sub-function attribution measured separately by cProfile",
        ),
    )
