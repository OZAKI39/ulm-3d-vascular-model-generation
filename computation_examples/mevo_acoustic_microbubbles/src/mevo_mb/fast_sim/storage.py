"""Actual-frame HDF5 dataset interface. No synthetic trajectory interpolation."""

import json
from time import perf_counter
from pathlib import Path
import h5py
import numpy as np
from .metrics import trajectory_metrics, DEFINITIONS


def save_run(path, result, physics, identity, config):
    path = Path(path)
    start = perf_counter()
    positions = result["positions_m"]
    metrics, events = trajectory_metrics(positions, physics.radius_m, result["time_s"])
    metrics_s = perf_counter() - start
    start = perf_counter()
    with h5py.File(path, "x") as h:
        h.attrs.update(
            schema="fast_pzt_f1_v1",
            complete=result["summary"]["status"] == "COMPLETE",
            physics_scope="EXPLORATORY",
            model_json=json.dumps(physics.metadata()),
            identity_json=json.dumps(identity),
            config_json=json.dumps(config),
            cluster_definitions_json=json.dumps(DEFINITIONS),
        )
        h["time_s"] = result["time_s"]
        h["control_index"] = result["control_index"]
        if "edge_event_ids" in result:
            h["edge_events/batch_i_j_linked"] = result["edge_event_ids"]
            h["edge_events/time_s"] = result["edge_event_time_s"]
            h["edge_events"].attrs["scope"] = (
                "Every accepted substep adjacency change; sorted by actual physical time; no interpolation"
            )
        h["particle_id"] = np.arange(positions.shape[2])
        h["radii_m"] = np.full(positions.shape[2], physics.radius_m)
        h.create_dataset("positions_m", data=positions, compression="gzip", compression_opts=1)
        for key, values in metrics.items():
            h.create_dataset("metrics/" + key, data=values, compression="gzip", compression_opts=1)
        h.create_dataset("sampled_lineage_json", data=json.dumps(events), dtype=h5py.string_dtype())
        h.attrs["controls_json"] = json.dumps(result["summary"]["controls"])
    result["summary"].update(
        metrics_wall_s=metrics_s,
        hdf5_wall_s=perf_counter() - start,
        dataset_wall_s=result["summary"]["wall_s"] + metrics_s + perf_counter() - start,
        cluster_definitions=DEFINITIONS,
        sampled_lineage_event_count=len(events),
        initial_cluster_count=metrics["cluster_count"][0].tolist(),
        final_cluster_count=metrics["cluster_count"][-1].tolist(),
        initial_largest_size=metrics["largest_size"][0].tolist(),
        final_largest_size=metrics["largest_size"][-1].tolist(),
        initial_global_rg_m=metrics["global_rg_m"][0].tolist(),
        final_global_rg_m=metrics["global_rg_m"][-1].tolist(),
    )
    return metrics, events
