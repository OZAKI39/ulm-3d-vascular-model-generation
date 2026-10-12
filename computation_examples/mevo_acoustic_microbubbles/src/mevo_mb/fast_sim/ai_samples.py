"""Read inherited HDF5 plus F6 logical groups, preserving paired split identity."""

from enum import Enum
import json
import h5py
import numpy as np


class EvidenceLevel(str, Enum):
    EXPLORATORY = "EXPLORATORY"
    MODEL_CHECKED = "MODEL_CHECKED"
    PHANTOM_VALIDATED = "PHANTOM_VALIDATED"


def read_group_sample(path, frame, history_frames=3):
    if history_frames < 1:
        raise ValueError("Positive history window required")
    with h5py.File(path) as h:
        if not (0 <= frame < len(h["time_s"])):
            raise IndexError(frame)
        first = max(0, frame - history_frames + 1)
        environment = json.loads(h.attrs["config_json"])["environment"]
        flow = {
            key: h["flow/" + key][first : frame + 1]
            for key in ("velocity_m_s", "gradient_s_inv", "valid", "flow_enabled")
        }
        enabled = flow["flow_enabled"]
        flow["applied_velocity_m_s"] = np.where(enabled[:, None, None, None], flow["velocity_m_s"], 0.0)
        flow["kind"] = environment["name"]
        return dict(
            source=json.loads(h["source"].attrs["identity_json"]),
            medium=json.loads(h["medium"].attrs["json"]),
            split_group=h.attrs["split_group"],
            partition=h.attrs["partition"],
            evidence_level=EvidenceLevel(h.attrs["evidence_level"]).value,
            bubbles=dict(particle_id=h["particle_id"][:], radii_m=h["radii_m"][:]),
            environment=environment,
            trajectory=dict(
                time_s=h["time_s"][first : frame + 1],
                positions_m=h["positions_m"][first : frame + 1],
                particle_id=h["particle_id"][:],
            ),
            morphology={
                key: h["metrics/" + key][first : frame + 1]
                for key in (
                    "labels",
                    "cluster_sizes",
                    "global_centroid_m",
                    "global_rg_m",
                    "global_shape_tensor_m2",
                )
            },
            flow=flow,
            wall_gap_m=h["wall_gap_m"][first : frame + 1],
            active=h["active"][first : frame + 1],
            control_index=h["control_index"][first : frame + 1],
            control=json.loads(h["control"].attrs["events_json"]),
            quality=json.loads(h["quality"].attrs["summary_json"]),
        )


def validate_partitions(paths):
    assignments = {}
    for path in paths:
        with h5py.File(path) as h:
            group = h.attrs["split_group"]
            part = h.attrs["partition"]
            if group in assignments and assignments[group] != part:
                raise ValueError("Paired scenario leaks across partitions")
            assignments[group] = part
    return assignments
