"""Causal F8 sample plus observed F9 displacement ledger; no training approval."""

import h5py
import numpy as np
from .f8_samples import read_navigation_sample


def read_audit_sample(path, frame, history_frames=3):
    sample = read_navigation_sample(path, frame, history_frames)
    now = float(sample["trajectory"]["time_s"][-1])
    first = float(sample["trajectory"]["time_s"][0])
    with h5py.File(path) as h:
        if "motion_ledger" not in h:
            sample["motion_ledger"] = None
        else:
            group = h["motion_ledger"]
            meta = group["step_meta"][:]
            ends = meta[:, 0] + meta[:, 1]
            stop = int(np.searchsorted(ends, now + 1e-13, side="right"))
            begin = int(np.searchsorted(meta[:, 0], first - 1e-13))
            events = group["wall_contact_events"][:]
            events = events[(events[:, 1] >= first - 1e-13) & (events[:, 1] + events[:, 2] <= now + 1e-13)]
            sample["motion_ledger"] = dict(
                step_meta=meta[begin:stop],
                displacement_m=group["displacement_m"][begin:stop],
                observed_cumulative_m=group["displacement_m"][:stop].sum(0),
                wall_contact_events=events,
                restriction="Observed accepted steps only. Full-run per-ID totals and paired timestep labels are audit outputs, never causal features.",
            )
    sample["split_restriction"] = (
        "All F5/F6/F7/F8/F9 original-Y 300ms initial-state controls, timesteps and audit reproductions share one family; no frame-wise or replay-wise train/test split."
    )
    sample["training_ready"] = False
    return sample
