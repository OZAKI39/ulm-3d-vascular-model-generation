"""Partition independent records by exact common control schedule; no async claim."""

from collections import OrderedDict
import hashlib
import json
import numpy as np


def control_families(records, max_batch=64):
    if not 1 <= max_batch <= 64:
        raise ValueError("Only bounded batches <=64 are admitted")
    seen = set()
    families = OrderedDict()
    for record in records:
        sample = record["sample_id"]
        if sample in seen:
            raise ValueError("Duplicate sample_id: " + sample)
        seen.add(sample)
        if not record.get("split_group"):
            raise ValueError("Explicit common train/test split family required")
        x = np.asarray(record["positions_m"], float)
        if x.shape != (32, 3) or not np.isfinite(x).all():
            raise ValueError("Full32 finite fixed-ID initial state required")
        payload = json.dumps(record["events"], sort_keys=True, separators=(",", ":"), allow_nan=False)
        key = hashlib.sha256(payload.encode()).hexdigest()
        families.setdefault(key, []).append(record)
    result = []
    for key, items in families.items():
        for begin in range(0, len(items), max_batch):
            members = items[begin : begin + max_batch]
            result.append(
                dict(
                    control_family_sha256=key,
                    events=members[0]["events"],
                    sample_ids=[r["sample_id"] for r in members],
                    positions_m=np.stack([r["positions_m"] for r in members]),
                    split_groups=[r["split_group"] for r in members],
                    training_ready=False,
                )
            )
    return result
