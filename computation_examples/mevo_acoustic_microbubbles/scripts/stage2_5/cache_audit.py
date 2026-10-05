"""Verify every published node key and retained-cycle hash without rerunning APECSS."""

from collections import Counter
import hashlib
import json
from pathlib import Path
import time

from mevo_mb.bubbles.adaptive_table import node_key
from mevo_mb.bubbles.apecss import RadialCase
from common import atomic_json, provenance


def main():
    start = time.perf_counter()
    root = Path("results/stage2_5/response_table/node_cache")
    binary = hashlib.sha256(Path("third_party/build/apecss/ultrasound_apecss").read_bytes()).hexdigest()
    counts = Counter()
    problems = []
    attempts = Counter()
    size = 0
    for p in sorted(root.glob("*/NODE.json")):
        m = json.loads(p.read_text())
        counts[m["status"]] += 1
        expected = node_key(RadialCase(**m["case"]), binary)
        if expected != p.parent.name or expected != m["key"]:
            problems.append(dict(node=str(p), error="key mismatch"))
        attempts[len(list(p.parent.glob("attempt_*")))] += 1
        if m["status"] == "COMPLETE":
            payload = p.parent / m["payload"]
            if (
                not payload.exists()
                or hashlib.sha256(payload.read_bytes()).hexdigest() != m["payload_sha256"]
            ):
                problems.append(dict(node=str(p), error="payload missing or changed"))
            else:
                size += payload.stat().st_size
    result = dict(
        status="PASS" if not problems else "FAIL",
        node_status_counts=dict(counts),
        attempts_per_key=dict(attempts),
        checked_payload_bytes=size,
        issues=problems,
        runtime_s=time.perf_counter() - start,
        provenance=provenance({}, ["APECSS v1.8; cache identity and retained-cycle file integrity"]),
    )
    atomic_json(Path("results/stage2_5/audit/CACHE_AUDIT.json"), result)
    print({k: v for k, v in result.items() if k != "provenance"}, flush=True)
    if problems:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
