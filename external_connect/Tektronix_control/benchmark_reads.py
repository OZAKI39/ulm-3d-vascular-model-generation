"""Compare sequential and batched status reads without changing outputs/settings.

Error-queue queries consume reported errors. No reset, clear, apply, trigger or
output command is sent. Close other instrument-control programs before running.
This measures the queries before ON, not the physical output switching delay.
"""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
from statistics import median
from time import perf_counter

from afg_control import BASE, Rig, load_config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=BASE / "config.json")
    parser.add_argument("--report", type=Path, default=BASE / "latency_report.json")
    parser.add_argument("--samples", type=int, choices=range(1, 11), default=3)
    args = parser.parse_args()
    cfg = load_config(args.config)
    report = {"timestamp_utc": datetime.now(timezone.utc).isoformat(),
              "scope": "Read-only; error queries consume queue entries. No output/settings changes. "
                       "Pre-ON query timing only; no measurement of physical switching.",
              "channels": {}}
    with Rig(cfg) as rig:
        report["devices"] = dict(rig.identities)
        for name in cfg["channels"]:
            session, _ = rig._target(name)
            cfg["query_batch_mode"] = "sequential"
            rig._errors(session)
            expected = rig.read_channel(name)
            row = {"initial_state": expected}
            for mode in ("sequential", "legacy_lines"):
                cfg["query_batch_mode"] = mode
                times = []
                reads = []
                for _ in range(args.samples):
                    read_started = perf_counter()
                    actual = rig.read_channel(name)
                    reads.append((perf_counter() - read_started) * 1000)
                    if actual != expected:
                        raise RuntimeError(f"{name}: readback changed during benchmark")
                    started = perf_counter()
                    actual = rig.read_channel(name, check_errors=True)
                    if actual != expected:
                        raise RuntimeError(f"{name}: checked readback changed during benchmark")
                    times.append((perf_counter() - started) * 1000)
                row[mode] = {"channel_read_median_ms": median(reads),
                             "pre_on_queries_median_ms": median(times),
                             "channel_read_samples_ms": reads, "pre_on_queries_samples_ms": times}
            report["channels"][name] = row
            print(f"{name}: pre-ON queries "
                  f"{row['sequential']['pre_on_queries_median_ms']:.1f} -> "
                  f"{row['legacy_lines']['pre_on_queries_median_ms']:.1f} ms", flush=True)
        # Confirm the response stream is aligned after all batches.
        for alias, session in rig.sessions.items():
            identity = rig._query(session, "*IDN?")
            if rig.identities[alias]["serial"] not in identity:
                raise RuntimeError(f"{alias}: unexpected identity after benchmark")
            rig._errors(session)
    args.report.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(f"Report: {args.report}")


if __name__ == "__main__":
    main()
