"""Wait for authorized M0 diagnostics, then analyze and close a failed gate."""

import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import time

from scripts.stage4m.common import BASE, write, sha


def clarify_startup_attempts():
    for path in (BASE / "dense_probe").glob("*/*/IDENTITY.json"):
        identity = json.loads(path.read_text())
        run = path.parent
        if not (run / "PROJECTION_RUNTIME.json").exists():
            write(
                run / "STARTUP_STATUS.json",
                dict(
                    status="INFRASTRUCTURE_STARTUP_FAILED",
                    accepted_steps=0,
                    reason="Shared compiled contact cache read before compiler completed; recorded supervisor stderr: shared library file too short.",
                    scientific_trajectory_started=False,
                    recovery="Independent per-run cache; unchanged compiled contact source and compiler flags.",
                ),
            )
        if identity["policy"] == "A1":
            write(
                run / "POLICY_SCOPE.json",
                dict(
                    policy="A1",
                    kind="hard subfloor ablation",
                    formula="w=0 if h<=0.001R, otherwise w=1",
                    allowed_duration_s=0.005,
                    production_candidate=False,
                    original_generic_banner_erratum="The generic smooth-removal sentence in the archived runner applies to A2. This A1 run used a hard subfloor switch, as specified by its config and archived force source.",
                ),
            )


def artifact_index():
    paths = []
    excluded = {"ARTIFACT_INDEX.json", "SYNC_RECEIPT.json"}
    for p in sorted(BASE.rglob("*")):
        if not p.is_file() or p.name in excluded or "ordered_contact_cache" in p.parts or "jobs" in p.parts:
            continue
        paths.append(dict(path=str(p.relative_to(BASE)), bytes=p.stat().st_size, sha256=sha(p)))
    write(
        BASE / "final_decision/ARTIFACT_INDEX.json",
        dict(
            created_utc=datetime.now(timezone.utc).isoformat(),
            artifact_count=len(paths),
            total_bytes=sum(p["bytes"] for p in paths),
            files=paths,
            excluded="Rebuildable compiler caches, live supervisor logs, this self-referential index and sync receipts. Scientific inputs are immutable external references identified by each run.",
        ),
    )


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--legacy-root", type=Path, required=True)
    p.add_argument("--stage4l-root", type=Path, required=True)
    a = p.parse_args()
    start = time.monotonic()
    while True:
        done = {
            (policy, case): list((BASE / "dense_probe" / policy).glob(f"*_{case}_{policy}_5ms/SUMMARY.json"))
            for policy in ("A1", "A2")
            for case in ("B2", "B3")
        }
        if all(len(v) == 1 for v in done.values()):
            break
        if any(len(v) > 1 for v in done.values()):
            raise RuntimeError("Multiple completed attempts require explicit result identity selection")
        if time.monotonic() - start > 3 * 3600:
            raise RuntimeError("M0 completion files not available within three hours")
        time.sleep(30)
    clarify_startup_attempts()
    from scripts.stage4m.audit_inputs import audit
    from scripts.stage4m.analyze import analyze_all
    from scripts.stage4m.final_report import main as report

    audit(a.legacy_root, a.stage4l_root)
    _, gate = analyze_all(a.legacy_root, a.stage4l_root)
    if gate["continue_formation"]:
        raise RuntimeError("M0 passed; do not finalize before authorized original-initial formation stage")
    from scripts.stage4m.validate_artifacts import validate

    validate()
    report()
    artifact_index()
    print("STAGE4M_FINAL_REPORT_READY", flush=True)


if __name__ == "__main__":
    main()
