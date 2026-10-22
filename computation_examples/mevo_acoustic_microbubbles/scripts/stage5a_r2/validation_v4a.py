"""Repair only V4's formal entry routing; kernel, settings and holdouts unchanged."""

import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from scripts.stage5a_r2.seal import ROOT, BASE
from scripts.stage5a_r2.validation_v4 import verify as verify_v4
from scripts.stage5a_r2.freeze import critical_hashes
from scripts.stage4l.dense_probe import sha, write


def hashes():
    return {
        f"scripts/stage5a_r2/{p}": sha(ROOT / "scripts/stage5a_r2" / p)
        for p in ["validation_v4a.py", "campaign_v4a.py", "run_v4a.py"]
    }


def verify():
    verify_v4()
    p = json.loads((BASE / "final/VALIDATION_PROTOCOL_V4A.json").read_text())
    g = json.loads((BASE / "final/TEST_GATE_V4A.json").read_text())
    assert p["source_sha256"] == critical_hashes() == g["source_sha256"]
    assert p["orchestration_sha256"] == hashes() == g["orchestration_sha256"]
    assert g["status"] == "PASS" and g["protocol_sha256"] == sha(BASE / "final/VALIDATION_PROTOCOL_V4A.json")
    assert g["junit_sha256"] == sha(ROOT / g["junit_path"])
    assert all(sha(ROOT / k) == v for k, v in g["test_source_sha256"].items())
    assert p["environment_sha256"] == sha(BASE / "environment/IDENTITY.json") == g["environment_sha256"]
    return p


def main():
    previous = verify_v4()
    target = BASE / "final/VALIDATION_PROTOCOL_V4A.json"
    assert not target.exists()
    p = json.loads(json.dumps(previous))
    p.update(
        version="STAGE5A_R2_S1_V6_STRICT_V4A_ENTRY_FIX",
        time_utc=datetime.now(timezone.utc).isoformat(),
        source_commit=(ROOT / "SOURCE_COMMIT_R2.txt").read_text().strip(),
        orchestration_sha256=hashes(),
        entry_fix="V4 used the V3 formal whitelist and both finer jobs were rejected before initialization. V4A validates the complete V4 contract in a separate versioned entry and calls the identical frozen loop. No change to any numerical parameter, kernel or holdout.",
        prior_protocol_sha256=sha(BASE / "final/VALIDATION_PROTOCOL_V4.json"),
    )
    write(target, p)
    oldgate = json.loads((BASE / "final/TEST_GATE_V4.json").read_text())
    tests = list(oldgate["test_source_sha256"]) + ["tests/test_stage5a_r2_entry.py"]
    xml = BASE / "environment/FORMAL_TESTS_V4A.xml"
    log = BASE / "environment/FORMAL_TESTS_V4A.log"
    with log.open("x") as f:
        result = subprocess.run(
            [sys.executable, "-m", "pytest", *tests, "-q", f"--junitxml={xml}"],
            cwd=ROOT,
            stdout=f,
            stderr=subprocess.STDOUT,
            timeout=240,
        )
    cases = list(ET.parse(xml).getroot().iter("testcase"))
    fail = sum(len(c.findall("failure")) + len(c.findall("error")) for c in cases)
    skip = sum(len(c.findall("skipped")) for c in cases)
    gate = dict(
        status="PASS" if result.returncode == 0 and fail == 0 else "FAIL",
        tests=len(cases),
        passed=len(cases) - fail - skip,
        failures=fail,
        skipped_or_xfail=skip,
        source_sha256=critical_hashes(),
        orchestration_sha256=hashes(),
        environment_sha256=sha(BASE / "environment/IDENTITY.json"),
        protocol_sha256=sha(target),
        junit_path=str(xml.relative_to(ROOT)),
        junit_sha256=sha(xml),
        test_source_sha256={n: sha(ROOT / n) for n in tests},
    )
    write(BASE / "final/TEST_GATE_V4A.json", gate)
    assert gate["status"] == "PASS"
    holdout = json.loads((BASE / "final/HOLDOUT_V4.json").read_text())
    assert holdout["status"] == "PASS"
    holdout.update(
        reused_without_reexecution=True,
        source_protocol="V4",
        source_sha256=sha(BASE / "final/HOLDOUT_V4.json"),
        reason="Identical kernel/parameters. Only entry routing fixed and separately tested; no holdout-driven retuning.",
    )
    write(BASE / "final/HOLDOUT_V4A.json", holdout)
    print(
        "V4A entry checks and regression:",
        gate["passed"],
        "PASS;",
        skip,
        "old xfail;",
        len(holdout["rows"]),
        "V4 probes reused",
    )


if __name__ == "__main__":
    main()
