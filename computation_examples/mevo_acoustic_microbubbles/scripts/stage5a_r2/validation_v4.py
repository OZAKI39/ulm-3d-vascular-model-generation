"""Freeze stricter parameters with the IDENTICAL V6 numerical kernel.

V3's coarse B2 Q1 discrepancy remains a failure. This is a bounded parameter
refinement, not a seventh integrator algorithm or a relaxed science gate.
"""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
import subprocess
import sys
import xml.etree.ElementTree as ET
import h5py
import numpy as np
from scripts.stage5a_r2.seal import ROOT, BASE
from scripts.stage5a_r2.freeze import verify as verify_v3, critical_hashes
from scripts.stage4l.dense_probe import sha, write


def orchestration_hashes():
    return {
        f"scripts/stage5a_r2/{n}": sha(ROOT / "scripts/stage5a_r2" / n)
        for n in ["validation_v4.py", "holdout_v4.py", "campaign_v4.py"]
    }


def verify():
    verify_v3()
    p = json.loads((BASE / "final/VALIDATION_PROTOCOL_V4.json").read_text())
    g = json.loads((BASE / "final/TEST_GATE_V4.json").read_text())
    assert p["source_sha256"] == critical_hashes() == g["source_sha256"]
    assert p["orchestration_sha256"] == orchestration_hashes() == g["orchestration_sha256"]
    assert p["environment_sha256"] == sha(BASE / "environment/IDENTITY.json") == g["environment_sha256"]
    assert g["status"] == "PASS" and g["protocol_sha256"] == sha(BASE / "final/VALIDATION_PROTOCOL_V4.json")
    assert g["junit_sha256"] == sha(ROOT / g["junit_path"])
    assert all(sha(ROOT / k) == v for k, v in g["test_source_sha256"].items())
    return p


def main():
    old = verify_v3()
    assert not (BASE / "final/VALIDATION_PROTOCOL_V4.json").exists()
    previous = json.loads((BASE / "fixtures/HOLDOUT_SELECTION_V3.json").read_text())
    selections = []
    for case in ["B2", "B3"]:
        source = next(q for q in previous if q["case"] == case)
        p = Path(source["source"])
        assert sha(p) == source["source_sha256"]
        with h5py.File(p) as h:
            n = len(h["time_s"])
            for fraction in [0.07, 0.31, 0.66]:
                k = int(fraction * (n - 1))
                raw = b"".join(
                    np.asarray(
                        h[key][k] if key in ["time_s", "positions_m", "active"] else h[key][:]
                    ).tobytes()
                    for key in ["time_s", "positions_m", "active", "radii_m", "particle_id"]
                )
                selections.append(
                    dict(
                        case=case,
                        source=str(p),
                        source_sha256=source["source_sha256"],
                        index=k,
                        time_s=float(h["time_s"][k]),
                        state_sha256=hashlib.sha256(raw).hexdigest(),
                        status="SEALED_NOT_EVALUATED",
                    )
                )
    selection = BASE / "fixtures/HOLDOUT_SELECTION_V4.json"
    assert not selection.exists()
    write(selection, selections)
    protocol = json.loads(json.dumps(old))
    protocol.update(
        version="STAGE5A_R2_S1_V6_STRICTER_CONFIGURATION_VALIDATION_V4",
        time_utc=datetime.now(timezone.utc).isoformat(),
        source_commit=(ROOT / "SOURCE_COMMIT_R2.txt").read_text().strip(),
        orchestration_sha256=orchestration_hashes(),
        holdout_selection_sha256=sha(selection),
        numerical_kernel_identical_to_V3=True,
        change_reason="V3 ordinary B2 Q1 error 2.202% fails 2%; V3 refined 1.754% passes. Preserve failure; select pre-existing tighter configuration and add another refinement.",
        prior_protocol_sha256=sha(BASE / "final/VALIDATION_PROTOCOL_V3.json"),
        budget_policy="Carry all V1 and V3 experimental costs forward; no new category allocation or per-run budget increase. At most one stricter-configuration campaign; no further retuning after V4.",
    )
    protocol["integrator"].update(
        max_dt_s=5e-6,
        atol_R=0.002 * 0.25,
        position_rtol=0.01 * 0.25,
        gap_rtol=0.02 * 0.25,
        driver_tolerance_factor=0.25,
        refinement=dict(
            max_dt_s=2.5e-6, all_time_error_tolerances_factor=0.25, driver_tolerance_factor=0.0625
        ),
    )
    protocol["integrator"]["recovery"].update(
        max_macro_dt_s=5e-8,
        refined_macro_dt_s=5e-8,
        refinement_note="Unchanged V6 binary recovery cap is 50 ns for both stricter configurations. Error budgets tighten by four; report actual S0 grids if used, without claiming this cap was halved.",
    )
    target = BASE / "final/VALIDATION_PROTOCOL_V4.json"
    write(target, protocol)
    tests = list(json.loads((BASE / "final/TEST_GATE_V3.json").read_text())["test_source_sha256"])
    xml = BASE / "environment/FORMAL_TESTS_V4.xml"
    log = BASE / "environment/FORMAL_TESTS_V4.log"
    with log.open("x") as stream:
        result = subprocess.run(
            [sys.executable, "-m", "pytest", *tests, "-q", f"--junitxml={xml}"],
            cwd=ROOT,
            stdout=stream,
            stderr=subprocess.STDOUT,
            timeout=240,
        )
    cases = list(ET.parse(xml).getroot().iter("testcase"))
    failures = sum(len(c.findall("failure")) + len(c.findall("error")) for c in cases)
    skipped = sum(len(c.findall("skipped")) for c in cases)
    gate = dict(
        status="PASS" if result.returncode == 0 and not failures else "FAIL",
        tests=len(cases),
        passed=len(cases) - failures - skipped,
        failures=failures,
        skipped_or_xfail=skipped,
        source_sha256=critical_hashes(),
        orchestration_sha256=orchestration_hashes(),
        environment_sha256=sha(BASE / "environment/IDENTITY.json"),
        protocol_sha256=sha(target),
        junit_path=str(xml.relative_to(ROOT)),
        junit_sha256=sha(xml),
        test_source_sha256={p: sha(ROOT / p) for p in tests},
    )
    write(BASE / "final/TEST_GATE_V4.json", gate)
    assert gate["status"] == "PASS"
    print("V4 identical-kernel freeze and fresh tests:", gate["passed"], "passed,", skipped, "old xfail")


if __name__ == "__main__":
    main()
