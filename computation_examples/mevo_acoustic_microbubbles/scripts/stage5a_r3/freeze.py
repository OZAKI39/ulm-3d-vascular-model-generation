"""Source-bound R3 scoped motion protocol; no automatic reference promotion."""

import argparse
from datetime import datetime, timezone
import hashlib
import subprocess
import sys
import xml.etree.ElementTree as ET
import h5py
import numpy as np
from scripts.stage5a_r3.common import ROOT, BASE, R2, read, write, sha, identity

PROTOCOL = BASE / "final/VALIDATION_PROTOCOL_R3_V1.json"
EXECUTABLE_SCRIPTS = [
    "common.py",
    "run.py",
    "resources.py",
    "freeze.py",
    "holdout.py",
    "campaign.py",
    "reference.py",
    "r2_reference.py",
]


def sources():
    paths = list((ROOT / "src").rglob("*.py"))
    paths += [ROOT / "scripts/stage5a_r3" / name for name in EXECUTABLE_SCRIPTS]
    paths += [
        ROOT / "scripts/stage4l/dense_probe.py",
        ROOT / "scripts/stage5a_r2/run.py",
        ROOT / "scripts/stage5a_r2/resources.py",
    ]
    return {str(p.relative_to(ROOT)): sha(p) for p in paths}


def verify():
    value = read(PROTOCOL)
    assert value["execution_source_sha256"] == sources(), "R3_SOURCE_CHANGED_AFTER_FREEZE"
    assert value["environment_sha256"] == sha(BASE / "environment/IDENTITY.json")
    assert value["holdout_selection_sha256"] == sha(BASE / "fixtures/HOLDOUT_R3_V1.json")
    assert value["test_gate_sha256"] == sha(BASE / "final/TEST_GATE_R3_V1.json")
    return value


def main():
    if PROTOCOL.exists():
        raise FileExistsError("Preserve existing protocol; never overwrite a frozen campaign")
    paths = [
        "tests/test_stage5a_r2.py",
        "tests/test_stage5a_r2_entry.py",
        "tests/test_stage5a_r3_geometry.py",
        "tests/test_stage5a_r3_qp.py",
        "tests/test_stage4l_lubrication.py",
        "tests/test_stage2_translation.py",
    ]
    junit = BASE / "environment/JUNIT_SOURCE_BOUND_R3_V1.xml"
    before = sources()
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", *paths, "--junitxml=" + str(junit)], cwd=ROOT
    )
    assert result.returncode == 0, "SOURCE_BOUND_TESTS_FAILED"
    assert before == sources()
    suites = ET.parse(junit).getroot()
    summary = {
        key: sum(int(s.get(key, 0)) for s in suites.iter("testsuite"))
        for key in ["tests", "failures", "errors", "skipped"]
    }
    assert summary["failures"] == summary["errors"] == 0
    write(
        BASE / "final/TEST_GATE_R3_V1.json",
        dict(
            status="PASS",
            junit=str(junit),
            junit_sha256=sha(junit),
            summary=summary,
            source_sha256=before,
            test_source_sha256={p: sha(ROOT / p) for p in paths},
            identity=identity(),
        ),
    )
    choices = []
    campaign = read(R2 / "final/CAMPAIGN_V4A/RESULT.json")["ledger"]
    for case in ["B2", "B3"]:
        row = next(
            r
            for r in campaign
            if r["case"] == case
            and r["duration"] == 0.005
            and r["kind"] == "nominal"
            and r.get("repeat") == 0
        )
        path = __import__("pathlib").Path(row["run"]) / "trajectories.h5"
        digest = sha(path)
        with h5py.File(path) as h:
            n = len(h["time_s"])
            # Unused internal indices selected by a fixed hash recipe. Earlier
            # all-ID descriptive audits saw trajectories; these numerical states
            # have not been used for R3 solver/layout development.
            selected = []
            for salt in range(3):
                seed = int(
                    hashlib.sha256(("R3-V1-unopened-QP-" + case + str(salt) + digest).encode()).hexdigest()[
                        :12
                    ],
                    16,
                )
                k = 100 + seed % (n - 200)
                while k in selected:
                    k += 1
                selected.append(k)
                raw = b"".join(
                    np.asarray(
                        h[key][k] if key in ["positions_m", "active", "time_s"] else h[key][:]
                    ).tobytes()
                    for key in ["positions_m", "active", "time_s", "particle_id", "radii_m"]
                )
                choices.append(
                    dict(
                        case=case,
                        source=str(path),
                        source_sha256=digest,
                        index=k,
                        time_s=float(h["time_s"][k]),
                        state_sha256=hashlib.sha256(raw).hexdigest(),
                        scope="New numerical solver holdout; trajectory-level descriptive lineage was previously inspected",
                    )
                )
    write(BASE / "fixtures/HOLDOUT_R3_V1.json", choices)
    protocol = dict(
        version="R3_V1",
        frozen_utc=datetime.now(timezone.utc).isoformat(),
        execution_source_sha256=before,
        environment_sha256=sha(BASE / "environment/IDENTITY.json"),
        test_gate_sha256=sha(BASE / "final/TEST_GATE_R3_V1.json"),
        holdout_selection_sha256=sha(BASE / "fixtures/HOLDOUT_R3_V1.json"),
        identity=identity(),
        method=dict(
            name="S1 endpoint resistance with frozen start forcing and finite-feature start contact; separately audited bounded S0 recovery",
            optimized=True,
            components=False,
            component_workers=1,
            pair_backend="cupy",
            nonlinear_strategy="quasi_newton",
            initial_dt_s=1e-8,
            output_clock_s=1e-5,
            nominal=dict(max_dt_s=5e-6, tolerance_factor=0.25, S0_cap_s=50e-9),
            refined=dict(max_dt_s=2.5e-6, tolerance_factor=0.0625, S0_cap_s=25e-9),
        ),
        numerical=dict(
            kkt=1e-9,
            resistance_residual=1e-11,
            contact_tolerance_m=1e-11,
            projection_m=0.0,
            S0_burst_macros=8,
            S0_max_bursts=1024,
        ),
        accuracy=dict(
            all_active_RMS_over_R=0.01,
            all_active_max_over_R=0.05,
            thresholds_source="STAGE5A_R3_PLAN.md preregistered before R3 outcomes",
            independent_global_bound="NOT_ESTABLISHED",
            gate_without_bound="NOT_EVALUABLE",
            Rg_endpoint=0.02,
            Rg_series=0.02,
            pair_quartiles=0.02,
            centroid="max(0.1R, 0.02 reference travel)",
            baseline_disagreement_blocks_promotion=True,
        ),
        budget=dict(
            new_total_s=14400,
            new_1ms_run_s=1200,
            new_5ms_run_s=2400,
            reference_total_s=28800,
            reference_run_s=4500,
            reference_1ms_run_s=1200,
            r2_protected_total_s=5400,
        ),
        formal_policy="One trajectory process at a time; complete duration only for speedup. Numerical continuation may be evaluated while all-active reference accuracy remains NOT_EVALUABLE. No promotion without all gates.",
    )
    write(PROTOCOL, protocol)
    print(dict(protocol=str(PROTOCOL), tests=summary, holdouts=len(choices)), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--verify", action="store_true")
    if p.parse_args().verify:
        print(verify()["version"])
    else:
        main()
