"""Resource-only V2 protocol: carry budgets, retain equivalent B2 executions."""

import argparse
import ast
from datetime import datetime, timezone
import hashlib
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
import h5py
import numpy as np
from scripts.stage5a_r3.common import ROOT, BASE, R2, read, write, sha, identity
from scripts.stage5a_r3.freeze import sources as original_sources
from scripts.stage5a_r3.campaign import numerical

PROTOCOL = BASE / "final/VALIDATION_PROTOCOL_R3_V2.json"


def sources():
    values = original_sources()
    values.update(
        {
            f"scripts/stage5a_r3/{name}": sha(ROOT / "scripts/stage5a_r3" / name)
            for name in ["freeze_v2.py", "holdout_v2.py", "campaign_v2.py"]
        }
    )
    return values


def verify():
    p = read(PROTOCOL)
    assert p["execution_source_sha256"] == sources(), "V2_SOURCE_CHANGED_AFTER_FREEZE"
    assert p["environment_sha256"] == sha(BASE / "environment/IDENTITY.json")
    assert p["holdout_selection_sha256"] == sha(BASE / "fixtures/HOLDOUT_R3_V2.json")
    assert p["test_gate_sha256"] == sha(BASE / "final/TEST_GATE_R3_V2.json")
    assert p["prior_campaign_sha256"] == sha(BASE / "final/CAMPAIGN_R3_V1/CLOSED_RESULT.json")
    return p


def main():
    assert not PROTOCOL.exists(), "Never replace a frozen protocol"
    old = read(BASE / "final/VALIDATION_PROTOCOL_R3_V1.json")
    closed = read(BASE / "final/CAMPAIGN_R3_V1/CLOSED_RESULT.json")
    assert closed["closed_early"]
    benchmark = read(BASE / "benchmarks/LAYOUT_FULL_TRIAL.json")
    assert benchmark["status"] == "PASS"
    before = sources()
    for name, digest in old["execution_source_sha256"].items():
        if name.startswith("src/"):
            assert before[name] == digest, name
    # Changing the protocol-validation guard did not change numerical execution.
    previous = Path(
        "/workspace/mevo_acoustic_microbubbles_stage5a_r3/.execution/fa2fba5608037ed45143c8ede501cadf415e9180/scripts/stage5a_r3/run.py"
    )

    def body(path):
        node = next(
            n for n in ast.parse(path.read_text()).body if isinstance(n, ast.FunctionDef) and n.name == "run"
        )
        assert isinstance(node.body[0], ast.If)
        node.body = node.body[1:]
        return ast.dump(node, include_attributes=False)

    assert body(previous) == body(ROOT / "scripts/stage5a_r3/run.py"), "NUMERICAL_DRIVER_CHANGED"
    paths = [
        "tests/test_stage5a_r2.py",
        "tests/test_stage5a_r2_entry.py",
        "tests/test_stage5a_r3_geometry.py",
        "tests/test_stage5a_r3_qp.py",
        "tests/test_stage4l_lubrication.py",
        "tests/test_stage2_translation.py",
        "tests/test_stage5a_r3_accounting.py",
        "tests/test_stage5a_r3_reporting.py",
    ]
    junit = BASE / "environment/JUNIT_SOURCE_BOUND_R3_V2.xml"
    result = subprocess.run(
        [sys.executable, "-m", "pytest", "-q", *paths, "--junitxml=" + str(junit)], cwd=ROOT
    )
    assert result.returncode == 0 and before == sources()
    suites = ET.parse(junit).getroot()
    counts = {
        key: sum(int(x.get(key, 0)) for x in suites.iter("testsuite"))
        for key in ["tests", "failures", "errors", "skipped"]
    }
    assert counts["failures"] == counts["errors"] == counts["skipped"] == 0
    write(
        BASE / "final/TEST_GATE_R3_V2.json",
        dict(
            status="PASS",
            summary=counts,
            junit=str(junit),
            junit_sha256=sha(junit),
            source_sha256=before,
            test_source_sha256={p: sha(ROOT / p) for p in paths},
            identity=identity(),
        ),
    )
    previous_indices = {(r["case"], r["index"]) for r in read(BASE / "fixtures/HOLDOUT_R3_V1.json")}
    selected = []
    ledger = read(R2 / "final/CAMPAIGN_V4A/RESULT.json")["ledger"]
    for case in ["B2", "B3"]:
        row = next(
            r
            for r in ledger
            if r["case"] == case
            and r["kind"] == "nominal"
            and r["duration"] == 0.005
            and r.get("repeat") == 0
        )
        path = Path(row["run"]) / "trajectories.h5"
        digest = sha(path)
        with h5py.File(path) as h:
            for salt in range(3):
                seed = int(
                    hashlib.sha256(("R3-V2-unopened-QP-" + case + str(salt) + digest).encode()).hexdigest()[
                        :12
                    ],
                    16,
                )
                k = 100 + seed % (len(h["time_s"]) - 200)
                while (case, k) in previous_indices:
                    k += 1
                previous_indices.add((case, k))
                raw = b"".join(
                    np.asarray(
                        h[key][k] if key in ["positions_m", "active", "time_s"] else h[key][:]
                    ).tobytes()
                    for key in ["positions_m", "active", "time_s", "particle_id", "radii_m"]
                )
                selected.append(
                    dict(
                        case=case,
                        source=str(path),
                        source_sha256=digest,
                        index=k,
                        time_s=float(h["time_s"][k]),
                        state_sha256=hashlib.sha256(raw).hexdigest(),
                        scope="New numerical-state holdout indices after layout freeze; complete source trajectories were previously inspected descriptively",
                    )
                )
    write(BASE / "fixtures/HOLDOUT_R3_V2.json", selected)
    reuse = []
    for r in closed["ledger"]:
        if r["case"] == "B2":
            assert numerical(r["summary"])
            i = read(Path(r["run"]) / "IDENTITY.json")
            assert (
                i["optimized"]
                and not i["components"]
                and i["component_workers"] == 1
                and i["pair_backend"] == "cupy"
            )
            for name, digest in i["source_sha256"].items():
                if name.startswith("src/"):
                    assert before[name] == digest
            reuse.append(
                dict(
                    case="B2",
                    kind=r["kind"],
                    duration_s=r["duration_s"],
                    run=r["run"],
                    identity_sha256=sha(Path(r["run"]) / "IDENTITY.json"),
                    summary_sha256=sha(Path(r["run"]) / "SUMMARY.json"),
                    source_commit=i["source_commit"],
                    reason="Identical scientific sources, driver numerical AST, inputs, parameters and B2 layout. Actual V1 execution retained; no reexecution claim.",
                )
            )
    method = dict(old["method"])
    method.pop("components", None)
    method["layouts"] = dict(
        B2=dict(components=False, component_workers=1), B3=dict(components=True, component_workers=1)
    )
    protocol = dict(
        version="R3_V2",
        frozen_utc=datetime.now(timezone.utc).isoformat(),
        execution_source_sha256=before,
        environment_sha256=sha(BASE / "environment/IDENTITY.json"),
        test_gate_sha256=sha(BASE / "final/TEST_GATE_R3_V2.json"),
        holdout_selection_sha256=sha(BASE / "fixtures/HOLDOUT_R3_V2.json"),
        prior_campaign_sha256=sha(BASE / "final/CAMPAIGN_R3_V1/CLOSED_RESULT.json"),
        layout_benchmark_sha256=sha(BASE / "benchmarks/LAYOUT_FULL_TRIAL.json"),
        method=method,
        numerical=old["numerical"],
        accuracy=old["accuracy"],
        budget=old["budget"],
        initial_spent_s=closed["spent_s"],
        equivalent_reused_B2_runs=reuse,
        driver_numerical_AST_unchanged=True,
        scientific_source_unchanged=True,
        identity=identity(),
        change="B3 global QP to component serial only; B2 global retained because initial compression dominates its measured total.",
        motivation="At real B3 1.5/2.5 ms states complete 5 us trials take 1.62/1.32 s global versus 0.35/0.29 s serial; same-ID paths differ at floating-point scale. Short single-substep benchmarks did not cover complete trials at these costly intermediate states.",
        formal_policy="One trajectory at a time. Existing per-run/group/overall budgets carry forward without reset; no tolerance, model, radius, geometry, S0 limit or accuracy-threshold change.",
    )
    write(PROTOCOL, protocol)
    print(
        dict(version="R3_V2", tests=counts, reused=len(reuse), initial_spent_s=closed["spent_s"]), flush=True
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--verify", action="store_true")
    if p.parse_args().verify:
        print(verify()["version"])
    else:
        main()
