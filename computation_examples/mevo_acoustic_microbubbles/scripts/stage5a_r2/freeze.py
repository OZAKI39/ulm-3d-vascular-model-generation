"""Hash-bound formal protocol and regression gate. Never reuse an old XML gate."""

from datetime import datetime, timezone
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
from scripts.stage4l.dense_probe import write, sha
from scripts.stage5a_r2.seal import ROOT, BASE


def critical_hashes():
    paths = list((ROOT / "src").rglob("*.py"))
    paths += [
        ROOT / "scripts/stage5a_r2" / p
        for p in ["run.py", "resources.py", "python.sh", "compare.py", "freeze.py", "frozen_reference.py"]
    ]
    return {str(p.relative_to(ROOT)): sha(p) for p in sorted(paths)}


def verify():
    protocol = json.loads((BASE / "final/VALIDATION_PROTOCOL_V1.json").read_text())
    gate = json.loads((BASE / "final/TEST_GATE_V1.json").read_text())
    assert protocol["source_sha256"] == critical_hashes(), "formal source differs from frozen protocol"
    assert gate["status"] == "PASS" and gate["source_sha256"] == critical_hashes()
    assert (
        protocol["environment_sha256"]
        == sha(BASE / "environment/IDENTITY.json")
        == gate["environment_sha256"]
    )
    assert gate["junit_sha256"] == sha(ROOT / gate["junit_path"])
    return protocol


def main():
    protocol_path = BASE / "final/VALIDATION_PROTOCOL_V1.json"
    assert not protocol_path.exists(), "frozen validation protocol cannot be overwritten"
    protocol = dict(
        version="STAGE5A_R2_S1_V4_VALIDATION_V1",
        time_utc=datetime.now(timezone.utc).isoformat(),
        source_commit=(ROOT / "SOURCE_COMMIT_R2.txt").read_text().strip(),
        source_sha256=critical_hashes(),
        environment_sha256=sha(BASE / "environment/IDENTITY.json"),
        plan_sha256=sha(ROOT / "docs/STAGE5A_R2_MOTION_PLAN.md"),
        solver=dict(
            name="ProxQP",
            version="0.7.2",
            eps_abs=1e-12,
            eps_rel=0.0,
            max_iter=1000,
            max_iter_in=1500,
            dense_dimension_cap=300,
            dense_rows_cap=1200,
            global_layout="sparse",
            fallback=[
                "OSQP 1.1.3 eps_abs=1e-12 max_iter=300000 time_limit=5s",
                "ProxQP sparse eps_abs=1e-13",
            ],
        ),
        integrator=dict(
            name="S1 semi-implicit Euler, adopted two half steps",
            resistance="A(x_end)",
            force_and_contact_normal="frozen at each substep start",
            nonlinear_strategy="quasi_newton",
            anderson_history=6,
            max_iterations=32,
            final_original_equation_kkt=1e-9,
            position_map_trigger=1e-7,
            max_dt_s=1e-5,
            initial_dt_s=1e-8,
            atol_R=0.002,
            position_rtol=0.01,
            gap_rtol=0.02,
            refinement=dict(max_dt_s=5e-6, all_time_error_tolerances_factor=0.25),
            output_clock_s=1e-5,
            coordinate_interpolation=False,
        ),
        numerical=dict(kkt=1e-9, linear_residual=1e-11, contact_tolerance_m=1e-11, projection=False),
        scientific=dict(
            Rg_endpoint_relative=0.02,
            Rg_series_relative_L2=0.02,
            pair_quartiles_relative=0.02,
            centroid_budget="max(0.1R, 0.02 reference centroid path length)",
        ),
        resources=dict(
            pair_backend="cupy",
            response_backend="numpy",
            BLAS_threads=1,
            scientific_workers=1,
            GPU_budget_GiB=4,
            RAM_budget_GiB=16,
            per_run_s=1800,
            formal_new_aggregate_s=10800,
            reference_aggregate_s=10800,
            repeat_target=3,
        ),
        state_policy="Solved + original independent certificate. Inaccurate candidate retained; choose a separately certified alternative native solve, never accept by status alone.",
        holdout_selection_sha256=sha(BASE / "fixtures/HOLDOUT_SELECTION.json"),
    )
    write(protocol_path, protocol)
    tests = [
        "tests/test_stage5a_motion.py",
        "tests/test_stage5a_osqp.py",
        "tests/test_stage5a_evidence.py",
        "tests/test_stage5a_r2.py",
        "tests/test_stage4l_lubrication.py",
        "tests/test_stage4l_metrics.py",
        "tests/test_stage4m_ablation.py",
    ]
    xml = BASE / "environment/FORMAL_TESTS_V1.xml"
    log = BASE / "environment/FORMAL_TESTS_V1.log"
    with log.open("x") as f:
        result = subprocess.run(
            [sys.executable, "-m", "pytest", *tests, "-q", f"--junitxml={xml}"],
            stdout=f,
            stderr=subprocess.STDOUT,
            timeout=240,
        )
    suite = ET.parse(xml).getroot()
    cases = list(suite.iter("testcase"))
    failures = sum(len(x.findall("failure")) + len(x.findall("error")) for x in cases)
    skipped = sum(len(x.findall("skipped")) for x in cases)
    gate = dict(
        status="PASS" if result.returncode == 0 and failures == 0 else "FAIL",
        tests=len(cases),
        passed=len(cases) - failures - skipped,
        failures=failures,
        skipped_or_xfail=skipped,
        source_sha256=critical_hashes(),
        environment_sha256=sha(BASE / "environment/IDENTITY.json"),
        protocol_sha256=sha(protocol_path),
        junit_path=str(xml.relative_to(ROOT)),
        junit_sha256=sha(xml),
        test_source_sha256={p: sha(ROOT / p) for p in tests},
    )
    write(BASE / "final/TEST_GATE_V1.json", gate)
    print(json.dumps({k: v for k, v in gate.items() if "sha256" not in k}), flush=True)
    assert gate["status"] == "PASS"


if __name__ == "__main__":
    main()
