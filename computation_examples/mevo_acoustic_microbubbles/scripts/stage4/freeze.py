"""Stage 4A: preserve Stage 3 identities and run existing regression checks."""
import json
from pathlib import Path
import subprocess
from scripts.stage4.common import ROOT, new_run, sha, write, finish


def main():
    expected = "dfe79a640fa54c164ce523dcdd8e129217148b5a"
    actual = subprocess.check_output(["git", "rev-parse", "stage3-minimal-validated^{}"], text=True).strip()
    assert actual == expected, (actual, expected)
    acceptance = json.loads(Path("PHYSICS_ACCEPTANCE.json").read_text())
    assert acceptance["status"] == "MINIMAL_SIMULATOR_VALIDATED_WITH_SCOPE"
    for item in acceptance["items"].values():
        for evidence in item["evidence"]:
            assert sha(evidence["path"]) == evidence["sha256"], evidence["path"]
    out = new_run("audit", {"baseline_tag": "stage3-minimal-validated", "baseline_commit": actual})
    protected = {str(p.relative_to(ROOT)): sha(p) for p in (ROOT / "results/stage3").rglob("*")
                 if p.is_file() and "__pycache__" not in p.parts}
    for item in acceptance["items"].values():
        for evidence in item["evidence"]:
            protected[evidence["path"]] = sha(evidence["path"])
    write(out / "STAGE3_BASELINE_FILES.json", protected)
    rows = []
    for label, cmd in [
        ("production_regression", [".venv/bin/pytest", "-q", "--junitxml=" + str(out / "pytest.xml")]),
        ("archive_reference", [".venv/bin/pytest", "-q", "-m", "deprecated_reference"]),
        ("ruff", [".venv/bin/ruff", "check", "src", "tests", "scripts", "validation/history_force_sensitivity"]),
    ]:
        process = subprocess.run(cmd, text=True, capture_output=True)
        (out / (label + ".log")).write_text(process.stdout + process.stderr)
        rows.append(dict(label=label, returncode=process.returncode, output=(process.stdout + process.stderr)[-700:]))
    changed = [name for name, digest in protected.items() if sha(ROOT / name) != digest]
    finish(out, dict(status="PASS" if all(r["returncode"] == 0 for r in rows) and not changed else "FAIL",
                     baseline_commit=actual, baseline_tag="stage3-minimal-validated", checks=rows,
                     protected_files=len(protected), changed=changed,
                     accepted_run_identity_check="All PHYSICS_ACCEPTANCE evidence SHA-256 matched"))


if __name__ == "__main__":
    main()
