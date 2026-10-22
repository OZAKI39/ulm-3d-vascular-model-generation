"""Persist production CI and independent archived-reference checks."""

import subprocess
from scripts.stage3.common import new_run, finish

out = new_run("software_checks", {}, ["Stage1/2/2.5 and Stage3 regression"])
rows = []
for name, args in [
    ("production_tests", [".venv/bin/pytest", "-q", "--junitxml=" + str(out / "pytest.xml")]),
    ("archived_tests", [".venv/bin/pytest", "-q", "-m", "deprecated_reference"]),
    ("ruff", [".venv/bin/ruff", "check", "src", "tests", "scripts", "validation/history_force_sensitivity"]),
]:
    p = subprocess.run(args, capture_output=True, text=True)
    (out / (name + ".log")).write_text(p.stdout + p.stderr)
    rows.append(
        dict(name=name, command=args, returncode=p.returncode, output_tail=(p.stdout + p.stderr)[-500:])
    )
finish(out, dict(status="PASS" if all(r["returncode"] == 0 for r in rows) else "FAIL", checks=rows))
