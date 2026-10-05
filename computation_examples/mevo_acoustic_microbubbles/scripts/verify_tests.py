"""Persist complete pytest/ruff evidence, including actual CPU and GPU solvers."""
from pathlib import Path
import subprocess
import sys

from mevo_mb.io.artifacts import write_json


def main():
    outcomes = {}
    for name, command in [("TEST_RESULTS", [sys.executable, "-m", "pytest", "-q", "--junitxml=docs/pytest.xml"]),
                          ("RUFF_RESULTS", [sys.executable, "-m", "ruff", "check", "src", "tests", "scripts"])]:
        with (Path("docs")/(name+".txt")).open("w") as stream:
            result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT, check=False)
        outcomes[name] = {"command": command, "returncode": result.returncode}
    write_json(Path("docs/TEST_STATE.json"), {"status": "PASS" if all(v["returncode"] == 0 for v in outcomes.values()) else "FAIL", "checks": outcomes})
    if any(v["returncode"] for v in outcomes.values()):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
