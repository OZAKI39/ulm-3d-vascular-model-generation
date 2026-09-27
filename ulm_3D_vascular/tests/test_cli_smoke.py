"""Help-only smoke tests for the two formal production entry points."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _help(script: str) -> str:
    completed = subprocess.run(
        [sys.executable, str(PROJECT_ROOT / script), "--help"],
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
    )
    assert completed.returncode == 0, completed.stderr
    return completed.stdout


def test_cfd_preprocess_help_smoke():
    output = _help("cfd_preprocess.py")
    assert "config" in output


def test_cfd_surface_prepare_help_has_no_recovery_mode():
    output = _help("cfd_surface_prepare.py")
    assert "config" in output
    assert "resume" not in output.lower()
    assert "recovery" not in output.lower()
