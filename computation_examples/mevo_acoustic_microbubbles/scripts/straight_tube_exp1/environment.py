"""Audited optional import guard for k-wave-python 0.6.2 Python-only execution.

Upstream unconditionally downloads two C++ binaries on import. They are not
used by the Python/CuPy solver. Preserve the original and patch only this
environment-controlled import guard; never modify numerical solver files.
"""

import importlib.metadata
import hashlib
from pathlib import Path
from .common import BASE, write


def install_guard():
    distribution = importlib.metadata.distribution("k-wave-python")
    assert distribution.version == "0.6.2"
    path = Path(distribution.locate_file("kwave/__init__.py"))
    original = path.with_name("__init__.py.f10l_original")
    old = "if not binaries_present():\n    install_binaries()"
    new = 'if os.environ.get("KWAVE_F10L_PYTHON_ONLY") != "1" and not binaries_present():\n    install_binaries()'
    if not original.exists():
        raw = path.read_bytes()
        assert raw.decode().count(old) == 1
        original.write_bytes(raw)
        path.write_text(raw.decode().replace(old, new))
    assert path.read_text() == original.read_text().replace(old, new)
    out = BASE / "environment"
    out.mkdir(parents=True, exist_ok=True)
    (out / "kwave_init_original.py").write_bytes(original.read_bytes())
    (out / "kwave_init_python_guard.py").write_bytes(path.read_bytes())
    write(
        out / "KWAVE_IMPORT_GUARD.json",
        dict(
            package_version=distribution.version,
            original_sha256=hashlib.sha256(original.read_bytes()).hexdigest(),
            installed_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
            exact_old=old,
            exact_new=new,
            solver_files_changed=False,
            scope="Only when KWAVE_F10L_PYTHON_ONLY=1; backend must be python, never cpp",
        ),
    )


if __name__ == "__main__":
    install_guard()
