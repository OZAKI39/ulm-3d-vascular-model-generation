"""Rebuild the academic CFD deck without invoking any scientific solver.

The PowerPoint is authored with @oai/artifact-tool; Python only orchestrates
the reproducible visual and deck build steps.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
OUTPUT = PROJECT / "references/cfd.pptx"
RUNTIME_PACKAGE = (
    Path.home()
    / ".cache/codex-runtimes/codex-primary-runtime/dependencies/node/node_modules"
    / "@oai/artifact-tool"
)


def ensure_artifact_tool(workspace: Path) -> None:
    package_json = RUNTIME_PACKAGE / "package.json"
    if not package_json.is_file():
        raise FileNotFoundError(f"Bundled @oai/artifact-tool not found: {package_json}")
    package = json.loads(package_json.read_text(encoding="utf-8"))
    if package.get("name") != "@oai/artifact-tool":
        raise RuntimeError(f"Unexpected package at {RUNTIME_PACKAGE}")

    (workspace / "node_modules/@oai").mkdir(parents=True, exist_ok=True)
    target = workspace / "node_modules/@oai/artifact-tool"
    try:
        os.symlink(RUNTIME_PACKAGE, target, target_is_directory=True)
    except OSError:
        link_path = str(target).replace("'", "''")
        source_path = str(RUNTIME_PACKAGE).replace("'", "''")
        subprocess.run(
            [
                "powershell",
                "-NoProfile",
                "-Command",
                f"New-Item -ItemType Junction -Path '{link_path}' -Target '{source_path}' | Out-Null",
            ],
            check=True,
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-visuals", action="store_true")
    parser.add_argument(
        "--output",
        type=Path,
        default=OUTPUT,
        help="PPTX destination; defaults to references/cfd.pptx",
    )
    args = parser.parse_args()
    output = args.output.resolve()

    if not args.skip_visuals:
        subprocess.run(
            [sys.executable, str(HERE / "generate_cfd_ppt_visuals.py")],
            cwd=PROJECT,
            check=True,
        )

    workspace = Path(tempfile.mkdtemp(prefix="cfd_ppt_artifact_tool_"))
    (workspace / "package.json").write_text(
        json.dumps({"private": True, "type": "module"}, indent=2) + "\n",
        encoding="utf-8",
    )
    ensure_artifact_tool(workspace)
    build_module = workspace / "build_cfd_ppt.mjs"
    shutil.copy2(HERE / "build_cfd_ppt.mjs", build_module)
    subprocess.run(
        [
            "node",
            str(build_module),
            str(PROJECT),
            str(output),
            str(HERE / "preview"),
        ],
        cwd=workspace,
        check=True,
    )
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
