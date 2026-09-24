#!/usr/bin/env python3
"""Run a command without a shell, retaining command, evidence, stdout and stderr."""
import argparse
import json
import os
import platform
import resource
import shlex
import subprocess
import sys
import time
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fem3d.audit import command, environment, sha256, timestamp, write_json

p = argparse.ArgumentParser()
p.add_argument("--label", required=True)
p.add_argument("--stage", choices=("0", "1", "1.5", "1.6", "1.7", "1.8", "2", "3"), default="0")
p.add_argument("--ranks", type=int, default=1)
p.add_argument("--input", action="append", default=[])
p.add_argument("--remote-dir")
p.add_argument("argv", nargs=argparse.REMAINDER)
a = p.parse_args()
argv = a.argv[1:] if a.argv[:1] == ["--"] else a.argv
if not argv:
    p.error("missing command")
root = Path(__file__).resolve().parents[1]
stamp = time.strftime("%Y%m%dT%H%M%S", time.gmtime()) + f"_{time.time_ns() % 1000000000:09d}"
stage = {"1.5": "stage01_5", "1.6": "stage01_6", "1.7": "stage01_7", "1.8": "stage01_8"}[a.stage] if a.stage in ("1.5", "1.6", "1.7", "1.8") else f"stage{int(a.stage):02d}"
directory = root / "logs" / stage / f"{stamp}_{a.label}"
directory.mkdir(parents=True)
inputs = [Path(v).resolve() for v in a.input]
for folder in ("src", "scripts", "tests", "configs", "remote", "inputs"):
    inputs += [x for x in (root/folder).rglob("*") if x.is_file() and
               x.suffix in (".py", ".sh", ".json", ".yml", ".yaml", ".toml", ".txt", ".stl", ".vtp", ".csv", ".npz")
               and not any(s.startswith(".") or s == "__pycache__" for s in x.relative_to(root).parts)]
inputs += [root/"pyproject.toml"]
contract = root / "reports/stage00/source_contract.json"
if contract.exists():
    inputs += [contract]
baseline = root/"reports/stage00/reference_baseline.json"
if baseline.exists():
    inputs += [baseline]
config = {str(x.relative_to(root)): x.read_text() for x in (root/"configs").glob("*") if x.is_file()}
info = {"command": argv, "command_display": shlex.join(argv), "hostname": platform.node(),
        "timestamp": timestamp(), "cwd": str(Path.cwd()), "mpi_ranks": a.ranks,
        "remote_work_dir": a.remote_dir, "configuration_snapshot": config,
        "input_sha256": {str(x): sha256(x) for x in sorted(set(inputs))},
        "git": command(["git", "-C", str(root), "rev-parse", "HEAD"]),
        "environment": environment(), "stdout_log": str(directory/"stdout.log"),
        "stderr_log": str(directory/"stderr.log")}
manifest_path = root/"remote/source_manifest.json"
if manifest_path.exists():
    manifest = json.loads(manifest_path.read_text())
    info["wsl_source_manifest_sha256"] = sha256(manifest_path)
    info["wsl_source_git"] = manifest.get("wsl_git")
env = os.environ.copy()
env.update(PYTHONDONTWRITEBYTECODE="1", FEM3D_RUN_ID=directory.name,
           FEM3D_STDOUT_LOG=str(directory/"stdout.log"), FEM3D_STDERR_LOG=str(directory/"stderr.log"),
           FEM3D_RUN_METADATA=str(directory/"metadata.json"), FEM3D_STAGE=str(a.stage))
start = time.perf_counter()
info["returncode"] = None
write_json(directory/"metadata.json", info)
try:
    with (directory/"stdout.log").open("w") as out, (directory/"stderr.log").open("w") as err:
        result = subprocess.run(argv, env=env, stdout=out, stderr=err)
    info["returncode"] = result.returncode
except Exception as exc:
    info["error"] = repr(exc)
    info["returncode"] = 127
finally:
    info["wall_time_s"] = time.perf_counter() - start
    info["finished_utc"] = timestamp()
    info["peak_child_rss_kib"] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    info["memory_scope"] = "OS wait4 child high-water mark; not an aggregate across MPI ranks; includes environment probe children"
    write_json(directory/"metadata.json", info)
print(json.dumps({"log_dir": str(directory), "returncode": info["returncode"], "wall_time_s": info["wall_time_s"]}))
for name in ("stdout.log", "stderr.log"):
    path = directory/name
    if path.exists():
        print(path.read_text()[-16000:])
sys.exit(info["returncode"])
