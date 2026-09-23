#!/usr/bin/env python3
import argparse
import json
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fem3d.audit import environment, write_json

p = argparse.ArgumentParser()
p.add_argument("--output", required=True)
args = p.parse_args()
data = environment()
write_json(Path(args.output).with_suffix(".json"), data)
lines = [f"Host: {data['hostname']}", f"UTC: {data['timestamp']}",
         f"Probe Python: {data['python_executable']}", f"Python: {data['python_version']}"]
for entry in data["commands"] + list(data["packages"].values()) + [data["petsc"], data["mpi_library"]]:
    lines += ["\n$ " + " ".join(entry["command"]), f"exit={entry['returncode']}",
              entry.get("stdout", ""), entry.get("stderr", entry.get("error", ""))]
Path(args.output).write_text("\n".join(lines))
print("\n".join(lines))
