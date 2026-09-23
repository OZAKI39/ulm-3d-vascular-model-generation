#!/usr/bin/env python3
import argparse
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fem3d.audit import reference_snapshot, write_json

p = argparse.ArgumentParser()
p.add_argument("output")
p.add_argument("--compare")
args = p.parse_args()
roots = ["/home/lzy/projects/ulm_microbubble_traj_gen_2D", "/home/lzy/projects/ulm_3D_vascular"]
snapshot = reference_snapshot(roots)
write_json(args.output, snapshot)
print(f"Recorded {len(snapshot['files'])} entries")
if args.compare:
    import json
    old = json.loads(Path(args.compare).read_text())["files"]
    new = snapshot["files"]
    changes = [k for k in old.keys() | new.keys() if old.get(k) != new.get(k)]
    write_json(Path(args.output).with_name("reference_integrity.json"),
               {"status": "FAIL" if changes else "PASS", "baseline": args.compare,
                "final": args.output, "entries": len(new), "changes": changes})
    print(f"Changes: {len(changes)}")
    sys.exit(bool(changes))
