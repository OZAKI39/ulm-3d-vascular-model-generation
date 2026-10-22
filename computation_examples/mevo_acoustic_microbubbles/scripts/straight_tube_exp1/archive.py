"""Return all F10L scientific outputs and logs with independent remote hashes."""

import subprocess
from .common import BASE, read, write, sha
from .remote import SERVER


def main():
    code = """from pathlib import Path
import json,hashlib,subprocess
root=Path(ROOT)
status=subprocess.run(['supervisorctl','status'],capture_output=True,text=True)
assert status.returncode in (0,3),status.stderr
states=status.stdout
assert 'mevo_f10l_' in states,'Missing F10L supervisor records'
assert not any(line.startswith('mevo_f10l_') and 'RUNNING' in line for line in states.splitlines()),'Wait for F10L tasks to exit'
files={str(p.relative_to(root)):dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 for p in sorted(root.rglob('*')) if p.is_file() and p.name!='REMOTE_MANIFEST.json'}
(root/'REMOTE_MANIFEST.json').write_text(json.dumps(dict(files=files,supervisor_status=states,
 gpu_state=subprocess.check_output(['nvidia-smi','--query-gpu=name,memory.used,utilization.gpu','--format=csv,noheader'],text=True)),indent=2))
""".replace("ROOT", repr(SERVER + "/results/straight_tube_exp1/f10l"))
    subprocess.run(["ssh", "mevo-new", "python3 -"], input=code, text=True, check=True)
    subprocess.run(
        ["rsync", "-a", "mevo-new:" + SERVER + "/results/straight_tube_exp1/f10l/", str(BASE) + "/"],
        check=True,
    )
    manifest = read(BASE / "REMOTE_MANIFEST.json")
    mismatches = []
    for rel, record in manifest["files"].items():
        path = BASE / rel
        if not path.exists() or path.stat().st_size != record["bytes"] or sha(path) != record["sha256"]:
            mismatches.append(rel)
    write(
        BASE / "analysis/TRANSFER_VERIFICATION.json",
        dict(
            status="PASS" if not mismatches else "FAIL",
            files=len(manifest["files"]),
            bytes=sum(v["bytes"] for v in manifest["files"].values()),
            mismatches=mismatches,
            remote_manifest_sha256=sha(BASE / "REMOTE_MANIFEST.json"),
        ),
    )
    if mismatches:
        raise RuntimeError(mismatches)
    print("Verified remote files:", len(manifest["files"]))


if __name__ == "__main__":
    main()
