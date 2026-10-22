"""Download immutable scientific outputs and verify every remote byte hash."""

import json
import subprocess
from .common import BASE, sha, write


def main():
    remote = "/workspace/mevo_acoustic_microbubbles/f10d/results/fast_pzt/f10d"
    code = r"""
import hashlib,json,subprocess,sys,platform
from pathlib import Path
import numpy,h5py,warp
root=Path('/workspace/mevo_acoustic_microbubbles/f10d/results/fast_pzt/f10d')
status=subprocess.run(['supervisorctl','status'],text=True,capture_output=True).stdout
assert not any('mevo_fast_pzt_f10d_' in line and 'RUNNING' in line for line in status.splitlines())
folder=root/'supervisor';folder.mkdir(exist_ok=True)
for p in Path('/etc/supervisor/conf.d').glob('mevo_fast_pzt_f10d_*.conf'):
 (folder/p.name).write_bytes(p.read_bytes())
(root/'ENVIRONMENT.json').write_text(json.dumps(dict(python=sys.version,platform=platform.platform(),numpy=numpy.__version__,h5py=h5py.__version__,warp=warp.__version__,supervisor=status,gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,driver_version,memory.total,utilization.gpu,memory.used','--format=csv'],text=True),external_compute_processes=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,process_name,used_gpu_memory','--format=csv'],text=True),guide_sha256=hashlib.sha256(Path('/etc/vast-agents-guide.md').read_bytes()).hexdigest()),indent=2)+'\n')
files={str(p.relative_to(root)):dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(root.rglob('*')) if p.is_file() and p.name!='REMOTE_MANIFEST.json'}
print(json.dumps(dict(files=files)))
"""
    response = subprocess.check_output(
        ["ssh", "mevo-new", "/workspace/mevo_acoustic_microbubbles/.venv/bin/python -"], input=code, text=True
    )
    manifest = json.loads(response)
    subprocess.run(["rsync", "-a", "mevo-new:" + remote + "/", str(BASE) + "/"], check=True)
    for rel, meta in manifest["files"].items():
        assert (BASE / rel).stat().st_size == meta["bytes"] and sha(BASE / rel) == meta["sha256"], rel
    write(BASE / "REMOTE_MANIFEST.json", manifest)
    write(
        BASE / "ARCHIVE_VERIFICATION.json",
        dict(
            status="PASS",
            remote_files_verified=len(manifest["files"]),
            bytes=sum(m["bytes"] for m in manifest["files"].values()),
            historical_input_files_reused=True,
        ),
    )
    print("Verified", len(manifest["files"]), "remote files")


if __name__ == "__main__":
    main()
