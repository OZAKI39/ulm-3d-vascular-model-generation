"""Fetch small irreplaceable outputs and verify every remote SHA-256."""

import json
import subprocess
from .common import ROOT, BASE, write, sha, js


def main():
    code = """from pathlib import Path
import hashlib,json,subprocess,sys,importlib.metadata
root=Path('/workspace/mevo_acoustic_microbubbles/f10b/results/fast_pzt/f10b')
files={str(p.relative_to(root)):dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(root.rglob('*')) if p.is_file()}
inputs=Path('/workspace/mevo_acoustic_microbubbles/project/results/fast_pzt/f5')
manifest=json.loads((inputs/'INPUT_MANIFEST.json').read_text())
for rel,meta in manifest.items(): assert hashlib.sha256((inputs/rel).read_bytes()).hexdigest()==meta['sha256'],rel
env=dict(python=sys.version,packages={k:importlib.metadata.version(k) for k in ['warp-lang','numpy','scipy','h5py','pyvista']},gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,driver_version,memory.used,utilization.gpu','--format=csv'],text=True),cpu=subprocess.check_output(['lscpu'],text=True),AGENTS_sha256=hashlib.sha256(Path('/etc/vast-agents-guide.md').read_bytes()).hexdigest())
jobs=subprocess.run(['supervisorctl','status'],capture_output=True,text=True).stdout
configs={p.name:p.read_text() for p in Path('/etc/supervisor/conf.d').glob('mevo_fast_pzt_f10b_*.conf')}
print(json.dumps(dict(files=files,environment=env,verified_input_count=len(manifest),jobs=jobs,supervisor_configs=configs)))
"""
    raw = subprocess.check_output(
        ["ssh", "mevo-new", "/workspace/mevo_acoustic_microbubbles/.venv/bin/python -"], input=code, text=True
    )
    manifest = json.loads(raw)
    subprocess.run(
        [
            "rsync",
            "-az",
            "mevo-new:/workspace/mevo_acoustic_microbubbles/f10b/results/fast_pzt/f10b/",
            str(BASE) + "/",
        ],
        check=True,
    )
    for rel, meta in manifest["files"].items():
        assert sha(BASE / rel) == meta["sha256"], rel
    write(BASE / "REMOTE_MANIFEST.json", manifest)
    original = js(ROOT / "results/fast_pzt/f5/INPUT_MANIFEST.json")
    for rel, meta in original.items():
        assert sha(ROOT / "results/fast_pzt/f5" / rel) == meta["sha256"], rel
    write(
        BASE / "ARCHIVE_VERIFICATION.json",
        dict(
            status="PASS",
            remote_files=len(manifest["files"]),
            bytes=sum(m["bytes"] for m in manifest["files"].values()),
            verified_original_inputs=len(original),
            scope="All remote F10B outputs including failures/logs; historical inputs reused and rehashed, no cache or redundant large dependency copy",
        ),
    )
    print("Archive PASS", len(manifest["files"]), "files")


if __name__ == "__main__":
    main()
