"""Seal original failed inputs, holdout selection and isolated environment."""
import hashlib
import importlib.metadata as md
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import h5py
import numpy as np
from scripts.stage4l.dense_probe import sha, write

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / 'results/stage5/motion_solver_r2'


def main():
    old = ROOT / 'results/stage5/motion_solver'
    manifest = json.loads((old / 'final/ARTIFACT_MANIFEST.json').read_text())
    entries = {r['path']: r for r in manifest['files']}
    rows = []
    for p in sorted((old / 'diagnostics').glob('*/*.npz')):
        rel = str(p.relative_to(ROOT))
        assert entries[rel]['sha256'] == sha(p), rel
        out = BASE / 'fixtures' / p.parent.name / p.name
        out.parent.mkdir(parents=True, exist_ok=True)
        assert not out.exists()
        shutil.copyfile(p, out)
        out.chmod(0o444)
        rows.append(dict(source=rel, path=str(out.relative_to(ROOT)), sha256=sha(p)))
    write(BASE / 'fixtures/MANIFEST.json', rows)
    selected = json.loads((old / 'audit_v2/BASELINE_SELECTION.json').read_text())
    manifest_rows = json.loads((old / 'audit_v2/BASELINE_MANIFEST_ROWS.json').read_text())
    holdout = []
    for row in selected:
        if not row['key'].endswith('5ms'):
            continue
        path = Path(row['run']) / 'trajectories.h5'
        digest = sha(path)
        # Old manifest content is retained verbatim; hash must occur in it.
        assert digest in json.dumps(manifest_rows), path
        with h5py.File(path) as h:
            n = len(h['time_s'])
            for index in [int(.37*(n-1)), int(.73*(n-1)), n-2]:
                raw = b''.join(np.asarray(h[k][index] if k in ['time_s','positions_m','active'] else h[k][:]).tobytes()
                               for k in ['time_s','positions_m','active','radii_m','particle_id'])
                holdout.append(dict(case=row['key'][:2], source=str(path), source_sha256=digest,
                                    index=index, time_s=float(h['time_s'][index]),
                                    state_sha256=hashlib.sha256(raw).hexdigest(), status='SEALED_NOT_EVALUATED'))
    write(BASE / 'fixtures/HOLDOUT_SELECTION.json', holdout)
    env = BASE / 'environment'
    lock = subprocess.check_output(['uv','pip','freeze','--python',str(ROOT/'.venv/bin/python')], text=True)
    (ROOT / 'requirements-stage5a-r2.lock.txt').write_text(lock)
    packages = []
    for d in md.distributions():
        packages.append(dict(name=d.metadata['Name'], version=d.version,
                             license=d.metadata.get('License-Expression') or d.metadata.get('License'),
                             license_files=d.metadata.get_all('License-File')))
    write(env / 'PACKAGES_LICENSES.json', packages)
    resource = dict(python=platform.python_version(), executable=os.sys.executable,
                    affinity=sorted(os.sched_getaffinity(0)), cpu_max=Path('/sys/fs/cgroup/cpu.max').read_text().strip(),
                    memory_max=Path('/sys/fs/cgroup/memory.max').read_text().strip(),
                    gpu=subprocess.check_output(['nvidia-smi','--query-gpu=name,driver_version,memory.total','--format=csv,noheader'],text=True),
                    disk_free_bytes=shutil.disk_usage(ROOT).free,
                    compiler=subprocess.check_output(['gcc','--version'],text=True).splitlines()[0],
                    builds='PyPI binary wheels; no user compilation flags; float64; no fast-math',
                    lock_sha256=sha(ROOT/'requirements-stage5a-r2.lock.txt'),
                    package_inventory_sha256=sha(env/'PACKAGES_LICENSES.json'))
    write(env / 'IDENTITY.json', resource)
    print(json.dumps(dict(fixtures=len(rows), holdout=len(holdout), environment=resource)))

if __name__ == '__main__':
    main()
