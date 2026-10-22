"""Read-only reuse of installed wheels into an independent, existing venv.

Adapted from F10U copy_local_torch.py. Copies files, never links or edits source.
Console entry points outside site-packages are intentionally not copied.
"""
import argparse
import hashlib
import importlib.metadata as md
import json
from pathlib import Path
import shutil
from packaging.requirements import Requirement

p = argparse.ArgumentParser()
p.add_argument('--source', type=Path, required=True)
p.add_argument('--destination', type=Path, required=True)
p.add_argument('--roots', nargs='+', required=True)
p.add_argument('--record', type=Path, required=True)
a = p.parse_args()
norm = lambda n: n.lower().replace('_', '-')
available = {norm(d.metadata['Name']): d for d in md.distributions(path=[str(a.source)])}
installed = {norm(d.metadata['Name']) for d in md.distributions(path=[str(a.destination)])}
queue = list(a.roots); seen = set(); rows = []
while queue:
    name = norm(queue.pop())
    if name in seen or name in installed:
        continue
    seen.add(name); d = available[name]
    for text in d.requires or []:
        req = Requirement(text)
        if req.marker and not req.marker.evaluate():
            continue
        queue.append(req.name)
    if name == 'cuda-toolkit':
        queue.extend(k for k in available if k.startswith('nvidia-'))
    count = size = 0
    for f in d.files or []:
        if '..' in f.parts or '__pycache__' in f.parts:
            continue
        src, dst = a.source / f, a.destination / f
        if not src.is_file():
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst); count += 1; size += src.stat().st_size
    rows.append(dict(name=name, version=d.version, files=count, bytes=size,
                     source_RECORD_sha256=hashlib.sha256((d.read_text('RECORD') or '').encode()).hexdigest()))
    print(name, d.version, size, flush=True)
a.record.parent.mkdir(parents=True, exist_ok=True)
a.record.write_text(json.dumps(dict(source=str(a.source), destination=str(a.destination),
    mode='independent files; source read-only; no symlinks', distributions=rows), indent=2)+'\n')
