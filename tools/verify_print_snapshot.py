#!/usr/bin/env python3
"""Read-only verification of the published printing snapshot (standard library)."""
from pathlib import Path
import hashlib
import json
import sys

ROOT = Path(__file__).resolve().parents[1]

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()

def main():
    metadata = ROOT / 'reports/github_print_sync'
    errors = []
    count = 0
    for row in (metadata / 'SHA256SUMS').read_text().splitlines():
        expected, relative = row.split('  ', 1)
        path = ROOT / relative
        if not path.is_file() or digest(path) != expected:
            errors.append(relative)
        count += 1
    manifest = json.loads((metadata / 'source_file_manifest.json').read_text())
    for item in manifest['entries']:
        if item['kind'] == 'symlink':
            path = ROOT / item['path']
            if not path.is_symlink() or str(path.readlink()) != item['target'] or not path.exists():
                errors.append(item['path'])
    print(json.dumps({'checked_files': count, 'errors': errors,
                     'status': 'PASS' if not errors else 'FAIL',
                     'scope': 'Snapshot bytes and portable symlinks only; not O3 smoothness or print qualification.'},
                    ensure_ascii=False, indent=2))
    return bool(errors)

if __name__ == '__main__':
    sys.exit(main())
