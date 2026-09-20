#!/usr/bin/env python3
"""Audit immutable SV1/SV1.1 evidence and the entire read-only FEM tree."""
import gzip
import json
import os
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'src'))
from sv_validation.sv12 import REPORT, validate_frozen
from sv_validation.provenance import inventory, compare_inventory, git_state, sha256, write_json, now

history = json.loads(gzip.decompress((REPORT/'history_baseline.json.gz').read_bytes()))
changes = []
for key, expected in history['files'].items():
    path = ROOT/key
    if not path.exists() and not path.is_symlink():
        changes.append({'path': key, 'reason': 'missing'})
        continue
    stat = path.lstat()
    if path.is_symlink():
        actual = {'symlink': os.readlink(path), 'mode': stat.st_mode}
    elif path.is_file():
        actual = {'sha256': sha256(path), 'size': stat.st_size,
                  'mtime_ns': stat.st_mtime_ns, 'mode': stat.st_mode}
    else:
        actual = {'directory': True, 'mode': stat.st_mode}
    if actual != expected:
        changes.append({'path': key, 'before': expected, 'after': actual})
for folder in history['scopes']:
    for key in inventory(ROOT/folder)['files']:
        full = folder+'/'+key
        if full not in history['files']:
            changes.append({'path': full, 'reason': 'added historical artifact'})
old = json.loads(gzip.decompress((REPORT/'old_fem_baseline.json.gz').read_bytes()))
current = inventory(old['root'])
comparison = compare_inventory(old, current)
git_after = git_state(old['root'])
git_unchanged = old['git'] == git_after
try:
    frozen = validate_frozen()
    frozen_error = None
except Exception as exc:
    frozen, frozen_error = False, str(exc)
result = {
    'status': 'PASS' if not changes and comparison['status'] == 'PASS' and git_unchanged and frozen else 'FAIL',
    'timestamp': now(), 'history_changes': changes,
    'historical_entries_checked': len(history['files']),
    'old_fem_entries_checked': len(current['files']),
    'old_fem_comparison': comparison, 'old_fem_git_unchanged': git_unchanged,
    'old_fem_git_after': git_after, 'frozen_inputs_unchanged': frozen,
    'frozen_input_error': frozen_error,
    'atime_excluded': True, 'comparison': 'SHA256, size, mtime_ns, mode, symlink target; all preexisting evidence'}
write_json(REPORT/'preservation_audit.json', result)
print('SV1.2 history and old FEM preservation:', result['status'], flush=True)
raise SystemExit(0 if result['status'] == 'PASS' else 1)
