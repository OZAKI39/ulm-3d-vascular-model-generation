"""Validate audit deliverables and compare read-only observations; never run project code."""
from pathlib import Path
import collections
import datetime
import hashlib
import json
import os
import re
import subprocess

A = Path(__file__).resolve().parent
S = Path('/home/lzy/projects/ulm_3D_vascular')
T = Path('/home/lzy/projects/hemocell_starter')

def load(name):
    return json.loads((A / name).read_text())

def dump(name, value):
    (A / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')

def sha(data):
    return hashlib.sha256(data).hexdigest()

def snapshot_source():
    rows = []
    for b, ds, fs in os.walk(S, followlinks=False):
        ds.sort()
        fs.sort()
        for n in ['.'] + ds + fs:
            if n == '.' and Path(b) != S:
                continue
            p = Path(b) if n == '.' else Path(b) / n
            st = p.lstat()
            row = [str(p.relative_to(S)), st.st_mode, st.st_size, st.st_mtime_ns,
                   st.st_ctime_ns, st.st_ino, st.st_nlink]
            if p.is_symlink():
                row.append(os.readlink(p))
            rows.append(row)
    return rows

before = load('source_stat_before.json')
after = snapshot_source()
b = {r[0]: r for r in before}
a = {r[0]: r for r in after}
added = sorted(a.keys() - b.keys())
removed = sorted(b.keys() - a.keys())
changed = sorted(k for k in a.keys() & b.keys() if a[k] != b[k])
source_unchanged = not (added or removed or changed)
hash_errors = []
hashes = load('authored_text_source_hashes.json')
for item in hashes:
    p = S / item['path']
    if not p.is_file() or sha(p.read_bytes()) != item['sha256']:
        hash_errors.append(item['path'])

# Target-side operations stay within the user's explicit read-command allowlist.
target_checks = []
target_commands = [
    ('status', ['git', '--no-optional-locks', '-C', str(T), 'status', '--short']),
    ('head', ['git', '-C', str(T), 'rev-parse', 'HEAD']),
    ('structure', ['find', str(T), '-maxdepth', '2', '-printf', '%p\t%y\t%s\t%T@\t%C@\t%i\n']),
]
for name, argv in target_commands:
    result = subprocess.run(argv, capture_output=True, env={**os.environ, 'GIT_OPTIONAL_LOCKS': '0'})
    old = (A / f'target_{name}_before.txt').read_bytes()
    (A / f'target_{name}_after.txt').write_bytes(result.stdout)
    target_checks.append({'name': name, 'argv': argv, 'exit_code': result.returncode,
                          'stderr': result.stderr.decode(errors='replace'),
                          'before_sha256': sha(old), 'after_sha256': sha(result.stdout),
                          'unchanged': result.returncode == 0 and result.stdout == old})
target_unchanged = all(r['unchanged'] for r in target_checks)
verification = {
    'checked_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'classification': 'CONFIRMED_FROM_DATA',
    'source_stat': {'scope': 'Entire source tree, including hidden and nested Git directories; lstat, no symlink traversal',
                    'fields': ['relative_path', 'mode', 'size', 'mtime_ns', 'ctime_ns', 'inode', 'nlink', 'symlink_target_if_any'],
                    'before_entries': len(before), 'after_entries': len(after),
                    'added': added, 'removed': removed, 'changed': changed, 'unchanged': source_unchanged},
    'source_authored_content': {'sha256_count': len(hashes), 'mismatches': hash_errors, 'unchanged': not hash_errors},
    'target_checks': target_checks,
    'limitations': ['Read access times (atime) are excluded.',
                    'Full 28 GiB source payload was not rehashed; full-tree metadata and 180 authored-file SHA256 baselines were compared.',
                    'Target verification covers Git status/HEAD and depth-2 metadata, not a full recursive content hash.'],
    'source_project_modified_by_audit': False,
    'hemocell_project_modified_by_audit': False,
    'source_observation_unchanged': source_unchanged and not hash_errors,
    'target_observation_unchanged': target_unchanged,
    'solver_or_project_script_executed': False,
    'source_or_target_tests_executed': False,
    'build_or_dependency_install_performed': False,
    'migration_or_geometry_copy_performed': False,
    'git_commit_or_push_performed': False,
    'all_audit_outputs_root': str(A),
    'result': 'PASS' if source_unchanged and not hash_errors and target_unchanged else 'REVIEW_REQUIRED',
}
dump('READ_ONLY_VERIFICATION.json', verification)

errors = []
inventory = load('migration_inventory.json')
required = {'source_path', 'component', 'purpose', 'entrypoint', 'solver_dependency', 'input', 'output',
            'units', 'coordinate_system', 'dependencies', 'experimental_relevance', 'migration_class',
            'modify_required', 'target_candidate', 'migration_status', 'evidence', 'unknowns'}
classes = {'KEEP_AS_IS', 'ADAPT', 'REWRITE', 'REFERENCE_ONLY', 'DATA_ONLY', 'DO_NOT_MIGRATE'}
confidence = {'CONFIRMED_FROM_CODE', 'CONFIRMED_FROM_CONFIG', 'CONFIRMED_FROM_DATA', 'INFERRED', 'UNVERIFIED'}
line_cache = {}
evidence_count = 0
for i, item in enumerate(inventory):
    if not required.issubset(item):
        errors.append(f'Inventory {i}: missing keys {sorted(required - item.keys())}')
    p = Path(item['source_path'])
    if not p.is_relative_to(S) or not p.exists():
        errors.append(f'Inventory {i}: invalid source path')
    if item['migration_class'] not in classes or item['solver_dependency'] not in {'independent', 'specific', 'mixed', 'unknown'}:
        errors.append(f'Inventory {i}: invalid category')
    if item['migration_status'] != 'DISCOVERED':
        errors.append(f'Inventory {i}: incorrect scope status')
    if not isinstance(item['entrypoint'], bool) or not isinstance(item['modify_required'], bool):
        errors.append(f'Inventory {i}: incorrect boolean type')
    for e in item['evidence']:
        evidence_count += 1
        ep = Path(e['file'])
        if not ep.exists() or e['classification'] not in confidence:
            errors.append(f'Inventory {i}: missing/invalid evidence {ep}')
        if 'line' in e:
            if str(ep) not in line_cache:
                line_cache[str(ep)] = ep.read_text().splitlines()
            if not 1 <= e['line'] <= max(1, len(line_cache[str(ep)])):
                errors.append(f'Inventory {i}: evidence line out of bounds')
headings = re.findall(r'^# (\d+)\. (.+)$', (A / 'MIGRATION_CONTEXT.md').read_text(), re.M)
if [int(n) for n, _ in headings] != list(range(1, 23)):
    errors.append('Context must contain exactly the requested 22 numbered sections')
for row in load('units_inventory.json'):
    if not {'quantity', 'value', 'unit', 'source_file', 'source_line_or_context', 'meaning', 'confidence'}.issubset(row) or row['confidence'] not in confidence:
        errors.append('Invalid units record')
for row in load('physics_parameters.json'):
    if row['likely_origin'] not in {'EXPERIMENTAL', 'OLD_SOLVER_PARAMETER', 'GEOMETRY_PARAMETER', 'PLACEHOLDER', 'UNKNOWN'} or row['status'] not in confidence:
        errors.append('Invalid physics classification')
risks = load('risk_register.json')
if {r['risk'] for r in risks} != {'UNIT_RISK', 'COORDINATE_RISK', 'BOUNDARY_LABEL_RISK', 'GEOMETRY_TOPOLOGY_RISK', 'SOLVER_COUPLING_RISK', 'EXPERIMENTAL_PROVENANCE_RISK', 'DEPENDENCY_RISK', 'DATA_DUPLICATION_RISK'}:
    errors.append('Required risk coverage incomplete')
for row in risks:
    if row['level'] not in {'LOW', 'MEDIUM', 'HIGH', 'UNVERIFIED'}:
        errors.append('Invalid risk level')
for row in load('geometry_feature_inventory.json'):
    if row['path'] and not Path(row['path']).is_file():
        errors.append('Geometry feature evidence missing')
json_files = sorted(A.glob('*.json'))
for p in json_files:
    json.loads(p.read_text())
if verification['result'] != 'PASS':
    errors.append('Read-only comparison requires review')
result = {'result': 'PASS' if not errors else 'FAIL', 'errors': errors,
          'inventory_entries': len(inventory), 'inventory_evidence_records': evidence_count,
          'context_numbered_sections': len(headings), 'units_records': len(load('units_inventory.json')),
          'physics_records': len(load('physics_parameters.json')), 'required_risks': len(risks),
          'json_files_parsed': len(json_files),
          'migration_class_counts': dict(collections.Counter(x['migration_class'] for x in inventory)),
          'scope': 'Report structure, enums, source path existence, evidence line bounds, read-only before/after comparison; no project tests or solver execution'}
dump('REPORT_VALIDATION.json', result)
report_files = sorted(p for p in A.iterdir() if p.is_file() and p.name != 'SHA256SUMS')
(A / 'SHA256SUMS').write_text(''.join(f'{sha(p.read_bytes())}  {p.name}\n' for p in report_files))
print(json.dumps({'read_only': verification['result'], 'report_validation': result,
                  'source_stat_entries': len(after), 'source_authored_hashes': len(hashes),
                  'target_checks_unchanged': target_unchanged, 'audit_files_hashed': len(report_files)}, ensure_ascii=False, indent=2))
if errors:
    raise SystemExit(1)
