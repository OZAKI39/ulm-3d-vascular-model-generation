#!/usr/bin/env python3
"""Read-only archive integrity checks; no solver, build or network activity."""
from pathlib import Path
import csv, hashlib, json
B = Path(__file__).resolve().parents[1]
prefix = 'review/hemocell_handoff_20260915/'
files = {str(p.relative_to(B)): p for p in B.rglob('*') if p.is_file()}
assert not any(p.is_symlink() for p in B.rglob('*'))
def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()
checks = {}
for line in (B / 'SHA256SUMS').read_text().splitlines():
    digest, rel = line.split('  ', 1)
    assert rel not in checks and rel in files
    assert sha(files[rel]) == digest, rel
    checks[rel] = digest
assert set(checks) == set(files) - {'SHA256SUMS'}
rows = list(csv.DictReader((B / 'SOURCE_MANIFEST.tsv').open(), delimiter='\t'))
assert {r['repo_relative_path'][len(prefix):] for r in rows} == set(files)
for row in rows:
    assert row['repo_relative_path'].startswith(prefix)
    rel = row['repo_relative_path'][len(prefix):]
    assert int(row['size_bytes']) == files[rel].stat().st_size
    if rel not in ('SHA256SUMS', 'SOURCE_MANIFEST.tsv'):
        assert row['sha256'] == sha(files[rel]), rel
state = json.loads((B / 'CURRENT_STATE.json').read_text())
assert state['rbc_stage1'] == 'FAIL_GEOMETRY_GATE'
assert state['rbc_timesteps'] == 0 and state['rbc_wall_interaction'] == 'ABSENT'
assert state['next_stage'] == 'RBC_GEOMETRY_AND_WALL_COMPATIBILITY_STAGE'
contract = json.loads((B / 'rbc_stage1/RBC_STAGE1_CONTRACT.json').read_text())
assert all(x['actual_steps'] == 0 and x['positions'] is None for x in contract['cases'].values())
assert contract['model']['actual_discrete_volume_um3'] == 45.046330078320075
medium = json.loads((B / 'new_medium/NEW_MEDIUM_NUMERICS_CONTRACT.json').read_text())
assert medium['script_sha256'] == sha(B / 'source/new_medium/scripts/prepare_numerics.py')
for row in csv.DictReader((B / 'provenance/ORIGINAL_FILE_HASHES.tsv').open(), delimiter='\t'):
    rel = row['repo_relative_path'][len(prefix):]
    assert row['source_sha256'] == row['handoff_sha256'] == sha(files[rel])
print(json.dumps({'status':'PASS','file_count':len(files),'checksums_checked':len(checks),
                  'total_MiB':sum(p.stat().st_size for p in files.values()) / 2**20,
                  'RBC_TIMESTEPS':0,'RBC_GEOMETRIC_FIT':'FAIL','WALL_INTERACTION_IMPLEMENTATION':'ABSENT'}, indent=2))
