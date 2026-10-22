"""Read-only audit of the specific mouse snapshot, its origins and Git ownership.

No deletion, copying, checkout, network access or scientific computation.
"""
from pathlib import Path
import hashlib
import json
import os
import subprocess

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'docs/cleanup_20261023'
SNAPSHOT = ROOT.parent / 'temp_storage/mouse_flow_github_sync_20261022'


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def git(path, *args):
    result = subprocess.run(['git', '-C', str(path), *args], capture_output=True, text=True)
    return {'returncode': result.returncode, 'stdout': result.stdout, 'stderr': result.stderr}


def save(name, value):
    (OUT / name).write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n')


def main():
    OUT.mkdir(exist_ok=True)
    sources = json.loads((ROOT / 'docs/audit/published_INCLUDED.json').read_text())
    records = []
    for row in sources:
        path = SNAPSHOT / row['path']
        original = Path(row['source'])
        records.append({**row, 'snapshot_exists': path.is_file(),
                        'snapshot_sha256': sha(path) if path.is_file() else None,
                        'original_exists': original.is_file(),
                        'original_sha256': sha(original) if original.is_file() else None})
    save('snapshot_content_comparison.json', records)
    mismatches = [r for r in records if r['snapshot_sha256'] != r['sha256']
                  or r['original_sha256'] != r['sha256']]
    save('snapshot_comparison_summary.json', {
        'files': len(records), 'logical_bytes': sum(r['size'] for r in records),
        'mismatches': mismatches,
        'status': 'REVIEW' if mismatches else 'PASS'})
    old_sources = json.loads((ROOT / 'docs/audit/source_map.json').read_text())
    source_review = []
    for row in old_sources:
        p = Path(row['source'])
        source_review.append({**row, 'original_current_sha256': sha(p) if p.is_file() else None,
                             'new_current_sha256': sha(ROOT / row['path'])})
    save('current_source_vs_original.json', source_review)
    candidates = [SNAPSHOT, ROOT.parent / 'temp_storage/ulm_particle_3d_particle0',
                  ROOT.parent / 'temp_storage/ulm-3d-vascular-model-generation',
                  ROOT.parent / 'computation_examples/_shared/active_network_development',
                  ROOT.parent / 'computation_examples/_shared/particle_workflow',
                  ROOT.parent / 'computation_examples/_shared/legacy_flow_workflow',
                  ROOT.parent / 'computation_examples/_shared/geometry_network',
                  ROOT.parent / 'computation_examples/_archives/github_sync_20260928_latest']
    states = []
    for p in candidates:
        states.append({'path': str(p), **{label: git(p, *args) for label, args in {
            'root': ('rev-parse', '--show-toplevel'),
            'git_common_dir': ('rev-parse', '--path-format=absolute', '--git-common-dir'),
            'head': ('rev-parse', 'HEAD'), 'branch': ('branch', '--show-current'),
            'status': ('status', '--porcelain=v1', '--untracked-files=all'),
            'ignored': ('ls-files', '--others', '--ignored', '--exclude-standard'),
            'upstream': ('rev-parse', '--abbrev-ref', '@{upstream}'),
            'local_unique_vs_cached_remotes': ('log', '--oneline', 'HEAD', '--not', '--remotes'),
        }.items()}})
    save('git_before.json', states)
    save('registered_worktrees_before.json', git(SNAPSHOT, 'worktree', 'list', '--porcelain'))
    links = []
    for start in [ROOT, SNAPSHOT, ROOT.parent / 'temp_storage']:
        maxdepth = 1 if start.name == 'temp_storage' else 100
        for base, dirs, files in os.walk(start, followlinks=False):
            dirs[:] = [d for d in dirs if d not in ('.venv', '.git', '__pycache__')]
            depth = len(Path(base).relative_to(start).parts)
            for name in dirs + files:
                p = Path(base) / name
                if p.is_symlink():
                    links.append({'path': str(p), 'target': os.readlink(p),
                                  'realpath': str(p.resolve()), 'exists': p.exists()})
            if depth >= maxdepth - 1:
                dirs[:] = []
    save('symlinks_before.json', links)
    print(json.dumps({'snapshot_files': len(records), 'mismatches': len(mismatches),
                      'source_original_changes_since_migration': sum(
                          r['original_current_sha256'] != r['original_sha256'] for r in source_review),
                      'git_candidates': len(states), 'links': len(links)}, indent=2), flush=True)


if __name__ == '__main__':
    main()
