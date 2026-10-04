"""Delete exact predeclared manifest entries after full preflight validation."""
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
from audit_inventory import AUDIT, PROJECT, ROOT, dump_json, extract_references, hash_file


def main():
    assert not (AUDIT/'DELETION_LOG.csv').exists(), 'Do not replay a destructive manifest.'
    baseline = json.loads((AUDIT/'BASELINE_TESTS.json').read_text())
    assert baseline['counts'] == {'passed': 412, 'failed': 3, 'skipped': 0}
    assert json.loads((AUDIT/'SCIENTIFIC_SMOKE_BEFORE.json').read_text())['loadable']
    assert json.loads((AUDIT/'CLEANUP_CONFIG.json').read_text())['AGGRESSIVE_TRACKED_CLEANUP'] is False
    rows = json.loads((AUDIT/'CLASSIFIED_INVENTORY.json').read_text())
    mapping = {r['path']: r for r in rows}
    with (AUDIT/'DELETE_MANIFEST.csv').open() as f:
        manifest = list(csv.DictReader(f))
    tracked = set(subprocess.check_output(['git','ls-files','-z'],cwd=ROOT).decode().split('\0'))
    refs, unresolved, _ = extract_references(rows)
    assert not [r for r in unresolved if r['kind']=='markdown_link'], 'Pre-existing missing linked artifacts must be resolved first.'
    file_refs = {r['target'] for r in refs if r['within_repository'] and not r['is_directory']}
    dir_refs = {r['target'] for r in refs if r['within_repository'] and r['is_directory'] and r['kind'] in ('markdown_link','json_path') and len(Path(r['target']).parts)>2}
    for item in manifest:
        path = item['path']; p = ROOT / path
        assert mapping[path]['candidate_category'] == 'DELETE_SAFE'
        assert path not in tracked and item['git_tracked'] == 'False'
        assert not item['referenced_by_review'] and path not in file_refs
        assert not any(str(parent) in dir_refs for parent in Path(path).parents)
        assert p.is_relative_to(PROJECT) and not p.is_relative_to(AUDIT)
        assert '.git' not in p.parts and 'references' not in p.parts and p.suffix != '.pdf'
        assert all(not parent.is_symlink() for parent in p.parents if parent != ROOT and parent.is_relative_to(ROOT))
        assert not (path.startswith(('clot_clearance/pd_clot/','clot_clearance/tests/','clot_clearance/scripts/','clot_clearance/configs/')) and p.suffix not in ('.pyc','.pyo','.nbc','.nbi'))
        assert p.lstat().st_size == int(item['size_bytes']), path
        if item['category'] == 'PYTEST_SOURCE_DUPLICATE':
            canonical = item['regenerable_from']
            assert canonical in mapping and mapping[canonical]['candidate_category'] != 'DELETE_SAFE'
            assert mapping[canonical]['hash'] == item['hash'] == hash_file(mapping[canonical])
    with ThreadPoolExecutor(max_workers=4) as pool:
        current = list(pool.map(hash_file, [mapping[item['path']] for item in manifest]))
    assert all(sha == item['hash'] for sha,item in zip(current,manifest)), 'Candidate changed after inventory; stop.'
    summary = json.loads((AUDIT/'PRE_DELETE_SUMMARY.json').read_text())
    assert len(manifest) == summary['DELETE_SAFE']['count']
    dump_json('DELETION_PREFLIGHT.json', {'passed': True, 'manifest_sha256': hashlib.sha256((AUDIT/'DELETE_MANIFEST.csv').read_bytes()).hexdigest(),
                                        'files': len(manifest), 'bytes': sum(int(r['size_bytes']) for r in manifest),
                                        'tracked_files': 0, 'referenced_files': 0, 'production_source_files': 0,
                                        'all_current_candidate_hashes_match_inventory': True})
    # Flush one record per unlink. No wildcard/rmtree operations are used.
    with (AUDIT/'DELETION_LOG.csv').open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['path','status','bytes_removed','category','hash'])
        writer.writeheader(); f.flush()
        for item in manifest:
            p = ROOT / item['path']
            p.unlink()
            writer.writerow(dict(path=item['path'],status='DELETED',bytes_removed=item['size_bytes'],category=item['category'],hash=item['hash']))
            f.flush()
    # Remove only newly empty ancestors of deleted files; preserve all milestone roots.
    parents = {p for item in manifest for p in (ROOT/item['path']).parents if p.is_relative_to(PROJECT) and p != PROJECT}
    removed = []
    for p in sorted(parents,key=lambda p:len(p.parts),reverse=True):
        if p.is_symlink() or not p.is_dir():
            continue
        rel = p.relative_to(PROJECT)
        if len(rel.parts)<=2 and rel.parts[0] in ('verification','results'):
            continue
        if p in [ROOT/r['target'] for r in refs if r['within_repository'] and r['is_directory']]:
            continue
        try:
            p.rmdir(); removed.append(str(p.relative_to(ROOT)))
        except OSError:
            pass
    dump_json('EMPTY_DIRECTORIES_REMOVED.json', removed)
    print(json.dumps({'files_removed':len(manifest),'bytes_removed':sum(int(r['size_bytes']) for r in manifest),'empty_directories_removed':len(removed)}),flush=True)


if __name__ == '__main__':
    main()
