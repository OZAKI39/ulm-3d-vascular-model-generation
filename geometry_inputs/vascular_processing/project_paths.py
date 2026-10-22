"""Explicit compatibility with byte-preserved provenance from the old directory."""
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def preserved_workspace_path(relative):
    """Locate a historical file, including CFD files deliberately left behind."""
    path=ROOT/relative
    if path.exists():return path
    plan=json.loads((ROOT/'migration/20260926/migration_plan.json').read_text())
    original=Path(plan['source'])/relative
    if not original.exists():raise FileNotFoundError(original)
    return original


def frozen_hash(mapping,path):
    """Find a saved digest by the physical file, including recorded old aliases."""
    if str(path) in mapping:return mapping[str(path)]
    resolved=Path(path).resolve()
    matches={value for key,value in mapping.items() if Path(key).resolve()==resolved}
    if len(matches)!=1:raise KeyError('Missing or conflicting frozen hash for '+str(path))
    return matches.pop()


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        while block:=stream.read(4*1024*1024):h.update(block)
    return h.hexdigest()


def verify_migrated_snapshot(snapshot,root=ROOT):
    """Verify history without rewriting original expected hashes.

    Only edits explicitly recorded by the migration may differ. Both the
    original archived bytes and the current edited bytes must match the ledger.
    A migration edit is reported separately; it is never called unchanged.
    """
    root=Path(root).resolve();ledger=root/'migration/20260926/path_compatibility_edits.json'
    edits=json.loads(ledger.read_text()) if ledger.exists() else {}
    changed=[];authorized=[];unchanged=0
    for name,expected in snapshot.items():
        path=Path(name);actual=digest(path) if path.is_file() else None
        if actual==expected:unchanged+=1;continue
        try:relative=str(path.resolve().relative_to(root))
        except ValueError:relative=None
        entry=edits.get(relative,{})
        original=root/entry.get('original_copy','__missing__')
        if (entry.get('before_sha256')==expected and entry.get('after_sha256')==actual
                and actual is not None and original.is_file() and digest(original)==expected):
            authorized.append(name)
        else:changed.append(name)
    return dict(passed=not changed,unchanged_file_count=unchanged,
        migration_only_changes=authorized,unexpected_changes=changed,
        all_content_unchanged=not changed and not authorized)
