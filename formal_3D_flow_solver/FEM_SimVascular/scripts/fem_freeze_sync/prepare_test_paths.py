"""Restore a legacy Stage Q source alias inside this checkout only; no CFD.

Original tests are preserved byte-for-byte. Their external/ source path is
materialized as an ignored local source copy. Only upstream .gitattributes is
restored from archival bytes: publishing does not enable upstream LFS rules.
"""
from pathlib import Path
import shutil
import hashlib
import json

ROOT = Path(__file__).resolve().parents[2]
source = ROOT / 'vendor/svMultiPhysics_stage_q'
alias = ROOT / 'external/sv13q/svMultiPhysics-reuse'
assert source.is_dir()
assert not alias.is_symlink(), 'Refusing to modify a symbolic link'
expected = json.loads((ROOT/'reports/sv1_3q/source_patch.json').read_text())['after']
attributes = json.loads((ROOT/'sync_metadata/upstream_attribute_map.json').read_text())
pointers = {r['upstream_relative_path']:r for r in json.loads((ROOT/'sync_metadata/upstream_lfs_pointer_inventory.json').read_text())['files']}
for name, digest in expected.items():
    original = ROOT/attributes[name] if name in attributes else source/name
    data = pointers[name]['original_pointer_text'].encode() if name in pointers else original.read_bytes()
    assert hashlib.sha256(data).hexdigest() == digest, name
    dest = alias/name
    if dest.exists():
        assert hashlib.sha256(dest.read_bytes()).hexdigest() == digest, 'Refusing to overwrite: '+name
    else:
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_bytes(data)
print('Legacy Stage Q tests resolve to the committed source snapshot inside this checkout.')
