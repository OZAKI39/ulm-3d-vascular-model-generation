import hashlib,json
from pathlib import Path
def test_protected_science_files_unchanged(report):
    before=json.loads((report/'data/particle_science_protection_manifest.json').read_text())
    assert len(before['files'])>100
    for name,spec in before['files'].items():
        p=Path(name);assert p.is_file() and p.stat().st_size==spec['size']
        assert hashlib.sha256(p.read_bytes()).hexdigest()==spec['sha256']
    snapshot=json.loads((report/'data/code_snapshot_manifest.json').read_text())
    for name,h in snapshot.items():assert hashlib.sha256((report/'server_bundle'/name).read_bytes()).hexdigest()==h
