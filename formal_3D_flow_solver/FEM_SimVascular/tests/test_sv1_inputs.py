import json
from pathlib import Path
from sv_validation.provenance import sha256

def test_formal_copies_match_readonly_sources(root):
    manifest=json.loads((root/"inputs/MANIFEST.json").read_text())
    assert len(manifest["files"]) == 5
    for item in manifest["files"].values():
        assert sha256(Path(item["source_path"])) == item["source_sha256"] == item["copied_sha256"] == sha256(root/item["copied_path"])

def test_active_cleanup_complete(root, report):
    audit=report("cleanup_summary")
    assert audit["status"] == "PASS"
    assert audit["active_text_scan"]["status"] == "PASS"
    assert not (root/"scripts/remote").exists()
    assert sorted(p.name for p in (root/"reports/sv0").iterdir()) == ["REPORT.md"]
