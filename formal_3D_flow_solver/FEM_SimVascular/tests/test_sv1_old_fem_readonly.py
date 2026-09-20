from sv_validation.provenance import sha256

def test_complete_old_project_audit(root, report):
    audit=report("old_fem_readonly_audit")
    assert audit["status"] == "PASS" and audit["git_state_unchanged"]
    assert audit["entries"] == len(report("old_fem_baseline")["files"])
    assert sha256(root/"reports/sv1/old_fem_baseline.json.gz") == audit["baseline_sha256"]
    assert sha256(root/"reports/sv1/old_fem_final.json.gz") == audit["final_sha256"]
