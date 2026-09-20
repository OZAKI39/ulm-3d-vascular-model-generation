from pathlib import Path
from sv_validation.provenance import sha256
from sv_validation.validation import parse_result_vtu

def test_native_executable_and_official_fluid_smoke(root, report):
    n=report("native_solver"); s=report("official_smoke")
    assert n["status"] == "PASS" and n["linked_libraries_complete"]
    assert sha256(Path(n["executable"])) == n["executable_sha256"]
    assert n["commit"] == (root/"external/svMultiPhysics_commit.txt").read_text().strip()
    assert not n["mesh_library_recompiled"]
    assert s["exit_code"] == 0 and s["status"] == "PASS"
    for f in s["result_files"]:
        assert sha256(root/f["path"]) == f["sha256"]
        parse_result_vtu(root/f["path"])
