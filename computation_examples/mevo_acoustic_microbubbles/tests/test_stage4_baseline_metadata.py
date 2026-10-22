import json
from pathlib import Path
import subprocess


def test_frozen_metadata_detects_modified_and_missing_files(tmp_path, monkeypatch):
    # CLI scripts are outside the installed src package.
    monkeypatch.syspath_prepend(str(Path(__file__).resolve().parents[1]))
    from scripts.stage4 import final_audit

    repo = tmp_path / "repo"
    repo.mkdir()
    names = ["PHYSICS_ACCEPTANCE.json", "docs/STAGE3_VALIDATION_REPORT.md",
             "docs/HISTORY_FORCE_SENSITIVITY.md", "docs/KNOWN_LIMITATIONS.md",
             "docs/METHODS_MODEL_DESCRIPTION.md"]
    for name in names:
        path = repo / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("Frozen test fixture: " + name + "\n")
    subprocess.run(["git", "init", "-q", str(repo)], check=True)
    subprocess.run(["git", "add", "."], cwd=repo, check=True)
    subprocess.run(["git", "-c", "user.name=Stage4 Test", "-c", "user.email=stage4-test@example.invalid",
                    "-c", "commit.gpgsign=false", "commit", "-qm", "Frozen metadata fixture"], cwd=repo, check=True)
    commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip()
    monkeypatch.setattr(final_audit, "ROOT", repo)
    out = tmp_path / "audit"
    out.mkdir()
    original = final_audit.stage3_metadata(commit, out)
    assert original["status"] == "PASS"
    assert len(original["files"]) == len(names)

    (repo / names[0]).write_text("Modified acceptance fixture\n")
    (repo / names[2]).unlink()
    changed = final_audit.stage3_metadata(commit, out)
    assert changed["status"] == "FAIL"
    assert {row["path"] for row in changed["files"] if not row["unchanged"]} == {names[0], names[2]}
    assert json.loads((out / "STAGE3_BASELINE_METADATA.json").read_text()) == changed
