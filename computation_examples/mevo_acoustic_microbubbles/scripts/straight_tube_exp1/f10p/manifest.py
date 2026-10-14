"""Final file-level integrity manifest, excluding only this manifest itself."""

from .common import BASE, ROOT, CONFIG, read, write, sha, identity


def main():
    assert read(BASE / "analysis/FINAL_VERIFICATION.json")["status"] == "PASS"
    files = [f for f in BASE.rglob("*") if f.is_file() and f.name != "ARTIFACT_MANIFEST.json"]
    files += list((ROOT / "scripts/straight_tube_exp1/f10p").glob("*.py"))
    files += list(CONFIG.parent.glob("*.json"))
    files += list((ROOT / "docs").glob("*F10P*.md"))
    files += [ROOT / "src/mevo_mb/fast_sim" / name for name in ["hard_sphere.py", "water_hard_sphere.py"]]
    records = {
        str(f.relative_to(ROOT)): dict(sha256=sha(f), bytes=f.stat().st_size)
        for f in sorted(set(files))
        if "__pycache__" not in str(f)
    }
    write(
        BASE / "ARTIFACT_MANIFEST.json",
        dict(
            identity=identity(),
            files=records,
            file_count=len(records),
            bytes=sum(r["bytes"] for r in records.values()),
            scope="F10P deliverables; inherited inputs retain F10N original hashes. Scientific failure preserved.",
        ),
    )
    print(len(records), "files", sum(r["bytes"] for r in records.values()), "bytes")


if __name__ == "__main__":
    main()
