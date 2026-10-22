"""Generate and verify a complete local delivery checksum manifest."""

from .common import BASE, ROOT, write, sha


def main():
    files = [
        p
        for p in BASE.rglob("*")
        if p.is_file() and p.name not in ["ARTIFACT_MANIFEST.json", "ARTIFACT_CHECK.json"]
    ]
    for folder, glob in [
        ("scripts/straight_tube_exp1/f10o", "*.py"),
        ("configs/straight_tube_exp1/f10o", "*"),
        ("docs", "*F10O*.md"),
    ]:
        files.extend(p for p in (ROOT / folder).glob(glob) if p.is_file())
    files.append(ROOT / "src/mevo_mb/fast_sim/water_b2.py")
    data = {str(p.relative_to(ROOT)): dict(sha256=sha(p), bytes=p.stat().st_size) for p in sorted(set(files))}
    write(BASE / "ARTIFACT_MANIFEST.json", dict(files=data, bytes=sum(r["bytes"] for r in data.values())))
    for rel, r in data.items():
        assert sha(ROOT / rel) == r["sha256"]
    write(
        BASE / "analysis/ARTIFACT_CHECK.json",
        dict(
            status="PASS",
            files=len(data),
            bytes=sum(r["bytes"] for r in data.values()),
            manifest_sha256=sha(BASE / "ARTIFACT_MANIFEST.json"),
        ),
    )
    print("Verified", len(data), "delivery files")


if __name__ == "__main__":
    main()
