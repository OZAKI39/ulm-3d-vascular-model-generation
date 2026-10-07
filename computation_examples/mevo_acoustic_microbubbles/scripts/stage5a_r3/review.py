"""Build a compact, reconstructible review bundle without changing raw results."""

import gzip
import subprocess
import zipfile
from scripts.stage5a_r3.common import ROOT, BASE, R2, B2_RUNS, write, sha
from scripts.stage5a_r3.current import campaign
from scripts.stage5a_r2.motion_visuals import wall_path

REQUIRED = [
    "R3_BASE_IDENTITY",
    "B2_FAILURE_REPRODUCTION",
    "B2_EVENT_ROOT_CAUSE",
    "WALL_CERTIFICATE_AUDIT",
    "REFERENCE_LINEAGE_AUDIT",
    "ALL_ACTIVE_TRAJECTORY_AUDIT",
    "S0_S1_USAGE_AND_REFINEMENT",
    "DYNAMIC_COST_PROFILE",
    "RESOURCE_LAYOUT_BENCHMARK",
    "MOTION_SOLVER_R3_GATE",
]


def main():
    assert campaign()["complete"]
    from scripts.stage5a_r3.artifacts import compress_ledgers, compact_comparisons

    compress_ledgers()
    compact_comparisons()
    from scripts.stage5a_r3.recompute_review import main as verify_compact

    verify_compact()
    review = BASE / "review"
    review.mkdir(exist_ok=True)
    rows = []

    def copy(source, relative, compress=False):
        out = review / relative
        out.parent.mkdir(parents=True, exist_ok=True)
        raw = source.read_bytes()
        out.write_bytes(gzip.compress(raw, mtime=0) if compress else raw)
        rows.append(
            dict(
                source=str(source),
                source_sha256=sha(source),
                review_path=str(out.relative_to(ROOT)),
                review_sha256=sha(out),
                bytes=out.stat().st_size,
                transformation="Lossless gzip, mtime=0" if compress else "Byte-identical copy",
            )
        )

    for kind, run_name in B2_RUNS.items():
        folder = R2 / "runs" / run_name
        for name in ["ACCEPTED.jsonl", "REJECTED.jsonl", "IDENTITY.json", "SUMMARY.json"]:
            p = folder / name
            compress = p.suffix == ".jsonl"
            copy(p, f"original_B2/{kind}/{name}" + (".gz" if compress else ""), compress)
        fixture = BASE / "diagnosis/reproduction_v1" / kind
        for p in sorted(fixture.iterdir()):
            if p.is_file() and p.suffix in [".npz", ".json"]:
                copy(p, f"B2_local_fixture/{kind}/{p.name}")
        local_run = next((BASE / "runs").glob(f"*local_v1_{kind}"))
        for name in ["ACCEPTED.jsonl", "REJECTED.jsonl", "IDENTITY.json", "SUMMARY.json", "trajectories.h5"]:
            copy(local_run / name, f"repaired_B2_terminal/{kind}/{name}")
    copy(wall_path(), "geometry/WALL.vtp")
    for p in (BASE / "fixtures").glob("B2_early_*.npz"):
        copy(p, "B2_early_actual_states/" + p.name)
    # Full source changes, not a fragment that hides alterations outside the solver.
    patch = subprocess.check_output(
        ["git", "diff", "07db8ed35f3a460a4ba315e42b15d257f8518b0a", "HEAD", "--", "."], cwd=ROOT
    )
    (review / "SOURCE_DIFF_FROM_R2.patch").write_bytes(patch)
    notes = """# Stage 5A-R3 lightweight scientific review

This directory contains real full-system B2 states, candidate QP/geometry arrays,
the actual open WALL mesh, and complete losslessly compressed original accepted
and rejected ledgers for both R2 failed runs. No particles are removed from the
fixtures. Local repaired trajectories start at the actual R2 terminal state;
they must not be mistaken for the separately recorded original-3-s formal runs.

Read the project docs/STAGE5A_R3_ROOT_CAUSE_ANALYSIS.md and final gate first.
The finite-feature geometric repair recovers B2 continuation. Independent
all-active global trajectory accuracy is not established; no promotion follows.
The early +30..50 us fixtures use distinct actual nominal/refined states with
explicit cold caches. A local difference does not identify either as truth.

REVIEW_SOURCES.json supplies hashes and exact copy/compression transformations.
SOURCE_DIFF_FROM_R2.patch records all project source changes through its named
parent commit. The containing Git archive also includes later packaging scripts.
The branch contains the required audit JSON, PNG/PDF figures and actual-state
NPZ data. Large trajectories, videos and every formal failed-candidate array are
indexed with actual paths/hashes by final/ARTIFACT_MANIFEST.json. No old failure
or gate is overwritten. V1 and V2 execution identities remain separate.

To read a retained original ledger: gzip -dc original_B2/nominal/REJECTED.jsonl.gz
For executable reproduction commands, see scripts/stage5a_r3/README.md.
"""
    (review / "README.md").write_text(notes)
    write(
        review / "REVIEW_SOURCES.json",
        dict(
            files=rows,
            source_diff_sha256=sha(review / "SOURCE_DIFF_FROM_R2.patch"),
            source_diff_through_commit=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
            ).strip(),
        ),
    )
    # ZIP paths preserve the project-relative tree. The directory is also checked
    # into Git; the standalone ZIP can remain a local/server downloadable asset.
    selected = set(
        p
        for p in review.rglob("*")
        if p.is_file() and p.suffix != ".zip" and p.name != "PACKAGE_IDENTITY.json"
    )
    selected.update(BASE / "final" / (name + ".json") for name in REQUIRED)
    selected.update((BASE / "final").glob("*PRECISION*.json"))
    selected.update((BASE / "environment").glob("JUNIT*.xml"))
    selected.update((BASE / "figures").glob("G*.png"))
    selected.update((BASE / "figures").glob("G*.pdf"))
    selected.update((BASE / "figures").glob("G*_SOURCES.json"))
    selected.update((BASE / "figures/data").glob("G*.npz"))
    selected.update((ROOT / "docs").glob("STAGE5A_R3*.md"))
    selected.update((ROOT / "scripts/stage5a_r3").glob("*.py"))
    selected.add(ROOT / "requirements-stage5a-r3.lock.txt")
    selected.add(ROOT / "scripts/stage5a_r3/README.md")
    output = BASE / "Stage5A_R3_review.zip"
    with zipfile.ZipFile(output, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for p in sorted(selected):
            assert p.is_file(), p
            z.write(p, str(p.relative_to(ROOT)))
    with zipfile.ZipFile(output) as z:
        assert z.testzip() is None
    write(
        review / "PACKAGE_IDENTITY.json",
        dict(
            path=str(output),
            bytes=output.stat().st_size,
            sha256=sha(output),
            member_count=len(selected),
            policy="Package excludes its own identity and final artifact manifest to avoid recursive hashes",
        ),
    )
    print(dict(package=str(output), MiB=output.stat().st_size / 2**20, copied_evidence=len(rows)))


if __name__ == "__main__":
    main()
