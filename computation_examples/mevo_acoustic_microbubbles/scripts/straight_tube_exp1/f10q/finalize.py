"""Finalize verified engineering metadata, portable source and artifact manifest."""

import subprocess
import tarfile
import io
from .common import ROOT, BASE, read, write, sha, identity


def main():
    verification = read(BASE / "VERIFICATION.json")
    assert verification["status"] == "PASS"
    g = read(BASE / "analysis/FINAL_GATE.json")
    g["ENGINEERING_PASS"] = "PASS"
    g["engineering_verification_sha256"] = sha(BASE / "VERIFICATION.json")
    write(BASE / "analysis/FINAL_GATE.json", g)
    costs = read(BASE / "analysis/COSTS.json")
    costs["additional_readonly_verification_s"] = sum(
        read(BASE / f)["wall_s"] for f in ["PRELIMINARY_VERIFICATION.json", "VERIFICATION.json"]
    )
    costs["numerical_plus_verification_s"] = (
        costs["total_numerical_s"] + costs["additional_readonly_verification_s"]
    )
    assert costs["numerical_plus_verification_s"] < 600
    write(BASE / "analysis/COSTS.json", costs)
    # These fields do not enter any plotted panel; preserve original render hashes too.
    p = BASE / "figures/FIGURE_PROVENANCE.json"
    v = read(p)
    v["render_time_input_sha256"] = v["input_sha256"].copy()
    for relative in list(v["input_sha256"]):
        v["input_sha256"][relative] = sha(ROOT / relative)
    v["post_render_metadata_update"] = (
        "Only ENGINEERING_PASS/verification hash and additional audit timing appended after successful final checks. All plotted controls, trajectories, metrics and route timings unchanged; images remain the visually reviewed files."
    )
    write(p, v)
    commit = identity()["execution_commit"]
    prefix = "computation_examples/mevo_acoustic_microbubbles/"
    paths = ["src", "configs/straight_tube_exp1", "scripts/straight_tube_exp1", "scripts/fast_pzt/python.sh"]
    paths += [str(p.relative_to(ROOT)) for p in (ROOT / "docs").glob("*F10Q*.md")]
    repo = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=ROOT, text=True).strip()
    raw = subprocess.check_output(["git", "archive", commit, *[prefix + p for p in paths]], cwd=repo)
    buf = io.BytesIO()
    with tarfile.open(fileobj=io.BytesIO(raw)) as source, tarfile.open(fileobj=buf, mode="w:gz") as target:
        for member in source:
            if not member.name.startswith(prefix):
                continue
            member.name = member.name[len(prefix) :]
            if member.name:
                target.addfile(member, source.extractfile(member) if member.isfile() else None)
    path = BASE / "execution" / ("delivery_source_" + commit + ".tar.gz")
    assert not path.exists()
    path.write_bytes(buf.getvalue())
    write(
        BASE / "DELIVERY.json",
        dict(
            branch="dev/fast-pzt-f10q-sequential-qp-xy-navigation-20261014",
            source_commit=commit,
            source_archive=str(path.relative_to(BASE)),
            source_archive_sha256=sha(path),
            ENGINEERING_PASS="PASS",
            numeric_seconds=costs["total_numerical_s"],
            numeric_and_verification_seconds=costs["numerical_plus_verification_s"],
            remote_returned=read(BASE / "REMOTE_MANIFEST.json")["bytes"],
            all_jobs_completed=True,
            labels=g["labels"],
            training_ready=False,
        ),
    )
    files = [
        p
        for p in BASE.rglob("*")
        if p.is_file() and "__pycache__" not in str(p) and p.name != "MANIFEST.json"
    ]
    files += [
        p
        for folder in ["scripts/straight_tube_exp1/f10q", "configs/straight_tube_exp1/f10q"]
        for p in (ROOT / folder).glob("*")
        if p.is_file()
    ]
    files += list((ROOT / "docs").glob("*F10Q*.md")) + [
        ROOT / "src/mevo_mb/fast_sim/sequential_hard_sphere.py",
        ROOT / "src/mevo_mb/fast_sim/water_sequential.py",
    ]
    table = {str(p.relative_to(ROOT)): dict(bytes=p.stat().st_size, sha256=sha(p)) for p in sorted(files)}
    write(
        BASE / "MANIFEST.json",
        dict(
            source_commit=commit,
            files=table,
            count=len(table),
            bytes=sum(r["bytes"] for r in table.values()),
            scope="All F10Q results, logs, inputs, immutable source packages, figures, reports and new source. Self excluded; protected historical inputs remain at original hashed paths.",
        ),
    )
    print(
        "Final manifest",
        len(table),
        "files",
        sum(r["bytes"] for r in table.values()),
        "bytes; total including audits",
        costs["numerical_plus_verification_s"],
    )


if __name__ == "__main__":
    main()
