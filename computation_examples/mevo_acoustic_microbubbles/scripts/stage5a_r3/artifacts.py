"""Strict, size-bounded R3 inventory. Run only after all scientific jobs finish."""

import argparse
import gzip
import json
import subprocess
from pathlib import Path
from scripts.stage5a_r3.common import ROOT, BASE, read, write, sha
from scripts.stage5a_r3.current import campaign, protocol_path

EXCLUDED = {
    "ARTIFACT_MANIFEST.json",
    "SERVER_ARTIFACT_MANIFEST.json",
    "GIT_ARTIFACTS.txt",
    "JSON_AUDIT.json",
    "SERVER_JSON_AUDIT.json",
    "DELIVERY_IDENTITY.json",
}
RAW_LEDGER_NAMES = {
    "COMPLETED_TRIALS.json",
    "RESOURCE_SAMPLES.json",
    "ACCEPTED.jsonl",
    "REJECTED.jsonl",
    "TRANSLATION_SOLVES.jsonl",
}


def strict(path):
    def invalid(value):
        raise ValueError("Nonstandard JSON value " + value)

    with path.open() as f:
        if path.suffix == ".jsonl":
            for line in f:
                if line.strip():
                    json.loads(line, parse_constant=invalid)
        else:
            json.load(f, parse_constant=invalid)


def compress_ledgers():
    """Lossless deterministic compression; original scientific ledgers untouched."""
    for folder in [BASE / "runs", BASE / "benchmarks/frozen_reference", BASE / "benchmarks/r2_protected"]:
        for p in folder.glob("*/*"):
            if p.is_file() and p.name in RAW_LEDGER_NAMES:
                out = BASE / "review/compact_ledgers" / p.parent.name / (p.name + ".gz")
                out.parent.mkdir(parents=True, exist_ok=True)
                raw = p.read_bytes()
                encoded = gzip.compress(raw, mtime=0)
                if not out.exists() or out.read_bytes() != encoded:
                    out.write_bytes(encoded)


def external_inputs():
    references = {}
    for p in BASE.rglob("IDENTITY.json"):
        for name, digest in read(p).get("input_sha256", {}).items():
            references.setdefault(name, {}).setdefault(digest, []).append(str(p))
    rows = []
    for name, hashes in sorted(references.items()):
        p = Path(name)
        assert p.is_file(), f"Missing actual input: {p}"
        digest = sha(p)
        assert set(hashes) == {digest}, f"Input hash mismatch: {p}"
        rows.append(
            dict(
                actual_server_path=name,
                bytes=p.stat().st_size,
                sha256=digest,
                referenced_by=hashes[digest],
                git_selected=False,
            )
        )
    return rows


def selected(path):
    relative = path.relative_to(BASE)
    if relative.parts[0] == "review" and path.suffix == ".h5" and path.stat().st_size <= 4 * 2**20:
        return True
    if path.stat().st_size > 4 * 2**20 or path.suffix in [".h5", ".mp4", ".zip", ".pstats"]:
        return False
    if path.name in RAW_LEDGER_NAMES or path.name == "SOURCE.zip":
        return False
    if path.name.startswith("REJECTED_CANDIDATE_") or path.name.startswith("FAILED_CANDIDATE_"):
        return False
    if path.name.endswith(".log") or path.name.startswith("DYNAMIC_INSTRUMENTED_"):
        return False
    if any(part.startswith("diagnostic_render_") for part in relative.parts):
        return False
    # Legacy large HDF5 remains external; its compact input/identity records remain.
    return True


def main(server=False):
    assert campaign()["complete"], "Do not hash or package mutable campaign outputs"
    if not server:
        compress_ledgers()
    rows, invalid, checked = [], [], 0
    for p in sorted(BASE.rglob("*")):
        if not p.is_file() or p.name in EXCLUDED or "__pycache__" in p.parts:
            continue
        if p.suffix in [".so", ".pyc"] or any("compiled" in part for part in p.parts):
            continue
        if p.suffix in [".json", ".jsonl"]:
            checked += 1
            try:
                strict(p)
            except Exception as exc:
                invalid.append(dict(path=str(p), error=repr(exc)))
        rows.append(
            dict(
                relative_path=str(p.relative_to(ROOT)),
                actual_path=str(p),
                bytes=p.stat().st_size,
                sha256=sha(p),
                git_selected=selected(p),
            )
        )
    audit_name = "SERVER_JSON_AUDIT.json" if server else "JSON_AUDIT.json"
    write(
        BASE / "final" / audit_name,
        dict(status="FAIL" if invalid else "PASS", files_checked=checked, invalid=invalid),
    )
    assert not invalid, invalid[:3]
    if server:
        write(
            BASE / "final/SERVER_ARTIFACT_MANIFEST.json",
            dict(
                root=str(ROOT),
                files=rows,
                external_input_files=external_inputs(),
                execution_snapshot_source=sha(__file__),
            ),
        )
        print(dict(server_files=len(rows), checked_JSON=checked))
        return
    remote_path = BASE / "final/SERVER_ARTIFACT_MANIFEST.json"
    remote = read(remote_path)
    indexed = {r["relative_path"]: r for r in rows}
    server_only = []
    for r in remote["files"]:
        if r["relative_path"] in indexed:
            q = indexed[r["relative_path"]]
            q.update(
                server_path=r["actual_path"],
                server_sha256=r["sha256"],
                server_copy_identical=q["sha256"] == r["sha256"],
            )
        else:
            server_only.append(dict(r, git_selected=False, local_copy_exists=False))
    write(
        BASE / "final/ARTIFACT_MANIFEST.json",
        dict(
            root=str(ROOT),
            files=rows,
            server_only_files=server_only,
            external_input_files=remote["external_input_files"],
            excluded_from_self_hash=sorted(EXCLUDED),
            remote_manifest_sha256=sha(remote_path),
            policy="Raw trajectories, complete rejected-candidate arrays, source ZIPs and videos stay at verified actual paths. Lightweight fixtures, audit JSON, figures/data and losslessly compressed full accepted/rejected ledgers are selected. No scientific record is downsampled in place.",
        ),
    )
    files = [r["relative_path"] for r in rows if r["git_selected"]]
    files += [
        str((BASE / "final" / name).relative_to(ROOT))
        for name in ["ARTIFACT_MANIFEST.json", "JSON_AUDIT.json", "DELIVERY_IDENTITY.json"]
    ]
    (BASE / "final/GIT_ARTIFACTS.txt").write_text("\n".join(files) + "\n")
    write(
        BASE / "final/DELIVERY_IDENTITY.json",
        dict(
            branch=subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
            actual_R2_archive_base="07db8ed35f3a460a4ba315e42b15d257f8518b0a",
            formal_execution_commits={
                "R3_V1": "fa2fba5608037ed45143c8ede501cadf415e9180",
                "R3_V2": "b5079ad19be81ffe9e3d4f5b7c548717b3d43d67",
            },
            early_precision_execution="5219d92db731df1f736dd4d59bc1a98c98609e8f",
            abbreviated_diagnostic_execution_alias={"933f370c3": "933f370c3ad4506848bfb867faff1e6c2c906ba4"},
            archive_parent_commit=subprocess.check_output(
                ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
            ).strip(),
            final_archive_commit="Resolve the commit containing this file; no recursive self-commit hash",
            formal_protocol_sha256=sha(protocol_path()),
            environment_sha256=sha(BASE / "environment/IDENTITY.json"),
            no_production_promotion=True,
        ),
    )
    print(
        dict(
            local_files=len(rows),
            server_only_files=len(server_only),
            git_files=len(files),
            selected_MiB=sum(r["bytes"] for r in rows if r["git_selected"]) / 2**20,
        )
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--server", action="store_true")
    main(parser.parse_args().server)
