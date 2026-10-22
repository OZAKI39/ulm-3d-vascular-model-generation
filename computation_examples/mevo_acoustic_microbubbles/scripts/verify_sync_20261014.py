"""Verify the curated F10Q synchronization without running scientific calculations."""

import argparse
import ast
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tarfile

import h5py
import numpy as np

from verify_sync_20261013 import scan_credentials
from export_sync_20261014 import PARENT, PREFIX, SOURCE

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT / "sync/20261014"


def js(path):
    return json.loads(path.read_text())


def data(path):
    return os.readlink(path).encode() if path.is_symlink() else path.read_bytes()


def sha(path):
    return hashlib.sha256(data(path)).hexdigest()


def git(*args):
    return subprocess.check_output(["git", *args], cwd=REPO)


def write(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seal", action="store_true", help="Seal the staged package after all files are staged")
    args = ap.parse_args()
    checks = {}
    copied, excluded = js(OUT / "COPIED_FILES.json"), js(OUT / "EXCLUDED_FILES.json")
    policy = js(OUT / "EXPORT_POLICY.json")
    assert policy["complete_scientific_delivery"] == SOURCE
    for name, entry in copied.items():
        assert sha(ROOT / name) == entry["sha256"], name
        assert len(data(ROOT / name)) == entry["bytes"], name
    for name in excluded:
        assert not (ROOT / name).exists(), name
    assert len(copied) == policy["copied_count"]
    assert len(excluded) == policy["excluded_count"]
    checks["exact_scientific_delta_files"] = len(copied)
    checks["explicit_omissions"] = len(excluded)
    code = js(OUT / "SCIENTIFIC_CODE_MANIFEST.json")
    for name, entry in code.items():
        assert sha(ROOT / name) == entry["sha256"], name
    checks["exact_current_code_config_test_dependency_files"] = len(code)
    p = js(OUT / "PROGRESS_PROVENANCE.json")
    assert sha(ROOT / p["destination"]) == p["sha256"]
    checks["current_progress_copy_unchanged"] = True

    q = ROOT / "results/straight_tube_exp1/f10q"
    remote = js(q / "REMOTE_MANIFEST.json")["files"]
    for name, entry in remote.items():
        assert sha(q / name) == entry["sha256"], name
    checks["all_F10Q_remote_returned_files"] = len(remote)
    original = js(q / "MANIFEST.json")["files"]
    omitted_full = []
    for name, entry in original.items():
        if name in excluded:
            assert excluded[name]["sha256"] == entry["sha256"]
            omitted_full.append(name)
        else:
            assert sha(ROOT / name) == entry["sha256"], name
    checks["full_F10Q_manifest_files_verified"] = len(original) - len(omitted_full)
    checks["full_F10Q_manifest_omitted_preview_files"] = omitted_full
    reg = js(q / "PREREGISTRATION.json")
    for name, digest in reg["input_sha256"].items():
        assert sha(ROOT / "results/straight_tube_exp1/f10n" / name) == digest, name
    checks["all_F10Q_preregistered_inputs_present"] = len(reg["input_sha256"])
    figures = js(q / "figures/FIGURE_PROVENANCE.json")
    for name, digest in figures["input_sha256"].items():
        assert sha(ROOT / name) == digest, name
    for name, digest in figures["figures_sha256"].items():
        assert sha(q / "figures" / name) == digest, name
    checks["F10Q_four_PNG_PDF_and_all_figure_inputs"] = True
    summary = ROOT / "results/straight_tube_exp1/f10q_parameter_summary_20261014"
    for name, digest in js(summary / "MANIFEST.json").items():
        assert sha(summary / name) == digest, name
    for name, digest in js(summary / "SUMMARY.json")["input_sha256"].items():
        assert sha(ROOT / name) == digest, name
    checks["all_step_PZT_force_summary_and_inputs"] = True
    gate = js(q / "analysis/FINAL_GATE.json")
    assert gate["ENGINEERING_PASS"] == "PASS" and not gate["training_ready"]
    assert gate["groups"]["NAV_X"]["TARGET_REACHED"]["X"]
    assert not gate["groups"]["NAV_Y"]["TARGET_REACHED"]["Y"]
    assert gate["groups"]["NAV_Y"]["GROUP_SHAPE_NUMERICALLY_ADMITTED"] == "FAIL"
    assert gate["groups"]["STRONG"]["BASIC_NUMERICS_PASS"] == "FAIL"
    assert gate["EXPERIMENTALLY_CALIBRATED"] == "NOT_VERIFIED"
    checks["successes_failures_and_uncalibrated_scope_preserved"] = True

    if args.seal:
        old = {}
        for line in git("ls-tree", "-r", PARENT).decode().splitlines():
            meta, name = line.split("\t", 1)
            old[name] = meta.split()[2]
        current = {}
        for line in git("ls-files", "--stage").decode().splitlines():
            meta, name = line.split("\t", 1)
            mode, oid, stage = meta.split()
            assert stage == "0"
            current[name] = oid
        assert set(old).issubset(current), "Historical files must not be deleted"
        changed = {n for n in current if current[n] != old.get(n)}
        assert all(n == "MEVO_SNAPSHOT.md" or n.startswith(PREFIX) for n in changed)
        for name in changed:
            assert git("hash-object", "--", name).decode().strip() == current[name], "Index differs from working file: " + name
        omit = {PREFIX + "sync/20261014/PACKAGE_MANIFEST.json", PREFIX + "sync/20261014/VERIFICATION.json"}
        files = {n: dict(bytes=len(data(REPO / n)), sha256=sha(REPO / n)) for n in sorted(changed - omit)}
    else:
        files = js(OUT / "PACKAGE_MANIFEST.json")["files"]

    legacy_literals = {PREFIX + "results/straight_tube_exp1/f10o/cleanup/F10N_VERIFY.json"}
    checks["preserved_legacy_python_dict_text_with_json_extension"] = sorted(legacy_literals)
    totals = dict(json_files=0, hdf5_files=0, hdf5_datasets=0, npz_files=0, npz_arrays=0,
                  source_archives=0, archive_members=0, python_files_parsed=0, text_files_scanned=0)
    for name, entry in files.items():
        path = REPO / name
        assert sha(path) == entry["sha256"] and len(data(path)) == entry["bytes"], name
        assert entry["bytes"] < policy["maximum_new_file_bytes"], name
        assert not any(part in {".venv", "__pycache__", ".ssh", ".cache"} for part in Path(name).parts)
        assert path.suffix not in {".mp4", ".avi", ".mov", ".whl", ".pyc"}
        if path.suffix == ".json":
            if name in legacy_literals:
                assert isinstance(ast.literal_eval(path.read_text()), dict)
            else:
                js(path)
                totals["json_files"] += 1
        if path.suffix == ".py":
            ast.parse(path.read_bytes(), filename=name)
            totals["python_files_parsed"] += 1
        if path.suffix == ".h5":
            with h5py.File(path) as handle:
                def read_dataset(_name, obj):
                    if isinstance(obj, h5py.Dataset):
                        _ = obj[()]
                        totals["hdf5_datasets"] += 1
                handle.visititems(read_dataset)
            totals["hdf5_files"] += 1
        if path.suffix == ".npz":
            with np.load(path, allow_pickle=False) as arrays:
                for key in arrays.files:
                    _ = arrays[key]
                    totals["npz_arrays"] += 1
            totals["npz_files"] += 1
        if path.name.endswith(".tar.gz"):
            with tarfile.open(path, "r:gz") as archive:
                for member in archive:
                    assert not Path(member.name).is_absolute() and ".." not in Path(member.name).parts
                    if member.isfile():
                        raw = archive.extractfile(member).read()
                        scan_credentials(name + ":" + member.name, raw)
                        totals["archive_members"] += 1
                    else:
                        assert member.isdir(), "Unexpected non-source archive member"
            totals["source_archives"] += 1
        if path.suffix in {".py", ".json", ".jsonl", ".txt", ".log", ".md", ".yaml", ".yml", ".toml", ".sh", ".diff", ".csv", ".conf"}:
            scan_credentials(name, data(path))
            totals["text_files_scanned"] += 1
    checks.update(totals)
    checks.update(package_files=len(files), package_bytes=sum(e["bytes"] for e in files.values()),
                  maximum_new_file_bytes=max(e["bytes"] for e in files.values()),
                  scientific_computation_started=False, unrelated_project_changes=False)
    if args.seal:
        write(OUT / "PACKAGE_MANIFEST.json", dict(parent=PARENT, source_delivery=SOURCE,
              scope="Actual new/updated Git payload; excludes this manifest and verification to avoid circular hashes. Published parent otherwise inherited.", files=files))
        write(OUT / "VERIFICATION.json", dict(status="PASS", date="2026-10-14", checks=checks,
              package_manifest_sha256=sha(OUT / "PACKAGE_MANIFEST.json"),
              note="Packaging/data-readability verification only; no new scientific simulation or physical calibration."))
    else:
        assert sha(OUT / "PACKAGE_MANIFEST.json") == js(OUT / "VERIFICATION.json")["package_manifest_sha256"]
    print(json.dumps(dict(status="PASS", checks=checks), indent=2))


if __name__ == "__main__":
    main()
