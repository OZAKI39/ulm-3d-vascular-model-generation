"""Verify the curated F10S synchronization without running scientific calculations."""

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
from export_sync_20261015 import PARENT, PREFIX, SOURCE

ROOT = Path(__file__).resolve().parents[1]
REPO = ROOT.parents[1]
OUT = ROOT / "sync/20261015"


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

    for stage in ["f10r", "f10s"]:
        base = ROOT / "results/straight_tube_exp1" / stage
        remote = js(base / "REMOTE_MANIFEST.json")["files"]
        for name, entry in remote.items():
            assert sha(base / name) == entry["sha256"], name
        checks[stage + "_all_returned_server_files"] = len(remote)
        manifest = js(base / "DELIVERY_MANIFEST.json")["files"]
        omissions = []
        for name, entry in manifest.items():
            if name in excluded:
                assert excluded[name]["sha256"] == entry["sha256"]
                omissions.append(name)
            else:
                assert sha(ROOT / name) == entry["sha256"], name
        checks[stage + "_full_delivery_files_verified"] = len(manifest) - len(omissions)
        checks[stage + "_omitted_preview_files"] = omissions
        reg = js(base / "PREREGISTRATION.json")
        for name, digest in reg["input_sha256"].items():
            assert sha(ROOT / name) == digest, name
        checks[stage + "_all_registered_inputs_present"] = len(reg["input_sha256"])
        figs = js(base / "figures/FIGURE_PROVENANCE.json")
        for name, digest in figs["input_sha256"].items():
            assert sha(ROOT / name) == digest, name
        for name, digest in figs.get("figures_sha256", figs.get("output_sha256", {})).items():
            assert sha(base / "figures" / name) == digest, name
        checks[stage + "_four_PNG_PDF_and_all_plot_inputs"] = True
    gate = js(ROOT / "results/straight_tube_exp1/f10s/analysis/FINAL_GATE.json")
    assert gate["groups"]["X_PAIR_3MS_DT1"]["TARGET_REACHED"] == "PASS"
    assert gate["groups"]["X_PAIR_3MS_DT1"]["STRICT_ENDPOINT_REACHED"] == "PASS"
    for group in gate["groups"].values():
        for key in ["BASIC_NUMERICS_PASS", "GROUP_POSITION_NUMERICALLY_ADMITTED", "GROUP_SHAPE_NUMERICALLY_ADMITTED", "GROUP_TOPOLOGY_NUMERICALLY_ADMITTED", "PARTICLE_ID_NUMERICALLY_ADMITTED"]:
            assert group[key] == "PASS"
    assert gate["PDMS_FIELD_ADMITTED"] == "FAIL" and gate["PDMS_FIELD"] == "NOT_RUN"
    assert gate["EXPERIMENTALLY_CALIBRATED"] == "NOT_VERIFIED" and not gate["training_ready"]
    old_gate = js(ROOT / "results/straight_tube_exp1/f10q/analysis/FINAL_GATE.json")
    assert old_gate["groups"]["NAV_Y"]["GROUP_SHAPE_NUMERICALLY_ADMITTED"] == "FAIL"
    assert old_gate["groups"]["STRONG"]["BASIC_NUMERICS_PASS"] == "FAIL"
    checks["current_successes_historical_failures_and_uncalibrated_scope_preserved"] = True
    # Read compact headers for all newly retained trajectories; no solver invocation.
    for stage in ["f10r", "f10s"]:
        paths = list((ROOT / "results/straight_tube_exp1" / stage).glob("*/runs/*/trajectory.h5"))
        for path in paths:
            with h5py.File(path) as h:
                np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
                assert h.attrs["flow"] == "STRICT_ZERO"
                assert h["positions_m"].shape[1:] == (32,3)
        checks[stage + "_original_32ID_zero_flow_trajectories"] = len(paths)

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
        omit = {PREFIX + "sync/20261015/PACKAGE_MANIFEST.json", PREFIX + "sync/20261015/VERIFICATION.json"}
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
        write(OUT / "VERIFICATION.json", dict(status="PASS", date="2026-10-15", checks=checks,
              package_manifest_sha256=sha(OUT / "PACKAGE_MANIFEST.json"),
              note="Packaging/data-readability verification only; no new scientific simulation or physical calibration."))
    else:
        assert sha(OUT / "PACKAGE_MANIFEST.json") == js(OUT / "VERIFICATION.json")["package_manifest_sha256"]
    print(json.dumps(dict(status="PASS", checks=checks), indent=2))


if __name__ == "__main__":
    main()
