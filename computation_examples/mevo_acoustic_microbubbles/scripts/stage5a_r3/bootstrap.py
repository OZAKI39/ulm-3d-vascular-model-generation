"""Verify the archived baseline and inventory this isolated environment."""

import importlib.metadata as md
import os
import platform
import shutil
import subprocess
from pathlib import Path
from scripts.stage5a_r3.common import ROOT, BASE, R2ROOT, R2, read, write, sha, B2_RUNS


def main():
    for folder in [
        "environment",
        "fixtures",
        "diagnosis",
        "benchmarks",
        "runs",
        "comparisons",
        "figures",
        "final",
    ]:
        (BASE / folder).mkdir(parents=True, exist_ok=True)
    base = read(BASE / "final/R3_BASE_IDENTITY.json")
    for name, digest in base["base_evidence_sha256"].items():
        assert sha(R2 / name) == digest, name
    test = read(R2 / "final/TEST_GATE_V4A.json")
    for field in ["source_sha256", "orchestration_sha256", "test_source_sha256"]:
        for name, digest in test[field].items():
            assert sha(R2ROOT / name) == digest, name
            assert sha(ROOT / name) == digest, name
    manifest = read(R2 / "final/ARTIFACT_MANIFEST.json")
    rows = {r["relative_path"]: r for r in manifest["files"]}
    verified = []
    for name in B2_RUNS.values():
        for file in [
            "CHECKPOINT.npz",
            "trajectories.h5",
            "ACCEPTED.jsonl",
            "REJECTED.jsonl",
            "IDENTITY.json",
            "INPUT_SNAPSHOT.npz",
        ]:
            rel = f"results/stage5/motion_solver_r2/runs/{name}/{file}"
            p = R2ROOT / rel
            assert sha(p) == rows[rel]["sha256"], str(p)
            verified.append(dict(path=str(p), bytes=p.stat().st_size, sha256=sha(p)))
    write(BASE / "fixtures/VERIFIED_R2_INPUTS.json", verified)
    lock = subprocess.check_output(
        ["uv", "pip", "freeze", "--python", str(ROOT / ".venv/bin/python")], text=True
    )
    (ROOT / "requirements-stage5a-r3.lock.txt").write_text(lock)
    packages = [
        dict(
            name=d.metadata["Name"],
            version=d.version,
            license=d.metadata.get("License-Expression") or d.metadata.get("License"),
            license_files=d.metadata.get_all("License-File"),
        )
        for d in md.distributions()
    ]
    write(BASE / "environment/PACKAGES_LICENSES.json", packages)
    data = dict(
        python=platform.python_version(),
        executable=os.sys.executable,
        affinity=sorted(os.sched_getaffinity(0)),
        cpu_max=Path("/sys/fs/cgroup/cpu.max").read_text().strip(),
        memory_max=Path("/sys/fs/cgroup/memory.max").read_text().strip(),
        gpu=subprocess.check_output(
            [
                "nvidia-smi",
                "--query-gpu=name,driver_version,memory.total,memory.used",
                "--format=csv,noheader",
            ],
            text=True,
        ),
        compiler=subprocess.check_output(["gcc", "--version"], text=True).splitlines()[0],
        disk_free_bytes=shutil.disk_usage(ROOT).free,
        dependencies_sha256=sha(ROOT / "requirements-stage5a-r3.lock.txt"),
        packages_sha256=sha(BASE / "environment/PACKAGES_LICENSES.json"),
        build="PyPI binary wheels; Numba float64 without fastmath; exact-rational diagnostics use Python fractions",
        existing_processes=subprocess.check_output(
            ["ps", "-eo", "pid,ppid,pcpu,pmem,comm", "--sort=-pcpu"], text=True
        ).splitlines()[:20],
    )
    write(BASE / "environment/IDENTITY.json", data)
    print("Verified R2 archive and frozen code; isolated R3 environment:", data["python"], flush=True)


if __name__ == "__main__":
    main()
