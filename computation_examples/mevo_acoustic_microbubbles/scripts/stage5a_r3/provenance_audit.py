"""Preserve a correction for early diagnostic live-deployment source snapshots."""

import hashlib
import subprocess
import zipfile
from pathlib import Path
from scripts.stage5a_r3.common import ROOT, BASE, read, write, sha


def main():
    repo = Path(subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=ROOT, text=True).strip())
    commit = subprocess.check_output(["git", "rev-parse", "bed890d80"], cwd=ROOT, text=True).strip()
    prefix = "computation_examples/mevo_acoustic_microbubbles/"
    paths = subprocess.check_output(
        ["git", "ls-tree", "-r", "--name-only", commit, "--", prefix + "src", prefix + "scripts/stage5a_r3"],
        cwd=repo,
        text=True,
    ).splitlines()
    content = {
        p[len(prefix) :]: subprocess.check_output(["git", "show", commit + ":" + p], cwd=repo)
        for p in paths
        if p.endswith(".py")
    }
    folder = BASE / "environment/source_corrections"
    folder.mkdir(parents=True, exist_ok=True)
    bundle = folder / "LOCAL_V1_LOADED_SOURCE.zip"
    with zipfile.ZipFile(bundle, "w", zipfile.ZIP_DEFLATED) as z:
        for name, data in content.items():
            z.writestr(name, data)
    rows = []
    driver = "scripts/stage5a_r3/run.py"
    for path in sorted((BASE / "runs").glob("*local_v1_B*L*/IDENTITY.json")):
        data = read(path)
        mismatches = {
            k: dict(loaded_sha256=hashlib.sha256(v).hexdigest(), snapshot_sha256=data["source_sha256"].get(k))
            for k, v in content.items()
            if k.startswith("src/") and data["source_sha256"].get(k) != hashlib.sha256(v).hexdigest()
        }
        rows.append(
            dict(
                run=str(path.parent),
                original_identity_sha256=sha(path),
                original_snapshot_driver_sha256=data["source_sha256"][driver],
                loaded_driver_sha256=hashlib.sha256(content[driver]).hexdigest(),
                numerical_kernel_mismatches=mismatches,
                driver_snapshot_matches=data["source_sha256"][driver]
                == hashlib.sha256(content[driver]).hexdigest(),
            )
        )
    assert rows and all(not r["numerical_kernel_mismatches"] for r in rows)
    report = dict(
        scope="Diagnostic local V1 campaign only; no formal campaign source is affected",
        status="CORRECTED_PROVENANCE_WITH_ORIGINALS_RETAINED",
        loaded_commit=commit,
        loaded_source_zip=str(bundle),
        bundle_sha256=sha(bundle),
        rows=rows,
        explanation="The local_campaign module imported its run function and numerical modules before forking two workers. During later levels, the live deployment received new driver/resource-accounting code. Snapshots read disk at each run start, so some describe the later on-disk driver rather than the already-loaded driver. No reload occurred; every original numerical kernel hash agrees. The supplied campaign-start source corrects provenance; original IDENTITY/SOURCE.zip files are unchanged. Subsequent layout, attribution, holdout and formal campaign use immutable execution snapshots.",
        proof_basis="Parent import order in bed890d80 local_campaign.py, first-run source snapshot, supervisor launch and run timestamps; no runtime reload in the campaign",
    )
    write(folder / "LOCAL_V1_PROVENANCE_CORRECTION.json", report)
    print(
        dict(runs=len(rows), driver_snapshot_mismatches=sum(not r["driver_snapshot_matches"] for r in rows)),
        flush=True,
    )


if __name__ == "__main__":
    main()
