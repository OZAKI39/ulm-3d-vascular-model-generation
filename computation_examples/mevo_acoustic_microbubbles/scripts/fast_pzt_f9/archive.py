"""Pull closed F9 results, verify server bytes and recorded Git source identities."""

import json
import shlex
import subprocess
import hashlib
from datetime import datetime, timezone
from .common import ROOT, BASE, sha, js, write

REMOTE = "/workspace/mevo_acoustic_microbubbles/f9/results/fast_pzt/f9"
SSH = ["ssh", "-o", "ConnectTimeout=20", "mevo-new"]


def main():
    # Never declare an in-progress HDF5 copy to be a closed scientific archive.
    code = """import json,hashlib
from pathlib import Path
root=Path(REMOTE)
for label in ['y0_5_v1','y1_10_v1','y1_5_v1','y2_10_v1']:
 assert json.loads((root/label/'COMPLETE.json').read_text())['status']=='COMPLETE'
rows={}
for p in sorted(root.rglob('*')):
 if not p.is_file() or 'dependencies' in p.relative_to(root).parts: continue
 rows[str(p.relative_to(root))]={'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
print(json.dumps(rows))
""".replace("REMOTE", repr(REMOTE))
    remote = json.loads(subprocess.check_output(SSH + ["python3 -c " + shlex.quote(code)], text=True))
    subprocess.run(
        ["rsync", "-a", "--exclude=dependencies", "mevo-new:" + REMOTE + "/", str(BASE) + "/"], check=True
    )
    for name, row in remote.items():
        p = BASE / name
        assert p.stat().st_size == row["bytes"] and sha(p) == row["sha256"], name
    cache = {}
    identities = []
    for p in sorted(list(BASE.rglob("IDENTITY.json")) + list(BASE.rglob("CHECKS.json"))):
        record = js(p)
        record = record.get("identity", record)
        if "execution_commit" not in record:
            continue
        commit = record["execution_commit"]
        for name, expected in record["source_sha256"].items():
            key = (commit, name)
            if key not in cache:
                blob = subprocess.check_output(
                    ["git", "show", commit + ":computation_examples/mevo_acoustic_microbubbles/" + name],
                    cwd=ROOT,
                )
                cache[key] = hashlib.sha256(blob).hexdigest()
            assert cache[key] == expected, (p, key)
        identities.append(dict(path=str(p.relative_to(BASE)), execution_commit=commit))
    report = dict(
        status="PASS",
        timestamp_UTC=datetime.now(timezone.utc).isoformat(),
        server=REMOTE,
        local=str(BASE),
        remote_file_count=len(remote),
        remote_bytes=sum(v["bytes"] for v in remote.values()),
        remote_files=remote,
        verified_execution_identities=identities,
        distinct_source_blobs=len(cache),
        omitted="Runtime/build caches excluded. All produced F9 scientific trajectories, diagnostic states and logs archived. No new historical dependency trees copied.",
    )
    write(BASE / "ARCHIVE_VERIFICATION.json", report)
    write(
        BASE / "ARTIFACT_MANIFEST.json",
        dict(
            timestamp_UTC=report["timestamp_UTC"],
            scope="F9 result files excluding this self-referential manifest; dependency builds excluded",
            files={
                str(p.relative_to(BASE)): dict(bytes=p.stat().st_size, sha256=sha(p))
                for p in sorted(BASE.rglob("*"))
                if p.is_file()
                and p.name != "ARTIFACT_MANIFEST.json"
                and "dependencies" not in p.relative_to(BASE).parts
            },
        ),
    )
    print(
        json.dumps(
            {k: v for k, v in report.items() if k not in ["remote_files", "verified_execution_identities"]},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
