"""Read-only audit of inherited archive hashes on both local disk and the server."""

import hashlib
from pathlib import Path
import subprocess
from .common import BASE, ROOT, js, write


def main():
    manifests = {
        name: js(ROOT / "results/fast_pzt" / name / "ARTIFACT_MANIFEST.json")["artifacts"]
        for name in ["f1", "f2", "f3"]
    }
    local = {}
    for name, rows in manifests.items():
        missing = []
        mismatch = []
        for row in rows:
            p = Path(row["local_path"])
            if not p.exists():
                missing.append(str(p))
            elif hashlib.sha256(p.read_bytes()).hexdigest() != row["sha256"]:
                mismatch.append(str(p))
        local[name] = dict(checked=len(rows), missing=missing, mismatched=mismatch)
    script = """import json,hashlib,pathlib,platform,subprocess
manifests=DATA
result={}
for name,rows in manifests.items():
 missing=[];mismatch=[]
 for row in rows:
  p=pathlib.Path(row['server_path'])
  if not p.exists():missing.append(str(p))
  elif hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:mismatch.append(str(p))
 result[name]=dict(checked=len(rows),missing=missing,mismatched=mismatch)
print(json.dumps(result))
""".replace("DATA", repr(manifests))
    raw = subprocess.check_output(["ssh", "vast4090", "python3 -"], input=script, text=True)
    import json

    server = json.loads(raw.strip().splitlines()[-1])
    write(
        BASE / "PRESERVATION.json",
        dict(
            local=local,
            server=server,
            scope="All inherited manifest entries, including files excluded from Git; scientific files and logs.",
        ),
    )
    assert all(
        not row["missing"] and not row["mismatched"] for side in [local, server] for row in side.values()
    )
    print("Inherited F1/F2/F3 preservation PASS", sum(len(v) for v in manifests.values()), "files per host")


if __name__ == "__main__":
    main()
