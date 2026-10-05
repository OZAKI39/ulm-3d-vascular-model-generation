"""Hash the delivered source, Stage 3 artifacts and the seven referenced radial nodes."""

import hashlib
import json
from pathlib import Path
import subprocess

root = Path.cwd()
acceptance = json.loads(Path("PHYSICS_ACCEPTANCE.json").read_text())
for item in acceptance["items"].values():
    for evidence in item["evidence"]:
        path = Path(evidence["path"])
        assert hashlib.file_digest(path.open("rb"), "sha256").hexdigest() == evidence["sha256"], path
paths = {Path(p) for p in subprocess.check_output(["git", "ls-files"], text=True).splitlines()}
paths.update(p for p in Path("results/stage3").rglob("*") if p.is_file() and "__pycache__" not in p.parts)
for item in acceptance["items"].values():
    paths.update(Path(e["path"]) for e in item["evidence"])
keys = {r["cache_key"] for r in json.loads(Path("results/stage3/history/STATE.json").read_text())["rows"]}
keys.update(
    r["cache_key"]
    for r in json.loads(Path("results/stage3/secondary/STATE.json").read_text())["pressure_rows"]
)
for key in keys:
    paths.update(
        p for p in (Path("results/stage2_5/response_table/node_cache") / key).rglob("*") if p.is_file()
    )
paths.update(p for p in Path("third_party/downloads/stage3").rglob("*") if p.is_file())
skip = {"results/stage3/audit/STAGE3_MANIFEST.json", "results/stage3/audit/LOCAL_SYNC_VERIFICATION.json"}
files = {
    str(p): dict(bytes=p.stat().st_size, sha256=hashlib.file_digest(p.open("rb"), "sha256").hexdigest())
    for p in sorted(paths)
    if str(p) not in skip
}
result = dict(
    git_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
    status=acceptance["status"],
    file_count=len(files),
    total_bytes=sum(v["bytes"] for v in files.values()),
    files=files,
)
Path("results/stage3/audit/STAGE3_MANIFEST.json").write_text(json.dumps(result, indent=2) + "\n")
print({k: v for k, v in result.items() if k != "files"})
