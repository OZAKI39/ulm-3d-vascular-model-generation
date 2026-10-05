"""Stage 4 evidence, prerequisites and immutable per-run identities."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess
import zipfile

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def require(kind):
    state = json.loads((ROOT / "results/stage4" / kind / "STATE.json").read_text())
    if state["status"] not in ("PASS", "PASS_WITH_SCOPE"):
        raise RuntimeError(f"Stage 4 prerequisite {kind}: {state['status']}")
    return state


def new_run(kind, config, prerequisites=()):
    for prerequisite in prerequisites:
        require(prerequisite)
    out = ROOT / "results/stage4" / kind / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out.mkdir(parents=True)
    sources = {str(p.relative_to(ROOT)): sha(p)
               for folder in ("src", "scripts/stage4", "scripts/stage4i5", "configs/stage4", "validation/history_branch_robustness")
               for p in (ROOT / folder).rglob("*") if p.is_file() and "__pycache__" not in p.parts}
    write(out / "IDENTITY.json", dict(
        created_UTC=datetime.now(timezone.utc).isoformat(), config=config,
        git_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        baseline_tag="stage3-minimal-validated", source_sha256=sources,
        prerequisites={k: require(k) for k in prerequisites},
    ))
    with zipfile.ZipFile(out / "source_snapshot.zip", "x", zipfile.ZIP_DEFLATED) as archive:
        for name in sources:
            archive.write(ROOT / name, name)
    return out


def finish(out, summary):
    result = dict(summary, run=str(Path(out).relative_to(ROOT)))
    write(out / "SUMMARY.json", result)
    write(out.parent / "STATE.json", result)
    print(json.dumps(result, allow_nan=False), flush=True)
    return result
