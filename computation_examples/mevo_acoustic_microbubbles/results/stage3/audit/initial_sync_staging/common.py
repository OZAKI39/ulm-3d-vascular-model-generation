"""Evidence utilities; no numerical physics in this module."""

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import subprocess


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")
    tmp.replace(path)


def provenance(config, papers):
    root = Path.cwd()
    return dict(
        created_UTC=datetime.now(timezone.utc).isoformat(),
        config=config,
        papers=papers,
        git_commit=subprocess.check_output(
            ["git", "-c", f"safe.directory={root}", "rev-parse", "HEAD"], text=True
        ).strip(),
        source_sha256={
            str(p): hashlib.sha256(p.read_bytes()).hexdigest()
            for folder in ("src", "references", "validation/history_force_sensitivity", "scripts/stage3", "configs/stage3")
            for p in Path(folder).rglob("*")
            if p.is_file() and "__pycache__" not in p.parts
        },
        environment="results/stage3/audit/SERVER_DOCTOR.json",
        random_seed=None,
    )


def new_run(kind, config, papers):
    out = Path("results/stage3") / kind / datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    out.mkdir(parents=True)
    atomic_json(out / "CONFIG.json", provenance(config, papers))
    return out


def resumable_run(kind, config, papers):
    """Resume only a matching deterministic experiment; per-case artifacts persist."""
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--resume", type=Path)
    args = parser.parse_args()
    if args.resume is None:
        return new_run(kind, config, papers)
    previous = json.loads((args.resume / "CONFIG.json").read_text())
    if previous["config"] != json.loads(json.dumps(config)):
        raise ValueError("Resume experiment configuration mismatch")
    current = provenance(config, papers)
    for name, digest in previous["source_sha256"].items():
        if name.startswith(("src/", "references/")) and current["source_sha256"].get(name) != digest:
            raise ValueError("Resume source changed: " + name + "; retain old run and start a new one")
    return args.resume


def finish(out, summary):
    summary = dict(summary, run=str(out))
    atomic_json(out / "SUMMARY.json", summary)
    atomic_json(out.parent / "STATE.json", summary)
    print(json.dumps(summary, allow_nan=False), flush=True)


def watermark(fig, text="PARTIALLY VALIDATED"):
    engine = fig.get_layout_engine()
    if engine is not None and type(engine).__name__ == "ConstrainedLayoutEngine":
        engine.set(rect=(0, 0.07, 1, 0.93))
    else:
        fig.subplots_adjust(bottom=0.12)
    fig.text(0.5, 0.012, text, ha="center", va="bottom", fontsize=9, color="darkred", weight="bold")
