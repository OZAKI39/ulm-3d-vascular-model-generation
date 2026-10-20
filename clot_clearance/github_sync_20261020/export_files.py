"""Create a curated export without changing the source archive (stdlib only)."""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import zipfile

LIMIT = 32 * 1024**2
META = "github_sync_20261020"
STATE_KEYS = {
    "x", "d", "atom_ids", "time", "step", "velocity", "scalars",
    "counters", "native_bonds", "Uc", "K", "native_U",
    "external_force", "external_drive_force", "external_drag_force",
    "external_wall_force", "wall_potential_J", "field_amplitude_Pa",
}
TEXT = {".csv", ".json", ".jsonl", ".log", ".txt", ".native", ".data", ".svg"}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for chunk in iter(lambda: f.read(2**20), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_sha(path):
    h = hashlib.sha256()
    with gzip.open(path, "rb") as f:
        for chunk in iter(lambda: f.read(2**20), b""):
            h.update(chunk)
    return h.hexdigest()


def gzip_write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as f:
        f.write(data)


def excluded(rel):
    p = Path(rel)
    if any(part in {".git", "__pycache__", ".pytest_cache", ".venv", "venv", "CMakeFiles"} for part in p.parts):
        return "generated cache/build metadata"
    if p.suffix in {".pyc", ".pyo", ".o", ".so", ".a", ".restart", ".bin", ".db", ".pkl", ".rsp"} or ".so." in p.name:
        return "rebuildable native binary, restart or redundant dense stream"
    if p.name in {"GIT_STAGE_PATHS.txt", ".DS_Store"} or ".rendering." in p.name:
        return "temporary staging/render artifact"
    if any(part in {"render_smoke", "camera_probes", "camera_probe"} for part in p.parts):
        return "superseded camera/render trial"
    if p.name.startswith(".env") or p.suffix in {".pem", ".key"}:
        return "environment/credential file"
    return None


def export(root, dest, plan, label):
    records = []
    for index, rel in enumerate(plan):
        src = root / rel
        if not src.is_file() or src.is_symlink():
            records.append({"source": rel, "action": "exclude", "reason": "not a regular file"})
            continue
        stat = src.stat()
        row = {"source": rel, "source_bytes": stat.st_size, "source_sha256": sha(src)}
        reason = excluded(rel)
        out = dest / rel
        if reason:
            row.update(action="exclude", reason=reason)
        elif out.is_file() and sha(out) == row["source_sha256"]:
            row.update(action="already_present", export=rel, export_sha256=row["source_sha256"])
        elif out.with_name(out.name + ".gz").is_file() and gzip_sha(out.with_name(out.name + ".gz")) == row["source_sha256"]:
            row.update(action="represented_by_existing_gzip", export=rel + ".gz",
                       export_sha256=sha(out.with_name(out.name + ".gz")))
        elif src.suffix == ".npz" and stat.st_size > 512 * 1024:
            with zipfile.ZipFile(src) as z:
                names = set(z.namelist())
                is_state = {"x.npy", "d.npy", "time.npy"}.issubset(names)
                latest_key = ("y_clot_irregular_current_method_v4/" in rel
                              and not src.name.startswith("frame_"))
                if is_state and not latest_key:
                    out = dest / META / "compact_states" / rel
                    out.parent.mkdir(parents=True, exist_ok=True)
                    fields = {}
                    with zipfile.ZipFile(out, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as w:
                        for name in sorted(names):
                            if name.removesuffix(".npy") in STATE_KEYS:
                                data = z.read(name)
                                zi = zipfile.ZipInfo(name, date_time=(1980, 1, 1, 0, 0, 0))
                                zi.compress_type = zipfile.ZIP_DEFLATED
                                w.writestr(zi, data)
                                fields[name] = hashlib.sha256(data).hexdigest()
                    with zipfile.ZipFile(out) as w:
                        assert all(hashlib.sha256(w.read(k)).hexdigest() == v for k, v in fields.items())
                    row.update(action="compact_state", fields_sha256=fields,
                               omitted_fields=sorted(names - set(fields)),
                               reason="Selected arrays preserved byte-for-byte; not a complete native checkpoint")
                else:
                    out.parent.mkdir(parents=True, exist_ok=True)
                    shutil.copy2(src, out)
                    row["action"] = "copy"
        elif src.suffix in TEXT and stat.st_size > 1024**2:
            out = out.with_name(out.name + ".gz")
            out.parent.mkdir(parents=True, exist_ok=True)
            with src.open("rb") as source, out.open("wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as target:
                shutil.copyfileobj(source, target, 2**20)
            assert gzip_sha(out) == row["source_sha256"]
            row["action"] = "gzip"
        elif stat.st_size >= LIMIT:
            row.update(action="exclude", reason="dense/duplicate artifact above 32 MiB; retained at source")
        else:
            out.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, out)
            row["action"] = "copy"
        if row["action"] in {"copy", "gzip", "compact_state"}:
            if out.stat().st_size >= LIMIT:
                row.update(action="exclude", reason="export still exceeds 32 MiB; retained at source")
                out.unlink()
            else:
                row.update(export=str(out.relative_to(dest)), export_bytes=out.stat().st_size,
                           export_sha256=sha(out))
        if src.stat().st_size != stat.st_size or src.stat().st_mtime_ns != stat.st_mtime_ns:
            raise RuntimeError(f"Source changed during export: {rel}")
        records.append(row)
        if index % 250 == 0:
            print(f"{label}: {index}/{len(plan)}", flush=True)
    path = dest / META / f"EXPORT_RECORDS.{label}.json.gz"
    gzip_write(path, json.dumps(records, ensure_ascii=False, indent=2).encode() + b"\n")
    from collections import Counter
    print(json.dumps({"label": label, "actions": dict(Counter(r["action"] for r in records)),
                      "export_MiB": sum(r.get("export_bytes", 0) for r in records) / 2**20}))


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--destination", type=Path, required=True)
    p.add_argument("--plan", type=Path, required=True)
    p.add_argument("--label", required=True)
    a = p.parse_args()
    export(a.source, a.destination, json.loads(a.plan.read_text()), a.label)
