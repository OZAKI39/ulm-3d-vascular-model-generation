"""Verify published bytes and optional original scientific hash contracts."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
META = Path(__file__).resolve().parent
SV = "formal_3D_flow_solver/FEM_SimVascular"
REPORT = ROOT / "particle_3d/reports/network_derived_flow_mb_validation_v1"


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def checked_path(relative):
    path = (ROOT / relative).resolve()
    if not path.is_relative_to(ROOT):
        raise ValueError(f"Path outside snapshot: {relative}")
    return path


def check(path, expected):
    if not path.is_file() or digest(path) != expected:
        raise ValueError(f"Missing file or SHA256 mismatch: {path}")


def scientific_inputs():
    path_map = json.loads((META / "original_path_map.json").read_text())
    prefixes = [
        (f"/home/lzy/projects/ulm_flow_mean_2p0_mmps/{SV}/frozen_reference", f"{SV}/upstream_stage_q_reference"),
        (f"/home/lzy/projects/ulm_flow_mean_2p0_mmps/{SV}", SV),
        ("/home/lzy/projects/ulm_particle_3d_particle0", ""),
        ("/home/lzy/projects/ulm_3D_vascular", "vascular_network"),
        ("/home/lzy/projects/formal_3D_flow_solver", "formal_3D_flow_solver"),
    ]

    def mapped(original):
        if original in path_map:
            return checked_path(path_map[original])
        for old, new in prefixes:
            if original == old or original.startswith(old + "/"):
                return checked_path(new + original[len(old):] if new else original[len(old):].lstrip("/"))
        raise ValueError(f"Unmapped original path: {original}")

    protected = json.loads((ROOT / "vascular_network/reports/a_network_1d0d_boundary_v1/data/protected_input_hashes.json").read_text())
    for original, expected in protected.items():
        check(mapped(original), expected)
    science = json.loads((REPORT / "data/particle_science_protection_manifest.json").read_text())["files"]
    sources = caches = 0
    for original, spec in science.items():
        if Path(original).suffix in {".pyc", ".pyo"}:
            caches += 1
            continue
        path = mapped(original)
        check(path, spec["sha256"])
        if path.stat().st_size != spec["size"]:
            raise ValueError(f"Size mismatch: {path}")
        sources += 1
    snapshot = json.loads((REPORT / "data/code_snapshot_manifest.json").read_text())
    for name, expected in snapshot.items():
        path = checked_path(str((REPORT / "server_bundle" / name).relative_to(ROOT)))
        check(path, expected)
    contract = json.loads((REPORT / "data/old_new_flow_contract.json").read_text())
    for label, original in contract["paths"].items():
        check(mapped(original), contract["sha256"][label])
        check(REPORT / "server_bundle/inputs" / f"{label}.vtu", contract["sha256"][label])
    print(json.dumps({"original_network_inputs": len(protected), "original_science_sources": sources,
                      "server_bundle_sources": len(snapshot), "omitted_generated_bytecode": caches,
                      "old_new_fields": "PASS"}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scientific-inputs", action="store_true")
    args = parser.parse_args()
    count = 0
    for line in (META / "SNAPSHOT_SHA256.txt").read_text().splitlines():
        expected, relative = line.split("  ", 1)
        check(checked_path(relative), expected)
        count += 1
    print(f"PASS: {count} published file hashes (manifest itself excluded).")
    if args.scientific_inputs:
        scientific_inputs()


if __name__ == "__main__":
    main()
