"""Verify the new delivery only; never reinterpret historical archive omissions."""

from pathlib import Path
import ast
import json
import hashlib
import numpy as np
import h5py
from scripts.stage4l.dense_probe import sha, write
from scipy.spatial.distance import pdist
from mevo_mb.coupling.geometry import geometry_from_config

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "results/stage5/motion_solver"


def main():
    checks = []
    for p in BASE.rglob("*.json"):
        json.loads(
            p.read_text(),
            parse_constant=lambda value: (_ for _ in ()).throw(ValueError(f"Nonfinite JSON {p}: {value}")),
        )
        checks.append(dict(check="strict_json", path=str(p.relative_to(ROOT)), status="PASS"))
    for p in BASE.rglob("*.jsonl"):
        for line in p.open():
            json.loads(line, parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))
        checks.append(dict(check="strict_jsonl", path=str(p.relative_to(ROOT)), status="PASS"))
    for p in (BASE / "runs").glob("*/trajectories.h5"):
        with h5py.File(p) as h:
            x = h["positions_m"][:]
            times = h["time_s"][:]
            assert np.isfinite(x).all() and np.isfinite(times).all() and np.all(np.diff(times) > 0)
            last = np.load(p.parent / "LAST_ACCEPTED_STATE.npz")
            assert float(last["time_s"]) == times[-1]
            np.testing.assert_array_equal(last["positions_m"], x[-1])
            summary = json.load(open(p.parent / "SUMMARY.json"))
            assert bool(h.attrs["complete"]) == (summary["status"] == "COMPLETE")
            assert len(times) == summary["accepted_steps"] + 1
            config = json.loads((p.parent / "CONFIG.json").read_text())
            tolerance = config["particles"]["contact_tolerance_m"]
            geometry = geometry_from_config(config)
            radii = h["radii_m"][:]
            gaps, wallgaps = [], []
            for index, frame in enumerate(x):
                active = h["active"][index]
                points, radius = frame[active], radii[active]
                assert np.all(radius == radius[0])
                gap = float(np.min(pdist(points)) - 2 * radius[0])
                wallgap, normal, valid = geometry.wall_gap_and_normal(points, radius)
                assert np.isfinite(wallgap).all() and np.isfinite(normal).all() and np.all(valid)
                assert min(gap, float(wallgap.min())) >= -tolerance
                gaps.append(gap)
                wallgaps.append(float(wallgap.min()))
            ledger = [json.loads(line) for line in (p.parent / "ACCEPTED.jsonl").open()]
            assert len(ledger) == len(times) - 1
            for index, row in enumerate(ledger):
                assert np.isfinite(row["kkt"]) and row["kkt"] <= 1e-9
                assert np.isfinite(row["linear_residual"]) and row["linear_residual"] <= 1e-11
                assert row["projection_calls"] == row["multi_radius_corrections"] == 0
                assert abs(row["minimum_pair_gap_m"] - gaps[index + 1]) < 1e-18
                assert abs(row["minimum_wall_gap_m"] - wallgaps[index + 1]) < 1e-18
            checks.append(
                dict(
                    check="accepted_ledger_and_independent_saved_geometry",
                    path=str(p.relative_to(ROOT)),
                    status="PASS",
                    minimum_pair_gap_m=min(gaps),
                    minimum_wall_gap_m=min(wallgaps),
                    tolerance_m=tolerance,
                )
            )
        checks.append(
            dict(
                check="actual_accepted_trajectory_and_last_state",
                path=str(p.relative_to(ROOT)),
                status="PASS",
                states=len(times),
            )
        )
    for p in [
        ROOT / "src/mevo_mb/bubbles/constrained_motion.py",
        *(ROOT / "scripts/stage5a").glob("*.py"),
        *(ROOT / "tests").glob("test_stage5a*.py"),
    ]:
        ast.parse(p.read_text())
        checks.append(dict(check="syntax", path=str(p.relative_to(ROOT)), status="PASS"))
    old = json.load(open(ROOT / "tests/data/stage5a_frozen_hashes.json"))
    content = (ROOT / "src/mevo_mb/bubbles/translation_backends.py").read_text()
    for node in ast.parse(content).body:
        if isinstance(node, ast.ClassDef):
            assert (
                hashlib.sha256(ast.get_source_segment(content, node).encode()).hexdigest() == old[node.name]
            )
    for name in ["forces.py", "contact.py", "secondary_ablation.py"]:
        assert sha(ROOT / "src/mevo_mb/bubbles" / name) == old[name]
    checks.append(dict(check="legacy_classes_and_physics_bytes", status="PASS"))
    figures = list((BASE / "figures").glob("FIG_M*.png"))
    pdfs = list((BASE / "figures").glob("FIG_M*.pdf"))
    assert len(figures) == len(pdfs) == 8
    from PIL import Image

    for p in figures:
        with Image.open(p) as im:
            im.verify()
    assert all(p.read_bytes().startswith(b"%PDF") for p in pdfs)
    gate = json.load(open(BASE / "final/MOTION_SOLVER_GATE.json"))
    assert gate["status"] == "MOTION_SOLVER_NUMERICAL_GATE_FAILED" and not gate["promotion"]
    assert not any("_5ms_" in p.name for p in (BASE / "runs").iterdir())
    checks.append(dict(check="eight_png_pdf_pairs_and_closed_5ms_gate", status="PASS"))
    sources = {
        str(p.relative_to(ROOT)): sha(p)
        for p in [
            ROOT / "docs/STAGE5A_MOTION_SOLVER_PLAN.md",
            ROOT / "docs/STAGE5A_MOTION_SOLVER_DEVELOPMENT.md",
            ROOT / "src/mevo_mb/bubbles/constrained_motion.py",
            ROOT / "src/mevo_mb/bubbles/osqp_motion.py",
            ROOT / "docs/STAGE5A_OSQP_SUPPLEMENT.md",
            ROOT / "requirements-stage5a-osqp.lock.txt",
            *(ROOT / "scripts/stage5a").glob("*.py"),
        ]
    }
    html = '<!doctype html><html lang="en"><meta charset="utf-8"><title>Stage5A motion solver evidence</title><style>body{font:16px system-ui;max-width:1100px;margin:40px auto;padding:0 20px;background:#f7f9fc;color:#17263c}img{max-width:100%;background:white}section{margin:40px 0;padding:20px;background:white;border-radius:10px}.status{padding:20px;border-left:5px solid #c73939;background:#fff0ef}a{color:#006b9e}</style><h1>Stage5A motion solver evidence</h1><p class="status"><b>FAILED_WITH_EVIDENCE — numerical gate failed.</b><br>OSQP 1.1.3 is installed and the frozen-state tests pass. New 1 ms attempts accepted no motion; 5 ms and longer runs were not launched. Frozen Stage4L remains unchanged. Failed-attempt wall times are not completion speedups.</p><p><a href="../../../docs/STAGE5A_MOTION_SOLVER_DEVELOPMENT.md">Development report</a> · <a href="final/MOTION_SOLVER_GATE.json">Machine-readable gate</a> · <a href="final/ARTIFACT_MANIFEST.json">Artifact hashes</a></p>'
    for p in sorted(figures):
        html += f'<section><h2>{p.stem.replace("_", " ")}</h2><a href="figures/{p.stem}.pdf">Vector PDF</a><br><img src="figures/{p.name}" alt="{p.stem}"></section>'
    html += "</html>"
    (BASE / "OPEN_RESULTS.html").write_text(html)
    manifest = [
        dict(path=str(p.relative_to(ROOT)), size_bytes=p.stat().st_size, sha256=sha(p))
        for p in sorted(BASE.rglob("*"))
        if p.is_file()
        and p.name not in ["DELIVERY_VERIFICATION.json", "ARTIFACT_MANIFEST.json"]
        and "ordered_cache" not in p.parts
        and "__pycache__" not in p.parts
    ]
    write(
        BASE / "final/ARTIFACT_MANIFEST.json",
        dict(
            files=manifest,
            source_sha256=sources,
            large_original_data="read-only by recorded server paths and hashes; not duplicated",
            excluded_rebuildable_cache_directories=["ordered_cache", "__pycache__"],
            self_hash_policy="manifest and verifier report excluded to avoid recursive hashes",
        ),
    )
    write(
        BASE / "final/DELIVERY_VERIFICATION.json",
        dict(
            status="PASS",
            scope="Artifact integrity, not a numerical promotion",
            checks=checks,
            source_sha256=sources,
            scientific_gate=gate["status"],
        ),
    )
    print("DELIVERY_VERIFICATION_PASS", len(checks), flush=True)


if __name__ == "__main__":
    main()
