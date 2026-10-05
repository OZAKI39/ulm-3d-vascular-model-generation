"""Command line entry points; every numerical task gets a new result directory."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from mevo_mb.io.artifacts import new_run, write_json
from mevo_mb.utils.config import load_config


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    doc = sub.add_parser("doctor")
    doc.add_argument("--output", type=Path)
    acoustic = sub.add_parser("acoustic").add_subparsers(dest="action", required=True)
    for action in ["precompute", "validate-superposition", "smoke"]:
        p = acoustic.add_parser(action)
        p.add_argument("--config", required=True, type=Path)
        p.add_argument("--results", default=Path("results"), type=Path)
        p.add_argument("--basis", type=Path)
        p.add_argument("--device", choices=["cpu", "gpu"])
    particle = sub.add_parser("bubbles").add_subparsers(dest="action", required=True).add_parser("run")
    particle.add_argument("--config", required=True, type=Path)
    particle.add_argument("--basis", required=True, type=Path)
    particle.add_argument("--pattern", default="baseline")
    particle.add_argument("--results", default=Path("results"), type=Path)
    visual = sub.add_parser("visualize")
    visual.add_argument("--run", required=True, type=Path)
    demo = sub.add_parser("demo")
    demo.add_argument("--reuse-acoustics", action="store_true")
    demo.add_argument("--results", default=Path("results"), type=Path)
    args = parser.parse_args()
    if args.command == "doctor":
        from mevo_mb.utils.doctor import doctor
        result = doctor()
        if args.output:
            write_json(args.output, result)
        print(json.dumps(result, indent=2))
        return
    if args.command == "visualize":
        from mevo_mb.acoustics.basis import AcousticBasisDatabase
        from mevo_mb.visualization.plots import plot_experiment
        control = json.loads((args.run/"acoustic/CONTROL.json").read_text())
        import numpy as np
        db = AcousticBasisDatabase.load_basis(control["basis_path"])
        field = db.combine(np.array(control["source_control_real_Pa"])+1j*np.array(control["source_control_imag_Pa"]))
        plot_experiment(args.run, field, load_config(args.run/"config.yaml"))
        return
    if args.command == "demo":
        from mevo_mb.coupling.campaign import run_demo
        run_demo(args.results, args.reuse_acoustics)
        return
    c = load_config(args.config)
    if getattr(args, "device", None):
        c["acoustics"]["device"] = args.device
        if args.device == "cpu":
            c["acoustics"]["backend"] = "python"
    run = new_run(args.results, f"{c['name']}_{args.action}", c)
    print(f"RUN_DIRECTORY={run}", flush=True)
    if args.command == "bubbles":
        from mevo_mb.coupling.experiment import run_particles
        run_particles(c, args.basis, run, args.pattern)
        return
    if args.action == "precompute":
        from mevo_mb.acoustics.kwave_backend import precompute
        print(precompute(c, run), flush=True)
    elif args.action == "validate-superposition":
        from mevo_mb.acoustics.validation import validate_superposition
        if args.basis is None:
            parser.error("--basis is required")
        validate_superposition(c, args.basis, run)
    else:
        import numpy as np
        from mevo_mb.acoustics.kwave_backend import KWaveCW
        p, _, meta = KWaveCW(c, run/"acoustic").run(np.array([e["amplitude_Pa"] for e in c["pzt"]]), "smoke")
        if not np.isfinite(p).all() or np.max(np.abs(p)) == 0:
            raise RuntimeError("Acoustic smoke did not produce finite nonzero pressure")
        write_json(run/"metadata.json", {"status": "PASS", **meta})


if __name__ == "__main__":
    main()
