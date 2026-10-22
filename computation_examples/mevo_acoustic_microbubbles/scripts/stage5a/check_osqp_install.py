"""Record isolated installation and fixed-QP replay; never advance a trajectory."""

from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import importlib.metadata as metadata
import json
import time

import numpy as np
import scipy
from scipy import sparse
import osqp

from mevo_mb.bubbles.constrained_motion import kkt_residual

ROOT = Path(__file__).resolve().parents[2]
SETTINGS = dict(
    verbose=False,
    eps_abs=1e-12,
    eps_rel=0.0,
    max_iter=100000,
    polishing=True,
    polish_refine_iter=10,
    adaptive_rho_interval=50,
    check_termination=25,
    warm_starting=True,
)


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def solve(a, b, g, lower, primal=None, dual=None):
    start = time.perf_counter()
    solver = osqp.OSQP(algebra="builtin")
    solver.setup(
        P=sparse.triu(a, format="csc"), q=-b, A=g.tocsc(), l=lower, u=np.full(len(lower), np.inf), **SETTINGS
    )
    # OSQP's lower-bound multiplier sign is opposite to our G U >= lower lambda.
    if primal is not None:
        solver.warm_start(x=primal, y=-dual)
    result = solver.solve(raise_error=False)
    u, lam = result.x, -result.y
    finite = bool(np.isfinite(u).all() and np.isfinite(lam).all())
    kkt, terms = kkt_residual(a, b, g, lower, u, lam) if finite else (None, {})
    row = dict(
        status=result.info.status,
        status_val=int(result.info.status_val),
        finite=finite,
        kkt=kkt,
        kkt_terms=terms,
        iterations=int(result.info.iter),
        setup_s=float(result.info.setup_time),
        solve_s=float(result.info.solve_time),
        polish_s=float(result.info.polish_time),
        wall_s=time.perf_counter() - start,
        polish_status=int(result.info.status_polish),
        pass_kkt=bool(finite and kkt <= 1e-9 and result.info.status_val == 1),
    )
    return u, lam, row


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    overlay = Path(osqp.__file__).resolve().parent.parent
    packages = dict(
        sorted((d.metadata["Name"], d.version) for d in metadata.distributions(path=[str(overlay)]))
    )
    source = ROOT / "tests/data/stage5a_failed_qp"
    with np.load(source / "FAILED_QP_STATE.npz") as h:
        d = {key: h[key] for key in h.files}
    a, g = sparse.load_npz(source / "A.npz"), sparse.load_npz(source / "G.npz")
    identity = dict(
        time_utc=datetime.now(timezone.utc).isoformat(),
        osqp=osqp.__version__,
        numpy=np.__version__,
        scipy=scipy.__version__,
        isolated_overlay=str(overlay),
        installed_packages=packages,
        algebra="builtin CPU",
        settings=SETTINGS,
        user_authorization="Install OSQP explicitly authorized in chat, 2026-10-06",
        input_sha256={
            p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in source.iterdir() if p.is_file()
        },
        scope="Installation and frozen-state QP replay; no scientific trajectory run",
    )
    write(args.output / "IDENTITY.json", identity)
    installed = {
        str(p.relative_to(overlay)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(overlay.rglob("*"))
        if p.is_file() and "__pycache__" not in p.parts
    }
    write(args.output / "INSTALLED_FILE_HASHES.json", installed)
    u, _, smoke = solve(
        sparse.eye(2, format="csc"), np.array([-1.0, 2.0]), sparse.eye(2, format="csc"), np.zeros(2)
    )
    smoke["analytic_max_error"] = float(np.max(np.abs(u - [0.0, 2.0])))
    write(args.output / "SMOKE.json", smoke)
    rows = []
    for name, primal, dual in [("cold", None, None), ("accepted_warm", d["warm_primal"], d["warm_dual"])]:
        u, lam, row = solve(a, d["b"], g, d["lower"], primal, dual)
        row["mode"] = name
        rows.append(row)
        np.savez_compressed(args.output / (name + ".npz"), u=u, multipliers=lam)
        write(args.output / "FAILED_QP_REPLAY.json", rows)
        print(json.dumps(row, allow_nan=False), flush=True)
    write(
        args.output / "SUMMARY.json",
        dict(
            installation="PASS",
            smoke=smoke,
            frozen_failure_replay="PASS" if all(r["pass_kkt"] for r in rows) else "FAILED_WITH_EVIDENCE",
            scientific_promotion=False,
            rows=rows,
        ),
    )


if __name__ == "__main__":
    main()
