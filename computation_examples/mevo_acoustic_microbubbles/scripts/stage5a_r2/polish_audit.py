"""Pre/post original saddle polish from immutable failed OSQP candidates."""

import json
import importlib.metadata as md
from pathlib import Path
import numpy as np
from scipy import sparse
from scipy.sparse.linalg import splu
from scripts.stage5a_r2.seal import BASE, ROOT
from scripts.stage4l.dense_probe import sha, write
from mevo_mb.bubbles.motion_qp_r2 import certificate


def main():
    out = BASE / "diagnosis/polish_parts"
    out.mkdir(exist_ok=True)
    replay = json.loads((BASE / "diagnosis/reconstruction_v1/REPLAY.json").read_text())
    rows = []
    for p in sorted((BASE / "fixtures").glob("*_osqp")):
        a = sparse.load_npz(p / "A.npz")
        g = sparse.load_npz(p / "G.npz")
        d = np.load(p / "QP_INPUT.npz")
        b = d["b"]
        lower = d["lower"]
        u = d["primal_candidate"]
        lam = d["multipliers_candidate"]
        context = next(
            r["reconstruction"] for r in replay if r["fixture"] == p.name and "reconstruction" in r
        )
        options = dict(velocity_scale=context["velocity_scale_m_s"], force_scale=context["force_scale_N"])
        before = certificate(a, b, g, lower, u, lam, **options)
        ix = np.flatnonzero(lam > 1e-10 * max(1.0, np.max(abs(lam), initial=0)))
        ca = g[ix]
        saddle = sparse.bmat([[a, -ca.T], [-ca, None]], format="csc")
        rhs = np.r_[b, -lower[ix]]
        answer = splu(saddle).solve(rhs)
        v = answer[: a.shape[0]]
        mult = np.zeros(g.shape[0])
        mult[ix] = answer[a.shape[0] :]
        after = certificate(a, b, g, lower, v, mult, **options)
        target = out / (p.name + ".npz")
        np.savez_compressed(
            target, primal_before=u, dual_before=lam, primal_after=v, dual_after=mult, selected_rows=ix
        )
        rows.append(
            dict(
                fixture=p.name,
                before=before,
                after=after,
                selected_rows=ix.tolist(),
                saddle_relative_residual=float(
                    np.linalg.norm(saddle @ answer - rhs) / max(1, np.linalg.norm(rhs))
                ),
                matrix_sha256={n: sha(p / n) for n in ["A.npz", "G.npz", "QP_INPUT.npz"]},
                candidates_sha256=sha(target),
            )
        )
    write(
        out / "SUMMARY.json",
        dict(
            rows=rows,
            scipy=md.version("scipy"),
            numpy=md.version("numpy"),
            original_polish_code_sha256=sha(ROOT / "src/mevo_mb/bubbles/osqp_motion.py"),
            script_sha256=sha(Path(__file__)),
            execution_scope="Supplementary local float64 SciPy replay of the original one-shot polish equations on the immutable saved candidates; raw earlier server replays stay unchanged. Near-singular saddle pivoting may vary by library.",
            original_problem_unchanged=True,
        ),
    )
    print([(r["fixture"], r["before"]["kkt"], r["after"]["kkt"]) for r in rows])


if __name__ == "__main__":
    main()
