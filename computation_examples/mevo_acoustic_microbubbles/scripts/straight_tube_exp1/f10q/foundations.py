"""Small curved-contact references, independent nonlinear KKT root equations."""

import time
import numpy as np
from scipy.optimize import root
from itertools import combinations
from mevo_mb.fast_sim.f10g_contact import resistance
from mevo_mb.fast_sim.sequential_hard_sphere import sequential_step
from scripts.straight_tube_exp1.f10p.foundations import run as old_foundations
from scripts.straight_tube_exp1.f10o.analyze import load
from .common import ROOT, BASE, CONFIG, read, write, identity


def nonlinear_reference(x, F, R, zeta, dt, A):
    """Independent nonlinear KKT system in dimensionless absolute displacement.
    Active subsets seeded by free overlap then checked at true finite sphere geometry.
    No call to production geometry/QP/SQP. <=4 bodies, exhaustive contact subsets.
    """
    pairs = list(combinations(range(len(x)), 2))
    u = x / R
    b = dt * F.ravel() / (zeta * R)
    free = np.linalg.solve(A, b)
    # Enumerate only plausible near pairs, not 6 distant redundant constraints.
    near = [k for k, (i, j) in enumerate(pairs) if np.linalg.norm(u[j] - u[i]) < 2.1]

    def data(s):
        y = u + s.reshape(u.shape)
        rels = np.array([y[j] - y[i] for i, j in pairs])
        dist = np.linalg.norm(rels, axis=1)
        B = np.zeros((len(pairs), u.size))
        for k, (i, j) in enumerate(pairs):
            B[k, 3 * i : 3 * i + 3] = -rels[k] / dist[k]
            B[k, 3 * j : 3 * j + 3] = rels[k] / dist[k]
        return dist - 2, B

    best = None
    for mask in range(1 << len(near)):
        ids = [near[k] for k in range(len(near)) if mask >> k & 1]

        def fun(v):
            s = v[: u.size]
            lam = v[u.size :]
            g, B = data(s)
            return np.r_[A @ s - b - B[ids].T @ lam, g[ids]]

        result = root(fun, np.r_[free, np.zeros(len(ids))], tol=1e-11)
        s = result.x[: u.size]
        lam = result.x[u.size :]
        g, _ = data(s)
        if np.max(abs(fun(result.x))) < 1e-10 and g.min() >= -1e-10 and lam.min(initial=0) >= -1e-10:
            obj = 0.5 * s @ A @ s - b @ s
            if best is None or obj < best[0]:
                best = (obj, x + R * s.reshape(x.shape))
    if best is None:
        raise RuntimeError("No independently verified nonlinear KKT reference")
    return best[1]


def main():
    start = time.perf_counter()
    p = read(CONFIG)
    p.update(read(CONFIG.parent / "algorithm_revision_1.json"))
    R = p["radius_m"]
    zeta = 6 * np.pi * 0.001 * R
    dt = 1e-6
    inherited = old_foundations()
    pair = np.array([[-R, 0, 0], [R, 0, 0]])
    square = np.array([[-R, -R, 0], [R, -R, 0], [-R, R, 0], [R, R, 0]])
    cases = [
        ("tangential_slide", pair, np.array([[0, -1e-11, 0], [0, 1e-11, 0]])),
        ("compress_slide_normal_change", pair, np.array([[1e-9, -1e-11, 0], [-1e-9, 1e-11, 0]])),
        ("release", pair, np.array([[-1e-10, 0, 0], [1e-10, 0, 0]])),
        (
            "four_multiple",
            square,
            -square / R * 1e-10 + np.array([[0, 1e-12, 0], [0, -1e-12, 0], [0, 1e-12, 0], [0, -1e-12, 0]]),
        ),
        (
            "active_switch",
            square,
            np.array([[1e-10, 1e-10, 0], [-1e-10, 1e-10, 0], [-1e-10, -1e-10, 0], [1e-10, -1e-10, 0]]),
        ),
    ]
    records = []
    for name, x, F in cases:
        A = resistance(x, R)["A"]
        r = sequential_step(x, F, R, zeta, dt, A, p)
        ref = nonlinear_reference(x, F, R, zeta, dt, A)
        err = float(np.linalg.norm(ref - r["positions_m"], axis=1).max())
        records.append(
            dict(
                name=name,
                input=dict(x_m=x, F_N=F, A=A),
                result=r,
                reference_positions_m=ref,
                error_m=err,
                status="PASS" if r["converged"] and err < 1e-13 else "FAIL",
            )
        )
    folder = ROOT / "results/straight_tube_exp1/f10p/validation_v1/runs/STRONG_DT1"
    _, x = load(folder)
    x = x[-1]
    saved = read(folder / "LAST_STEP.json")
    F = np.array(saved["fields"]["total_force_N"])
    A = resistance(x, R)["A"]
    r = sequential_step(x, F, R, zeta, dt, A, p)
    write(
        BASE / "local/FOUNDATIONS_v2.json",
        dict(
            identity=identity(),
            inherited_F10P=inherited,
            cases=records,
            frozen_failure=r,
            status="PASS"
            if inherited["status"] == "PASS"
            and all(v["status"] == "PASS" for v in records)
            and r["converged"]
            else "FAIL",
            wall_s=time.perf_counter() - start,
        ),
    )
    print("cases", [(v["name"], v["status"], v["error_m"], len(v["result"]["iterations"])) for v in records])
    print(
        "frozen",
        r["gates"],
        [(v["index"], v["true_force_residual_N"], v["true_complementarity_m"]) for v in r["iterations"]],
    )


if __name__ == "__main__":
    main()
