"""Re-audit archived candidates in SI units; preserve the original diagnostic."""

import json
from collections import defaultdict
import numpy as np
from scipy import sparse
from scripts.stage5a_r2.seal import BASE
from scripts.stage5a_r2.validation_figures import local
from scripts.stage4l.dense_probe import write, sha
from mevo_mb.bubbles.motion_qp_r2 import certificate


def main():
    folder = BASE / "benchmarks/changing_qp"
    original = folder / "RESULTS.json"
    manifest_path = BASE / "fixtures/development_sequence/MANIFEST.json"
    manifest = json.loads(manifest_path.read_text())
    rows = json.loads(original.read_text())["rows"]
    reports = []
    groups = defaultdict(dict)
    for row in rows:
        item = next(i for i in manifest if i["case"] == row["case"] and i["frame"] == row["frame"])
        p = str(local(item["prefix"]))
        a = sparse.load_npz(p + "_A.npz")
        g = sparse.load_npz(p + "_G.npz")
        d = np.load(p + "_data.npz")
        sol = np.load(local(row["candidate_file"]))
        b = d["b"] * (1 + row["index"] * 1e-5 if row["mode"] == "perturbed" else 1)
        vs = item["velocity_scale_m_s"]
        fs = vs * item["resistance_scale_N_s_m"]
        cert = certificate(
            a, b, g, d["lower"], sol["u"], sol["lambda_"], velocity_scale=vs, force_scale=fs, dt=1e-7
        )
        key = (row["case"], row["mode"], row["index"])
        groups[key][row["solver"]] = (sol["u"].copy(), sol["reaction"].copy(), cert["objective"])
        reports.append(
            dict(
                case=row["case"],
                mode=row["mode"],
                index=row["index"],
                solver=row["solver"],
                native_status=row["diagnostics"].get("status"),
                certificate=cert,
                velocity_scale_m_s=vs,
                force_scale_N=fs,
                candidate_sha256=sha(local(row["candidate_file"])),
                original_wall_s=row["wall_s"],
            )
        )
    agreements = []
    for key, solvers in groups.items():
        u0, r0, o0 = solvers["prox"]
        for name, (u, r, o) in solvers.items():
            agreements.append(
                dict(
                    case=key[0],
                    mode=key[1],
                    index=key[2],
                    solver=name,
                    reference="prox",
                    velocity_relative_max=float(np.max(abs(u - u0)) / max(np.max(abs(u0)), 1e-30)),
                    reaction_relative_L2=float(np.linalg.norm(r - r0) / max(np.linalg.norm(r0), 1e-30)),
                    objective_relative=float(abs(o - o0) / max(abs(o0), 1e-30)),
                )
            )
    write(
        folder / "AUDITED_SI.json",
        dict(
            source_sha256=sha(original),
            manifest_sha256=sha(manifest_path),
            rows=reports,
            cross_solver_comparison=agreements,
            correction="Original diagnostic was generated when the generic candidate checker defaulted force_scale to 1. Its fields labelled force_residual_l2_N/velocity_feasibility_m_s were normalized, NOT SI. This supplementary audit uses the archived resistance and velocity scales and recomputes every original matrix residual. Original RESULTS.json is preserved unchanged.",
            timing_scope="Six simultaneous diagnostic CPU workers; measured wall times describe diagnostic throughput, not isolated formal trajectory latency.",
            multiplier_policy="Compare U, objective and G^T lambda; individual multipliers are not required to be unique.",
        ),
    )
    print("Re-audited", len(reports), "saved candidates; original evidence unchanged.")


if __name__ == "__main__":
    main()
