"""Freeze implicit pair-backend evidence and predeclare deterministic short replays."""
import json
import shutil

import h5py
import numpy as np

from scripts.stage4i5.common import ROOT, BASE, GOLD, sha, write


def main():
    out = BASE / "protocol"
    out.mkdir(exist_ok=False)
    state_path = ROOT / "results/stage2/performance/STATE.json"
    state = json.loads(state_path.read_text())
    historical_path = ROOT / state["run"] / "SUMMARY.json"
    historical = json.loads(historical_path.read_text())
    assert state["status"] == "PASS" and state["rows"] == historical["rows"]
    expected = "37dac3ffde66fab18a9617c44e9e049e00a415ff5e9ce27460a8ed146415b424"
    assert sha(state_path) == expected
    for source, name in [(state_path, "PAIR_BACKEND_STATE.json"), (historical_path, "PAIR_BACKEND_HISTORICAL_SUMMARY.json")]:
        shutil.copy2(source, out / name)
    wins = []
    for n in sorted({r["N"] for r in state["rows"]}):
        by = {r["backend"]: r for r in state["rows"] if r["N"] == n}
        if by["cupy"]["synthetic_step_s"] < by["numpy"]["synthetic_step_s"]:
            wins.append(n)
    crossover = min(wins)
    write(out / "SUPPLEMENTARY_PROVENANCE.json", dict(
        gold_manifest_sha256=sha(GOLD / "GOLD_BASELINE_MANIFEST.json"),
        implicit_input=str(state_path.relative_to(ROOT)), sha256=sha(state_path),
        historical_summary_sha256=sha(historical_path), crossover=crossover,
        local_and_server_sha256_agree=True,
        limitation="Identified after sealing, not independently hashed at original run start. Historical benchmark rows agree. Gold manifest is unchanged. Auto uses GPU at N>=256 unless explicit resource fallback; no original lifetime GPU timing exists."))
    segments = []
    for label in ("B2", "B3"):
        case = GOLD / "cases" / label
        rows = json.loads((case / "metrics.json").read_text())["rows"]
        # Dense selection uses saved geometry/contact evidence only, never candidate performance.
        dense = max(range(len(rows)-1), key=lambda i: (
            (rows[i]["pair_regime"]["contact_fraction"] + rows[i+1]["pair_regime"]["contact_fraction"])/2,
            -i))
        with h5py.File(case / "particles/trajectories.h5") as h:
            for name, index, last in [("S1", 0, 2), ("S2", dense, dense+1)]:
                t0, t1 = float(h["time_s"][index]), float(h["time_s"][last])
                mask = h["active"][index]
                sid = label+"_"+name
                np.savez_compressed(out / (sid+".npz"), time_s=t0, end_time_s=t1,
                    particle_id=h["particle_id"][:], active=mask,
                    positions_m=h["positions_m"][index], radii_m=h["radii_m"][:],
                    **{key: h[key][index] for key in h if key.endswith(("_N", "_m_s", "_Pa", "_m3"))})
                segments.append(dict(id=sid, case=label, segment=name, start_frame=index,
                    end_frame=last, start_time_s=t0, end_time_s=t1, active_count=int(mask.sum()),
                    snapshot_sha256=sha(out / (sid+".npz")),
                    saved_contact_fraction=[rows[i]["pair_regime"]["contact_fraction"] for i in range(index,last+1)],
                    selection="Initial low-contact 0.04 s" if name=="S1" else "Maximum mean saved contact-pair fraction over adjacent frames; earliest tie, approximately 0.02 s"))
    write(out / "PROTOCOL.json", dict(segments=segments, pair_crossover=crossover,
        gold_manifest_sha256=sha(GOLD / "GOLD_BASELINE_MANIFEST.json"),
        baseline="P0 original RK2, original NNLS, original per-particle exact interior certificate",
        benchmark="Sequential same-server workers, BLAS/OMP=1; initialization and diagnostic overhead reported separately; cProfile replay not used as speed denominator",
        gates=dict(P1_max_position_error_m=1e-12, P1_branch_agreement=1., P1_exact_contact_events=True,
                   P1_exact_dt_and_flags=True, contact_missed_pairs=0, wall_penetration_tolerance_m=1e-12,
                   P2_max_position_error_m=1e-9, P2_branch_agreement=1., P2_contact_equivalence=True,
                   P3_max_position_error_m=1.78e-7, P3_cluster_relative_tolerance=.01,
                   P3_validity_absolute_tolerance=.001, P3_branch_agreement=1.),
        note="P1 tolerances are gates, not changes to physical tolerances. No long P0 rerun. Short replay agreement does not prove unsaved full-run identity or local physical accuracy."))
    for path in out.iterdir():
        path.chmod(0o444)
    out.chmod(0o555)
    print(json.dumps(dict(protocol=str(out.relative_to(ROOT)), segments=segments)), flush=True)


if __name__ == "__main__":
    main()
