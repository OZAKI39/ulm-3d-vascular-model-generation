"""Exact signed ledger, paired timestep sensitivity and separate credibility flags."""

import argparse
import json
import h5py
import numpy as np
from .common import BASE, ROOT, js, sha, write
from .validation import validate_ledger

TERMS = [
    "flow",
    "primary",
    "secondary",
    "mobility",
    "pair_contact",
    "wall",
    "boundary",
    "accepted",
    "residual",
]


def load_case(path):
    summary = js(path / "SUMMARY.json")
    cfg = js(path / "CONFIG.json")
    gate = js(path / "GATE_CROSSING.json")
    with h5py.File(path / "trajectory.h5") as h:
        d = dict(
            t=h["time_s"][:],
            x=h["positions_m"][:, 0],
            c=h["metrics/global_centroid_m"][:, 0],
            rg=h["metrics/global_rg_m"][:, 0],
            fraction=h["metrics/largest_size"][:, 0] / 32,
            error=h["trajectory/target_error_m"][:],
            wall_norm=h["wall_stats"][:, 0, :, 3].sum(1) / 32,
            controls=json.loads(h["control"].attrs["events_json"]),
        )
        assert np.all(h["active"][:] == 1)
        np.testing.assert_allclose(d["c"], d["x"].mean(1), atol=1e-17, rtol=0)
        np.testing.assert_allclose(
            d["error"],
            np.linalg.norm(d["c"][:, 1:] - np.asarray(cfg["environment"]["target_m"])[1:], axis=1),
            atol=1e-17,
            rtol=0,
        )
        assert d["controls"] == summary["controls"]
        assert np.isclose(np.sum(h["trajectory/physical_dt_s"][:]), 0.3, atol=1e-12, rtol=0)
        if "motion_ledger" in h:
            ledger = {
                k: h["motion_ledger/" + k][:]
                for k in ["displacement_m", "step_meta", "wall_contact_events", "per_ID_cumulative_m"]
            }
            ledger["rejected_intervals"] = []
            d["ledger_check"] = validate_ledger(ledger, h["positions_m"][:], d["t"])
            d["ledger"] = ledger
    assert summary["status"] == "COMPLETE" and summary["numerical_quality"] == "PASS"
    lo = np.asarray(gate["accepted_before_m"])
    hi = np.asarray(gate["accepted_after_m"])
    cross = lo + gate["fraction"] * (hi - lo)
    np.testing.assert_array_equal(cross, np.asarray(gate["derived_positions_m"]))
    assert abs(cross.mean(0)[0] - cfg["environment"]["x_gate_m"]) < 1e-14
    d.update(summary=summary, cfg=cfg, gate=gate, path=path)
    return d


def ledger_report(d):
    ledger = d["ledger"]
    dx = ledger["displacement_m"]
    meta = ledger["step_meta"]
    ev = ledger["wall_contact_events"]
    ends = meta[:, 0] + meta[:, 1]
    cum = np.vstack([np.zeros((1, 9, 3)), np.cumsum(dx, axis=0)])
    take = np.searchsorted(np.r_[0, ends], d["t"] - 1e-12)
    d["ledger_at_samples"] = cum[take]
    group = d["gate"]
    k = np.argmin(np.abs(meta[:, 0] - group["bracket_s"][0]))
    assert abs(meta[k, 0] - group["bracket_s"][0]) < 1e-12
    gate_components = cum[k] + group["fraction"] * dx[k]
    np.testing.assert_allclose(
        gate_components[7], np.asarray(group["derived_centroid_m"]) - d["c"][0], atol=5e-17, rtol=0
    )
    # Fractional attribution at the gate is a bookkeeping interpolation of whole-step terms.
    # It is not a physical chronology of the within-step projection operations.
    d["ledger_gate_components"] = gate_components
    by_id = np.zeros((32, 3))
    by_id_norm = np.zeros(32)
    if len(ev):
        ids = ev[:, 3].astype(int)
        np.add.at(by_id, ids, ev[:, 11:14])
        np.add.at(by_id_norm, ids, ev[:, 14])
    np.testing.assert_allclose(by_id, ledger["per_ID_cumulative_m"][:, 5], atol=2e-18, rtol=1e-9)
    faces = {}
    for f in np.unique(ev[:, 19]).astype(int):
        e = ev[ev[:, 19] == f]
        faces[str(f)] = dict(
            count=len(e),
            signed_sum_m=e[:, 11:14].sum(0),
            length_sum_m=e[:, 14].sum(),
            first_step_start_s=e[0, 1],
            last_step_start_s=e[-1, 1],
            IDs=np.unique(e[:, 3]).astype(int),
        )
    significant = ev[(ev[:, 14] > 1e-15) & (ev[:, 4] > 0)]
    cancellation = np.sum(significant[:, 11:14] * significant[:, 20:23], axis=1) < 0
    return dict(
        check=d["ledger_check"],
        final_mean_displacement_by_term_m={t: dx[:, i].sum(0) for i, t in enumerate(TERMS)},
        gate_interpolated_mean_by_term_m={t: gate_components[i] for i, t in enumerate(TERMS)},
        gate_ledger_scope="Fractional whole-step bookkeeping to the accepted-step linear gate intersection; not substep physical event timing.",
        mean_pair_path_m=meta[:, 4].sum(),
        mean_wall_path_m=meta[:, 5].sum(),
        per_ID_signed_wall_m=by_id,
        per_ID_wall_path_m=by_id_norm,
        top_Z_IDs=np.argsort(-abs(by_id[:, 2]))[:5],
        first_wall_step_bracket_s=None if not len(ev) else [ev[0, 1], ev[0, 1] + ev[0, 2]],
        wall_face_indices=faces,
        round_counts={str(i): int(np.sum(ev[:, 4] == i)) for i in [0, 3, 6, 9, 12]},
        opposite_to_preceding_PBD_fraction=None if not len(significant) else float(cancellation.mean()),
        alternation_filter="Only nonzero wall vectors >1e-15m after PBD rounds3/6/9/12; raw events all retained.",
        max_boundary_step_m=float(np.linalg.norm(dx[:, 6], axis=1).max(initial=0)),
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preliminary", action="store_true")
    a = ap.parse_args()
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    paths = {f"{c}_10": ROOT / "results/fast_pzt/f8/y_formal_v1" / c for c in ["Y0", "Y1", "Y2"]}
    paths.update({f"{c}_5": ROOT / "results/fast_pzt/f8/y_half_v1" / c for c in ["Y1", "Y2"]})
    for c in ["Y0_5", "Y1_10", "Y1_5", "Y2_10", "Y2_5"]:
        p = BASE / (c.lower() + "_v1") / "run"
        if (p.parent / "COMPLETE.json").exists():
            paths[c] = p
    d = {k: load_case(p) for k, p in paths.items()}
    proto = js(ROOT / "configs/fast_pzt/f9/protocol.json")
    thresholds = proto["diagnostics"]
    rows = {}
    plot = {}
    provenance = {}
    for name, v in d.items():
        s = v["summary"]
        rows[name] = dict(
            path=str(v["path"].relative_to(ROOT)),
            numerical_quality=s["numerical_quality"],
            dt_s=v["cfg"]["dt_s"],
            gate=v["gate"],
            final_centroid_m=v["c"][-1],
            displacement_m=v["c"][-1] - v["c"][0],
            final_fraction=v["fraction"][-1],
            integration_wall_s=s["wall_s"],
            gpu_stream_elapsed_s=s["gpu_stream_elapsed_s"],
            process_body_wall_s=s["process_body_wall_s"],
            output_and_postprocessing_s=s.get("output_and_postprocessing_s"),
            audit=s.get("audit", dict(enabled=False)),
            quality=s["quality"],
        )
        if "ledger" in v:
            rows[name]["ledger"] = ledger_report(v)
            rows[name]["WALL_CORRECTION_SIGNIFICANT"] = bool(
                abs(v["ledger"]["displacement_m"][:, 5, 2].sum())
                > thresholds["WALL_CORRECTION_SIGNIFICANT"]["threshold_m"]
            )
            rows[name]["CONTACT_SENSITIVE_direct_group_component"] = bool(
                abs(v["ledger"]["displacement_m"][:, 4, 2].sum())
                > thresholds["CONTACT_SENSITIVE"]["threshold_m"]
            )
            rows[name]["contact_indirect_arrangement_effect"] = (
                "NOT_ISOLATED; small signed group sum does not imply PBD is physically irrelevant"
            )
            plot[name + "_ledger_cumulative"] = v["ledger_at_samples"]
        for k in ["t", "c", "rg", "fraction", "error", "wall_norm"]:
            plot[name + "_" + k] = v[k]
        provenance[name] = dict(
            execution_commit=js(v["path"] / "IDENTITY.json")["execution_commit"],
            source_files={str(p.relative_to(ROOT)): sha(p) for p in v["path"].iterdir() if p.is_file()},
            sample_role="audit_reproduction"
            if name != "Y0_5" and "ledger" in v
            else "new_half_step_numerical_control"
            if name == "Y0_5"
            else "inherited_F8",
        )
    paired = {}
    for condition in ["Y0", "Y1", "Y2"]:
        if condition + "_5" not in d:
            continue
        coarse, fine = d[condition + "_10"], d[condition + "_5"]
        np.testing.assert_allclose(coarse["t"], fine["t"], atol=1e-14, rtol=0)
        delta = fine["c"] - coarse["c"]
        transverse = np.linalg.norm(delta[:, 1:], axis=1)
        gp = np.asarray(fine["gate"]["derived_centroid_m"]) - np.asarray(coarse["gate"]["derived_centroid_m"])
        plot[condition + "_step_centroid_delta"] = delta
        row = dict(
            gate_YZ_difference_m=float(np.linalg.norm(gp[1:])),
            gate_signed_XYZ_difference_m=gp,
            gate_time_difference_s=fine["gate"]["time_s"] - coarse["gate"]["time_s"],
            max_same_time_centroid_difference_m=np.linalg.norm(delta, axis=1).max(),
            max_same_time_YZ_difference_m=transverse.max(),
            max_Rg_relative=np.max(abs(fine["rg"] - coarse["rg"]) / coarse["rg"]),
            max_fraction_difference=np.max(abs(fine["fraction"] - coarse["fraction"])),
            max_same_ID_RMS_difference_m=np.sqrt(
                np.mean(np.sum((fine["x"] - coarse["x"]) ** 2, axis=2), axis=1)
            ).max(),
            STEP_SENSITIVE=bool(np.linalg.norm(gp[1:]) > thresholds["STEP_SENSITIVE"]["threshold_m"]),
            first_same_time_YZ_above_s={
                str(nm): float(coarse["t"][hit[0]])
                if len(hit := np.flatnonzero(transverse > nm * 1e-9))
                else None
                for nm in [1, 10, 50, 150]
            },
        )
        if "ledger" in coarse and "ledger" in fine:
            diff = fine["ledger_at_samples"] - coarse["ledger_at_samples"]
            row["final_signed_component_difference_m"] = {t: diff[-1, i] for i, t in enumerate(TERMS)}
            row["wall_face_sets_equal"] = set(rows[condition + "_10"]["ledger"]["wall_face_indices"]) == set(
                rows[condition + "_5"]["ledger"]["wall_face_indices"]
            )
            plot[condition + "_step_ledger_delta"] = diff
            # Nonoverlapping 10ms windows, reporting actual differences, not confidence intervals.
            edges = np.arange(0, len(coarse["t"]), 20)
            if edges[-1] != len(coarse["t"]) - 1:
                edges = np.r_[edges, len(coarse["t"]) - 1]
            changes = np.diff(diff[edges, 5, 2])
            idx = np.argsort(-abs(changes))[:6]
            row["largest_10ms_wall_difference_changes"] = [
                dict(start_s=coarse["t"][edges[j]], end_s=coarse["t"][edges[j + 1]], change_m=changes[j])
                for j in idx
            ]
        paired[condition] = row
    effects = {}
    for comparator in ["Y0", "Y1"]:
        if comparator not in paired:
            continue
        c = d[comparator + "_10"]
        f = d["Y2_10"]
        increment = f["gate"]["derived_centroid_m"][2] - c["gate"]["derived_centroid_m"][2]
        bound = paired["Y2"]["gate_YZ_difference_m"] + paired[comparator]["gate_YZ_difference_m"]
        limit = max(5e-8, 3 * bound)
        effects[comparator] = dict(
            gate_target_direction_increment_m=increment,
            numerical_pair_bound_m=bound,
            required_increment_m=limit,
            TRANSVERSE_EFFECT_RESOLVED=bool(increment > limit),
            common_time_final_delta_m=f["c"][-1] - c["c"][-1],
            gate_error_improvement_m=c["gate"]["transverse_error_m"] - f["gate"]["transverse_error_m"],
        )
    report = dict(
        cases=rows,
        timestep=paired,
        effects=effects,
        provenance=provenance,
        definitions=thresholds,
        fixed_target_registration_sha256=sha(
            ROOT / "results/fast_pzt/f8/y_static_v1/TARGETS_BEFORE_TRAJECTORIES.json"
        ),
        limitations=[
            "Within-trajectory bookkeeping, not a no-wall counterfactual",
            "Two-step differences are deterministic sensitivities, not statistical confidence intervals",
            "Wall geometric constraint is not calibrated near-wall hydrodynamics",
            "Unrecorded ledger components remain missing, never reconstructed from coarse saved frames",
        ],
    )
    write(out / ("PRELIMINARY.json" if a.preliminary else "RESULTS.json"), report)
    np.savez_compressed(out / ("preliminary_plot_data.npz" if a.preliminary else "plot_data.npz"), **plot)
    if not a.preliminary:
        assert set(effects) == {"Y0", "Y1"}
        assert all("ledger" in rows[k] for k in ["Y0_5", "Y1_10", "Y1_5", "Y2_10"])
        assert not all(e["TRANSVERSE_EFFECT_RESOLVED"] for e in effects.values())
        write(
            BASE / "FINAL_GATE.json",
            dict(
                software="PASS",
                ledger_closure="PASS",
                audit_noninterference="PASS",
                Y0_half_step="COMPLETE",
                historical_F8_labels="UNCHANGED",
                numerical="NUMERICAL_PASS",
                scientific="WALL_SENSITIVE / NAVIGATION_NOT_RESOLVED",
                navigation="COHESION_ONLY",
                physics="PHYSICS_UNCALIBRATED",
                evidence_level="EXPLORATORY",
                audit_reproductions_not_independent_samples=True,
                training_ready=False,
                figures="PENDING",
                archive="PENDING",
                progress_document="PENDING",
            ),
        )
    print(
        json.dumps(
            {
                k: dict(
                    step_sensitive=v["STEP_SENSITIVE"], gate_YZ_difference_um=v["gate_YZ_difference_m"] * 1e6
                )
                for k, v in paired.items()
            },
            indent=2,
        )
    )
    for k, v in rows.items():
        if "ledger" in v:
            print(
                k,
                "mean signed displacement [um]:",
                {
                    t: np.round(np.asarray(z) * 1e6, 9).tolist()
                    for t, z in v["ledger"]["final_mean_displacement_by_term_m"].items()
                },
            )


if __name__ == "__main__":
    main()
