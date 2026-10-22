"""Actual saved trajectory plane crossings; no velocity extrapolation as evidence."""

import json
from time import perf_counter
import h5py
import numpy as np
from .common import ROOT, BASE, identity, js, sha, write
from .science import OriginalDomain


def crossing(t, c, xplane):
    ids = np.flatnonzero((c[:-1, 0] < xplane) & (c[1:, 0] >= xplane))
    if not len(ids):
        return dict(status="NOT_OBSERVED")
    i = int(ids[0])
    w = float((xplane - c[i, 0]) / (c[i + 1, 0] - c[i, 0]))
    return dict(
        status="OBSERVED_SAVED_BRACKET",
        left_index=i,
        fraction=w,
        time_bracket_s=t[i : i + 2],
        time_s=float(t[i] + w * (t[i + 1] - t[i])),
        centroid_m=c[i] + w * (c[i + 1] - c[i]),
        interpolation="linear within real saved50us bracket, no smoothing",
    )


def main():
    start = perf_counter()
    domain = OriginalDomain()
    rows = []
    arrays = {}
    for step in ["10", "5"]:
        for mode in ["OFF", "S_1000", "O_1000"]:
            p = ROOT / f"results/fast_pzt/f10a/extended_{step}_v1/{mode}/trajectory.h5"
            with h5py.File(p) as f:
                t = f["time_s"][:]
                x = f["positions_m"][:, 0]
                c = x.mean(1)
                c0 = c[0]
                env = json.loads(f.attrs["config_json"])["environment"]
                control = json.loads(f.attrs["controls_json"])
                active = f["active"][:]
                valid = f["flow/valid"][:]
                flow = f["flow/velocity_m_s"][:, 0].mean(1)
                pressure = f["control/pressure_Pa"][:]
                gap = f["wall_gap_m"][:]
                quality = f.attrs["numerical_quality"]
            xplane = float(c0[0] + 2e-5)
            cross = crossing(t, c, xplane)
            radius = 1.78e-6
            field = np.all((x >= domain.lower + radius) & (x <= domain.upper - radius), axis=(1, 2))
            initial = domain.sample(x[0])
            target = np.array(env["target_m"])
            target_support = domain.sample(target)
            target_proposal = c0 + np.array([2e-5, 0, 1e-6])
            proposal_support = domain.sample(target_proposal)
            row = dict(
                mode=mode,
                dt_s=int(step) / 1e6,
                source=str(p.relative_to(ROOT)),
                sha256=sha(p),
                historical_target_m=target,
                initial_centroid_m=c0,
                initial_positions_m=x[0],
                historical_task=env["task_type"],
                historical_target_is_downstream=False,
                plane_x_m=xplane,
                plane_definition="retrospective20um diagnostic, not historical navigation objective",
                crossing=cross,
                quality=str(quality),
                all_active=bool(active.all()),
                CFD_saved_valid=bool(valid.all()),
                field_saved_valid=bool(field.all()),
                minimum_saved_wall_gap_m=float(gap.min()),
                observed_field_support_duration_s=float(t[-1])
                if field.all()
                else float(t[np.flatnonzero(~field)[0]]),
                field_exit_time="RIGHT_CENSORED at10ms, no integrated exit or residence time known",
                initial_CFD_m_s=initial["CFD_velocity_m_s"].mean(0),
                historical_target_support=target_support,
                proposed_target_m=target_proposal,
                proposed_target_support=proposal_support,
                control=control,
                final_displacement_m=c[-1] - c0,
                pressure_range_Pa=[float(np.min(abs(pressure))), float(np.max(abs(pressure)))],
            )
            rows.append(row)
            arrays[f"{step}_{mode}_t"] = t
            arrays[f"{step}_{mode}_c"] = c
            arrays[f"{step}_{mode}_flow"] = flow
    pairs = []
    for step in ["10", "5"]:
        by = {r["mode"]: r for r in rows if r["dt_s"] == int(step) / 1e6}
        off, original, opt = [by[m] for m in ["OFF", "S_1000", "O_1000"]]
        co = np.array(opt["crossing"]["centroid_m"])
        cf = np.array(off["crossing"]["centroid_m"])
        cs = np.array(original["crossing"]["centroid_m"])
        t = arrays[f"{step}_OFF_t"]
        d = arrays[f"{step}_O_1000_c"] - arrays[f"{step}_OFF_c"]
        limit = opt["crossing"]["left_index"] + 2
        pair = dict(
            dt_s=int(step) / 1e6,
            O_minus_OFF_plane_Z_m=float(co[2] - cf[2]),
            O_minus_S_plane_Z_m=float(co[2] - cs[2]),
            O_minus_OFF_10ms_Z_m=float(d[-1, 2]),
            O_total_plane_Z_m=float(co[2] - np.array(opt["initial_centroid_m"])[2]),
            O_plane_target_YZ_error_m=float(np.linalg.norm((co - np.array(opt["proposed_target_m"]))[1:])),
            O_positive_increment_at_all_saved_times_after_initial=bool(np.all(d[1:, 2] > 0)),
            minimum_early_paired_Z_rate_m_s=float(np.min(np.diff(d[:limit, 2]) / np.diff(t[:limit]))),
            full_10ms_paired_Z_rate_range_m_s=[
                float(np.min(np.diff(d[:, 2]) / np.diff(t))),
                float(np.max(np.diff(d[:, 2]) / np.diff(t))),
            ],
        )
        pairs.append(pair)
    delta = abs(pairs[0]["O_minus_OFF_plane_Z_m"] - pairs[1]["O_minus_OFF_plane_Z_m"])
    result = dict(
        rows=rows,
        paired=pairs,
        dt_pair_plane_increment_difference_m=delta,
        original50nm_effect_gate_at_plane=all(r["O_minus_OFF_plane_Z_m"] >= 50e-9 for r in pairs),
        paired_signal_above_3x_dt_difference=all(r["O_minus_OFF_plane_Z_m"] > 3 * delta for r in pairs),
        numerical_scope="Centroid only. Two timesteps and saved-bracket interpolation are not a rigorous error bound; internal morphology remains sensitive.",
        control_time_scope="To20um plane: measured crossing. In stored field: observed for10ms, exit censored.",
        navigation_established=False,
        training_ready=False,
        identity=identity(),
        wall_s=perf_counter() - start,
    )
    f8p = ROOT / "results/fast_pzt/f8/y_static_v1/TARGETS_BEFORE_TRAJECTORIES.json"
    f8 = js(f8p)
    gates = {}
    for mode in ["Y0", "Y1", "Y2"]:
        gp = ROOT / f"results/fast_pzt/f8/y_formal_v1/{mode}/GATE_CROSSING.json"
        gate = js(gp)
        if gate["status"] == "CROSSED":
            before = np.array(gate["accepted_before_m"])
            after = np.array(gate["accepted_after_m"])
            cross = before + gate["fraction"] * (after - before)
            np.testing.assert_allclose(cross, gate["derived_positions_m"], rtol=0, atol=1e-18)
            assert abs(cross.mean(0)[0] - f8["environment"]["x_gate_m"]) < 1e-14
        gates[mode] = dict(source=str(gp.relative_to(ROOT)), sha256=sha(gp), gate=gate)
    result["F8_original_preregistered_near_wall_gate"] = dict(
        source=str(f8p.relative_to(ROOT)),
        sha256=sha(f8p),
        target_m=f8["target_m"],
        initial_centroid_m=f8["environment"]["initial_centroid_m"],
        design=f8["design"],
        actual_crossing_evidence=gates,
        scope="Historical near-wall different location/control budget. Not the1kPa interior target and not calibrated B2 evidence.",
    )
    write(BASE / "analysis/NAVIGATION.json", result)
    np.savez_compressed(BASE / "analysis/NAVIGATION_SERIES.npz", **arrays)
    selected = next(r for r in rows if r["mode"] == "O_1000" and r["dt_s"] == 1e-5)
    proposal = dict(
        status="PREREGISTERED_PROPOSAL_NOT_EXECUTED",
        source=selected["source"],
        source_sha256=selected["sha256"],
        frame=0,
        positions_m=selected["initial_positions_m"],
        target_m=selected["proposed_target_m"],
        plane_x_m=selected["plane_x_m"],
        target_radius_YZ_m=0.2e-6,
        task="First downstream-plane crossing atX0+20um; fixed targetY0,Z0+1um. Report target error and signed paired effect even on failure. Do not move target after observing result.",
        maximum_duration_s=0.003,
        stop_at="First centroid crossing, invalid state, or3ms. Outside-lumen/disconnected/lost-ID is failure.",
        controls="OFF/S/O existing F10A events, fixed for window; no new optimizer in this proposed paired baseline",
        control_events={r["mode"]: r["control"] for r in rows if r["dt_s"] == 1e-5},
        initial_median_pressure_Pa=1000,
        per_source_cap_Pa=20000,
        command_squared_budget_Pa2=float(
            np.sum(
                np.array(next(r for r in rows if r["mode"] == "S_1000")["control"][0]["amplitude_Pa"]) ** 2
            )
        ),
        phases_rad_bounds=[-float(np.pi), float(np.pi)],
        amplitude_budget_note="Exact squared sum of archived S commands. Not electric/acoustic power.",
        constraints="Same initial median1kPa; peak local pressure<=1.01*S; command squared sum<=S with1e-8 relative numerical tolerance.",
        dt_s=[1e-5, 5e-6],
        quality=dict(
            overlap_m=17.8e-9,
            single_correction_m=178e-9,
            centroid_dt_error_m=10e-9,
            minimum_resolved_acoustic_increment_m=50e-9,
            required_target_error_m=0.2e-6,
            largest_fraction=0.9,
            all_IDs_active=True,
        ),
        expectation="Current fixed1kPaO predicted insufficient by retrospective actual trajectories. Baseline falsification/authority mapping, not expected navigation success.",
        scope="Interior mother branch; avoid near-wall extension. Any longer-window/newtarget/controller proposal requires separate preregistration before running.",
        training_ready=False,
    )
    write(BASE / "analysis/NAVIGATION_PREREGISTRATION.json", proposal)
    print([(r["mode"], r["dt_s"], r["crossing"]["time_s"]) for r in rows])
    print(pairs)


if __name__ == "__main__":
    main()
