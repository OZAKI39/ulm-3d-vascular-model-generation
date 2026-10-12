"""Offline actual-data navigation gates, timestep comparison and compact figures."""

import argparse
import json
import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mevo_mb.fast_sim.ai_samples import read_group_sample, validate_partitions
from .common import BASE, F5, js, write, sha


def load_case(path):
    summary = js(path / "SUMMARY.json")
    with h5py.File(path / "trajectory.h5") as h:
        data = {
            k: h[k][:]
            for k in [
                "time_s",
                "positions_m",
                "trajectory/target_error_m",
                "metrics/global_centroid_m",
                "metrics/global_rg_m",
                "metrics/largest_size",
                "active",
                "wall_gap_m",
                "control/pressure_Pa",
            ]
        }
    return summary, data


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", default="navigation_v1")
    a = ap.parse_args()
    campaign = BASE / a.label
    out = campaign / "analysis"
    out.mkdir(exist_ok=False)
    protocol = js(campaign / "REGISTERED_PROTOCOL.json")
    paths = sorted(p.parent for p in campaign.glob("*/SUMMARY.json"))
    cases = {p.name: load_case(p) for p in paths}
    coarse, ca = cases["zero_flow_feedback"]
    fine, fi = cases["zero_flow_feedback_dt5"]
    common = np.intersect1d(ca["time_s"], fi["time_s"])
    ic = np.searchsorted(ca["time_s"], common)
    iff = np.searchsorted(fi["time_s"], common)
    ce = np.linalg.norm(
        ca["metrics/global_centroid_m"][ic, 0] - fi["metrics/global_centroid_m"][iff, 0], axis=1
    )
    rg = (
        abs(ca["metrics/global_rg_m"][ic, 0] - fi["metrics/global_rg_m"][iff, 0])
        / fi["metrics/global_rg_m"][iff, 0]
    )
    fractions = abs(ca["metrics/largest_size"][ic, 0] - fi["metrics/largest_size"][iff, 0]) / 32
    dt = dict(
        common_frames=len(common),
        common_duration_s=float(common[-1]),
        centroid_max_difference_m=float(ce.max()),
        target_error_max_difference_m=float(
            abs(ca["trajectory/target_error_m"][ic] - fi["trajectory/target_error_m"][iff]).max()
        ),
        Rg_max_relative_difference=float(rg.max()),
        largest_fraction_max_difference=float(fractions.max()),
        coarse_quality=coarse["numerical_quality"],
        fine_quality=fine["numerical_quality"],
    )
    dt["status"] = (
        "PASS"
        if coarse["numerical_quality"] == fine["numerical_quality"] == "PASS"
        and ce.max() <= protocol["dt_centroid_error_limit_m"]
        and rg.max() <= protocol["dt_rg_relative_limit"]
        and fractions.max() <= protocol["dt_largest_fraction_limit"]
        and abs(common[-1] - 1) < 1e-10
        else "FAIL"
    )
    write(out / "TIMESTEP_COMPARISON.json", dt)
    causal_threshold = max(protocol["min_resolvable_control_effect_m"], 3 * ce.max())
    rows = []
    for env in protocol["environments"]:
        prefix = env["name"] + "_"
        off, od = cases[prefix + "off"]
        old, old_d = cases[prefix + "unchanged"]
        for key, (s, d) in cases.items():
            if not key.startswith(prefix):
                continue
            end = d["metrics/global_centroid_m"][-1, 0]
            t = float(d["time_s"][-1])
            oi = int(np.searchsorted(od["time_s"], t))
            ui = int(np.searchsorted(old_d["time_s"], t))
            matched = (
                oi < len(od["time_s"])
                and ui < len(old_d["time_s"])
                and abs(od["time_s"][oi] - t) < 1e-10
                and abs(old_d["time_s"][ui] - t) < 1e-10
            )
            delta_off = np.linalg.norm(end - od["metrics/global_centroid_m"][oi, 0]) if matched else np.nan
            delta_old = np.linalg.norm(end - old_d["metrics/global_centroid_m"][ui, 0]) if matched else np.nan
            supported = bool(
                s["numerical_quality"] == "PASS"
                and s["navigation"]["arrived_and_held"]
                and delta_off > causal_threshold
                and delta_old > causal_threshold
                and (dt["status"] == "PASS" if env["name"] == "zero_flow" else True)
            )
            # Slow-flow has its admitted2.5us quality check, but not a full halved-dt pair.
            row = dict(
                case=key,
                environment=env["name"],
                completed_physical_s=t,
                numerical_quality=s["numerical_quality"],
                failure_reason=s.get("reason"),
                wall_model_enabled=env["wall_enabled"],
                background_flow_enabled=env["flow_enabled"],
                actual_physical_steps=s["actual_substep_count"],
                overlap_max_m=s["quality"]["max_accepted_penetration_m"],
                contact_correction_max_m=s["quality"]["max_contact_pass_correction_m"],
                correction_over_free_Stokes_path=s["quality"].get("correction_over_free_Stokes_path"),
                navigation=s["navigation"],
                dt_s=s["nominal_dt_s"],
                wall_s=s["wall_s"],
                process_body_wall_s=s["process_body_wall_s"],
                gpu_stream_elapsed_s=s["gpu_stream_elapsed_s"],
                difference_from_flow_only_m=float(delta_off),
                difference_from_unchanged_m=float(delta_old),
                minimum_resolvable_difference_m=causal_threshold,
                paired_common_endpoint=matched,
                supported_exploratory_navigation=supported,
                pressure_min_Pa=float(abs(d["control/pressure_Pa"]).min()),
                pressure_max_Pa=float(abs(d["control/pressure_Pa"]).max()),
                scope="EXPLORATORY; not physical calibration",
            )
            row["interpretation"] = (
                "Numerical failure; not evidence of physical inability to navigate"
                if s["numerical_quality"] != "PASS"
                else "Arrival and hold observed; resolution/causality gates evaluated separately"
                if s["navigation"]["arrived_and_held"]
                else "CFD advection dominates tested source controls; no claim of globally optimized unreachability"
                if env["name"] == "Y_CFD"
                else "Tested fixed control does not reach and hold the prescribed target in the registered horizon"
            )
            if "FAILURE_PROPOSAL.json" in [p.name for p in (campaign / key).iterdir()]:
                proposal = js(campaign / key / "FAILURE_PROPOSAL.json")
                row["attempted_physical_steps_including_failed_interval"] = (
                    s["actual_substep_count"] + proposal["attempted_steps"]
                )
                row["failed_interval_bracket_s"] = [t, t + proposal["dt_s"] * proposal["attempted_steps"]]
            rows.append(row)
    write(
        out / "NAVIGATION_RESULTS.json",
        dict(
            rows=rows,
            criterion="Arrival+50ms dwell+>=90% connected, complete numerical PASS, resolved differences from both OFF and unchanged; zero-flow additionally requires key10/5us check",
            threshold_m=causal_threshold,
            do_not_infer_universal_controllability=True,
        ),
    )
    # Existing A is read only: exact archived release suffix. No rerun.
    cs, cd = cases["F5_C_legacy"]
    with h5py.File(F5 / "science_v1/main/trajectory.h5") as h:
        np.testing.assert_array_equal(cd["positions_m"][0], h["positions_m"][600])
        a_end = h["metrics/global_centroid_m"][-1, 0]
        a_rg = h["metrics/global_rg_m"][-1, 0]
        a_largest = int(h["metrics/largest_size"][-1, 0])
    comparison = dict(
        C_quality=cs["numerical_quality"],
        C_complete_s=float(cd["time_s"][-1]),
        same_exact_F5_switch_state=True,
        A_source_sha256=sha(F5 / "science_v1/main/trajectory.h5"),
        A_centroid_final_m=a_end,
        C_centroid_final_m=cd["metrics/global_centroid_m"][-1, 0],
        A_rg_final_m=float(a_rg),
        C_rg_final_m=float(cd["metrics/global_rg_m"][-1, 0]),
        A_largest_final=a_largest,
        C_largest_final=int(cd["metrics/largest_size"][-1, 0]),
        interpretation="Graph cohesion versus flow-only release, not gas fusion/fragmentation",
    )
    write(out / "F5_C_VS_ARCHIVED_A.json", comparison)
    integrity = []
    for path in paths:
        with h5py.File(path / "trajectory.h5") as h:
            groups = ["source", "medium", "flow", "bubbles", "control", "trajectory", "morphology", "quality"]
            assert all(k in h for k in groups)
            assert h.attrs["evidence_level"] == "EXPLORATORY"
            for name in ["positions_m", "time_s", "trajectory/physical_dt_s", "control/pressure_Pa"]:
                assert np.isfinite(h[name][:]).all(), (path, name)
            t = h["time_s"][:]
            assert len(t) == 1 or np.all(np.diff(t) > 0)
            assert h["positions_m"].shape[0] == len(t)
            assert h["particle_id"][:].tolist() == list(range(32))
            assert np.array_equal(
                h["control_index"][:],
                np.searchsorted(
                    [e["time_s"] for e in json.loads(h["control"].attrs["events_json"])],
                    t + 1e-13,
                    side="right",
                )
                - 1,
            )
        sample = read_group_sample(path / "trajectory.h5", len(t) - 1)
        assert sample["flow"]["applied_velocity_m_s"].shape == sample["flow"]["velocity_m_s"].shape
        if not sample["environment"]["flow_enabled"]:
            assert np.all(sample["flow"]["applied_velocity_m_s"] == 0)
        integrity.append(
            dict(
                case=path.name,
                sha256=sha(path / "trajectory.h5"),
                frames=len(t),
                history_frames=len(sample["trajectory"]["time_s"]),
            )
        )
    partitions = validate_partitions([p / "trajectory.h5" for p in paths])
    write(
        out / "DATA_VALIDATION.json",
        dict(
            status="PASS",
            cases=integrity,
            partitions=partitions,
            validation_partition="empty: only two independent target/initial-state families; never split paired controls or frames",
        ),
    )
    # Curves are direct saved observations; plot failed prefixes with dotted lines.
    colors = {
        "off": "#787878",
        "unchanged": "#467aa7",
        "target_phase": "#b779aa",
        "optimized": "#8064a2",
        "feedback": "#2f295f",
    }
    labels = {
        "off": "Ultrasound off",
        "unchanged": "Unchanged control",
        "target_phase": "Target phase",
        "optimized": "Static optimization",
        "feedback": "Feedback",
    }
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
            "savefig.facecolor": "white",
        }
    )
    fig, axs = plt.subplots(2, 3, figsize=(12, 6), layout="constrained")
    plotdata = {}
    for col, env in enumerate(protocol["environments"]):
        for mode in [
            *env["events"],
            *(["feedback"] if env["name"] in protocol["feedback_policy"]["environments"] else []),
        ]:
            key = env["name"] + "_" + mode
            s, d = cases[key]
            t = d["time_s"]
            err = d["trajectory/target_error_m"] * 1e6
            style = "-" if s["numerical_quality"] == "PASS" else ":"
            axs[0, col].plot(t, err, style, color=colors[mode], label=labels[mode])
            axs[1, col].plot(t, d["metrics/largest_size"][:, 0] / 32, style, color=colors[mode])
            plotdata[key + "_time_s"] = t
            plotdata[key + "_error_um"] = err
            plotdata[key + "_largest_fraction"] = d["metrics/largest_size"][:, 0] / 32
        axs[0, col].axhline(2, color="#999999", ls="--", lw=0.8)
        axs[0, col].set(
            title={
                "zero_flow": "No background flow",
                "slow_affine": "Slow affine flow",
                "Y_CFD": "Original Y-vessel CFD",
            }[env["name"]],
            ylabel="Target error (µm)",
            xlabel="Time (s)",
        )
        axs[1, col].set(ylabel="Largest connected fraction", xlabel="Time (s)", ylim=(0, 1.05))
        axs[0, col].legend(fontsize=7)
    fig.savefig(out / "Fig_F6_navigation.png", dpi=250)
    fig.savefig(out / "Fig_F6_navigation.pdf")
    plt.close(fig)
    fig, axs = plt.subplots(1, 3, figsize=(12, 3.5), layout="constrained")
    for key, label in [
        ("zero_flow_feedback", "F5 local mobility"),
        ("zero_flow_feedback_dt5", "Half step"),
        ("zero_flow_feedback_RPY", "Free-space RPY"),
        ("zero_flow_feedback_legacy_gradient", "Legacy nodal gradient"),
    ]:
        s, d = cases[key]
        axs[0].plot(d["time_s"], d["trajectory/target_error_m"] * 1e6, label=label)
    axs[0].axhline(2, color="#999999", ls="--", lw=0.8)
    axs[0].set(xlabel="Time (s)", ylabel="Target error (µm)")
    axs[0].legend(fontsize=7)
    axs[1].plot(common, ce * 1e6, label="Centroid difference")
    axs[1].axhline(0.5, color="#999999", ls="--")
    axs[1].set(xlabel="Time (s)", ylabel="10 vs 5 µs difference (µm)")
    axs[2].plot(common, 100 * rg)
    axs[2].axhline(5, color="#999999", ls="--")
    axs[2].set(xlabel="Time (s)", ylabel="Rg difference (%)")
    fig.savefig(out / "Fig_F6_sensitivity.png", dpi=250)
    fig.savefig(out / "Fig_F6_sensitivity.pdf")
    plt.close(fig)
    np.savez_compressed(
        out / "PLOT_DATA.npz", **plotdata, dt_time_s=common, dt_centroid_difference_m=ce, dt_rg_relative=rg
    )
    admission = js(BASE / "admission_v1_checks/CHECKS.json")
    eff = js(BASE / "admission_v1_efficiency/EFFICIENCY.json")
    ref = js(BASE / "reference_v2/REFERENCE.json")
    campaign_report = js(campaign / "CAMPAIGN.json")
    write(
        out / "ANALYSIS_IDENTITY.json",
        dict(
            analysis_source_sha256=sha(__file__),
            ai_reader_sha256=sha(__import__("mevo_mb.fast_sim.ai_samples", fromlist=["__file__"]).__file__),
            simulation_identity=js(campaign / "IDENTITY.json"),
            archived_campaign_sha256=sha(campaign / "CAMPAIGN.json"),
            postprocessing_only=True,
        ),
    )
    gate = dict(
        software_admission=admission["status"],
        adaptive_and_batch=eff["status"],
        local_translation_reference=ref["status"],
        timestep=dt["status"],
        data="PASS",
        complete_valid_trajectories=campaign_report["valid_trajectories"],
        total_trajectories=campaign_report["total_trajectories"],
        all_formal_numerical_pass=all(s["numerical_quality"] == "PASS" for s, _ in cases.values()),
        supported_navigation_cases=[v["case"] for v in rows if v["supported_exploratory_navigation"]],
        physics_evidence="EXPLORATORY",
        translation_calibration="UNCALIBRATED",
        physical_validation=False,
    )
    gate["validation_acceptance"] = (
        "PASS"
        if all(
            gate[k] == "PASS"
            for k in [
                "software_admission",
                "adaptive_and_batch",
                "local_translation_reference",
                "timestep",
                "data",
            ]
        )
        and gate["all_formal_numerical_pass"]
        else "FAIL"
    )
    gate["navigation_hypothesis"] = (
        "SUPPORTED_ONLY_FOR_LISTED_CASES"
        if gate["supported_navigation_cases"]
        else "NOT_SUPPORTED_BY_THIS_CAMPAIGN"
    )
    write(out / "FINAL_GATE.json", gate)
    print(json.dumps(gate, indent=2))


if __name__ == "__main__":
    main()
