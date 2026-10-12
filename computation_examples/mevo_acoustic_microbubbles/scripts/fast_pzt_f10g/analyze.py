"""Actual-data admission, common-prefix morphology and diagnostic cost ledger."""

import gzip
import json
import numpy as np
from scipy import linalg
from scripts.fast_pzt_f10f.analyze import compare as prefix_compare
from mevo_mb.fast_sim.f10g_contact import resistance, compare, shape
from .common import BASE, ROOT, CASES, state, js, write, identity


def gz(p):
    with gzip.open(p, "rt") as f:
        return json.load(f)


def main():
    out = BASE / "analysis"
    out.mkdir(exist_ok=True)
    lubrication = js(BASE / "lubrication_v1/RESULTS.json")
    prototype = js(BASE / "prototype_v1/DENSE_STATIC.json")
    costs = js(BASE / "integration_v1/RESULTS.json")
    small = js(BASE / "prototype_v1/SMALL_SYSTEMS.json")
    morph = []
    trajectories = []
    spectra = []
    stiffness = []
    for case in CASES:
        c = prefix_compare(
            ROOT / f"results/fast_pzt/f10e/contact_v1/{case}_C0", BASE / f"integration_v1/{case}_B2_K2"
        )
        c["interpretation"] = (
            "Different contact and mobility models: descriptive agreement with old accepted prefix, NOT a reference truth or same-model temporal convergence test."
        )
        c["contact_path_interpretation"] = (
            "PBD iteration path versus implicit Newton iteration path; algorithmic workload diagnostics, not the same physical contact observable."
        )
        morph.append(c)
        steps = gz(BASE / f"integration_v1/{case}_B2_K2/SOFT_STEPS.json.gz")
        ledger = gz(BASE / f"integration_v1/{case}_B2_K2/SUBSTEP_LEDGER.json.gz")["steps"]
        assert all(r["accepted"] and not r["error_bits"] for r in steps)
        assert all(
            r["accepted"] and r["final_overlap_m"] <= 17.8e-9 and r["wall_correction_max_m"] == 0
            for r in ledger
        )
        maxlaw = 0.0
        maxnet = 0.0
        maxmean = 0.0
        maxmove = 0.0
        geometry_offset = 0.0
        frames = []
        geometric_noncontacts = []
        max_torque = 0.0
        for s in steps:
            z = s["result"]
            f = np.array(z["pair_force_N"])
            delta = np.array(z["geometric_overlap_m"])
            maxlaw = max(maxlaw, float(np.max(np.abs(f - z["kc_N_m"] * delta))))
            pos = np.array(z["positions_m"])
            contacts = np.array(z["pairs"])
            mask = (delta == 0) & (f > 0)
            for pair, force in zip(contacts[mask], f[mask]):
                geometric_noncontacts.append(
                    dict(
                        step=s["step"],
                        pair=pair,
                        force_N=force,
                        gap_m=float(np.linalg.norm(pos[pair[1]] - pos[pair[0]]) - 3.56e-6),
                    )
                )
            torque = np.sum(np.cross(pos - pos.mean(0), np.array(z["contact_force_N"])), axis=0)
            max_torque = max(max_torque, float(np.linalg.norm(torque)))
            geometry_offset = max(
                geometry_offset, float(np.max(np.abs(np.array(z["frozen_normal_overlap_m"]) - delta)))
            )
            maxnet = max(maxnet, float(np.linalg.norm(z["net_contact_force_N"])))
            maxmean = max(maxmean, float(np.linalg.norm(z["centroid_contact_shift_m"])))
            maxmove = max(maxmove, z["max_single_iteration_m"])
            frames.append(
                dict(
                    step=s["step"],
                    time_s=s["step"] * 1e-5,
                    shape=shape(np.array(s["actual_accepted_positions_m"]), 1.78e-6),
                )
            )
        trajectories.append(
            dict(
                case=case,
                accepted_steps=len(steps),
                end_s=len(steps) * 1e-5,
                max_overlap_m=max(r["final_overlap_m"] for r in ledger),
                max_single_iteration_m=maxmove,
                max_frozen_vs_geometric_overlap_m=geometry_offset,
                max_actual_geometry_spring_law_force_error_N=maxlaw,
                max_net_contact_force_N=maxnet,
                geometric_noncontact_reactions=geometric_noncontacts,
                max_contact_torque_N_m=max_torque,
                max_contact_centroid_shift_m=maxmean,
                force_law_scope="Force uses initial frozen normals and linearized overlap. Actual geometric overlap also checked independently every step.",
                total_step_wall_s=sum(s["full_step_s"] for s in steps),
                field_plus_input_readback_s=sum(s["field_and_read_s"] for s in steps),
                resistance_assembly_s=sum(s["resistance_assembly_s"] for s in steps),
                soft_initialization_s=sum(s["result"]["initialization_s"] for s in steps),
                soft_contact_s=sum(s["result"]["contact_s"] for s in steps),
                commit_plus_readback_s=sum(s["commit_and_read_s"] for s in steps),
                mean_newton_iterations=float(np.mean([len(s["result"]["iterations"]) for s in steps])),
                max_newton_iterations=max(len(s["result"]["iterations"]) for s in steps),
                exact_step_shapes=frames,
            )
        )
        _, t, _ = state(case)
        op = resistance(np.array(t["initial_positions_m"]), 1.78e-6)
        Dinv = np.zeros_like(op["A"])
        for i, block in enumerate(op["blocks"]):
            Dinv[3 * i : 3 * i + 3, 3 * i : 3 * i + 3] = np.linalg.inv(block)
        spectral_radius = float(np.max(np.abs(linalg.eigvals(np.eye(96) - 0.5 * Dinv @ op["A"]))))
        spectra.append(
            dict(
                case=case,
                relaxed_Jacobi_spectral_radius=spectral_radius,
                interpretation="Slow low-frequency collective modes remain; eight passes not an algebraic reference",
            )
        )
        for mode in ["B1", "B2"]:
            selected = sorted(
                [r for r in prototype["rows"] if r["case"] == case and r["mode"] == mode],
                key=lambda r: r["multiplier"],
            )
            stiffness.append(
                dict(
                    case=case,
                    mode=mode,
                    half_vs_double=compare(
                        np.array(selected[0]["result"]["positions_m"]),
                        np.array(selected[-1]["result"]["positions_m"]),
                    ),
                    rows=[
                        dict(
                            multiplier=s["multiplier"],
                            max_overlap_m=max(s["result"]["geometric_overlap_m"]),
                            max_correction_m=s["result"]["max_single_iteration_m"],
                            pass_=s["contact_budget_pass"],
                            shape=s["shape"],
                            self_solver_tolerance_check=s["tight_reference"],
                            timing=s["timing"],
                        )
                        for s in selected
                    ],
                )
            )
    processes = [dict(label=p.name, **js(p)) for p in BASE.glob("*_PROCESS.json")]
    write(
        out / "MORPHOLOGY.json",
        dict(
            comparisons=morph,
            stiffness=stiffness,
            solver_tolerance_convergence="PASS_FIXED_STATE_ONLY",
            temporal_convergence="NOT_EVALUATED",
            training_ready=False,
        ),
    )
    write(
        out / "LUBRICATION.json",
        dict(
            spectra=spectra,
            results=lubrication,
            interpretation="LUBRICATION_CONVERGENCE_CONTACT_RELEVANT on the same two frozen states; no Jacobi4 trajectory was run.",
        ),
    )
    write(out / "TRAJECTORY_AUDIT.json", dict(rows=trajectories, training_ready=False))
    write(
        out / "COST_AND_YIELD.json",
        dict(
            processes=processes,
            full_process_wall_s=sum(p["full_process_wall_s"] for p in processes),
            trajectory_costs=costs,
            short_numerical_completions=2,
            complete_10ms_samples=0,
            independent_new_seeds=0,
            scientific_failures=0,
            deliberate_engineering_rejections=2,
            full_matched_model_OFF_S_O_triplets=0,
            production_valid_sample_rate_per_hour=None,
            cost_per_calibrated_sample=None,
            scope="Do not use legacy recorder valid_samples_per_gpu_hour for these unequal short diagnostic windows; no production sample-rate claim.",
        ),
    )
    gate = dict(
        status="ENGINEERING_COMPLETE_BOUNDED_CONTACT_PASS",
        LUBRICATION_REFERENCE_PASS=True,
        LUBRICATION_CONVERGENCE_CHECKED=True,
        LUBRICATION_CONVERGENCE_NOT_DECISIVE=False,
        SOFT_CONTACT_TWO_BUBBLE_PASS=small["SOFT_CONTACT_TWO_BUBBLE_PASS"],
        SOFT_CONTACT_FOUR_BUBBLE_PASS=small["SOFT_CONTACT_FOUR_BUBBLE_PASS"],
        SOFT_CONTACT_32_STATIC_PASS="B2_K2_ONLY_BOTH_CASES",
        SOFT_CONTACT_SHORT_TRAJECTORY_PASS="B2_K2: O1102 5ms; preparation1104 2.35ms",
        CONTACT_NUMERICAL_PASS="BOUNDED_WINDOWS_ONLY",
        MORPHOLOGY_NOT_CONVERGED=True,
        PRODUCTION_NOT_ADMITTED=True,
        PHYSICS_UNCALIBRATED=True,
        GEOMETRIC_SPRING_LAW_EXACT=False,
        full_10ms_admitted=False,
        default_unchanged=True,
        training_ready=False,
        reasons=[
            "Two fixed-state and two short-window contact checks passed, without relaxed budgets.",
            "Between-model fixed-ID/tensor/topology differences exceed inherited descriptive thresholds; time-step convergence of new model not measured.",
            "CPU 96DOF reference is not a production GPU contact implementation; matched full end-to-end production cost is unmeasured.",
            "Frozen-normal spring is a numerical contact candidate; shell contact stiffness and near-field forces remain uncalibrated.",
            "Linearized contact may still exert small repulsion after true geometric separation: 1 / 16 pair-step events. Not an exact geometric unilateral spring.",
        ],
        source=identity(),
    )
    write(out / "FINAL_GATE.json", gate)
    for c in morph:
        print(
            c["reference"].split("/")[-1],
            c["accepted_common_end_s"],
            {
                k: v
                for k, v in c["values"].items()
                if k
                in [
                    "centroid_max_m",
                    "fixed_ID_max_m",
                    "Kabsch_max_m",
                    "Rg_max_relative",
                    "aspect_max_relative",
                    "event_count",
                ]
            },
        )
    for row in trajectories:
        print(
            row["case"],
            "overlap nm",
            row["max_overlap_m"] * 1e9,
            "pass nm",
            row["max_single_iteration_m"] * 1e9,
            "frozengeom nm",
            row["max_frozen_vs_geometric_overlap_m"] * 1e9,
        )


if __name__ == "__main__":
    main()
