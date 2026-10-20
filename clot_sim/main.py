"""Read this file first: initialization → predictor → native step → audits → commit.

The step equations and acceptance gates are unchanged from the accepted V4.
A short validation run stops after N accepted steps without compressing the load
history. Frozen 3.6 ms evidence is never rerun by the test command.
"""

from pathlib import Path
import argparse
import time
import resource
import traceback
import csv
import gzip
import json
import numpy as np
from .config import (
    BASE,
    MIN_H,
    END,
    RISE,
    K,
    G,
    V,
    INPUTS,
    LIB,
    LABEL,
    HET_LABEL,
    CONFIG_PATH,
)
from .io import save, dump, load_npz, sha, cv, verify_inputs
from .material import SofteningLPS
from .damage import predict
from .fragmentation import Lineage
from .contact import Contact
from .analysis import energy_residual
from .solver import NativeSolver


def run(stage, case, max_steps=None, checkpoint=None):
    start = time.perf_counter()
    cpu = time.process_time()
    case = Path(case).resolve()
    case.mkdir(parents=True, exist_ok=True)
    verify_inputs()
    continuation = stage == "continue"
    target = END if continuation else RISE
    if continuation:
        assert (
            checkpoint is not None
            or json.loads((case / "PILOT_RESULT.json").read_text())["status"]
            == "PILOT_PASS"
        )
        checkpoint = Path(checkpoint) if checkpoint else case / "checkpoints/pilot"
        side = json.loads((checkpoint / "ledger.json").read_text())
        saved = load_npz(checkpoint / "accepted.npz")
        restart = checkpoint / "native.restart"
    else:
        assert not (
            case / "history.csv"
        ).exists(), "Never overwrite an existing physical trajectory"
        side = dict(
            W=0.0,
            Wdrive=0.0,
            Wdrag=0.0,
            Wwall=0.0,
            Jc=0.0,
            Wnew=0.0,
            WL2=0.0,
            scale=1e-24,
            rejects=0,
            substeps=0,
            max_relative_R=0.0,
            max_N20_200=0,
            errors={},
            events={},
            wall_previous_s=0.0,
            cpu_previous_s=0.0,
            snap_next=25e-6,
            checkpoint_next=0.00025,
            registry={},
            splits=[],
        )
        restart = None
        dump(
            case / "REQUEST.json",
            dict(
                formal_case="H_Y_IRREGULAR_LPS_V4",
                one_physical_trajectory=True,
                pilot_end_s=RISE,
                end_s=END,
                material_label=HET_LABEL,
                load_label=LABEL,
                configuration_sha256=sha(CONFIG_PATH),
            ),
        )
    ns = NativeSolver(case, restart=restart)
    ref = ns.ref
    n = ns.n
    b = ref["bonds"]
    bt = ref["bond_types"]
    cut = ref["initial_cut"]
    model = SofteningLPS(ref, K=K, G=G)
    contact = Contact(ref)
    on = ns.field(13)
    off = ns.field(14)
    mass = ns.mass
    s = ns.dynamics()
    surface = ns.tool.evaluate(s["x"], s["velocity"], s["d"], s["time"])
    if continuation:
        for key in ["x", "velocity", "d", "kappa", "atom_ids"]:
            assert np.array_equal(s[key], saved[key]), key
        assert (
            abs(s["scalars"][:6] - saved["scalars"][:6]).max() < 1e-17
            and abs(s["time"] - saved["time"]) < 1e-14
        )
        ns.tool.restore_checkpoint(checkpoint / "surface.npz", s)
        nodework = load_npz(checkpoint / "node_work.npz")["node_work"]
        initial = load_npz(case / "keys/initial.npz")
    else:
        ns.tool.commit(surface)
        initial = s
        nodework = np.zeros(ns.N)
        save(case / "keys/initial.npz", s)
    lineage = Lineage(ref, side if continuation else None)
    lineage.update(surface["topology"], s)
    previous = s
    hist = (case / "history.csv").open("a" if continuation else "w")
    writer = None
    components = gzip.open(
        case / "CONNECTED_COMPONENT_HISTORY.jsonl.gz", "at" if continuation else "wt"
    )
    damage = gzip.open(case / "damage_events.jsonl.gz", "at" if continuation else "wt")

    def error(key, value, tolerance):
        side["errors"][key] = max(side["errors"].get(key, 0.0), float(value))
        assert value <= tolerance, (key, value, tolerance)

    def disk_checkpoint(label, state):
        dest = case / "checkpoints" / label
        dest.mkdir(parents=True, exist_ok=False)
        ns.l.command(f'write_restart {dest/"native.restart"}')
        save(dest / "accepted.npz", state)
        ns.tool.save_checkpoint(dest / "surface.npz", state)
        save(dest / "node_work.npz", dict(node_work=nodework))
        payload = dict(side, **lineage.data())
        payload["wall_previous_s"] = (
            side["wall_previous_s"] + time.perf_counter() - start
        )
        payload["cpu_previous_s"] = side["cpu_previous_s"] + time.process_time() - cpu
        dump(dest / "ledger.json", payload)
        dump(
            case / "LAST_LEGAL_CHECKPOINT.json",
            dict(
                path=str(dest.relative_to(case)),
                time_s=state["time"],
                step=state["step"],
                native_restart_sha256=sha(dest / "native.restart"),
            ),
        )

    if not continuation:
        disk_checkpoint("initial", s)
    accepted_count = 0
    failure = None
    try:
        while previous["time"] < target - 1e-14 and (
            max_steps is None or accepted_count < max_steps
        ):
            p = previous
            boundary = min(t for t in ns.tool.boundaries if t > p["time"] + 1e-14)
            startup = (
                BASE / 8
                if p["time"] < 2e-6 - 1e-15
                else (
                    BASE / 4
                    if p["time"] < 4e-6 - 1e-15
                    else BASE / 2 if p["time"] < 8e-6 - 1e-15 else BASE
                )
            )
            h = min(startup, boundary - p["time"])
            depth = 0
            while True:
                xp = p["x"].copy()
                xp[:n] += h * (
                    p["velocity"][:n] + h * p["force"][:n] / (2 * mass[:n, None])
                )
                kp, dp = predict(model, ref, xp, p["kappa"], on, off)
                if np.max(dp - p["d"]) <= 0.1:
                    break
                assert h > MIN_H * (1 + 1e-12), "MINIMUM_EVENT_SUBSTEP_EXHAUSTED"
                h /= 2
                depth += 1
                side["rejects"] += 1
            side["substeps"] += int(depth > 0)
            s = ns.step(h)
            dd = s["d"] - p["d"]
            sc = s["scalars"]
            nf = np.flatnonzero((s["d"] == 1) & (p["d"] < 1) & ~cut)
            error("predicted_x_m", abs(xp - s["x"]).max(), 2e-18)
            error("predicted_damage", abs(dp - s["d"]).max(), 1e-10)
            error("predicted_kappa", abs(kp - s["kappa"]).max(), 1e-10)
            error(
                "irreversibility",
                max(0.0, float((-dd).max()), float((p["kappa"] - s["kappa"]).max())),
                1e-10,
            )
            error(
                "velocity_VV_m_s",
                abs(
                    s["velocity"][:n]
                    - p["velocity"][:n]
                    - h * (p["force"][:n] + s["force"][:n]) / (2 * mass[:n, None])
                ).max(),
                1e-10,
            )
            error(
                "slip_velocity_m_s",
                abs(s["force_stage_velocity"] - s["velocity"]).max(),
                1e-10,
            )
            assert (
                dd.max() <= 0.1 + 1e-10
                and np.array_equal(s["x"][n:], ref["X"][n:])
                and not s["velocity"][n:].any()
            )
            assert (
                s["native_bonds"] == int((s["d"] < 1).sum())
                and not s["d"][bt == 3].any()
                and np.isfinite(s["force"]).all()
            )
            work = ns.tool.work(p, s)
            assert work[1] <= 1e-17
            W = side["W"] + float(work.sum())
            Jc = side["Jc"] + contact.energy_for_pairs(s["x"], b[nf])
            res = energy_residual(W, s, initial, Jc)
            scale = max(
                side["scale"],
                abs(W),
                abs(sc[0] - initial["scalars"][0]) + s["K"] + sc[1] + s["Uc"],
            )
            rel = abs(res) / scale
            assert rel <= 0.01, ("ENERGY_RESIDUAL", rel, res, scale)
            if s["step"] % 100 == 0 or np.any(dd):
                new = model.evaluate(s["x"], s["d"])
                old = model.evaluate(s["x"], p["d"]) if np.any(dd) else new
                uc, fc, pairs = contact.evaluate(s["x"], s["d"])
                error(
                    "native_force_N",
                    abs(new["force"] - s["material_force"]).max(),
                    1e-12,
                )
                error(
                    "force_port_N",
                    abs(
                        s["force"] - s["material_force"] - fc - s["external_force"]
                    ).max(),
                    1e-12,
                )
                error("native_U_J", abs(new["U"] - sc[0]), 1e-17)
                error(
                    "dissipation_increment_J",
                    abs(old["U"] - new["U"] - (sc[1] - p["scalars"][1])),
                    1e-17,
                )
                error("contact_U_J", abs(uc - s["Uc"]), 1e-17)
            surface = ns.tool.evaluate(s["x"], s["velocity"], s["d"], s["time"])
            error(
                "surface_replay_force_N",
                abs(
                    surface["force"] + surface["wall_force"] - s["external_force"]
                ).max(),
                1e-15,
            )
            assert s["failed_tension_max_N"] <= 1e-20
            assert s["wall_clearance"].min() >= -2e-6, "WALL_CENTER_PENETRATION_GATE"
            # All gates passed; only here commit diagnostics and the external work ledger.
            ns.tool.commit(surface)
            dx = s["x"] - p["x"]
            nodework += 0.5 * np.sum(
                (p["external_force"] + s["external_force"]) * dx, axis=1
            )
            side.update(
                W=W,
                Jc=Jc,
                scale=scale,
                Wdrive=side["Wdrive"] + float(work[0]),
                Wdrag=side["Wdrag"] + float(work[1]),
                Wwall=side["Wwall"] + float(work[2]),
                max_relative_R=max(side["max_relative_R"], float(rel)),
            )
            side["Wnew"] += 0.5 * float(
                np.sum(
                    (p["external_accessible_force"] + s["external_accessible_force"])
                    * dx
                )
            )
            side["WL2"] += 0.5 * float(
                np.sum((p["external_L2_force"] + s["external_L2_force"]) * dx)
            )
            error("node_work_sum_J", abs(nodework.sum() - W), 1e-17)
            top = surface["topology"]
            if len(nf):
                lineage.update(top, s)
            for key, condition in [
                ("first_softening", np.any((s["d"] > 0) & ~cut)),
                ("first_failure", len(nf) > 0),
                ("first_free_component", bool(lineage.current)),
            ]:
                if condition and key not in side["events"]:
                    side["events"][key] = dict(
                        time_s=s["time"],
                        step=s["step"],
                        evidence="keys/" + key + ".npz",
                    )
                    save(case / "keys" / (key + ".npz"), s)
            freeids = np.flatnonzero(top["free"][:n])
            sizes = np.bincount(top["labels"], minlength=len(top["attached"]))
            free_sizes = sizes[~top["attached"]]
            target_labels = np.flatnonzero(
                ~top["attached"] & (sizes >= 20) & (sizes <= 200)
            )
            target_ids = np.flatnonzero(np.isin(top["labels"][:n], target_labels))
            Ntarget = len(target_labels)
            if Ntarget and "first_target" not in side["events"]:
                side["events"]["first_target"] = dict(
                    time_s=s["time"],
                    step=s["step"],
                    evidence="keys/first_target_observed.npz",
                )
                save(case / "keys/first_target_observed.npz", s)
            if Ntarget > side["max_N20_200"]:
                side["max_N20_200"] = Ntarget
                save(case / "keys/max_target_count.npz", s)
            geom = ns.tool.metrics(surface, s["x"])
            row = dict(
                step=s["step"],
                time_s=s["time"],
                dt_s=h,
                phase=ns.tool.phase(s["time"]),
                W_J=W,
                W_drive_J=side["Wdrive"],
                W_drag_J=side["Wdrag"],
                W_wall_J=side["Wwall"],
                wall_potential_J=float(s["wall_potential_J"]),
                min_wall_clearance_m=float(s["wall_clearance"].min()),
                wall_contact_nodes=int((s["wall_clearance"] < 20e-6).sum()),
                U_J=float(sc[0]),
                K_J=s["K"],
                D_J=float(sc[1]),
                Uc_J=s["Uc"],
                Jc_J=Jc,
                R_J=float(res),
                relative_R=float(rel),
                numerical_damping_work_J=0.0,
                W_accessible_new_faces_J=side["Wnew"],
                W_L2_surface_J=side["WL2"],
                bulk_failed=int(((s["d"] == 1) & ~cut & (bt == 1)).sum()),
                interface_failed=int(((s["d"] == 1) & ~cut & (bt == 2)).sum()),
                softened=int(((s["d"] > 0) & ~cut).sum()),
                detached_nodes=len(freeids),
                free_components=len(free_sizes),
                largest_nodes=int(free_sizes.max(initial=0)),
                N20_200=Ntarget,
                detached_volume_m3=float(ref["V"][freeids].sum()),
                detached_mass_kg=float(mass[freeids].sum()),
                target_mass_kg=float(mass[target_ids].sum()),
                target_volume_m3=float(ref["V"][target_ids].sum()),
                target_equivalent_nominal_nodes=float(ref["V"][target_ids].sum() / V),
                max_dd=float(dd.max()),
                max_speed_m_s=float(np.linalg.norm(s["velocity"], axis=1).max()),
                rejects=side["rejects"],
                substeps=side["substeps"],
                slip_iterations=int(s["slip_iterations"]),
                **geom,
            )
            if writer is None:
                writer = csv.DictWriter(hist, fieldnames=list(row))
                if not continuation:
                    writer.writeheader()
            writer.writerow(row)
            changed = np.flatnonzero(dd > 0)
            if len(changed):
                damage.write(
                    json.dumps(
                        dict(
                            step=s["step"],
                            time_s=s["time"],
                            original_bond_ids=changed + 1,
                            d=s["d"][changed],
                            kappa=s["kappa"][changed],
                        ),
                        default=cv,
                    )
                    + "\n"
                )
            if (
                s["time"] >= side["snap_next"] - 1e-14
                or abs(s["time"] - boundary) < 1e-14
            ):
                save(case / "keys" / f'frame_{s["step"]:06d}.npz', s)
                components.write(
                    json.dumps(
                        dict(
                            time_s=s["time"],
                            step=s["step"],
                            components=lineage.records(s, top, nodework),
                        ),
                        default=cv,
                    )
                    + "\n"
                )
                side["snap_next"] = max(side["snap_next"] + 25e-6, s["time"] + 25e-6)
            previous = s
            accepted_count += 1
            if (
                s["time"] >= side["checkpoint_next"] - 1e-14
                and abs(s["time"] - target) > 1e-14
            ):
                side["checkpoint_next"] += 0.00025
                disk_checkpoint(f'step_{s["step"]:06d}', s)
            if s["step"] % 250 == 0:
                hist.flush()
                components.flush()
                damage.flush()
                dump(
                    case / "PROGRESS.json",
                    dict(
                        status="RUNNING",
                        stage=stage,
                        time_s=s["time"],
                        step=s["step"],
                        accepted_steps=s["step"],
                        rejected_predictions=side["rejects"],
                        max_relative_R=side["max_relative_R"],
                        N20_200=Ntarget,
                        largest_nodes=row["largest_nodes"],
                        wall_s=side["wall_previous_s"] + time.perf_counter() - start,
                    ),
                )
                print(
                    "STEP",
                    s["step"],
                    s["time"],
                    "N20",
                    Ntarget,
                    "largest",
                    row["largest_nodes"],
                    "R",
                    rel,
                    "wall",
                    time.perf_counter() - start,
                    flush=True,
                )
        name = (
            ("final" if continuation else "pilot")
            if previous["time"] >= target - 1e-14
            else f"short_validation_step_{previous['step']:06d}"
        )
        save(case / "keys" / (name + ".npz"), previous)
        disk_checkpoint(name, previous)
    except BaseException:
        failure = traceback.format_exc()
        save(case / "last_accepted.npz", previous)
        if s is not previous:
            save(case / "failed_native_trial.npz", s)
        (case / "FAILURE.txt").write_text(failure)
        print(failure, flush=True)
    finally:
        hist.close()
        components.close()
        damage.close()
        ns.close()
    result = dict(
        status=(
            "NUMERICAL_FAILURE"
            if failure
            else (
                "SHORT_VALIDATION_PASS"
                if previous["time"] < target - 1e-14
                else "COMPLETED" if continuation else "PILOT_PASS"
            )
        ),
        stage=stage,
        accepted_time_s=previous["time"],
        target_time_s=target,
        formal_cases=int(max_steps is None),
        accepted_steps=previous["step"],
        stage_accepted_steps=accepted_count,
        rejected_predictions=side["rejects"],
        substeps=side["substeps"],
        max_relative_R=side["max_relative_R"],
        max_N20_200=side["max_N20_200"],
        events=side["events"],
        errors=side["errors"],
        wall_s=side["wall_previous_s"] + time.perf_counter() - start,
        cpu_s=side["cpu_previous_s"] + time.process_time() - cpu,
        peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024,
        exception=failure,
        reference_sha256=sha(INPUTS / "reference.npz"),
        library_sha256=sha(LIB),
        model_label=HET_LABEL,
        load_label=LABEL,
    )
    dump(
        case
        / ("NATIVE_RESULT.json" if continuation or failure else "PILOT_RESULT.json"),
        result,
    )
    dump(case / "LINEAGE.json", lineage.data())
    print("STAGE_RESULT", json.dumps(result, default=cv), flush=True)
    if failure:
        raise RuntimeError(failure)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    inspect = commands.add_parser(
        "inspect", help="Verify frozen inputs and summarize the reference geometry"
    )
    run_parser = commands.add_parser(
        "run", help="Run unchanged V4 dynamics into a new output directory"
    )
    run_parser.add_argument("--output", type=Path, required=True)
    run_parser.add_argument("--stage", choices=["pilot", "continue"], default="pilot")
    run_parser.add_argument(
        "--steps",
        type=int,
        help="Stop after this many accepted steps; keep the original time/load schedule",
    )
    run_parser.add_argument(
        "--checkpoint",
        type=Path,
        help="Continue this checkpoint in the same output case",
    )
    args = parser.parse_args()
    if args.command == "inspect":
        from .analysis import reference_summary

        verify_inputs()
        print(json.dumps(reference_summary(), indent=2))
        return
    if args.steps is not None and args.steps < 1:
        parser.error("--steps must be positive")
    run(args.stage, args.output, args.steps, args.checkpoint)


if __name__ == "__main__":
    main()
