"""One admitted cross-model counterfactual, conditional dt/2; no old-batch resume."""

import argparse
from time import perf_counter
import h5py
import numpy as np
from mevo_mb.fast_sim.f10c_engine import FullParticleBatch
from mevo_mb.fast_sim.f10h_diagnostic import GeometricContactDiagnostic
from mevo_mb.fast_sim.f10d_metrics import process
from scripts.fast_pzt_f10c.runner import run
from scripts.fast_pzt_f9.common import objects
from .common import ROOT, BASE, identity, js, write, sha


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--label", required=True)
    ap.add_argument("--half", action="store_true")
    ap.add_argument("--recover", action="store_true")
    args = ap.parse_args()
    out = BASE / args.label
    out.mkdir(parents=True, exist_ok=False)
    started = perf_counter()
    source = identity()
    write(out / "IDENTITY.json", source)
    cfg = js(ROOT / "configs/fast_pzt/f10k/topology_case.json")
    if args.recover:
        parent = BASE / "topology_dt10_v1"
        summ = js(parent / "trial/SUMMARY.json")
        old = js(parent / "IDENTITY.json")
        with h5py.File(parent / "trial/trajectory.h5") as h:
            assert h["time_s"][-1] == cfg["duration_s"]
            assert h["substeps/accepted_end_s"].shape[0] == 300
        assert summ["status"] == "COMPLETE" and summ["numerical_quality"] == "NUMERICAL_PASS"
        assert js(parent / "RESTART_CHECK.json")["status"] == "PASS"
        counts = summ["topology"][0]["counts"]
        write(
            out / "RESULTS.json",
            dict(
                status=summ["status"],
                numerical_quality=summ["numerical_quality"],
                topology=summ["topology"],
                persistent_events=counts["merge_persistent"] + counts["split_persistent"],
                absolute_start_s=cfg["selected"]["absolute_start_s"],
                dt_s=1e-5,
                duration_s=cfg["duration_s"],
                timing=js(parent / "trial/TIMING.json"),
                scientific_execution_commit=old["execution_commit"],
                source_path=str((parent / "trial").relative_to(ROOT)),
                recovered_only_metadata=True,
                no_science_rerun=True,
                failure="CONFIG.json write requested missing identity protocol_sha256 AFTER all trajectory, summary and diagnostics had been durably written.",
                preserved_files={
                    str(p.relative_to(parent)): sha(p) for p in parent.rglob("*") if p.is_file()
                },
                training_ready=False,
            ),
        )
        write(
            out / "RECOVERED_CONFIG.json",
            dict(
                initial_positions_source=cfg["selected"],
                events=cfg["rebased_events"],
                duration_s=cfg["duration_s"],
                dt_s=1e-5,
                protocol_sha256=sha(ROOT / "configs/fast_pzt/f10k/protocol.json"),
                reconstruct_from_immutable_execution=old["execution_commit"],
            ),
        )
        print("Recovered metadata only; original scientific files unchanged", flush=True)
        return
    selected = cfg["selected"]
    parent = ROOT / selected["source"]
    assert sha(parent) == selected["sha256"] and not selected["reasons"]
    with h5py.File(parent) as h:
        x = h["positions_m"][selected["frame"], 0]
        successor = h["positions_m"][selected["frame"] + 1, 0]
        interval = float(h["time_s"][selected["frame"] + 1] - h["time_s"][selected["frame"]])
        np.testing.assert_array_equal(h["particle_id"][:], np.arange(32))
        source_active = h["trajectory/active"][selected["frame"], 0]
    assert source_active.all()
    events = cfg["rebased_events"]
    obj = objects()
    p, b, r, g, f = obj
    # The adapter uses a fresh local clock with a preserved absolute-time offset.
    # All selected controls are constant; averaged acoustics/steady CFD are autonomous.
    assert len(events) == 1 and events[0]["time_s"] == 0
    if not args.half:
        replay = FullParticleBatch(x, p, r, b, vessel=g, flow=f)
        replay.control(events[0])
        replay.advance(1e-5, round(interval / 1e-5), graphs=False)
        err = float(np.linalg.norm(replay.x.numpy() - successor, axis=1).max())
        same = np.array_equal(replay.x.numpy(), successor)
        spatial = g.query(x)
        cf = f.sample(x)
        record = dict(
            status="PASS" if err <= 1e-15 and replay.quality()["error_bits"] == 0 else "FAIL",
            reconstructed_M0_interval_s=interval,
            max_position_difference_m=err,
            bitwise_positions=same,
            quality=replay.quality(),
            CFD_valid=cf["valid"],
            wall_gap_m=spatial["distance_m"] - p.radius_m,
            absolute_start_s=selected["absolute_start_s"],
            local_start_s=0.0,
            exact_B2_restart_claim=False,
            control=events[0],
            statement="Same M0 saved interval re-evaluated with fresh caches. B2 starts a new Markov counterfactual, does not inherit M0 contacts or persistence. Reset statistics do not drive motion.",
            source_sha256=sha(parent),
            replayed_positions_m=replay.x.numpy(),
            archived_successor_m=successor,
        )
        write(out / "RESTART_CHECK.json", record)
        if record["status"] != "PASS":
            write(
                out / "RESULTS.json",
                dict(
                    status="B2_TOPOLOGY_CASE_NOT_ADMITTED",
                    reason="RESTART_CHECK_FAILED",
                    training_ready=False,
                ),
            )
            return
    else:
        previous = js(BASE / "topology_recovery_v1/RESULTS.json")
        assert previous["status"] == "COMPLETE" and previous["numerical_quality"] == "NUMERICAL_PASS"
        assert previous["persistent_events"] > 0, (
            "No positive event; dt/2 not authorized by preregistered value gate"
        )
        assert js(BASE / "topology_dt10_v1/RESTART_CHECK.json")["status"] == "PASS"
    dt = 5e-6 if args.half else 1e-5
    mode = "B2_GEOM_K2_DT5" if args.half else "B2_GEOM_K2_DT10"
    timing = run(
        out / "trial",
        x,
        events,
        obj,
        duration=cfg["duration_s"],
        physical_dt_s=dt,
        graphs=False,
        simulator_class=GeometricContactDiagnostic,
        numerical_mode=mode,
        substep_protocol=True,
        export_observer=True,
        kc_N_m=cfg["kc_N_m"],
        sample_ids=["F10K_counterfactual_F10D_1101_frame42"],
        postprocess=process,
        identity_provider=lambda: source,
        dataset_metadata=dict(
            initial_state_family=["F10D_seed_1101_preparation"],
            split_groups=["F10D_seed_1101"],
            split_group="F10D_seed_1101",
            labels=["EXPLORATORY", "PHYSICS_UNCALIBRATED", "CROSS_MODEL_INITIALIZATION"],
            independent_scientific_initial_states=0,
            initial_state_sources=[
                dict(parent_sha256=sha(parent), saved_frame=selected["frame"], exact_restart=False)
            ],
            original_physics_mode=mode,
            postprocessing_mode="F10D_optimized",
        ),
    )
    summ = js(out / "trial/SUMMARY.json")
    with h5py.File(out / "trial/trajectory.h5", "a") as h:
        h["absolute_time_s"] = h["time_s"][:] + selected["absolute_start_s"]
        h["edge_events/absolute_time_s"] = h["edge_events/time_s"][:] + selected["absolute_start_s"]
        h.attrs["absolute_start_s"] = selected["absolute_start_s"]
        h.attrs["restart_semantics"] = (
            "CROSS_MODEL_EXPLORATORY_COUNTERFACTUAL; new B2 graph history; no past M0 dwell transfer"
        )
    counts = (
        summ["topology"][0]["counts"] if isinstance(summ["topology"], list) else summ["topology"]["counts"]
    )
    write(
        out / "RESULTS.json",
        dict(
            status=summ["status"],
            numerical_quality=summ["numerical_quality"],
            topology=summ["topology"],
            persistent_events=counts["merge_persistent"] + counts["split_persistent"],
            absolute_start_s=selected["absolute_start_s"],
            dt_s=dt,
            duration_s=cfg["duration_s"],
            timing=timing,
            whole_job_s=perf_counter() - started,
            training_ready=False,
        ),
    )
    print(summ["status"], summ["numerical_quality"], counts, flush=True)


if __name__ == "__main__":
    main()
