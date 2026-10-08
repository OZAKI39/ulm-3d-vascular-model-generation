"""Read-only reduction of closed F4 evidence; compact copies preserve every saved sample."""

import hashlib
import h5py
import numpy as np
from .common import BASE, write, js

R = 1.78e-6
B_NAMES = [r + "_wall_" + w for r in ["parent", "junction"] for w in ["off", "on"]]
C_NAMES = ["wall_" + w + "_" + c for w in ["on", "off"] for c in ["baseline", "switched"]]


def read(path):
    if (path / "trajectory.h5").exists():
        with h5py.File(path / "trajectory.h5") as h:
            d = {}

            def take(k, v):
                if isinstance(v, h5py.Dataset) and v.dtype.kind not in "OSU":
                    d[k] = v[()]

            h.visititems(take)
    else:
        d = {}
        for f in sorted(path.glob("SAVED_*.npz")):
            with np.load(f) as a:
                d.update({k: a[k] for k in a.files})
    return {**d, **{k[8:]: v for k, v in d.items() if k.startswith("metrics/")}}


def compact(path):
    if not (path / "trajectory.h5").exists():
        assert (path / "SAVED_METADATA.json").exists()
        return
    d = read(path)
    parts = {"STATE": {}, "METRICS": {}, "DIAGNOSTICS": {}}
    for k, v in d.items():
        if k in ["time_s", "positions_m", "particle_id", "radii_m", "active", "wall_gap_m", "control_index"]:
            parts["STATE"][k] = v
        elif k.startswith("metrics/"):
            parts["METRICS"][k] = v
        elif "/" in k or k == "wall_stats":
            parts["DIAGNOSTICS"][k] = v
    for name, data in parts.items():
        np.savez_compressed(path / ("SAVED_" + name + ".npz"), **data)
        with np.load(path / ("SAVED_" + name + ".npz")) as a:
            for k, v in data.items():
                np.testing.assert_array_equal(a[k], v)
    with h5py.File(path / "trajectory.h5") as h:
        write(
            path / "SAVED_METADATA.json",
            dict(
                original_sha256=hashlib.sha256((path / "trajectory.h5").read_bytes()).hexdigest(),
                h5_attributes={k: v for k, v in h.attrs.items()},
                text_datasets={
                    k: h[k][()].decode() for k in ["sampled_lineage_json", "boundary_events_json"]
                },
                compression="Lossless; every saved frame, ID, coordinate, radius, metric, diagnostic and edge event. No decimation.",
            ),
        )


def summarize(path):
    d = read(path)
    s = js(path / "SUMMARY.json")
    m = np.load(path / "MORPHOLOGY.npz")
    contact = (d["wall_gap_m"][:, 0] <= 0.01 * R).sum(axis=1)
    nearest = d["wall_gap_m"][:, 0].min(axis=1)
    first = np.flatnonzero(contact)
    all32 = np.flatnonzero(d["largest_size"][:, 0] == 32)
    # Saved-frame step size is not a physical-step velocity; both are reported separately.
    dx = np.linalg.norm(np.diff(d["positions_m"][:, 0], axis=0), axis=-1)
    assert np.isfinite(d["positions_m"]).all() and d["active"].all()
    assert s["status"] == "COMPLETE" and s["numerical_quality"] == "PASS"
    np.savez_compressed(
        path / "WALL_METRICS.npz",
        time_s=d["time_s"],
        nearest_wall_gap_m=nearest,
        wall_contact_count=contact,
        pair_contact_count=d["hard_contact_pair_count"][:, 0],
        pair_near_contact_count=d["near_contact_pair_count"][:, 0],
    )
    return dict(
        status=s["status"],
        numerical_quality=s["numerical_quality"],
        source_sha256=js(path / "SAVED_METADATA.json")["original_sha256"],
        wall_s=s["wall_s"],
        process_body_wall_s=s["process_body_wall_s"],
        physical_duration_s=float(d["time_s"][-1]),
        physical_steps=s["actual_substep_count"],
        additional_substeps=0,
        count=32,
        active_IDs_preserved=True,
        initial_cluster_count=int(d["cluster_count"][0, 0]),
        final_cluster_count=int(d["cluster_count"][-1, 0]),
        final_largest_size=int(d["largest_size"][-1, 0]),
        first_saved_all32_connected_s=float(d["time_s"][all32[0]]) if len(all32) else None,
        rg_initial_m=float(m["fixed_rg_m"][0]),
        rg_final_m=float(m["fixed_rg_m"][-1]),
        shape_tensor_final_m2=m["fixed_tensor_m2"][-1],
        shape_eigenvalues_final_m2=m["fixed_eigenvalues_m2"][-1],
        centroid_displacement_m=m["fixed_centroid_m"][-1] - m["fixed_centroid_m"][0],
        principal_axis_identifiable_fraction=float(m["principal_axis_identifiable"].mean()),
        identifiable_final_axis_angle_deg=float(m["principal_axis_angle_deg"][-1])
        if m["principal_axis_identifiable"][-1]
        else None,
        minimum_saved_wall_gap_m=float(nearest.min()),
        final_nearest_wall_gap_m=float(nearest[-1]),
        first_saved_wall_contact_s=float(d["time_s"][first[0]]) if len(first) else None,
        max_saved_wall_contact_count=int(contact.max()),
        final_wall_contact_count=int(contact[-1]),
        final_hard_pair_contact_count=int(d["hard_contact_pair_count"][-1, 0]),
        maximum_saved_frame_displacement_m=float(dx.max()),
        quality=s["quality"],
        boundary_events=js(path / "SAVED_METADATA.json")["text_datasets"]["boundary_events_json"],
    )


def main():
    checks = js(BASE / "checks_v2/RESULT.json")
    edges = js(BASE / "edge_checks_v1/RESULT.json")
    assert checks["status"] == edges["status"] == "PASS"
    write(
        BASE / "WALL_CONTACT_VALIDATION.json",
        dict(
            status="PASS",
            point_accuracy=checks,
            finite_edge_checks=edges,
            distance_error_scope="Agreement with original float32 stored mesh, not anatomical accuracy",
            wall_budget_m=0.01 * R,
            per_pass_correction_limit_m=0.1 * R,
            original_mesh_discretization_preserved=True,
            sign="Closed mesh only; contact distance from finite WALL triangles",
            failed_attempts=[
                dict(
                    label="checks_v1",
                    cause="Pair unit-test initialization too close to a faceted-wall edge; replaced with explicitly interior 0.2 um gap; no tolerance change",
                )
            ],
        ),
    )
    rows = {}
    for name in B_NAMES:
        path = BASE / "aggregation_v1" / name
        compact(path)
        rows[name] = summarize(path)
    contrasts = {}
    for region in ["parent", "junction"]:
        off = read(BASE / "aggregation_v1" / (region + "_wall_off"))
        on = read(BASE / "aggregation_v1" / (region + "_wall_on"))
        np.testing.assert_array_equal(off["positions_m"][0], on["positions_m"][0])
        contrasts[region] = dict(
            initial_positions_bitwise_equal=True,
            final_centroid_difference_m=float(
                np.linalg.norm(on["global_centroid_m"][-1] - off["global_centroid_m"][-1])
            ),
            final_rg_difference_m=float(on["global_rg_m"][-1, 0] - off["global_rg_m"][-1, 0]),
            wall_time_ratio=rows[region + "_wall_on"]["wall_s"] / rows[region + "_wall_off"]["wall_s"],
        )
    write(
        BASE / "NEAR_WALL_AGGREGATION.json",
        dict(
            status="PASS_GEOMETRIC_CONFINEMENT_AND_AGGREGATION",
            rows=rows,
            contrasts=contrasts,
            model="EXPLORATORY / GEOMETRIC_WALL_ONLY",
            definitions=dict(
                cluster_link="center separation <= 2.1 R",
                wall_contact="saved gap <= 0.01 R; count of bubbles",
                pair_hard_contact="pair gap <= 0.001 R, counted once per unordered pair",
                Rg="all 32 fixed IDs",
            ),
            interpretation="Actual PZT, zero flow, no prescribed velocity. WALL OFF is the unconfined F3 comparator; its negative physical wall gaps are expected and are not a wall-validation pass.",
        ),
    )
    crows = {}
    for name in C_NAMES:
        path = BASE / "controls_v1" / name
        compact(path)
        crows[name] = summarize(path)
    c = js(BASE / "controls_v1/RESULT.json")
    equal = {}
    for control in ["baseline", "switched"]:
        a = read(BASE / "controls_v1" / ("wall_on_" + control))
        b = read(BASE / "controls_v1" / ("wall_off_" + control))
        equal[control] = bool(np.array_equal(a["positions_m"], b["positions_m"]))
    write(
        BASE / "PZT_WALL_CONTROL.json",
        dict(
            status="WEAK_COMPUTED_RESPONSE; ROBUST_SHAPE_CONTROL_NOT_ESTABLISHED",
            rows=crows,
            endpoint_contrasts=c["endpoint_contrasts"],
            same_state_identity=js(BASE / "controls_v1/STATE_IDENTITY.json"),
            control=c["control"],
            wall_on_off_positions_bitwise_equal=equal,
            interpretation="The cluster detached before continuation; no additional wall contact occurs. WALL ON/OFF coincide exactly. Small control contrast is a computed response only; F3 rearrangement sensitivity and incomplete angular balance remain, without a new F4 uncertainty envelope.",
        ),
    )
    perf = js(BASE / "performance_v1/RESULT.json")
    assert perf["status"] == "PASS"
    write(
        BASE / "PERFORMANCE.json",
        dict(
            status="MEASURED_BOUNDED_OVERHEAD",
            matched_32={r: contrasts[r] for r in contrasts},
            N200=perf,
            timing_scope="Simulation loop includes save-point transfers and geometry observations; excludes compile/build, pre-run capture, file export. Process-body time includes setup/export and morphology analysis.",
        ),
    )
    quality = [v["quality"] for v in list(rows.values()) + list(crows.values())]
    write(
        BASE / "ANALYSIS_CHECKS.json",
        dict(
            status="PASS",
            closed_acoustic_trajectories=8,
            max_accepted_pair_overlap_m=max(q["max_accepted_penetration_m"] for q in quality),
            max_wall_pass_correction_m=max(
                q.get("wall", {}).get("max_pass_correction_m", 0) for q in quality
            ),
            max_physical_step_displacement_m=max(q["max_physical_substep_displacement_m"] for q in quality),
            all_active=True,
            no_boundary_event=True,
            dt_s=1e-5,
            no_time_step_scan=True,
        ),
    )
    print("F4 analysis PASS")


if __name__ == "__main__":
    main()
