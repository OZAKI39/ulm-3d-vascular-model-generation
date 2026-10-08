"""Static preflight and preregistration, before any F5 scientific trajectory."""

import numpy as np
from mevo_mb.flow.fields import SimVascularFlowField
from mevo_mb.fast_sim.vessel_geometry import VesselGeometry, ORIGINAL_SHA256
from .common import BASE, INPUT, ROOT, write, js, sha, F4_ARCHIVE


def main():
    source = INPUT / "result_121.vtu"
    assert sha(source) == "131f7cb9b3e4ac9887bb0aa344d0eb2c0a2c8ac1d3b19fafb093ecd9a1045620"
    assert sha(INPUT / "f4/labelled_lumen_m.vtp") == ORIGINAL_SHA256
    state = js(INPUT / "CFD_STATE.json")
    assert state["status"] == "PASS_WITH_SCOPE" and state["snapshot_sha256"] == sha(source)
    assert state["coordinate_scale_to_m"] == state["velocity_scale_to_m_s"] == 1
    field = SimVascularFlowField(source, spatial_method="tetra_p1")
    geometry = VesselGeometry(INPUT / "f4/labelled_lumen_m.vtp")
    x = np.load(INPUT / "f4/parent_initial.npz")["positions_m"]
    u = field.velocity(x, 0)
    dist, _, _, inside = geometry.query(x)
    assert inside.all() and (dist >= 1.78e-6).all()
    bounds = np.array([[-0.0101982, 0.0101982], [-0.0017982, 0.0017982], [-0.0017982, 0.0017982]])
    assert (x >= bounds[:, 0]).all() and (x <= bounds[:, 1]).all()
    # Straight, frozen-velocity preview is only a domain feasibility indicator.
    # It is not used to advance, tune, or select a scientific trajectory.
    preview = x + 0.3 * u
    field.velocity(preview, 0)
    assert (preview >= bounds[:, 0]).all() and (preview <= bounds[:, 1]).all()
    points = np.asarray(field.grid.points, float)
    cells = field.sampler.cells
    centers = points[cells].mean(1)
    groups = {}
    masks = dict(
        parent=(centers[:, 0] < -0.003) & (np.linalg.norm(centers[:, 1:], axis=1) < 0.0008),
        bifurcation=(abs(centers[:, 0]) < 0.001) & (np.linalg.norm(centers[:, 1:], axis=1) < 0.0008),
        daughter=(centers[:, 0] > 0.004) & (abs(centers[:, 1]) > 0.002),
    )
    for label, mask in masks.items():
        candidates = np.flatnonzero(mask)
        groups[label] = centers[candidates[np.linspace(0, len(candidates) - 1, 12).astype(int)]]
    groups["near_wall"] = x[::4]
    np.savez_compressed(INPUT / "flow_validation_points.npz", **groups)
    controls = js(ROOT / "configs/fast_pzt/f4/controls.json")
    # Use archived actual parent controls, unchanged, rather than search controls.
    event = controls["regions"]["parent"]["baseline"]["event"]
    event = {k: event[k] for k in ["amplitude_Pa", "phase_rad", "enabled"]}
    on = dict(event, time_s=0.0, flow_enabled=True, secondary_enabled=True, mode="ACOUSTIC_ON_FLOW_ON")
    off = dict(
        event,
        time_s=0.3,
        enabled=[False] * 4,
        flow_enabled=True,
        secondary_enabled=False,
        mode="ULTRASOUND_OFF_FLOW_ON",
    )
    protocol = dict(
        base_f4=F4_ARCHIVE,
        positions_m=x,
        initial_source="Unmodified archived F4 parent cloud",
        duration_s=0.6,
        dt_s=1e-5,
        sample_interval_s=0.0005,
        wall_budget_s=240,
        events=[on, off],
        cases=[
            "main_ON_then_ULTRASOUND_OFF_FLOW_ON",
            "branch_A_ULTRASOUND_OFF_FLOW_ON",
            "branch_B_ULTRASOUND_OFF_FLOW_OFF",
        ],
        branch_A="Exact suffix of main run, saved from the 300 ms switch state; no duplicate rerun",
        optional_C="Not preselected; only add if A/B leave an unresolved direct acoustic question",
        preflight=dict(
            valid_initial_count=len(x),
            minimum_wall_gap_m=float(dist.min() - 1.78e-6),
            velocity_min_m_s=u.min(0),
            velocity_max_m_s=u.max(0),
            maximum_pair_velocity_difference_m_s=float(
                np.linalg.norm(u[:, None] - u[None, :], axis=-1).max()
            ),
            frozen_velocity_preview_inside_flow_and_acoustic=True,
        ),
        domain_policy=dict(
            acoustic_on="F4 FIELD_DOMAIN_EXIT, ID and crossing time retained",
            ultrasound_off="No acoustic lookup or field-domain exit",
            flow="Invalid tetra membership fails the run; no extrapolation; physical cap exits remain legal",
        ),
        physics_scope="EXPLORATORY / PRESCRIBED_STEADY_FLOW / GEOMETRIC_WALL_ONLY",
    )
    write(BASE / "ACOUSTIC_RELEASE_PROTOCOL.json", protocol)
    write(ROOT / "configs/fast_pzt/f5/protocol.json", protocol)
    write(
        BASE / "FLOW_SOURCE_IDENTITY.json",
        dict(
            source_server="/workspace/mevo_acoustic_microbubbles/results/stage4/simvascular/steady/20261002T141810826990Z/run/8-procs/result_121.vtu",
            source_local=str(source),
            sha256=sha(source),
            original_vessel_sha256=ORIGINAL_SHA256,
            state=state,
            CPU_import_validation=js(INPUT / "SIMVASCULAR_FLOW_IMPORT_VALIDATION.json"),
            node_count=field.grid.n_points,
            tetra_count=field.grid.n_cells,
            velocity_array="Velocity",
            units="m; m/s",
            velocity_range_m_s=[
                float(np.linalg.norm(field.grid["Velocity"], axis=1).min()),
                float(np.linalg.norm(field.grid["Velocity"], axis=1).max()),
            ],
            CFD_rerun=False,
            fluid_scope="Archived steady water-parameter numerical CFD; rho=1000 kg/m3, mu=0.001 Pa s; no mesh-independence or experimental validation",
        ),
    )
    print("Preflight and immutable scene registration complete")


if __name__ == "__main__":
    main()
