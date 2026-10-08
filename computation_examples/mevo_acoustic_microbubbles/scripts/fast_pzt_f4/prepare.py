"""Read original mesh and characterize field intersection before GPU development."""

import numpy as np
import h5py
from mevo_mb.fast_sim.vessel_geometry import VesselGeometry, initialize_particles_in_vessel
from .common import ROOT, BASE, js, write


def main():
    g = VesselGeometry(BASE / "inputs/labelled_lumen_m.vtp")
    planes = js(BASE / "inputs/PLANES.json")
    ident = g.prepare(BASE / "inputs/prepared_geometry.npz", planes)
    write(BASE / "GEOMETRY_IDENTITY.json", ident)
    with h5py.File(BASE / "inputs/basis.h5") as h:
        axes = [h["axes/" + n + "_m"][:] for n in ["x", "y", "z"]]
    bounds = np.array([[a.min(), a.max()] for a in axes])
    margin = 1.78e-6 + 2e-8
    valid = bounds + np.array([margin, -margin])
    regions = {}
    for name, target in [("parent", [-0.002, 0, -0.00129]), ("junction", [0.001, 0.0003, -0.001])]:
        x, info = initialize_particles_in_vessel(g, valid, target)
        np.savez_compressed(
            BASE / "inputs" / (name + "_initial.npz"),
            positions_m=x,
            particle_id=np.arange(len(x)),
            radii_m=np.full(len(x), 1.78e-6),
        )
        regions[name] = info
    coverage = {}
    for name, end in [
        ("parent", planes["INLET"]["origin_m"]),
        ("daughter_A", planes["OUTLET_A"]["origin_m"]),
        ("daughter_B", planes["OUTLET_B"]["origin_m"]),
    ]:
        t = np.linspace(0, 1, 2001)
        x = t[:, None] * np.array(end)
        d, _, _, inside = g.query(x)
        keep = inside & (d >= 1.78e-6) & (x >= valid[:, 0]).all(1) & (x <= valid[:, 1]).all(1)
        coverage[name] = dict(
            centerline_sample_count=int(keep.sum()),
            sample_spacing_m=float(np.linalg.norm(end) / 2000),
            centerline_covered_distance_from_junction_m=[
                float(t[keep].min() * np.linalg.norm(end)),
                float(t[keep].max() * np.linalg.norm(end)),
            ],
            full_branch_covered=bool(keep.all()),
            note="Centerline coverage only; full supported lumen is the mesh/field intersection",
        )
    write(
        BASE / "GEOMETRY_FIELD_OVERLAP.json",
        dict(
            vessel_bounds_m=ident["bounds_m"],
            acoustic_bounds_m=bounds,
            usable_center_bounds_m=valid,
            center_margin_m=margin,
            radius_m=1.78e-6,
            valid_region="Closed original lumen interior AND distance to actual WALL >= R AND center inside field bounds inset by R+20 nm",
            bounding_box_is_not_lumen=True,
            field_boundary="FIELD_DOMAIN_EXIT; not a physical wall; never reflect or extrapolate",
            terminals=planes,
            coverage=coverage,
            initialization=regions,
            whole_vessel_covered=False,
            wall_normal_convention="Outward; closed-lumen interior has positive signed WALL distance",
        ),
    )
    config = dict(
        dt_s=1e-5,
        physical_substeps=1,
        mobility_iterations=2,
        contact_iterations=12,
        wall_pass_after_pair_iterations=[3, 6, 9, 12],
        wall_epsilon_over_radius=0.01,
        pair_epsilon_over_radius=0.01,
        contact_pass_limit_over_radius=0.1,
        wall_path_safety="Far-wall Lipschitz bound; near-wall finite triangle swept distance via static GPU BVH, double precision refinement; reject unresolved path",
        execution_mode="stable_reduction",
        flow=dict(enabled=False),
        vessel=dict(
            enabled=True,
            geometry_source="original_stage4_y",
            boundary_mode="wall_only_open_terminals",
            wall_contact_enabled=True,
        ),
        physics_scope="EXPLORATORY",
        wall_physics="GEOMETRIC_WALL_ONLY",
        aggregation_duration_s=0.12,
        continuation_duration_s=0.12,
        sample_interval_s=0.0002,
        wall_budget_s=180,
        fixed_count=32,
        control_selection="F3 baseline and amplitude/enable; no post-trajectory tuning",
        uncertainty="Inherited F3 field-gradient, local mobility and initial-condition sensitivity remain",
    )
    write(ROOT / "configs/fast_pzt/f4/protocol.json", config)
    print(
        {
            n: dict(
                gap_um=v["minimum_wall_gap_m"] * 1e6,
                center=np.mean(
                    np.load(BASE / "inputs" / (n + "_initial.npz"))["positions_m"], axis=0
                ).tolist(),
            )
            for n, v in regions.items()
        }
    )


if __name__ == "__main__":
    main()
