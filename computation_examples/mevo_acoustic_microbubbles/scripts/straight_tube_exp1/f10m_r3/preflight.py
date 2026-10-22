"""CPU geometry, topology, exposed-wall and unresolved-medium admission checks."""

import time
import shutil
import copy
import numpy as np
import h5py
from scipy.ndimage import label
from scipy.spatial.distance import pdist
from mevo_mb.coupling.pdms_cavity import apparatus, PDMSCavity
from mevo_mb.coupling.straight_tube import ZeroFlow, disc_points
from mevo_mb.control.pzt import PZTArray
from mevo_mb.acoustics.spatial_medium import medium_arrays
from .common import BASE, ROOT, configs, read, write, sha, identity


def main():
    start = time.perf_counter()
    c = configs()
    out = BASE / "analysis"
    out.mkdir(parents=True, exist_ok=True)
    assert not (BASE / "CONFIG_SNAPSHOT.json").exists(), "Keep attempts; do not overwrite recorded inputs"
    cavity, elements = apparatus(c["geometry"])
    b = c["protocol"]["bubbles"]
    r = b["radius_m"]
    checks = {}
    checks["confirmed_dimensions_and_four_centered_faces"] = (
        len(elements) == 4 and cavity.R == 0.0008 and cavity.L == 0.01625
    )
    checks["both_X_sealed"] = cavity.lumen.ends == ("closed", "closed")
    checks["no_bottom_perforation"] = all(cavity.labels([x, 0, -0.0045]) == 1 for x in cavity.ports)
    checks["top_void_open_and_mirrored"] = (
        all(cavity.labels([x, 0, 0.0045]) in (3, 4) for x in cavity.ports) and sum(cavity.ports) == 0
    )
    checks["no_X_void_beyond_seals"] = all(cavity.implicit([s * 0.016251, 0, 0]) < 0 for s in [-1, 1])
    checks["sealed_end_surface_tags"] = np.array_equal(
        cavity.boundary_labels([[-cavity.L, 0, 0], [cavity.L, 0, 0]]), [1, 2]
    )
    checks["open_top_surface_tags"] = np.array_equal(
        cavity.boundary_labels([[x, 0, cavity.H] for x in cavity.ports]), [3, 4]
    )
    bad = copy.deepcopy(c["geometry"])
    bad["end_conditions"][0] = "open_fail_stop"
    try:
        PDMSCavity(bad)
        checks["reject_legacy_open_X_config"] = False
    except ValueError:
        checks["reject_legacy_open_X_config"] = True

    # Dense geometric topology grid, unrelated to any acoustic PPW claim.
    axes = [
        np.linspace(-cavity.L, cavity.L, 651),
        np.linspace(-0.00085, 0.00085, 35),
        np.linspace(-0.00085, 0.00455, 109),
    ]
    pts = np.stack(np.meshgrid(*axes, indexing="ij"), -1)
    labels = cavity.labels(pts)
    void = np.isin(labels, [2, 3, 4])
    components, n = label(void)
    ids = []
    for point in [[0, 0, 0], [-0.01085, 0, 0.0045], [0.01085, 0, 0.0045]]:
        ix = tuple(np.argmin(abs(a - v)) for a, v in zip(axes, point))
        ids.append(int(components[ix]))
    checks["true_union_one_six_neighbor_component"] = n == 1 and ids == [1, 1, 1]
    topidx = np.argmin(abs(axes[2] - cavity.H))
    checks["both_open_mouths_reach_exterior"] = all(
        void[np.argmin(abs(axes[0] - x)), 17, topidx]
        and not void[np.argmin(abs(axes[0] - x)), 17, topidx + 1]
        for x in cavity.ports
    )
    write(
        out / "CONNECTIVITY.json",
        dict(
            status="PASS" if checks["true_union_one_six_neighbor_component"] else "FAIL",
            components=n,
            main_and_ports_component_ids=ids,
            grid_spacing_m=5e-5,
            shape=void.shape,
            void_nodes=int(void.sum()),
            method="6-neighbor voxel flood fill of actual shared membership function; analytic intersections checked separately",
            phase_scope="Void connectivity, not evidence of liquid filling",
            X_seals="explicit surface tags; no void outside X planes",
        ),
    )
    del pts, labels, void, components

    gaps = []
    for e in elements:
        center = np.array(e["position_m"]) * 1e3
        for port in cavity.ports * 1e3:
            # Each axis-aligned distance lower bound has an attainable witness.
            if e["id"].startswith("side"):
                face = np.array([np.sign(port) * 3.5, center[1], 0.0])
                delta = face[:2] - [port, 0.0]
                hole = np.r_[np.array([port, 0.0]) + 0.75 * delta / np.linalg.norm(delta), 0.0]
                distance = np.linalg.norm(delta) - 0.75
            else:
                face = center.copy()
                hole = np.array([port + np.sign(center[0] - port) * 0.75, 0.0, 0.0])
                distance = abs(center[0] - port) - 0.75
            assert np.linalg.norm(face - center) <= 3.5 + 1e-12
            assert np.isclose(np.linalg.norm(face - hole), distance, atol=1e-12)
            gaps.append(
                dict(
                    pzt=e["id"],
                    port_x_mm=port,
                    minimum_face_to_hole_volume_gap_mm=float(distance),
                    face_witness_mm=face,
                    hole_witness_mm=hole,
                    method="Attained analytic lower bound for axis-aligned face / vertical cylinder",
                )
            )
    checks["no_PZT_face_hole_intersection"] = min(g["minimum_face_to_hole_volume_gap_mm"] for g in gaps) > 0
    checks["adjacent_end_face_hole_gap_4p65mm"] = np.isclose(
        min(g["minimum_face_to_hole_volume_gap_mm"] for g in gaps), 4.65, atol=1e-7
    )
    array = PZTArray.from_config({"pzt": elements, "acoustics": {"frequency_Hz": 3e6}})
    array.set_amplitudes([100, 200, 300, 400])
    array.set_phases([0, np.pi / 2, np.pi, 3 * np.pi / 2])
    array.set_enabled([True, True, False, True])
    checks["exactly_4_independent_source_controls"] = np.allclose(
        array.get_control_vector(), [100, 200j, 0, -400j]
    )
    try:
        array.set_amplitudes([1, 2, 3])
        checks["reject_legacy_three_controls"] = False
    except ValueError:
        checks["reject_legacy_three_controls"] = True

    # Independent-radius tests include the true re-entrant intersection curve.
    path = np.array([[x, 0, z] for x in cavity.ports for z in np.linspace(0.0006, 0.0012, 13)])
    checks["sphere_path_through_both_junctions_without_false_main_wall"] = all(
        cavity.query(q, r)["valid"] for q in path
    )
    probes = [[0, 0, 0], [0, 0.0007, 0], [-0.01625 + r, 0, 0], [0.01625 - r, 0, 0]]
    for x in cavity.ports:
        probes += [[x, 0, 0.0008], [x, 0, 0.001], [x + 0.00074, 0, 0.001], [x, 0.0004, 0.00068]]
        for t in np.linspace(0.15, 6.1, 13):
            edge = cavity.curve(t, x)
            inward = -np.array(
                [(edge[0] - x) / cavity.a, edge[1] / cavity.a + edge[1] / cavity.R, edge[2] / cavity.R]
            )
            inward /= np.linalg.norm(inward)
            probes.extend([edge + inward * r * 0.95, edge + inward * r * 2, edge - inward * r * 0.05])
    probes = np.array(probes)
    wall = [cavity.query(q, r) for q in probes]
    # Sphere-surface sampling independently tests union membership (not distance).
    i = np.arange(1024)
    z = 1 - 2 * (i + 0.5) / len(i)
    phi = i * np.pi * (3 - np.sqrt(5))
    sphere = np.stack([np.sqrt(1 - z * z) * np.cos(phi), np.sqrt(1 - z * z) * np.sin(phi), z], -1)
    checked = 0
    for q, row in zip(probes, wall):
        if row["valid"]:
            assert cavity.implicit(q + r * sphere).min() >= -1e-12
            checked += 1
    checks["finite_spheres_checked_against_1024_surface_directions"] = checked > 10
    projection = []
    test_points = [
        [-cavity.L + r - 0.05 * r, 0, 0],
        [cavity.L - r + 0.05 * r, 0, 0],
        [0, cavity.R - r + 0.05 * r, 0],
    ]
    test_points += [[x + cavity.a - r + 0.05 * r, 0, 0.001] for x in cavity.ports]
    test_points += [probes[8]]
    for q in test_points:
        corrected, status, cost = cavity.project(q, r)
        projection.append(dict(input=q, output=corrected, status=status, correction_m=cost))
    checks["finite_radius_projection_both_ends_main_and_port_walls"] = all(
        x["status"] == "PASS" for x in projection
    )
    rejected = cavity.project([cavity.L - r + 0.2 * r, 0, 0], r)
    checks["original_0p1R_correction_budget_retained"] = rejected[1] == "WALL_BUDGET_REJECTED"
    ceiling_point = np.array([cavity.ports[0], 0, cavity.ceiling - r / 2])
    ceiling_result = cavity.project(ceiling_point, r)
    checks["unknown_surface_scope_rejected_without_wall_projection"] = ceiling_result[
        1
    ] == "UNCONFIRMED_FREE_SURFACE_SCOPE" and np.array_equal(ceiling_result[0], ceiling_point)
    flow = ZeroFlow().query(probes)
    checks["zero_velocity_and_gradient"] = np.count_nonzero(flow[0]) == np.count_nonzero(flow[1]) == 0
    np.savez_compressed(
        BASE / "wall_reference.npz",
        positions_m=probes,
        signed_distance_m=[q["wall_distance_m"] for q in wall],
        normals=[q["inward_normal"] for q in wall],
        closest_m=[q["closest_wall_m"] for q in wall],
        projection_input=np.array(test_points),
        projection_output=np.array([p["output"] for p in projection]),
    )

    initial = np.array(
        [
            [0, 0, 0],
            [0, 20e-6, 0],
            [0, -20e-6, 0],
            [0, 0, 20e-6],
            [-cavity.L + r + 1e-9, 0, 0],
            [cavity.L - r - 1e-9, 0, 0],
            [0, cavity.R - r - 1e-9, 0],
            [0, -cavity.R + r + 1e-9, 0],
            *[[x, 0, z] for x in cavity.ports for z in [0.0006, 0.0008, 0.0009, 0.0012]],
        ]
    )
    assert len(initial) == b["count"] and pdist(initial).min() > 2 * r
    checks["initial_spheres_legal_and_nonoverlapping"] = all(cavity.query(q, r)["valid"] for q in initial)
    np.savez_compressed(
        BASE / "initial_state.npz", positions_m=initial, particle_id=np.arange(len(initial)), radius_m=r
    )
    (BASE / "inputs").mkdir(exist_ok=True)
    shutil.copy2(ROOT / "results/straight_tube_exp1/f10m/inputs/response.h5", BASE / "inputs/response.h5")
    assert sha(BASE / "inputs/response.h5") == b["response_sha256"]

    # Labels and unknown-phase NaNs on a diagnostic grid, never a solve grid.
    ma = [
        np.linspace(-0.01725, 0.01725, 346),
        np.linspace(-0.0055, 0.0055, 111),
        np.linspace(-0.0055, 0.0055, 111),
    ]
    reference = read(ROOT / "results/straight_tube_exp1/f10m/materials/REFERENCE.json")
    material = BASE / "materials/regions_and_unresolved_media.h5"
    material.parent.mkdir(exist_ok=True)
    with h5py.File(material, "w") as h:
        shape = tuple(len(a) for a in ma)
        for key, axis in zip("xyz", ma):
            h.create_dataset(key + "_m", data=axis)
        ds = {
            name: h.create_dataset(name, shape=shape, dtype=dt, compression="gzip", shuffle=True)
            for name, dt in [
                ("region_id", "u1"),
                ("boundary_id", "u1"),
                ("sound_speed_map", "f4"),
                ("density_map", "f4"),
                ("alpha_coeff_map", "f4"),
            ]
        }
        y, z = np.meshgrid(ma[1], ma[2], indexing="ij")
        for j, x in enumerate(ma[0]):
            q = np.stack([np.full_like(y, x), y, z], -1)
            reg = cavity.labels(q)
            bound = cavity.boundary_labels(q)
            ds["region_id"][j] = reg
            ds["boundary_id"][j] = bound
            for name, pdms_val, water_val in [
                ("sound_speed_map", reference["zero_frequency_sound_speed_m_s"], 1500),
                ("density_map", reference["density_kg_m3"], 1000),
                ("alpha_coeff_map", reference["alpha_coeff_dB_MHz_power_cm"], 0),
            ]:
                value = np.where(reg == 1, pdms_val, np.where(reg == 2, water_val, np.nan))
                value[np.isin(bound, [1, 2])] = np.nan
                ds[name][j] = value
        h.attrs["status"] = "INCOMPLETE_NOT_SOLVER_INPUT; ports/exterior/seal acoustic parameters unknown"
        h.attrs["region_ids"] = (
            "0 exterior unknown;1 PDMS reference;2 main cavity liquid DEMO;3/4 top void PHASE_UNKNOWN"
        )
        h.attrs["boundary_ids"] = (
            "1/2 X mechanical closed (material/thickness unknown);3/4 top geometrically open"
        )
        h.attrs["scope"] = (
            "Diagnostic geometry grid 100 um; not an acoustic discretization or water-filled-port claim"
        )
    try:
        medium_arrays(
            dict(maps_h5=str(material), maps_sha256=sha(material), alpha_power=reference["alpha_power"]),
            shape,
        )
        checks["solver_rejects_unresolved_media"] = False
    except ValueError:
        checks["solver_rejects_unresolved_media"] = True
    with h5py.File(material) as h:
        ix = np.argmin(abs(ma[0] + 0.01085))
        iy = np.argmin(abs(ma[1]))
        iz = np.argmin(abs(ma[2] - 0.003))
        checks["port_not_silently_water"] = h["region_id"][ix, iy, iz] == 3 and np.isnan(
            h["sound_speed_map"][ix, iy, iz]
        )

    grid = read(ROOT / "results/straight_tube_exp1/f10m/analysis/RESOURCE_PREFLIGHT.json")
    grid["estimated_four_basis_seconds"] = grid.pop("estimated_three_basis_seconds") * 4 / 3
    grid["scope"] = (
        "CONDITIONAL_LITERATURE_WATER_ENVELOPE_ONLY; actual minimum speed unknown until phases/seals known"
    )
    grid["source_count"] = 4
    grid["actual_apparatus_weights_generated"] = False
    grid["estimate_method"] += (
        " Four bases scale old three-basis estimate by 4/3; same source envelope. Air/solid seals may invalidate this candidate grid."
    )
    rim = np.vstack([disc_points(e) for e in elements])
    bounds = np.array(grid["axes_bounds_m"])
    checks["four_rims_inside_candidate_grid"] = bool(((rim >= bounds[:, 0]) & (rim <= bounds[:, 1])).all())
    reasons = [
        "TOP_FILL_HEIGHT_PHASE_AND_EXTERIOR_UNKNOWN",
        "X_SEAL_MATERIAL_THICKNESS_AND_COUPLING_UNKNOWN",
        "MEASUREMENT_REFERENCE_PROVISIONAL",
        "EXPERIMENTAL_FREQUENCY_AND_VOLTAGE_PRESSURE_TRANSFER_UNKNOWN",
        "NO_FREQUENCY_MATCHED_APECSS_FOR_CANDIDATE_ON",
        "CONDITIONAL_MEMORY_ESTIMATE_EXCEEDS_24_GiB",
    ]
    write(
        out / "ACOUSTIC_ADMISSION.json",
        dict(
            status="NOT_RUN",
            actual_basis_count=0,
            reasons=reasons,
            old_three_fields_reused=False,
            unknowns_never_rigid_or_water=True,
            training_ready=False,
        ),
    )
    write(out / "RESOURCE_PREFLIGHT.json", grid)
    write(
        BASE / "materials/IDENTITY.json",
        dict(
            sha256=sha(material),
            status="UNRESOLVED_NOT_SOLVER_INPUT",
            reference_sha256=sha(ROOT / "results/straight_tube_exp1/f10m/materials/REFERENCE.json"),
        ),
    )
    write(
        out / "GEOMETRY.json",
        dict(
            status="PASS",
            block_dimensions_m=c["geometry"]["block_dimensions_m"],
            main_radius_m=cavity.R,
            main_length_m=2 * cavity.L,
            side_path_m=0.0045 - cavity.R,
            top_ports=c["geometry"]["top_ports"],
            ends=["closed", "closed"],
            pzt=elements,
            face_hole_gaps=gaps,
            end_opening_coverage_fraction=[1.0, 1.0],
            end_seal_material=None,
            end_seal_acoustics="UNCONFIRMED",
            side_source_PDMS_gap_m=0,
            source_count=4,
        ),
    )
    write(BASE / "CONFIG_SNAPSHOT.json", c)
    write(BASE / "INPUT_IDENTITY.json", identity())
    write(
        out / "CPU_CHECKS.json",
        dict(
            status="PASS" if all(checks.values()) else "FAIL",
            checks=checks,
            wall_reference_points=len(probes),
            projection=projection,
            wall_s=time.perf_counter() - start,
            identity=identity(),
        ),
    )
    print({"checks": checks, "wall_s": time.perf_counter() - start}, flush=True)
    assert all(checks.values())


if __name__ == "__main__":
    main()
