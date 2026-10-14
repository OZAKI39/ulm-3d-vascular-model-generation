"""Meaningful CPU geometry, material identity and input rejection checks."""

import copy
import numpy as np
import h5py
from mevo_mb.coupling.pdms_block import apparatus, audit, material_labels
from mevo_mb.coupling.straight_tube import ZeroFlow
from .common import BASE, configs, write, read, sha, identity
from .preflight import reference_material


def rejected(fn):
    try:
        fn()
    except ValueError:
        return True
    return False


def main():
    c = configs()
    g = c["geometry"]
    m = c["materials"]
    lumen, e = apparatus(g)
    a = audit(g)
    checks = {}
    points = np.array(
        [
            [0, 0, 0],
            [0, 0.000799, 0],
            [0, 0.000801, 0],
            [0.016249, 0, 0],
            [0.016251, 0, 0],
            [0, 0.004499, 0],
            [0, 0.004501, 0],
            [0, 0, 0.004501],
        ]
    )
    checks["three_regions_at_interfaces"] = np.array_equal(
        material_labels(points, g), [2, 2, 1, 2, 0, 1, 0, 0]
    )
    checks["centering_and_side_PDMS_3p7mm"] = np.isclose(a["side_PDMS_path_m"], 0.0037) and a["centered"]
    checks["end_disk_covers_entire_opening"] = np.isclose(a["end_opening_coverage_fraction"], 1)
    checks["end_contact_area_excludes_hole"] = np.isclose(
        a["end_face_PDMS_contact_area_m2"], np.pi * (0.0035**2 - 0.0008**2)
    )
    checks["user_confirmed_x_plus_closed"] = lumen.ends == ("open_fail_stop", "closed")
    checks["correct_new_length_in_contains"] = (
        lumen.contains(np.array([[0.015, 0, 0]]), 1.78e-6)[0]
        and not lumen.contains(np.array([[0.016251, 0, 0]]), 1.78e-6)[0]
    )
    checks["finite_radius_wall"] = not lumen.contains(np.array([[0, 0.0008 - 1e-6, 0]]), 1.78e-6)[0]
    gg = copy.deepcopy(g)
    gg["lumen_center_local_m"] = [0, 0.0001, 0]
    checks["off_center_rejected"] = rejected(lambda: apparatus(gg))
    gg = copy.deepcopy(g)
    gg["outer_diameter_m"] = 0.0018
    checks["fictitious_tube_outer_diameter_rejected"] = rejected(lambda: apparatus(gg))
    gg = copy.deepcopy(g)
    gg["pzt"][0]["normal_local"] = [0, -1, 0]
    checks["outward_normal_rejected"] = rejected(lambda: apparatus(gg))
    gg = copy.deepcopy(g)
    gg["pzt"].append(copy.deepcopy(gg["pzt"][0]))
    checks["fourth_source_rejected"] = rejected(lambda: apparatus(gg))
    gg = copy.deepcopy(g)
    gg["pzt"][0]["active_diameter_m"] = 0.009
    checks["old_9mm_diameter_rejected"] = rejected(lambda: apparatus(gg))
    gg = copy.deepcopy(g)
    gg["rotation_local_to_world"] = [[0, -1, 0], [1, 0, 0], [0, 0, 1]]
    gg["center_m"] = [0.01, 0.02, 0.03]
    ll, _ = apparatus(gg)
    checks["rotation_translation_material_equivariance"] = np.array_equal(
        material_labels(ll.world(points), gg), material_labels(points, g)
    )
    checks["PDMS_2p5MHz_not_silently_extrapolated"] = rejected(lambda: reference_material(m, 2.5e6))
    mm = copy.deepcopy(m)
    mm["shared_alpha_power"] = 2
    checks["incompatible_material_power_rejected"] = rejected(lambda: reference_material(mm, 5e6))
    ref = reference_material(m, 5e6)
    roundtrip = ref["alpha_coeff_dB_MHz_power_cm"] * 5 ** ref["alpha_power"] * 100 * np.log(10) / 20
    checks["Np_m_to_dB_MHz_power_cm_roundtrip"] = np.isclose(roundtrip, 214.8, rtol=1e-12)
    checks["phase_speed_reference_exact"] = ref["phase_speed_m_s"] == 1029
    velocity, grad = ZeroFlow().query(np.repeat(points, 512, axis=0))
    checks["4096_CPU_zero_flow"] = np.count_nonzero(velocity) == np.count_nonzero(grad) == 0
    info = read(BASE / "materials/IDENTITY.json")
    path = BASE / "materials/material_maps.h5"
    checks["material_map_sha_identity"] = sha(path) == info["sha256"]
    with h5py.File(path) as h:
        axes = [h[v + "_m"][:] for v in "xyz"]
        ids = h["material_id"]
        index = tuple(np.argmin(abs(v)) for v in axes)
        checks["material_map_lumen_at_origin"] = ids[index] == 2
        checks["three_regions_in_saved_map"] = set(np.unique(ids[index[0]])) == {0, 1, 2}
        for field in ["sound_speed_map", "density_map", "alpha_coeff_map"]:
            checks[field + "_shape_matches"] = h[field].shape == tuple(info["grid_shape"])
    pp = np.array(
        [
            [0.01085, 0, 0.00449],
            [0.01085, 0.000749, 0.003],
            [0.01085, 0.000751, 0.003],
            [0.01085, 0, 0.0003],
            [0.01085, 0, 0.00451],
        ]
    )
    checks["top_port_material_union_and_open_surface"] = np.array_equal(
        material_labels(pp, g), [2, 2, 1, 2, 0]
    )
    checks["top_port_diameter_user_confirmed"] = g["top_port"]["diameter_m"] == 0.0015
    checks["top_port_x_derived_from_5p4mm"] = np.isclose(0.01625 - g["top_port"]["center_xy_m"][0], 0.0054)
    with h5py.File(path) as h:
        ix = np.argmin(abs(h["x_m"][:] - g["top_port"]["center_xy_m"][0]))
        iy = np.argmin(abs(h["y_m"][:]))
        iz = np.argmin(abs(h["z_m"][:] - 0.003))
        checks["top_port_present_in_saved_material_map"] = h["material_id"][ix, iy, iz] == 2
    write(
        BASE / "analysis/CPU_CHECKS.json",
        dict(status="PASS" if all(checks.values()) else "FAIL", checks=checks, identity=identity()),
    )
    assert all(checks.values()), checks
    print("PASS", len(checks), "CPU checks")


if __name__ == "__main__":
    main()
