"""Measured centered channel in a PDMS block; no fictitious thin tube wall."""

from dataclasses import dataclass
import numpy as np
from .straight_tube import StraightTube


@dataclass
class CylindricalLumen:
    radius_m: float
    length_m: float
    center: np.ndarray
    rotation: np.ndarray
    ends: tuple
    world = StraightTube.world
    local = StraightTube.local
    gaps = StraightTube.gaps
    contains = StraightTube.contains

    @classmethod
    def from_config(cls, g):
        if any(k in g for k in ("wall_thickness_m", "outer_diameter_m")):
            raise ValueError("PDMS cast channel has no independent tube wall or outer diameter")
        dims = np.asarray(g["block_dimensions_m"], float)
        r = float(g["lumen_radius_m"])
        if not np.array_equal(dims, [0.0325, 0.009, 0.009]) or r != 0.0008:
            raise ValueError("F10M measured block and lumen dimensions must be retained")
        if np.any(g["lumen_center_local_m"]):
            raise ValueError("Confirmed channel must be centered")
        rot = np.asarray(g["rotation_local_to_world"], float)
        if not np.allclose(rot.T @ rot, np.eye(3), atol=1e-12) or not np.isclose(np.linalg.det(rot), 1):
            raise ValueError("Invalid rotation")
        ends = tuple(g["end_conditions"])
        if len(ends) != 2 or any(e not in ("open_fail_stop", "closed") for e in ends):
            raise ValueError("Explicit port boundary required")
        return cls(r, dims[0], np.asarray(g["center_m"], float), rot, ends)


def apparatus(g):
    lumen = CylindricalLumen.from_config(g)
    elements = []
    wanted = {"side_y_minus": ([0, -0.0045, 0], [0, 1, 0]), "side_y_plus": ([0, 0.0045, 0], [0, -1, 0])}
    if [e["id"] for e in g["pzt"]] != ["side_y_minus", "side_y_plus", "end_x_plus"]:
        raise ValueError("Exactly two Y sources and one X end source, no Z sources")
    for e in g["pzt"]:
        c = np.asarray(e["center_local_m"], float)
        n = np.asarray(e["normal_local"], float)
        if (
            e["active_diameter_m"] != 0.007
            or not np.isfinite(c).all()
            or not np.isclose(np.linalg.norm(n), 1)
        ):
            raise ValueError("7 mm diameter and finite unit normal required")
        if e["id"] in wanted:
            cc, nn = wanted[e["id"]]
            if not np.allclose(c, cc, rtol=0, atol=1e-12) or not np.array_equal(n, nn):
                raise ValueError("Side sources must contact confirmed Y faces and face inward")
        elif not np.array_equal(n, [-1, 0, 0]) or c[0] < 0.01625 - 1e-12:
            raise ValueError("End source must face inward and not invade PDMS")
        elements.append(
            dict(
                id=e["id"],
                position_m=lumen.world(c).tolist(),
                normal=(n @ lumen.rotation.T).tolist(),
                geometry="disc",
                size_m=[0.007],
                amplitude_Pa=0.0,
                phase_rad=0.0,
                enabled=True,
            )
        )
    return lumen, elements


def circle_overlap(r, R, d):
    if d >= r + R:
        return 0.0
    if d <= abs(R - r):
        return np.pi * min(r, R) ** 2
    a = np.arccos((d * d + r * r - R * R) / (2 * d * r))
    b = np.arccos((d * d + R * R - r * r) / (2 * d * R))
    return r * r * a + R * R * b - 0.5 * np.sqrt((-d + r + R) * (d + r - R) * (d - r + R) * (d + r + R))


def audit(g):
    lumen, elements = apparatus(g)
    end = g["pzt"][2]
    c = np.asarray(end["center_local_m"])
    projected = circle_overlap(lumen.radius_m, 0.0035, np.linalg.norm(c[1:]))
    touching = np.isclose(c[0], lumen.length_m / 2, atol=1e-12, rtol=0)
    area = projected if touching else 0.0
    contact = np.pi * 0.0035**2 - area
    return dict(
        status="PASS",
        block_dimensions_m=g["block_dimensions_m"],
        lumen_radius_m=lumen.radius_m,
        lumen_length_m=lumen.length_m,
        centered=True,
        independent_tube_wall=False,
        side_PDMS_path_m=0.0045 - lumen.radius_m,
        source_count=3,
        pzt=elements,
        no_z_sources=True,
        end_conditions=lumen.ends,
        end_mounting_status=g["end_mounting_status"],
        end_disc_projected_opening_overlap_m2=projected,
        end_disc_contact_opening_overlap_m2=area,
        end_opening_coverage_fraction=area / (np.pi * lumen.radius_m**2),
        end_face_PDMS_contact_area_m2=contact,
        end_mounting_confirmed=g["end_mounting_status"] == "CONFIRMED",
        port_boundary_changed_automatically=False,
        side_gap_to_PDMS_m=0.0,
        opening_area_m2=np.pi * lumen.radius_m**2,
    )


def material_labels(points, g):
    """0 exterior; 1 PDMS; 2 channel liquid. Grid-node interfaces are explicit."""
    lumen = CylindricalLumen.from_config(g)
    q = lumen.local(points)
    block = (abs(q) <= np.asarray(g["block_dimensions_m"]) / 2 + 1e-14).all(axis=-1)
    radial = np.linalg.norm(q[..., 1:], axis=-1) <= lumen.radius_m + 1e-14
    port = np.zeros(radial.shape, bool)
    if "top_port" in g:
        top_port_audit(g)
        p = g["top_port"]
        x, y = p["center_xy_m"]
        r = p["diameter_m"] / 2
        port = (
            ((q[..., 0] - x) ** 2 + (q[..., 1] - y) ** 2 <= r * r + 1e-18)
            & (q[..., 2] >= 0)
            & (q[..., 2] <= 0.0045 + 1e-14)
        )
    liquid = block & (radial | port)
    return np.where(liquid, 2, np.where(block, 1, 0)).astype(np.uint8)


def top_port_audit(g):
    p = g["top_port"]
    x, y = p["center_xy_m"]
    r = p["diameter_m"] / 2
    if r <= 0 or r >= 0.0045 or p["axis_local"] != [0, 0, 1] or p["z_span_m"] != [0, 0.0045]:
        raise ValueError("Explicit positive Z-connected port geometry required")
    if abs(x) + r >= 0.01625 or abs(y) + r >= 0.0045 or abs(y) >= 0.0008 + r:
        raise ValueError("Top port leaves block or does not intersect main lumen")
    return dict(
        status="PASS",
        diameter_m=2 * r,
        center_top_local_m=[x, y, 0.0045],
        center_axis_local_m=[x, y, 0],
        axial_distance_to_end_PZT_m=0.01625 - x,
        boundary=p["post_injection_boundary"],
        x_plus="closed",
        intersection_with_lumen=True,
        reference_status=p["axial_reference_status"],
        model_depth_assumption=p["depth_model"],
        no_simulated_port_flow=True,
        bubble_junction_motion="NOT_ADMITTED; main-bore GPU stops before axial port neighborhood",
    )
