"""Stage 4B: smooth implicit Y, planar caps, stable IDs and explicit SI exports."""
import time

import numpy as np
import pyvista as pv
import yaml

from mevo_mb.coupling.bifurcation import BifurcationGeometry
from scripts.stage4.common import ROOT, new_run, write, sha, finish


def build_surface(g):
    step = g.config["surface_spacing_m"]
    daughter = next(iter(g.directions.values()))
    # Include all rounded ends before clipping. Only final caps have sharp rims.
    extent_y = (g.daughter_length + 3*g.daughter_radius)*daughter[1] + 2*g.parent_radius
    bounds = np.array([[-g.parent_length-5*g.parent_radius,
                       (g.daughter_length+3*g.daughter_radius)*daughter[0]+2*g.parent_radius],
                      [-extent_y, extent_y], [-2*g.parent_radius, 2*g.parent_radius]])
    lower = np.floor(bounds[:, 0]/step).astype(int)
    upper = np.ceil(bounds[:, 1]/step).astype(int)
    grid = pv.ImageData(dimensions=upper-lower+1, spacing=(step,)*3, origin=lower*step)
    points = grid.points
    values = np.empty(len(points))
    for start in range(0, len(points), 100000):
        values[start:start+100000] = g._tube_values(points[start:start+100000])[0]
    grid["smooth_lumen_level_m"] = values
    surface = grid.contour([0.0], scalars="smooth_lumen_level_m", method="flying_edges").triangulate().clean()
    assert surface.n_open_edges == 0
    surface = surface.compute_normals(consistent_normals=True, auto_orient_normals=True,
                                      split_vertices=False, inplace=False)
    for name in ["INLET", "OUTLET_A", "OUTLET_B"]:
        plane = g.planes()[name]
        surface = surface.clip_closed_surface(normal=-np.array(plane["normal"]),
                                               origin=plane["origin_m"], tolerance=1e-9)
        surface = surface.triangulate().clean(tolerance=1e-10, absolute=True)
        print(f"Clipped {name}: {surface.n_points} points, {surface.n_open_edges} open edges", flush=True)
    surface = surface.triangulate().clean(tolerance=1e-10, absolute=True)
    # Construct the symmetric control geometry from one half, avoiding tiny
    # clipping/triangulation differences being confounded with mirror controls.
    half = surface.clip(normal=(0, 1, 0), origin=(0, 0, 0), invert=False).triangulate().clean()
    half.points[np.abs(half.points[:, 1]) < 1e-10, 1] = 0.0
    reflected_points = np.asarray(half.points).copy()
    reflected_points[:, 1] *= -1
    reflected_faces = half.faces.reshape(-1, 4).copy()
    reflected_faces[:, 1:] = reflected_faces[:, 1:][:, ::-1]
    reflected = pv.PolyData(reflected_points, reflected_faces)
    surface = half.merge(reflected, merge_points=True, tolerance=1e-10).triangulate().clean()
    surface = surface.compute_normals(consistent_normals=True, auto_orient_normals=True,
                                      split_vertices=False, inplace=False)
    faces = surface.faces.reshape(-1, 4)[:, 1:]
    tri = np.asarray(surface.points, float)[faces]
    centers = tri.mean(axis=1)
    tags = np.ones(surface.n_cells, np.int32)
    for name, face_id in [("INLET", 2), ("OUTLET_A", 3), ("OUTLET_B", 4)]:
        plane = g.planes()[name]
        distance = np.abs((tri-np.array(plane["origin_m"])) @ np.array(plane["normal"]))
        tags[np.max(distance, axis=1) < 2e-8] = face_id
    surface.cell_data["ModelFaceID"] = tags
    surface.cell_data["BranchID"] = np.where(centers[:, 0] <= 0, 0,
                                            np.where(centers[:, 1] >= 0, 1, 2)).astype(np.int32)
    surface.field_data["coordinate_scale_to_m"] = [1.0]
    surface.field_data["geometry_hash"] = [g.identity]
    return surface


def main():
    config = yaml.safe_load((ROOT / "configs/stage4/bifurcation.yaml").read_text())
    g = BifurcationGeometry(config["geometry"])
    out = new_run("geometry", config, ["audit"])
    started = time.perf_counter()
    surface = build_surface(g)
    path = out / "labelled_lumen_m.vtp"
    surface.save(path)
    connected = surface.connectivity()
    tags = surface["ModelFaceID"]
    face_info = {}
    for name, face_id in [("WALL", 1), ("INLET", 2), ("OUTLET_A", 3), ("OUTLET_B", 4)]:
        face = surface.extract_cells(tags == face_id).extract_surface().triangulate()
        face.save(out / (name + ".vtp"))
        face_info[name] = dict(face_id=face_id, triangles=face.n_cells, area_m2=float(face.area))
    paths = []
    for branch_id, (direction, length) in enumerate([(np.array([-1,0,0]), g.parent_length),
                                                     *[(d, g.daughter_length) for d in g.directions.values()]]):
        points = np.linspace(0, length, 101)[:, None] * direction
        line = pv.PolyData(points)
        line.lines = np.r_[101, np.arange(101)]
        line.cell_data["BranchID"] = np.full(line.n_cells, branch_id, np.int32)
        paths.append(line)
    pv.merge(paths).save(out / "centerlines.vtp")
    write(out / "PLANES.json", g.planes())
    write(out / "BRANCH_IDS.json", {"PARENT":0, "BRANCH_A":1, "BRANCH_B":2,
                                   "UNASSIGNED":-1, "WALL_LOSS":-2})
    random = np.random.default_rng(310).uniform([-.014,-.008,-.002], [.012,.008,.002], (2000, 3))
    mirrored = random.copy()
    mirrored[:,1] *= -1
    symmetry = float(np.max(np.abs(g.wall_distance(random)-g.wall_distance(mirrored))))
    tests = dict(closed=surface.n_open_edges == 0, one_connected_lumen=len(np.unique(connected["RegionId"])) == 1,
                 four_face_ids=set(np.unique(tags)) == {1,2,3,4}, positive_volume=surface.volume > 0,
                 exact_implicit_mirror=symmetry < 1e-14,
                 cap_area_symmetry=abs(face_info["OUTLET_A"]["area_m2"]-face_info["OUTLET_B"]["area_m2"])
                 /face_info["OUTLET_A"]["area_m2"] < 1e-4)
    finish(out, dict(status="PASS" if all(tests.values()) else "FAIL", tests=tests,
                     geometry_hash=g.identity, surface=str(path.relative_to(ROOT)), surface_sha256=sha(path),
                     points=surface.n_points, triangles=surface.n_cells, volume_m3=float(surface.volume),
                     faces=face_info, implicit_symmetry_error_m=symmetry, elapsed_s=time.perf_counter()-started,
                     units="Coordinates metres, area m2, volume m3; explicit generated geometry contract",
                     scope="Smooth log-sum-exp tube union; flat terminal caps; positive-y half reflected for exact surface symmetry; no anatomical claim"))


if __name__ == "__main__":
    main()
