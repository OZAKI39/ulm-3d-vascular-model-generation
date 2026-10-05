"""Stage 4D mesh glue: reuse official SimVascular/TetGen and existing mesh diagnostics."""
import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree
import yaml

from scripts.stage4.common import ROOT, new_run, require, sha, write, finish


def compact(points, triangles):
    used, inverse = np.unique(triangles, return_inverse=True)
    face = pv.PolyData(points[used], np.column_stack([np.full(len(triangles),3), inverse.reshape(-1,3)]))
    face.point_data["GlobalNodeID"] = (used+1).astype(np.int32)
    return face


def audit(raw, output, workflow):
    sys.path.insert(0, str(workflow / "mesh_generate/src"))
    from vascular_validation.mesh_diagnostics import volume_validity, tetra_sicn, quality_summary
    grid = pv.read(raw / "volume.vtu")
    if set(grid.cells_dict) != {pv.CellType.TETRA}:
        raise ValueError("Expected linear tetrahedra from the existing mesher")
    points = np.array(grid.points, dtype=float)
    tetra = np.asarray(grid.cells_dict[pv.CellType.TETRA], np.int64)
    validity = volume_validity(points, tetra)
    quality = tetra_sicn(points, tetra)
    stats = quality_summary(quality)
    if not stats["q_min"] > 0:
        raise ValueError("Nonpositive signed tetrahedron quality")
    all_faces = np.concatenate([tetra[:, indices] for indices in ([0,1,2],[0,1,3],[0,2,3],[1,2,3])])
    owners = np.tile(np.arange(len(tetra)),4)
    _, first, counts = np.unique(np.sort(all_faces,axis=1),axis=0,return_index=True,return_counts=True)
    assert np.max(counts) <= 2
    boundary = all_faces[first[counts == 1]].copy()
    adjacent = owners[first[counts == 1]]
    xyz = points[boundary]
    direction = xyz.mean(axis=1)-points[tetra[adjacent]].mean(axis=1)
    wrong = np.einsum("ij,ij->i",np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),direction) < 0
    boundary[wrong] = boundary[wrong][:,[0,2,1]]
    lookup = {tuple(sorted(face)):i for i,face in enumerate(boundary)}
    tags = np.zeros(len(boundary),np.int32)
    tree = cKDTree(points)
    for face_id in (1,2,3,4):
        face = pv.read(raw / ("face_"+str(face_id)+".vtp"))
        distances, index = tree.query(face.points)
        assert distances.max() < 1e-10
        for triangle in index[face.faces.reshape(-1,4)[:,1:]]:
            j = lookup[tuple(sorted(triangle))]
            assert tags[j] == 0
            tags[j] = face_id
    assert np.all(tags > 0)
    (output / "mesh-surfaces").mkdir(parents=True)
    grid.point_data["GlobalNodeID"] = np.arange(1,len(points)+1,dtype=np.int32)
    grid.cell_data["GlobalElementID"] = np.arange(1,len(tetra)+1,dtype=np.int32)
    grid.cell_data["minSICN"] = quality
    grid.save(output / "mesh-complete.mesh.vtu")
    surface = compact(points,boundary)
    surface.cell_data["ModelFaceID"] = tags
    surface.cell_data["GlobalElementID"] = (adjacent+1).astype(np.int32)
    surface.cell_data["BranchID"] = np.where(surface.cell_centers().points[:,0] <= 0,0,
                                            np.where(surface.cell_centers().points[:,1] >= 0,1,2))
    surface.save(output / "mesh-complete.exterior.vtp")
    assert surface.n_open_edges == 0
    face_info = {}
    for name,face_id in [("WALL",1),("INLET",2),("OUTLET_A",3),("OUTLET_B",4)]:
        selected = tags == face_id
        face = compact(points,boundary[selected])
        face.cell_data["GlobalElementID"] = (adjacent[selected]+1).astype(np.int32)
        face.save(output / "mesh-surfaces" / (name+".vtp"))
        face_info[name] = dict(triangles=face.n_cells,area_m2=float(face.area))
    np.savez_compressed(output / "mesh_arrays.npz", points_m=points,tetra=tetra,
                        boundary_triangles=boundary,facet_tags=tags,adjacent_tetra=adjacent)
    return dict(validity=validity,quality=stats,faces=face_info,points=len(points),tetrahedra=len(tetra),
                diagnostic_source=str(workflow / "mesh_generate/src/vascular_validation/mesh_diagnostics.py"),
                diagnostic_source_sha256=sha(workflow / "mesh_generate/src/vascular_validation/mesh_diagnostics.py"))


def main():
    c = yaml.safe_load((ROOT / "configs/stage4/bifurcation.yaml").read_text())
    workflow = yaml.safe_load((ROOT / "configs/stage4/workflow.yaml").read_text())
    source = require("geometry")
    out = new_run("simvascular/mesh", dict(config=c, workflow=workflow), ["synthetic_flow"])
    request = dict(surface=str(ROOT / source["surface"]), output=str(out / "raw"),
                   edge_size_m=c["geometry"]["mesh_edge_size_m"])
    write(out / "MESH_REQUEST.json", request)
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1", QT_QPA_PLATFORM="offscreen",
               OMP_NUM_THREADS="1", MEVO_STAGE4_MESH_REQUEST=str(out / "MESH_REQUEST.json"),
               XDG_CONFIG_HOME=str(out / "app_config"), XDG_CACHE_HOME=str(out / "app_cache"))
    env.pop("PYTHONPATH",None)
    command = [workflow["simvascular_launcher"],"-python","--",str(ROOT / "scripts/stage4/mesh_embedded_sv.py")]
    with (out / "mesh.log").open("w") as log:
        process = subprocess.run(command,env=env,stdout=log,stderr=subprocess.STDOUT)
    write(out / "EXECUTION.json", dict(command=command, returncode=process.returncode,
                                       launcher_sha256=sha(workflow["simvascular_launcher"])))
    if process.returncode:
        finish(out,dict(status="FAIL",error="Official mesher failed; see preserved mesh.log"))
        raise RuntimeError("Meshing failed")
    report = audit(out / "raw",out / "SV_MESH",Path(workflow["local_meshing_workflow"]))
    mesh_path = out / "SV_MESH/mesh-complete.mesh.vtu"
    finish(out,dict(status="PASS",**report,geometry_hash=source["geometry_hash"],
                    source_surface_sha256=source["surface_sha256"], mesh=str(mesh_path.relative_to(ROOT)),
                    mesh_sha256=sha(mesh_path),coordinate_scale_to_m=1.,
                    source="Existing official SimVascular/TetGen workflow; no custom CFD mesher/solver"))


if __name__ == "__main__":
    main()
