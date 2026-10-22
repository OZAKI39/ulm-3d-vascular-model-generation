"""Independent P1 barycentric/raw-node checks of the unchanged VTKFlowField API."""
import time
import numpy as np
import pyvista as pv
from mevo_mb.flow.fields import VTKFlowField
from scripts.stage4.common import ROOT, new_run, require, write, finish, sha


def main():
    source = require("simvascular/steady")
    mesh = require("simvascular/mesh")
    out = new_run("simvascular/import", source, ["simvascular/steady"])
    path = ROOT / source["snapshot"]
    grid = pv.read(path)
    points = np.asarray(grid.points,float)
    field = VTKFlowField(path,velocity_array=source["velocity_array"],
                         coordinate_scale_to_m=source["coordinate_scale_to_m"],
                         velocity_scale_to_m_s=source["velocity_scale_to_m_s"],spatial_method="tetra_p1")
    raw = np.asarray(grid[source["velocity_array"]],float)
    tetra = grid.cells_dict[pv.CellType.TETRA]
    rng = np.random.default_rng(42026002)
    batches = {}
    nodes = rng.choice(grid.n_points,400,replace=False)
    batches["nodes"] = (points[nodes],raw[nodes])
    for name,pool,count,alpha in [("interior",np.arange(len(tetra)),200,1.),
                                  ("bifurcation",np.flatnonzero(np.abs(grid.cell_centers().points[:,0])<.001),200,1.)]:
        cells = rng.choice(pool,count,replace=False)
        weights = rng.dirichlet(np.ones(4)*alpha,count)
        batches[name] = (np.einsum("ni,nij->nj",weights,points[tetra[cells]]),
                         np.einsum("ni,nij->nj",weights,raw[tetra[cells]]))
    surface = pv.read((ROOT / mesh["mesh"]).parent / "mesh-surfaces/WALL.vtp")
    # Find boundary-face owners in the actual RESULT node order, independent of probe interpolation.
    from scipy.spatial import cKDTree
    _,global_nodes = cKDTree(grid.points).query(surface.points)
    wall_tri = global_nodes[surface.faces.reshape(-1,4)[:,1:]]
    owners = {}
    for index,cell in enumerate(tetra):
        for inds in ([0,1,2],[0,1,3],[0,2,3],[1,2,3]):
            owners[tuple(sorted(cell[inds]))] = index
    chosen = rng.choice(len(wall_tri),200,replace=False)
    near_x,near_u = [],[]
    for j in chosen:
        tri = wall_tri[j]
        owner = tetra[owners[tuple(sorted(tri))]]
        # 1% tetra-centroid + 99% boundary-face centroid: inside, genuinely near wall.
        near_x.append(.99*points[tri].mean(axis=0)+.01*points[owner].mean(axis=0))
        near_u.append(.99*raw[tri].mean(axis=0)+.01*raw[owner].mean(axis=0))
    batches["near_wall"] = (np.array(near_x),np.array(near_u))
    results = {}
    for name,(x,expected) in batches.items():
        start = time.perf_counter()
        actual = field.velocity(x,123.)
        error = np.linalg.norm(actual-expected,axis=1)
        results[name] = dict(count=len(x),max_absolute_error_m_s=float(error.max()),
                             relative_L2=float(np.linalg.norm(actual-expected)/max(np.linalg.norm(expected),1e-30)),
                             elapsed_s=time.perf_counter()-start,pass_test=bool(error.max()<1e-10))
        np.savez_compressed(out / (name+".npz"),positions_m=x,raw_P1_velocity_m_s=expected,imported_velocity_m_s=actual)
    query = np.concatenate([v[0] for v in batches.values()])[:256]
    timings = []
    for _ in range(30):
        start = time.perf_counter()
        field.velocity(query,0.)
        timings.append(time.perf_counter()-start)
    outside_rejected = False
    try:
        field.velocity([[1.,1.,1.]],0.)
    except ValueError:
        outside_rejected = True
    report = dict(status="PASS" if all(r["pass_test"] for r in results.values()) and outside_rejected else "FAIL",
                  comparisons=results,outside_rejected=outside_rejected,snapshot=source["snapshot"],
                  snapshot_sha256=sha(path),mesh_sha256=mesh["mesh_sha256"],velocity_array=source["velocity_array"],
                  coordinate_scale_to_m=1.,velocity_scale_to_m_s=1.,
                  spatial_method="tetra_p1",barycentric_roundoff_tolerance=field.fields[0].sampler.roundoff,
                  rejected_locator_candidates=field.fields[0].sampler.fallback_queries,
                  batch256_median_s=float(np.median(timings)),batch256_minimum_s=float(min(timings)),
                  interface="velocity(positions_m, time_s)",scope="Interpolation fidelity to native P1 data; no claim of local CFD physical accuracy")
    write(out / "SIMVASCULAR_FLOW_IMPORT_VALIDATION.json",report)
    finish(out,report)


if __name__ == "__main__":
    main()
