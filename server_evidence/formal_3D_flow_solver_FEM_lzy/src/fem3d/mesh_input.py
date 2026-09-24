"""Frozen tagged-surface adapter. This is the only geometry unit conversion."""
import json
from pathlib import Path
import numpy as np
from .audit import sha256

FACET_NAMES = {1: "WALL", 2: "OUTLET_03", 3: "OUTLET_01", 4: "INLET", 5: "OUTLET_02"}
PORT_TAGS = {"inlet": 4, "outlet_01": 3, "outlet_02": 5, "outlet_03": 2}


def read_contract(root):
    root = Path(root)
    lock = json.loads((root/"configs/stage01_source_lock.json").read_text())
    path = Path(lock["source_contract_path"])
    if sha256(path) != lock["source_contract_sha256"]:
        raise ValueError("FAIL: frozen Stage 0 source contract hash mismatch")
    contract = json.loads(path.read_text())
    if contract["status"] != "PASS" or contract["geometry_units"] != "um":
        raise ValueError("FAIL: source contract status/units invalid")
    for pathkey, hashkey, lockkey in (("geometry_path", "geometry_sha256", "source_surface_sha256"),
                                     ("meter_geometry_path", "meter_geometry_sha256", "meter_stl_sha256")):
        if sha256(contract[pathkey]) != contract[hashkey] or contract[hashkey] != lock[lockkey]:
            raise ValueError(f"FAIL: frozen {pathkey} hash mismatch")
    mapping = {p["name"]:p["surface_entity_id"] for p in contract["inlets"]+contract["outlets"]}
    if mapping != PORT_TAGS or contract["wall_entity_id"] != 1:
        raise ValueError("FAIL: boundary semantics do not match frozen contract")
    if len(contract["inlets"]) != 1 or len(contract["outlets"]) != 3:
        raise ValueError("FAIL: port counts changed")
    return contract, lock


def exact_weld(points_m, triangles, tolerance_m=0.0):
    """Merge only identical coordinates. No movement, averaging or tolerance repair."""
    if tolerance_m != 0:
        raise ValueError("Tolerance welding is not authorized; distinct points must stay distinct")
    points = np.asarray(points_m, dtype=np.float64)
    triangles = np.asarray(triangles, dtype=np.int64)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise ValueError("Invalid coordinates")
    if triangles.ndim != 2 or triangles.shape[1] != 3 or triangles.min()<0 or triangles.max()>=len(points):
        raise ValueError("Invalid triangles")
    _, first, inverse = np.unique(points, axis=0, return_index=True, return_inverse=True)
    unique = points[first].copy()  # retain an original coordinate exactly
    welded = inverse[triangles]
    cross = np.cross(unique[welded[:,1]]-unique[welded[:,0]], unique[welded[:,2]]-unique[welded[:,0]])
    if (np.linalg.norm(cross,axis=1) <= 0).any():
        raise ValueError("Degenerate input triangle; no geometry repair is permitted")
    distance = float(np.linalg.norm(unique[inverse]-points,axis=1).max())
    if distance != 0 or not np.array_equal(unique[welded], points[triangles]):
        raise ValueError("Welding changed geometry")
    return unique, welded, {
        "points_before_welding":len(points),"unique_points_after_welding":len(unique),
        "merged_record_count":len(points)-len(unique),"maximum_merge_distance_m":distance,
        "welding_tolerance_m":0.0,"triangle_count_before":len(triangles),"triangle_count_after":len(welded),
        "retained_vertices_moved":0,"method":"exact coordinate equality, deterministic lexicographic order"}


def surface_topology(points, triangles):
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    directed = np.concatenate([triangles[:,[0,1]],triangles[:,[1,2]],triangles[:,[2,0]]])
    edges, inv, counts = np.unique(np.sort(directed,axis=1),axis=0,return_inverse=True,return_counts=True)
    orientation = np.bincount(inv,weights=np.where(directed[:,0]<directed[:,1],1,-1))
    graph=coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(points),len(points)))
    components=connected_components(graph,directed=False,return_labels=False)
    data={"boundary_edges":int((counts==1).sum()),"nonmanifold_edges":int((counts>2).sum()),
          "inconsistent_oriented_edges":int((orientation!=0).sum()),"connected_surface_components":int(components)}
    if data != {"boundary_edges":0,"nonmanifold_edges":0,"inconsistent_oriented_edges":0,"connected_surface_components":1}:
        raise ValueError(f"Closed source topology failed: {data}")
    return data


def adapt_surface(root):
    import pyvista as pv
    contract, lock = read_contract(root)
    surface = pv.read(contract["geometry_path"])
    if "CellEntityIds" not in surface.cell_data:
        raise ValueError("Tagged VTP must contain CellEntityIds")
    packed = np.asarray(surface.faces).reshape(-1,4)
    if (packed[:,0] != 3).any() or len(packed) != 67262:
        raise ValueError("Source triangle topology changed")
    tags = np.asarray(surface.cell_data["CellEntityIds"], dtype=np.int32)
    if set(tags) != set(FACET_NAMES) or int((tags==1).sum()) != 67071:
        raise ValueError("Source tag partition changed")
    # Exactly one conversion, here. All downstream arrays and sizes are metres.
    points_m = np.asarray(surface.points,dtype=np.float64) * 1e-6
    unique, welded, welding = exact_weld(points_m,packed[:,1:])
    used = np.unique(welded)
    remap = np.full(len(unique),-1,dtype=np.int64)
    remap[used] = np.arange(len(used))
    points = unique[used]
    triangles = remap[welded]
    welding.update(unreferenced_unique_point_records=int(len(unique)-len(used)),
                   mesher_boundary_vertex_count=len(points),
                   unreferenced_record_policy="Not geometric surface vertices: omit only records unused by every frozen triangle; all triangles unchanged")
    topology = surface_topology(points,triangles)
    edge = np.concatenate([triangles[:,[0,1]],triangles[:,[1,2]],triangles[:,[2,0]]])
    lengths=np.linalg.norm(points[edge[:,0]]-points[edge[:,1]],axis=1)
    metadata={"input_units":"um","solver_units":"m","conversion_factor":1e-6,
              "area_conversion_factor":1e-12,"source_surface_sha256":contract["geometry_sha256"],
              "source_contract_sha256":lock["source_contract_sha256"],"welding":welding,"topology":topology,
              "boundary_edge_lengths_m":dict(zip(["minimum","P5","median","P95","maximum"],np.quantile(lengths,[0,.05,.5,.95,1]).tolist())),
              "facet_tag_meanings":FACET_NAMES,"cell_tag_meanings":{100:"FLUID"}}
    return points,triangles,tags,metadata,contract


def load_mesh_config(path, profile):
    config=json.loads(Path(path).read_text())
    if set(config["profiles"]) != {"development_coarse","development_medium"}:
        raise ValueError("Only two development profiles are authorized")
    selected=config["profiles"][profile]
    size=selected["bulk_target_size_m"]
    if not np.isfinite(size) or size<=0 or config["solver_units"]!="m":
        raise ValueError("Invalid SI mesh size")
    if selected["optional_interior_size_control"] is not None:
        raise ValueError("Optional control is reserved; only the audited constant interior target is supported")
    return config,selected
