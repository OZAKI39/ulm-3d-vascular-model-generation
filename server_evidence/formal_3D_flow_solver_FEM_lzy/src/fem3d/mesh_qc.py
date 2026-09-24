"""Geometry-only checks: no field spaces, variational forms or solves."""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from .mesh_input import FACET_NAMES


def quantiles(values):
    return dict(zip(("minimum", "P1", "P5", "median", "P95", "maximum"),
                    np.quantile(values, [0, .01, .05, .5, .95, 1]).tolist()))


def tetra_volumes(points, tetra):
    coords = points[tetra]
    return np.einsum("ij,ij->i", coords[:,1]-coords[:,0],
                     np.cross(coords[:,2]-coords[:,0], coords[:,3]-coords[:,0])) / 6


def boundary_partition(points, tetra, triangles, tags):
    """Extract exterior from cell adjacency, independently of surface labels."""
    local_faces = np.array([[1,2,3], [0,3,2], [0,1,3], [0,2,1]])
    faces = tetra[:, local_faces].reshape(-1, 3)
    keys, first, inverse, counts = np.unique(np.sort(faces, axis=1), axis=0,
        return_index=True, return_inverse=True, return_counts=True)
    if np.any(counts > 2):
        raise ValueError("Nonmanifold tetra face adjacency")
    sort = np.argsort(inverse, kind="stable")
    offsets = np.r_[0, np.cumsum(counts)]
    shared = np.flatnonzero(counts==2)
    pair = np.c_[sort[offsets[shared]]//4, sort[offsets[shared]+1]//4]
    graph = coo_matrix((np.ones(len(pair)*2),
        (np.r_[pair[:,0], pair[:,1]], np.r_[pair[:,1], pair[:,0]])), shape=(len(tetra),len(tetra)))
    components = int(connected_components(graph, directed=False, return_labels=False))
    indices = first[counts==1]
    exterior = faces[indices].copy()
    owners = indices//4
    coords = points[exterior]
    normals = np.cross(coords[:,1]-coords[:,0], coords[:,2]-coords[:,0])
    to_interior = points[tetra[owners]].mean(axis=1) - coords.mean(axis=1)
    flip = np.einsum("ij,ij->i", normals, to_interior)>0
    exterior[flip] = exterior[flip][:,[0,2,1]]
    _, all_inverse = np.unique(np.vstack([np.sort(exterior,axis=1), np.sort(triangles,axis=1)]),
                               axis=0, return_inverse=True)
    n = len(exterior)
    ext_keys, labelled_keys = all_inverse[:n], all_inverse[n:]
    coverage = np.bincount(labelled_keys, minlength=all_inverse.max()+1)
    unlabeled = int(np.count_nonzero(coverage[ext_keys]==0))
    multiple = int(np.count_nonzero(coverage[ext_keys]>1))
    nonexterior = int(np.count_nonzero(~np.isin(labelled_keys,ext_keys)))
    invalid_tags = int(np.count_nonzero(~np.isin(tags,list(FACET_NAMES))))
    if unlabeled or multiple or nonexterior or invalid_tags:
        raise ValueError(f"Boundary partition failed: unlabelled={unlabeled}, multiple={multiple}, nonexterior={nonexterior}, invalid={invalid_tags}")
    values = np.zeros(len(coverage),dtype=np.int32)
    values[labelled_keys] = tags
    exterior_tags = values[ext_keys]
    result = {"connected_fluid_components": components, "exterior_facet_count": n,
              "interior_facet_count": int(len(shared)), "nonmanifold_facets": 0,
              "unlabelled_exterior_facets": unlabeled, "multiply_labelled_exterior_facets": multiple,
              "tagged_nonexterior_facets": nonexterior,
              "boundary_counts": {name: int(np.count_nonzero(exterior_tags==tag)) for tag,name in FACET_NAMES.items()},
              "normal_orientation_method": "Each exterior triangle is oriented away from its unique adjacent tetrahedron centroid"}
    return exterior, exterior_tags, result


def triangle_geometry(points, triangles):
    c = points[triangles]
    vectors = np.cross(c[:,1]-c[:,0],c[:,2]-c[:,0])/2
    return np.linalg.norm(vectors,axis=1), c.mean(axis=1), vectors


def fidelity(points, exterior, tags, source_points, source_triangles, source_tags, tolerance):
    used = np.unique(exterior)
    displacement, match = cKDTree(source_points).query(points[used])
    max_displacement = float(displacement.max())
    if max_displacement > tolerance or len(np.unique(match)) != len(source_points) or len(used)!=len(source_points):
        raise ValueError(f"Frozen boundary vertices changed: maximum displacement={max_displacement}")
    correspondence = np.full(len(points),-1,dtype=np.int64)
    correspondence[used] = match
    source_keys = np.c_[np.sort(source_triangles,axis=1),source_tags]
    final_keys = np.c_[np.sort(correspondence[exterior],axis=1),tags]
    def sorted_rows(array):
        return array[np.lexsort(array[:,::-1].T)]
    if not np.array_equal(sorted_rows(source_keys),sorted_rows(final_keys)):
        raise ValueError("Frozen tagged boundary triangles changed")
    area, _, _ = triangle_geometry(points,exterior)
    source_area, _, _ = triangle_geometry(source_points,source_triangles)
    return {"matched_tagged_source_triangles": len(exterior),
            "source_triangle_count": len(source_triangles), "final_triangle_count": len(exterior),
            "maximum_boundary_displacement_m": max_displacement, "roundoff_tolerance_m": tolerance,
            "coordinates_exactly_equal": bool(np.array_equal(points[used],source_points[match])),
            "tagged_triangle_connectivity_equal": True, "surface_remeshing": False,
            "source_surface_area_m2": float(source_area.sum()), "final_surface_area_m2": float(area.sum()),
            "source_bounds_m": [source_points.min(axis=0).tolist(),source_points.max(axis=0).tolist()],
            "final_boundary_bounds_m": [points[used].min(axis=0).tolist(),points[used].max(axis=0).tolist()],
            "comparison_method": "Bijective boundary vertex correspondence followed by equality of every triangle's vertex set AND its original tag; nearest distances only measure fidelity, never assign tags"}


def port_metrics(points, exterior, tags, contract, config):
    result = {}
    tolerance = config["fidelity"]["roundoff_factor"]*np.finfo(float).eps*np.abs(points).max()
    for port in contract["inlets"]+contract["outlets"]:
        selected = exterior[tags==port["surface_entity_id"]]
        area, centroids, vectors = triangle_geometry(points,selected)
        total = float(area.sum())
        centroid = (centroids*area[:,None]).sum(axis=0)/total
        normal = vectors.sum(axis=0)
        normal /= np.linalg.norm(normal)
        expected_normal = np.array(port["outward_normal"])
        origin = np.array(port["plane_origin_m"])
        plane_error = float(np.abs((points[selected].reshape(-1,3)-origin)@expected_normal).max())
        relative_error = abs(total-port["area_m2"])/port["area_m2"]
        shift = float(np.linalg.norm(centroid-origin))
        dot = float(normal@expected_normal)
        passed = (len(selected)==port["triangle_count"] and relative_error<=config["fidelity"]["port_relative_area_tolerance"]
            and dot>config["fidelity"]["normal_dot_minimum"] and shift<=tolerance
            and plane_error<=port["max_plane_deviation_m"]+tolerance)
        result[port["name"]] = {"entity_id": port["surface_entity_id"], "facet_count": len(selected),
            "area_m2": total, "source_area_m2": port["area_m2"],
            "absolute_error_m2": abs(total-port["area_m2"]), "relative_error": relative_error,
            "centroid_m": centroid.tolist(), "centroid_displacement_m": shift,
            "outward_normal": normal.tolist(), "normal_dot_product": dot,
            "max_plane_deviation_m": plane_error, "source_max_plane_deviation_m": port["max_plane_deviation_m"],
            "status": "PASS" if passed else "FAIL"}
    return result


def audit_mesh(data, source, contract, config):
    points, tetra = data["points_m"], data["tetra"]
    if not np.isfinite(points).all():
        raise ValueError("Nonfinite tetra coordinates")
    volume = tetra_volumes(points,tetra)
    validity = {"zero_volume": int(np.count_nonzero(volume==0)),
                "negative_volume": int(np.count_nonzero(volume<0)),
                "nonfinite_volume": int(np.count_nonzero(~np.isfinite(volume)))}
    if any(validity.values()):
        raise ValueError(f"Invalid tetrahedra: {validity}")
    if not np.all(data["cell_tags"]==100):
        raise ValueError("FLUID cell tag lost")
    exterior, tags, partition = boundary_partition(points,tetra,data["boundary_triangles"],data["facet_tags"])
    if partition["connected_fluid_components"]!=1:
        raise ValueError("Fluid is not one connected component")
    tolerance = config["fidelity"]["roundoff_factor"]*np.finfo(float).eps*np.abs(source["points_m"]).max()
    boundary = fidelity(points,exterior,tags,source["points_m"],source["triangles"],source["facet_tags"],tolerance)
    ports = port_metrics(points,exterior,tags,contract,config)
    if any(p["status"]!="PASS" for p in ports.values()):
        raise ValueError(f"Port comparison failed: {ports}")
    area, centers, vectors = triangle_geometry(points,exterior)
    source_area, source_centers, source_vectors = triangle_geometry(source["points_m"],source["triangles"])
    origin = source["points_m"].mean(axis=0)
    source_volume = abs(float(np.einsum("ij,ij->",source_centers-origin,source_vectors)/3))
    relative_volume_error = abs(float(volume.sum())-source_volume)/source_volume
    if relative_volume_error>1e-10:
        raise ValueError("Cell volume does not equal frozen closed-surface enclosed volume")
    edges = np.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
    lengths = np.linalg.norm(points[tetra[:,edges[:,1]]]-points[tetra[:,edges[:,0]]],axis=2)
    ratio = lengths.max(axis=1)/lengths.min(axis=1)
    quality = data["min_sicn"]
    if not np.isfinite(quality).all() or not np.isfinite(ratio).all() or (quality<=0).any():
        raise ValueError("Nonfinite/nonpositive Gmsh signed quality or edge ratio")
    threshold = config["quality"]["advisory_min_sicn_below"]
    low = np.flatnonzero(quality<threshold)
    order = np.argsort(quality)[:20]
    cell_centers = points[tetra].mean(axis=1)
    distances, nearest = cKDTree(centers).query(cell_centers[order])
    worst = [{"cell_index": int(i), "gmsh_element_id": int(data["gmsh_element_ids"][i]),
        "min_sicn": float(quality[i]), "volume_m3": float(volume[i]), "edge_ratio": float(ratio[i]),
        "centroid_m": cell_centers[i].tolist(), "nearest_boundary_patch": FACET_NAMES[int(tags[j])],
        "distance_to_nearest_boundary_triangle_center_m": float(d)} for i,j,d in zip(order,nearest,distances)]
    low_nearest = cKDTree(centers).query(cell_centers[low])[1] if len(low) else np.array([],dtype=int)
    return {"hard_gate_status": "PASS", "quality_status": "ADVISORY / MANUAL_REVIEW" if len(low) else "PASS",
        "status": "CONDITIONAL PASS" if len(low) else "PASS", "tetrahedron_count": len(tetra),
        "vertex_count": len(points), "cell_tags": {"100":len(tetra)}, "validity": validity,
        "topology": partition, "boundary_fidelity": boundary, "ports": ports,
        "volume_closure": {"tetra_volume_m3": float(volume.sum()), "source_enclosed_volume_m3": source_volume,
            "relative_error": relative_volume_error, "tolerance":1e-10},
        "quality": {"volume_m3":quantiles(volume), "max_to_min_edge_ratio":quantiles(ratio),
            "gmsh_min_sicn":quantiles(quality), "metric": "Gmsh minSICN: signed inverse condition number; higher is better, ideal regular tetrahedron 1",
            "advisory_threshold":threshold, "advisory_count":len(low), "advisory_fraction":len(low)/len(tetra),
            "threshold_role":config["quality"]["advisory_threshold_role"],
            "low_quality_nearest_boundary_counts":{name:int(np.count_nonzero(tags[low_nearest]==tag)) for tag,name in FACET_NAMES.items()},
            "location_method":"Nearest boundary triangle center is approximate context only; junction/narrow-vessel interpretation requires image review",
            "worst_elements":worst}}
