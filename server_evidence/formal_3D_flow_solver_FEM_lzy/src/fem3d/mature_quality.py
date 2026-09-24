"""One Gmsh minSICN evaluator for external meshes; no mesh generation."""
import numpy as np
from scipy.spatial import cKDTree

from .mesh_input import FACET_NAMES
from .mesh_qc import boundary_partition, quantiles, tetra_volumes, triangle_geometry
from .planar_port import p2_proxy


def evaluate_gmsh(data, msh_path=None):
    import gmsh
    points, tetra = data['points_m'], data['tetra']
    if not np.isfinite(points).all():
        raise ValueError('Nonfinite coordinates cannot enter the quality evaluator')
    gmsh.initialize()
    try:
        gmsh.option.setNumber('General.Terminal', 0)
        gmsh.model.add('imported_mesh_quality_only')
        gmsh.model.addDiscreteEntity(3, 100)
        gmsh.model.mesh.addNodes(3, 100, np.arange(1, len(points)+1), points.ravel())
        ids = np.arange(1, len(tetra)+1)
        gmsh.model.mesh.addElementsByType(100, 4, ids, (tetra+1).ravel())
        gmsh.model.addPhysicalGroup(3, [100], 100, 'FLUID')
        offset = len(ids)
        for tag in np.unique(data['facet_tags']):
            tag = int(tag)
            triangles = data['boundary_triangles'][data['facet_tags'] == tag]
            gmsh.model.addDiscreteEntity(2, tag)
            gmsh.model.mesh.addElementsByType(tag, 2, np.arange(offset+1, offset+1+len(triangles)), (triangles+1).ravel())
            gmsh.model.addPhysicalGroup(2, [tag], tag, FACET_NAMES.get(tag, f'UNKNOWN_{tag}'))
            offset += len(triangles)
        quality = np.asarray(gmsh.model.mesh.getElementQualities(ids, 'minSICN'))
        if msh_path is not None:
            gmsh.option.setNumber('Mesh.MshFileVersion', 4.1)
            gmsh.option.setNumber('Mesh.Binary', 1)
            gmsh.write(str(msh_path))
        data.update(min_sicn=quality, gmsh_element_ids=ids, gmsh_node_ids=np.arange(1,len(points)+1))
        return {'gmsh_version': gmsh.__version__, 'mesh_generation_calls': 0, 'optimizer_calls': 0,
                'metric': 'minSICN', 'coordinates_modified': False}
    finally:
        gmsh.finalize()


def summarize(data, boundary):
    points, tetra, q = data['points_m'], data['tetra'], data['min_sicn']
    volumes = tetra_volumes(points, tetra)
    validity = {'zero_volume': int(np.count_nonzero(volumes == 0)),
                'negative_volume': int(np.count_nonzero(volumes < 0)),
                'nonfinite_volume': int(np.count_nonzero(~np.isfinite(volumes)))}
    exterior, tags, topology = boundary_partition(points, tetra, data['boundary_triangles'], data['facet_tags'])
    _, centers, vectors = triangle_geometry(points, exterior)
    low = np.flatnonzero(q < 0.1)
    cell_centers = points[tetra].mean(axis=1)
    tree = cKDTree(centers)
    distance, nearest = tree.query(cell_centers[low])
    counts = {name: int(np.count_nonzero(tags[nearest] == tag)) for tag, name in FACET_NAMES.items()}
    residual = [{'cell_index': int(i), 'min_sicn': float(q[i]), 'centroid_m': cell_centers[i].tolist(),
                 'volume_m3': float(volumes[i]), 'vertex_coordinates_m': points[tetra[i]].tolist(),
                 'nearest_boundary_triangle_id': int(j), 'nearest_boundary_patch': FACET_NAMES[int(tags[j])],
                 'distance_to_nearest_boundary_triangle_center_m': float(d)} for i, j, d in zip(low, nearest, distance)]
    origin = points.mean(axis=0)
    enclosed = abs(float(np.einsum('ij,ij->', centers-origin, vectors)/3))
    closure = abs(float(volumes.sum())-enclosed)/enclosed if enclosed else float('inf')
    return {'measurement_status': 'PASS', 'validity': validity, 'topology': topology,
            'cell_tags_valid': bool(np.all(data['cell_tags'] == 100)),
            'boundary_fidelity': boundary, 'volume_closure': {'relative_error': closure, 'volume_m3': float(volumes.sum())},
            'quality': {'min_sicn': quantiles(q), 'total_below_0_1': len(low),
                        'cap_adjacent_below_0_1': len(low)-counts['WALL'],
                        'low_quality_nearest_boundary_counts': counts,
                        'metric': 'Gmsh minSICN; same API as Stage 1/1.7',
                        'location_method': 'Nearest exterior triangle center, as in Stage 1.7',
                        'nonfinite_quality_count': int(np.count_nonzero(~np.isfinite(q))),
                        'nonpositive_quality_count': int(np.count_nonzero(q <= 0))},
            'residual_cells': residual, 'proxy': p2_proxy(tetra),
            'unused_stored_vertex_count': len(points)-len(np.unique(tetra)),
            'fem_space_created': False, 'fem_solved': False}
