"""Independent exterior extraction and physical boundary identity checks."""
import numpy as np
from scipy.spatial import cKDTree

from .mesh_qc import boundary_partition, fidelity


def compact_surface(data):
    triangles = data.get('boundary_triangles', data.get('triangles'))
    used = np.unique(triangles)
    return {'points_m': data['points_m'][used], 'triangles': np.searchsorted(used, triangles),
            'facet_tags': data['facet_tags']}


def verify_boundary(data, source, tolerance=1e-15, rim_coordinates=None):
    """Never round coordinates, guess markers or change a stored coordinate."""
    source = compact_surface(source)
    points, tri, tags = data['points_m'], data['boundary_triangles'], data['facet_tags']
    result = {'status': 'FAIL', 'tolerance_m': tolerance, 'coordinate_rounding': False,
              'geometry_pass': False, 'connectivity_pass': False, 'references_pass': False,
              'maximum_boundary_displacement_m': None, 'maximum_wall_displacement_m': None,
              'maximum_cap_displacement_m': None, 'maximum_rim_displacement_m': None}
    try:
        if not np.isfinite(points).all():
            raise ValueError('Non-finite coordinates')
        exterior, ext_tags, partition = boundary_partition(points, data['tetra'], tri, tags)
        used = np.unique(exterior)
        distances, matched = cKDTree(source['points_m']).query(points[used])
        result['maximum_boundary_displacement_m'] = float(distances.max())
        result['exterior_partition'] = partition
        bijective = len(used) == len(source['points_m']) and len(np.unique(matched)) == len(source['points_m'])
        result['geometry_pass'] = bool(bijective and distances.max() <= tolerance)
        if not bijective:
            raise ValueError('Boundary vertices are not in bijective correspondence')
        restored = np.empty_like(source['points_m'])
        restored[matched] = points[used]
        displacements = np.linalg.norm(restored - source['points_m'], axis=1)
        for name, mask in [('wall', source['facet_tags'] == 1), ('cap', source['facet_tags'] != 1)]:
            vertices = np.unique(source['triangles'][mask])
            result[f'maximum_{name}_displacement_m'] = float(displacements[vertices].max()) if len(vertices) else 0.
        if rim_coordinates is None:
            wall_vertices = np.unique(source['triangles'][source['facet_tags'] == 1])
            cap_vertices = np.unique(source['triangles'][source['facet_tags'] != 1])
            rim = np.intersect1d(wall_vertices, cap_vertices)
        else:
            d, rim = cKDTree(source['points_m']).query(rim_coordinates)
            if np.any(d != 0):
                raise ValueError('Contract rim coordinates missing from source boundary')
        result['maximum_rim_displacement_m'] = float(displacements[rim].max()) if len(rim) else 0.
        # Also test untagged identity, so lost references are distinguished from geometry.
        untagged = fidelity(points, exterior, np.ones(len(exterior), int), source['points_m'],
                            source['triangles'], np.ones(len(source['triangles']), int), tolerance)
        result['connectivity_pass'] = bool(untagged['tagged_triangle_connectivity_equal'])
        exact = fidelity(points, exterior, ext_tags, source['points_m'], source['triangles'], source['facet_tags'], tolerance)
        result.update(status='PASS', references_pass=True, exact_coordinates=exact['coordinates_exactly_equal'],
                      boundary_fidelity=exact, triangles=len(exterior), source_triangles=len(source['triangles']),
                      method='Bijective unrounded float64 vertex matching; unordered triangle triples plus original reference; independent tetra exterior extraction')
    except (ValueError, IndexError) as exc:
        result['error'] = str(exc)
    return result
