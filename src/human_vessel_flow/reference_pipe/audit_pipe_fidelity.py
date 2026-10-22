"""Original P1 cap trace quadrature; no field correction."""
import numpy as np

def boundary_profile(poly, radius, mean=.1):
    """Exact degree-four quadrature of the squared P1-trace error on cap triangles."""
    xyz = np.asarray(poly.points, float)
    tri = poly.faces.reshape(-1, 4)[:, 1:]
    vertices = xyz[tri]
    areas = np.linalg.norm(np.cross(vertices[:, 1] - vertices[:, 0], vertices[:, 2] - vertices[:, 0]), axis=1) / 2
    # Perimeter IDs from topology, never a radial tolerance that could erase valid nodes.
    edges, counts = np.unique(np.sort(np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1), axis=0, return_counts=True)
    rim = np.unique(edges[counts == 1])
    values = 2 * mean * (1 - (xyz[:, 0]**2 + xyz[:, 1]**2) / radius**2)
    rim_raw = float(np.max(abs(values[rim])))
    # The analytical rim value is zero; native Float32 coordinates can differ by one ULP.
    values[rim] = 0.
    bary, weights = [], []
    for a, b, w in [(.445948490915965, .108103018168070, .223381589678011),
                    (.091576213509771, .816847572980459, .109951743655322)]:
        bary.extend([[a, a, b], [a, b, a], [b, a, a]]); weights.extend([w] * 3)
    bary, weights = np.array(bary), np.array(weights)
    query = np.einsum('qi,tij->tqj', bary, vertices)
    exact = 2 * mean * (1 - np.sum(query[..., :2]**2, axis=-1) / radius**2)
    actual = np.einsum('qi,ti->tq', bary, values[tri])
    integrate = lambda x: float(np.einsum('t,q,tq->', areas, weights, x))
    q = float(np.dot(areas, values[tri].mean(axis=1)))
    result = dict(profile_relative_L2=float(np.sqrt(integrate((actual-exact)**2) / integrate(exact**2))),
                  Q_discrete_m3_s=q, Q_exact_m3_s=np.pi * radius**2 * mean,
                  Q_relative_error=abs(q / (np.pi * radius**2 * mean) - 1),
                  maximum_nodal_velocity_m_s=float(values.max()), rim_nodes=len(rim),
                  rim_velocity_m_s=0., rim_analytic_coordinate_roundoff_m_s=rim_raw,
                  profile_scaling_factor=1., node_count=poly.n_points, triangles=len(tri))
    result['PASS'] = bool(result['profile_relative_L2'] <= .01 and result['Q_relative_error'] <= .005 and rim_raw < 1e-7)
    return result, values
