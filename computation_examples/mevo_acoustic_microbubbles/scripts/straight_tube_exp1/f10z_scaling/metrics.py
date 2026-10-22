"""Fixed-ID cloud statistics normalized by actual N; no inferred trajectories."""
import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import connected_components
from .common import R

def validate(t, x, ids):
    t, x, ids = np.asarray(t), np.asarray(x), np.asarray(ids)
    if x.ndim != 3 or x.shape[0] != len(t) or x.shape[2] != 3 or x.shape[1] < 2 or not len(t):
        raise ValueError('Expected nonempty (time, N>=2, xyz) trajectory')
    if not np.isfinite(x).all() or not np.isfinite(t).all() or t[0] != 0 or np.any(np.diff(t) <= 0):
        raise ValueError('Finite states and strictly increasing physical times required')
    if not np.array_equal(ids, np.arange(x.shape[1])):
        raise ValueError('Fixed contiguous particle IDs must match actual N')
    return t, x

def series(t, x):
    t, x = validate(t, x, np.arange(x.shape[1]))
    center = x.mean(1)
    q = x - center[:, None]
    tensor = np.einsum('tni,tnj->tij', q, q) / x.shape[1]
    eigen, axes = np.linalg.eigh(tensor)
    identifiable = (eigen[:, -1] - eigen[:, -2]) > 1e-3 * eigen[:, -1]
    pair = np.array([pdist(a) for a in x])
    edge = pair <= 2.1 * R
    counts, largest = [], []
    for e in edge:
        n, labels = connected_components(squareform(e), directed=False)
        counts.append(n)
        largest.append(np.bincount(labels).max())
    return dict(time_s=t, centroid_m=center, centered_positions_m=q, tensor_m2=tensor,
        rg_m=np.sqrt(np.maximum(eigen.sum(1), 0)),
        aspect=np.sqrt(eigen[:, -1] / np.maximum(eigen[:, 0], 1e-30)),
        directional_xy=np.sqrt(tensor[:, 0, 0] / np.maximum(tensor[:, 1, 1], 1e-30)),
        principal_axis=axes[:, :, -1], axis_identifiable=identifiable,
        components=np.array(counts), largest_component=np.array(largest), edges=edge,
        min_distance_m=pair.min(1), contact_count=np.sum(pair <= 2 * R + 1e-9, axis=1))

def edge_events(t, edges):
    """Every fixed-pair transition, including transient edges; endpoint events censored."""
    out = {}
    for k in range(edges.shape[1]):
        idx = np.flatnonzero(edges[1:, k] != edges[:-1, k]) + 1
        if len(idx):
            out[str(k)] = [(float(t[i]), bool(edges[i, k])) for i in idx]
    return out
