"""Preserve all raw SWC components and parent edges, including cycles."""
from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


@dataclass
class SWCGraph:
    ids: np.ndarray
    xyz_um: np.ndarray
    radius_um: np.ndarray
    parent: np.ndarray
    type_code: np.ndarray
    edges: np.ndarray  # array indices, parent -> child
    component: np.ndarray
    degree: np.ndarray
    raw_xyz: np.ndarray

    @property
    def index(self):
        return {int(n): i for i, n in enumerate(self.ids)}

    @property
    def length_um(self):
        return np.linalg.norm(self.xyz_um[self.edges[:, 1]] - self.xyz_um[self.edges[:, 0]], axis=1)


def load_swc(path, xyz_scale_um):
    """xyz_scale_um must be supplied from provenance, never inferred from extent."""
    raw = np.loadtxt(Path(path), comments='#', ndmin=2)
    if raw.shape[1] != 7 or not np.all(np.isfinite(raw)):
        raise ValueError('Expected seven finite SWC columns')
    for col in (0, 1, 6):
        if not np.all(raw[:, col] == np.rint(raw[:, col])):
            raise ValueError('Noninteger SWC identifier')
    ids, parent = raw[:, 0].astype(np.int64), raw[:, 6].astype(np.int64)
    if len(np.unique(ids)) != len(ids):
        raise ValueError('Duplicate node IDs: refusing ambiguous graph')
    index = {int(n): i for i, n in enumerate(ids)}
    if any(int(p) not in index for p in parent if p != -1):
        raise ValueError('Missing parent: refusing silent edge removal')
    edges = np.asarray([(index[int(p)], i) for i, p in enumerate(parent) if p != -1], dtype=int).reshape(-1, 2)
    u, v = edges.T
    adj = coo_matrix((np.ones(2*len(edges)), (np.r_[u, v], np.r_[v, u])), shape=(len(ids), len(ids))).tocsr()
    _, comp = connected_components(adj, directed=False)
    degree = np.bincount(edges.ravel(), minlength=len(ids))
    scale = np.asarray(xyz_scale_um, dtype=float)
    if scale.shape != (3,) or np.any(scale <= 0) or not np.all(np.isfinite(scale)):
        raise ValueError('Invalid coordinate scale')
    return SWCGraph(ids, raw[:, 2:5]*scale, raw[:, 5], parent, raw[:, 1].astype(int), edges, comp, degree, raw[:, 2:5])


def audit_graph(g):
    n, e, c = len(g.ids), len(g.edges), len(np.unique(g.component))
    unique_xyz, counts = np.unique(g.xyz_um, axis=0, return_counts=True)
    degree, count = np.unique(g.degree, return_counts=True)
    return dict(nodes=n, edges=e, components=c, cycle_rank=e-n+c,
                structural_roots=g.ids[g.parent == -1].tolist(),
                isolated_nodes=g.ids[g.degree == 0].tolist(), duplicate_node_ids=[],
                self_parent_ids=g.ids[g.parent == g.ids].tolist(),
                repeated_coordinate_groups=int(np.sum(counts > 1)),
                repeated_coordinate_excess=int(np.sum(counts[counts > 1]-1)),
                repeated_coordinates_um=unique_xyz[counts > 1].tolist(),
                zero_length_edge_ids=np.flatnonzero(g.length_um == 0).tolist(),
                nonpositive_radius_ids=g.ids[g.radius_um <= 0].tolist(),
                degree_distribution=dict(zip(map(str, degree), map(int, count))),
                degree1_count=int(np.sum(g.degree == 1)), branch_points=int(np.sum(g.degree > 2)),
                short_edge_threshold_um=0.01, extremely_short_edge_ids=np.flatnonzero(g.length_um < 0.01).tolist(),
                small_radius_threshold_um=0.1, extremely_small_radius_ids=g.ids[g.radius_um < 0.1].tolist(),
                min_edge_length_um=float(g.length_um.min()), min_radius_um=float(g.radius_um.min()),
                max_radius_um=float(g.radius_um.max()), total_length_um=float(g.length_um.sum()),
                bbox_um=[g.xyz_um.min(0).tolist(), g.xyz_um.max(0).tolist()],
                swc_model_limit='One parent per node: cannot encode general anastomotic networks without duplicated nodes/additional edges. Rooted valid SWC is a forest; malformed directed cycles are retained and audited. No inference of biological absence of loops.')
