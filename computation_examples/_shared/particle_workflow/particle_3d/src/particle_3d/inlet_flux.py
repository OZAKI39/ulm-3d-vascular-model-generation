"""Exact positive-part P1 boundary flux and flux-density inlet sampling (SI)."""
import numpy as np
from .audit import tetra_connectivity


def positive_pieces(vertices, values):
    """Clip a triangle at q=0; retain exact linear q at new polygon vertices."""
    polygon = [(np.asarray(x, float), float(q)) for x, q in zip(vertices, values)]
    clipped = []
    for (a, qa), (b, qb) in zip(polygon, polygon[1:] + polygon[:1]):
        if qa >= 0:
            clipped.append((a, qa))
        if (qa < 0 < qb) or (qb < 0 < qa):
            clipped.append((a + qa / (qa - qb) * (b - a), 0.))
    pieces = []
    for k in range(1, len(clipped)-1):
        x = np.array([clipped[j][0] for j in [0, k, k+1]])
        q = np.array([clipped[j][1] for j in [0, k, k+1]])
        area = np.linalg.norm(np.cross(x[1]-x[0], x[2]-x[0]))/2
        weight = area * q.mean()
        if weight > 0:
            pieces.append((x, q, float(weight)))
    return pieces


class InletFluxSampler:
    def __init__(self, triangles, normal_velocity):
        self.triangles = np.asarray(triangles, float)
        self.q = np.asarray(normal_velocity, float)
        if self.triangles.shape != (*self.q.shape, 3) or self.q.ndim != 2 or self.q.shape[1] != 3:
            raise ValueError('Triangle vertices and nodal normal velocities required')
        if not np.isfinite(self.triangles).all() or not np.isfinite(self.q).all():
            raise ValueError('Finite inlet data required')
        pieces = [(i, x, q, w) for i, (xyz, vals) in enumerate(zip(self.triangles, self.q))
                  for x, q, w in positive_pieces(xyz, vals)]
        self.weights = np.zeros(len(self.triangles))
        for i, _, _, w in pieces:
            self.weights[i] += w
        self.Q_m3_s = float(self.weights.sum())
        self.pieces = pieces
        if not pieces:
            raise ValueError('No positive incoming volume flux')
        self.cdf = np.cumsum([p[3] for p in pieces]); self.cdf /= self.cdf[-1]
        self.cdf[-1] = 1.
        self.proposed = self.accepted = 0

    def sample(self, rng, n=1):
        # Piece selection followed by exact Dirichlet-mixture sampling:
        # q=sum(q_i*b_i), and b_i times uniform triangle density is Dirichlet(2,1,1).
        # This is equivalent to q/qmax rejection, but accepts every proposal.
        idx = np.searchsorted(self.cdf, rng.random(n), side='right')
        xyz = np.array([self.pieces[i][1] for i in idx])
        q = np.array([self.pieces[i][2] for i in idx])
        vertex = (rng.random(n)[:, None] * q.sum(axis=1)[:, None] >= np.cumsum(q, axis=1)).sum(axis=1)
        g = rng.exponential(size=(n, 3))
        g[np.arange(n), vertex] += rng.exponential(size=n)
        bary = g/g.sum(axis=1)[:, None]
        self.proposed += n; self.accepted += n
        return np.einsum('ni,nij->nj', bary, xyz), np.array([self.pieces[i][0] for i in idx])

    def expectation(self):
        # Dirichlet(1,1,1): E[b_i b_j]=(1+delta_ij)/12.
        return sum(w * ((x.sum(axis=0)*q.sum()+q@x)/(4*q.sum())) for _, x, q, w in self.pieces)/self.Q_m3_s


def frozen_boundary_flux(mesh, flow, boundaries):
    wanted = {}
    for role, surface in boundaries.items():
        if role == 'WALL': continue
        faces = surface.faces.reshape(-1, 4)[:, 1:]
        ids = np.asarray(surface.point_data['GlobalNodeID'], int)[faces]-1
        for row in ids: wanted[tuple(sorted(row))] = None
    tets = tetra_connectivity(mesh)
    for tet in tets:
        for omit in range(4):
            key = tuple(sorted(np.delete(tet, omit)))
            if key in wanted: wanted[key] = mesh.points[tet].mean(axis=0)
    result = {}; samplers = {}
    for role, surface in boundaries.items():
        if role == 'WALL': continue
        faces = surface.faces.reshape(-1, 4)[:, 1:]
        ids = np.asarray(surface.point_data['GlobalNodeID'], int)[faces]-1
        xyz = np.asarray(mesh.points, dtype=np.float64)[ids]
        normals = np.cross(xyz[:, 1]-xyz[:, 0], xyz[:, 2]-xyz[:, 0])
        area = np.linalg.norm(normals, axis=1)/2
        normals /= (2*area[:, None])
        centers = np.array([wanted[tuple(sorted(row))] for row in ids], dtype=np.float64)
        normals[np.einsum('ij,ij->i', normals, centers-xyz.mean(axis=1)) > 0] *= -1
        if role == 'INLET': normals *= -1
        q = np.einsum('ijk,ik->ij', np.asarray(flow.point_data['Velocity'])[ids], normals)
        sampler = InletFluxSampler(xyz, q); samplers[role] = sampler
        result[role] = dict(positive_Q_m3_s=sampler.Q_m3_s, signed_Q_m3_s=float(np.sum(area*q.mean(axis=1))),
                            area_m2=float(area.sum()), triangles=len(xyz), normal_direction='INWARD' if role=='INLET' else 'OUTWARD')
    signed_in = result['INLET']['signed_Q_m3_s']
    signed_out = sum(v['signed_Q_m3_s'] for k, v in result.items() if k.startswith('OUTLET'))
    result['mass_balance'] = dict(signed_in_minus_out_m3_s=signed_in-signed_out,
        relative_signed_residual=(signed_in-signed_out)/signed_in, role='FROZEN_FEM_DISCRETIZATION_DIAGNOSTIC_NO_CORRECTION')
    return result, samplers
