"""Conservative geometric audit, using a second method on suspect observations.

The archived winding/edge tests remain unchanged. Plane cuts and ray parity
provide independent checks; a near-surface band is never counted as leakage.
"""
import numpy as np
from py_scripts.single_rbc_benchmark.physics import order_vertices, unwrap, geometry
from py_scripts.single_rbc_benchmark.quality import points_inside, mesh_intersections


def checked_mesh(ids, positions, faces, length):
    f = np.asarray(faces, dtype=int)
    if f.ndim != 2 or f.shape[1] != 3 or f.min() < 0 or f.max() >= len(ids):
        raise ValueError('INVALID_CONNECTIVITY')
    if any(len(set(row)) != 3 for row in f):
        raise ValueError('REPEATED_FACE_VERTEX')
    v = unwrap(order_vertices(ids, positions, len(ids)), f, (length, length, length))
    original=order_vertices(ids,positions,len(ids))
    for i,j in ((0,1),(1,2),(2,0)):
        delta=original[f[:,j]]-original[f[:,i]]
        delta[:,:2]-=length*np.rint(delta[:,:2]/length)
        if not np.allclose(v[f[:,j]]-v[f[:,i]],delta,rtol=0,atol=1e-5):
            raise ValueError('PERIODIC_EDGE_CLOSURE_FAILED')
    g = geometry(v, f)
    if not g['closed'] or g['euler_characteristic'] != 2:
        raise ValueError('NOT_CLOSED_ORIENTED_SPHERE_TOPOLOGY')
    return v, g


def require_same_frame(membrane_phase,membrane_step,probe_phase,probe_step):
    if (membrane_phase,int(membrane_step)) != (probe_phase,int(probe_step)):
        raise ValueError('PROBE_MEMBRANE_TIME_MISMATCH')


def surface_distance(points, vertices, faces):
    """Euclidean point/triangle distance: projection if interior, else edges."""
    tri = np.asarray(vertices)[faces]
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    e, f = b-a, c-a
    normal = np.cross(e, f)
    nn = np.einsum('ij,ij->i', normal, normal)
    ee, ff, ef = (np.einsum('ij,ij->i', x, y) for x, y in ((e,e),(f,f),(e,f)))
    den = ee*ff-ef*ef
    if np.any(nn <= 1e-24) or np.any(den <= 1e-24):
        raise ValueError('DEGENERATE_TRIANGLE')
    result, nearest = [], []
    for p in np.asarray(points):
        ap = p-a
        pe = np.einsum('ij,ij->i', ap, e)
        pf = np.einsum('ij,ij->i', ap, f)
        u, v = (ff*pe-ef*pf)/den, (ee*pf-ef*pe)/den
        d2 = np.einsum('ij,ij->i', ap, normal)**2/nn
        d2[(u < 0) | (v < 0) | (u+v > 1)] = np.inf
        for x, y in ((a,b),(b,c),(c,a)):
            edge = y-x
            s = np.clip(np.einsum('ij,ij->i', p-x, edge)/np.einsum('ij,ij->i', edge, edge), 0, 1)
            q = p-x-s[:,None]*edge
            d2 = np.minimum(d2, np.einsum('ij,ij->i', q, q))
        k = int(np.argmin(d2))
        result.append(np.sqrt(max(0, d2[k])))
        nearest.append(k)
    return np.asarray(result), np.asarray(nearest)


def ray_parity(point, vertices, faces, tolerance=1e-9):
    """Three non-axis-aligned rays, solving independent 3x3 systems.

    Hits on edges/vertices make that ray ambiguous. Agreement of the remaining
    directions is required. This does not use the winding formula.
    """
    tri = np.asarray(vertices)[faces]
    votes = []
    for direction in ((1, .371, .193), (.217, 1, .593), (.443, .279, 1)):
        d = np.asarray(direction, float)
        d /= np.linalg.norm(d)
        mat = np.stack((np.broadcast_to(d, (len(tri), 3)), -(tri[:,1]-tri[:,0]), -(tri[:,2]-tri[:,0])), axis=2)
        ok = np.abs(np.linalg.det(mat)) > 1e-12
        sol = np.linalg.solve(mat[ok], (tri[ok,0]-point)[...,None])[...,0]
        t, u, v = sol.T
        possible = (t > tolerance) & (u >= -tolerance) & (v >= -tolerance) & (u+v <= 1+tolerance)
        edge = possible & ((u <= tolerance) | (v <= tolerance) | (u+v >= 1-tolerance))
        if not edge.any():
            votes.append(bool(np.count_nonzero(possible) % 2))
    return votes[0] if len(votes) >= 2 and len(set(votes)) == 1 else None


def membership(points, expected_inside, vertices, faces, band=1e-4):
    winding = points_inside(points, vertices, faces)
    distance, triangle = surface_distance(points, vertices, faces)
    uncertain = distance <= band
    mismatch = (winding != expected_inside) & ~uncertain
    verified = np.zeros(len(points), bool)
    ambiguous = np.zeros(len(points), bool)
    for i in np.flatnonzero(mismatch):
        second = ray_parity(points[i], vertices, faces)
        ambiguous[i] = second is None or second != winding[i]
        verified[i] = second is not None and second == winding[i]
    return dict(winding=winding, distance=distance, nearest_triangle=triangle,
                uncertain=uncertain, raw_mismatch=winding != expected_inside,
                confirmed_mismatch=verified, ambiguous=ambiguous)


def plane_cut_intersection(a, b, tolerance=1e-5):
    """Independent triangle/triangle plane-cut interval overlap.

    Coplanarity and near contact are returned as uncertain, not as penetration.
    A proper intersection needs vertices strictly on both sides of each plane
    and a common segment longer than the explicit spatial tolerance.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    na, nb = np.cross(a[1]-a[0], a[2]-a[0]), np.cross(b[1]-b[0], b[2]-b[0])
    if min(np.linalg.norm(na), np.linalg.norm(nb)) < 1e-12:
        return dict(status='DEGENERATE', length=None)
    na, nb = na/np.linalg.norm(na), nb/np.linalg.norm(nb)
    da, db = (a-b[0])@nb, (b-a[0])@na
    if da.min() > tolerance or da.max() < -tolerance or db.min() > tolerance or db.max() < -tolerance:
        return dict(status='SEPARATE', length=0.)
    direction = np.cross(na, nb)
    if np.linalg.norm(direction) < 1e-8:
        return dict(status='COPLANAR_OR_NEAR_PARALLEL_UNCERTAIN', length=None)
    direction /= np.linalg.norm(direction)

    def cut(t, distances):
        points = [t[i] for i in range(3) if abs(distances[i]) <= tolerance]
        for i, j in ((0,1),(1,2),(2,0)):
            if distances[i]*distances[j] < 0:
                points.append(t[i]+(t[j]-t[i])*distances[i]/(distances[i]-distances[j]))
        return np.asarray(points)

    ca, cb = cut(a, da), cut(b, db)
    if not len(ca) or not len(cb):
        return dict(status='SEPARATE', length=0.)
    aa, bb = ca@direction, cb@direction
    length = float(min(aa.max(), bb.max())-max(aa.min(), bb.min()))
    transverse = da.min() < -tolerance and da.max() > tolerance and db.min() < -tolerance and db.max() > tolerance
    status = 'PROPER_INTERSECTION' if length > tolerance and transverse else ('SEPARATE' if length < -tolerance else 'CONTACT_UNCERTAIN')
    return dict(status=status, length=length, a_cut=ca.tolist(), b_cut=cb.tolist())


def independent_intersections(vertices, faces, tolerance=1e-5):
    legacy = mesh_intersections(vertices, faces)
    # Independently repeat the broad phase; do not use the archived first-12 cap.
    tri = np.asarray(vertices)[faces]
    lo, hi = tri.min(axis=1), tri.max(axis=1)
    ii, jj = np.where(np.triu(np.all(lo[:,None] <= hi[None]+tolerance, axis=2)
                            & np.all(hi[:,None] >= lo[None]-tolerance, axis=2), 1))
    independent = ~np.any(faces[ii,:,None] == faces[jj,None,:], axis=(1,2))
    confirmed, uncertain = [], []
    for i, j in zip(ii[independent], jj[independent]):
        evidence = plane_cut_intersection(tri[i], tri[j], tolerance)
        if evidence['status'] == 'SEPARATE':
            continue
        item = dict(triangles=[int(i),int(j)], vertex_ids=[faces[i].tolist(),faces[j].tolist()],
                    coordinates=[tri[i].tolist(),tri[j].tolist()], **evidence)
        (confirmed if evidence['status'] == 'PROPER_INTERSECTION' else uncertain).append(item)
    return dict(legacy_count=legacy['nonadjacent_intersections'], confirmed_count=len(confirmed),
                uncertain_count=len(uncertain), first_confirmed=confirmed[:1], first_uncertain=uncertain[:1],
                tolerance=tolerance, scope='saved geometry only; all pairs sharing any vertex excluded')


def paired_statistics(samples, confidence_multiplier=1.96):
    """One finite sample mask for both the mean and uncertainty calculation."""
    x = np.asarray(samples, float)
    x = x[np.isfinite(x)]
    if not len(x):
        return dict(n=0, mean=None, ci=None)
    mean = float(x.mean())
    half = confidence_multiplier*float(x.std(ddof=1))/np.sqrt(len(x)) if len(x) > 1 else None
    return dict(n=len(x), mean=mean, ci=None if half is None else [mean-half,mean+half],
                assumption='independent samples; correlated time frames require block analysis')
