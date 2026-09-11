"""Unmodified CPU geometry functions extracted from the verified repair; provenance in rbc_geometry_sources.json."""
from pathlib import Path
import numpy as np
from cloud_geometry import order_vertices, geometry

def unwrap(v,faces,lengths,periodic=(True,True,False)):
    v=np.asarray(v,float);out=np.full_like(v,np.nan);seen=np.zeros(len(v),bool);out[0]=v[0];seen[0]=True
    adjacency=[set() for _ in v]
    for face in faces:
        for i,j in zip(face,np.roll(face,-1)):adjacency[i].add(j);adjacency[j].add(i)
    queue=[0];lengths=np.array(lengths,float);mask=np.array(periodic)
    while queue:
        i=queue.pop()
        for j in adjacency[i]:
            delta=v[j]-v[i];delta[mask]-=lengths[mask]*np.rint(delta[mask]/lengths[mask]);point=out[i]+delta
            if seen[j]:
                if not np.allclose(point,out[j],atol=1e-5):raise ValueError('INCONSISTENT_UNWRAPPING_OR_CELL_SPANS_BOX')
            else:out[j]=point;seen[j]=True;queue.append(j)
    if not seen.all():raise ValueError('DISCONNECTED_CELL')
    return out

def _segment_triangle(p,q,t,tol=1e-9):
    direction=q-p;e1=t[:,1]-t[:,0];e2=t[:,2]-t[:,0]
    h=np.cross(direction,e2);det=np.einsum('ij,ij->i',e1,h)
    inv=np.divide(1.,det,out=np.zeros_like(det),where=np.abs(det)>tol)
    s=p-t[:,0];u=np.einsum('ij,ij->i',s,h)*inv
    qq=np.cross(s,e1);v=np.einsum('ij,ij->i',direction,qq)*inv
    along=np.einsum('ij,ij->i',e2,qq)*inv
    return (np.abs(det)>tol)&(u>=-tol)&(v>=-tol)&(u+v<=1+tol)&(along>=-tol)&(along<=1+tol)

def mesh_intersections(vertices,faces):
    """AABB broad phase, nonadjacent triangle edges, coplanar separating axes."""
    tri=np.asarray(vertices)[faces];lo=tri.min(axis=1);hi=tri.max(axis=1)
    overlap=np.all(lo[:,None]<=hi[None,:]+1e-9,axis=2)&np.all(hi[:,None]>=lo[None,:]-1e-9,axis=2)
    i,j=np.where(np.triu(overlap,1))
    independent=~np.any(faces[i,:,None]==faces[j,None,:],axis=(1,2));i=i[independent];j=j[independent]
    a,b=tri[i],tri[j];hit=np.zeros(len(i),bool)
    for k in range(3):
        hit|=_segment_triangle(a[:,k],a[:,(k+1)%3],b)
        hit|=_segment_triangle(b[:,k],b[:,(k+1)%3],a)
    if len(i):
        normal=np.cross(a[:,1]-a[:,0],a[:,2]-a[:,0]);norm=np.linalg.norm(normal,axis=1)
        distance=np.max(np.abs(np.einsum('nkj,nj->nk',b-a[:,:1],normal)),axis=1)
        for k in np.where((distance<=1e-9*np.maximum(norm,1e-30))&~hit)[0]:
            axes=[x for x in range(3) if x!=np.argmax(np.abs(normal[k]))]
            aa=a[k][:,axes];bb=b[k][:,axes];separated=False
            for polygon in (aa,bb):
                for edge in np.roll(polygon,-1,axis=0)-polygon:
                    axis=np.array([-edge[1],edge[0]]);pa=aa@axis;pb=bb@axis
                    if pa.max()<pb.min()-1e-9 or pb.max()<pa.min()-1e-9:separated=True
            hit[k]=not separated
    return dict(nonadjacent_intersections=int(hit.sum()),candidate_pairs=len(i),first_pairs=np.column_stack([i[hit],j[hit]])[:12].tolist(),
                scope='saved frames; pairs sharing a vertex excluded; not a continuous collision proof')

def points_inside(points,vertices,faces):
    """Closed-surface solid-angle winding test, chunked to bound memory."""
    tri=np.asarray(vertices)[faces];answer=[]
    for start in range(0,len(points),32):
        q=tri[None]-np.asarray(points)[start:start+32,None,None,:]
        a,b,c=q[:,:,0],q[:,:,1],q[:,:,2];la=np.linalg.norm(a,axis=2);lb=np.linalg.norm(b,axis=2);lc=np.linalg.norm(c,axis=2)
        numerator=np.einsum('nki,nki->nk',a,np.cross(b,c))
        denominator=la*lb*lc+np.einsum('nki,nki->nk',a,b)*lc+np.einsum('nki,nki->nk',b,c)*la+np.einsum('nki,nki->nk',c,a)*lb
        angle=np.sum(2*np.arctan2(numerator,denominator),axis=1)
        answer.extend(np.abs(angle)>2*np.pi)
    return np.asarray(answer,bool)

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

def fluid_summary(directory,c,reader):
    p=Path(directory)/'profiles.csv'
    if not p.exists() or p.stat().st_size<120:return None
    a=reader(p);shear=a[a['phase']=='shear'];result={}
    if not len(shear):return result
    selected=shear[shear['strain']>=2.-1e-8]
    if not len(selected):selected=shear[shear['step']==shear['step'].max()]
    bins=np.unique(selected['bin']);mean={k:np.array([selected[selected['bin']==b][k].mean() for b in bins]) for k in ('z','ux','uy','uz','density')}
    z=mean['z'];expected=c['protocol']['shear_rate']*(z-c['geometry']['gap']/2)
    rate,intercept=np.polyfit(z,mean['ux'],1)
    mass=[]
    for phase,step in dict.fromkeys(zip(a['phase'].tolist(),a['step'].tolist())):
        q=a[(a['phase']==phase)&(a['step']==step)];mass.append(float(np.average(q['density'],weights=q['count'])) if not (Path(directory)/'moments.csv').exists() else float(q['count'].sum()))
    result.update(mean_profile={k:v.tolist() for k,v in mean.items()},window_strain=[float(selected['strain'].min()),float(selected['strain'].max())],
                  profile_relative_l2=float(np.linalg.norm(mean['ux']-expected)/max(np.linalg.norm(expected),1e-30)),
                  effective_shear=float(rate),effective_shear_relative_error=float(abs(rate/c['protocol']['shear_rate']-1)),
                  fitted_wall_velocities=[float(intercept),float(intercept+rate*c['geometry']['gap'])],
                  near_wall_density_relative_error=float(max(abs(mean['density'][[0,-1]]/(c['dpd']['number_density']*c['dpd']['mass'])-1))),
                  mass_relative_drift=float(max(abs(np.array(mass)/mass[0]-1))))
    mp=Path(directory)/'moments.csv'
    if mp.exists():
        m=reader(mp);ms=m[(m['phase']=='shear')&(m['strain']>=selected['strain'].min()-1e-8)]
        result.update(temperature_mean=float(ms['temperature'].mean()),temperature_relative_error=float(abs(ms['temperature'].mean()/c['dpd']['kBT']-1)),
                      temperature_observed_range=[float(ms['temperature'].min()),float(ms['temperature'].max())],
                      particle_count_range=[int(m['N'].min()),int(m['N'].max())],inner_count_range=[int(m['inner_N'].min()),int(m['inner_N'].max())],
                      wall_crossings_max=int(m['wall_crossings'].max()),max_particle_speed=float(m['max_speed'].max()))
    else:result.update(temperature_status='LBM_WITHOUT_THERMAL_FLUCTUATIONS',wall_crossings_status='NO_PARTICLE_TRACKING_IN_LBM')
    return result

def local_flow_summary(directory,reader,window=(2.,4.)):
    p=Path(directory)/'local_flow.csv'
    if not p.exists():return None
    a=reader(p)
    if not len(a):return None
    a=a[(a['phase']=='shear')&(a['strain']>=window[0]-1e-8)&(a['strain']<=window[1]+1e-8)]
    if not len(a):return None
    keys=list(dict.fromkeys(zip(a['ix'].tolist(),a['iy'].tolist(),a['iz'].tolist())))
    rows=[]
    for ix,iy,iz in keys:
        r=a[(a['ix']==ix)&(a['iy']==iy)&(a['iz']==iz)]
        blocks=[]
        for lo,hi in zip(np.linspace(window[0],window[1],5)[:-1],np.linspace(window[0],window[1],5)[1:]):
            q=r[(r['strain']>=lo)&(r['strain']<(hi if hi<window[1] else hi+1e-8))]
            if len(q):blocks.append([float(q[k].mean()) for k in ('ux','uy','uz')])
        rows.append(dict(cell=[int(ix),int(iy),int(iz)],position=[float(r[k][0]) for k in ('x','y','z')],velocity=[float(r[k].mean()) for k in ('ux','uy','uz')],
                         block_std=np.std(blocks,axis=0,ddof=1).tolist() if len(blocks)>1 else None))
    return dict(window=[float(a['strain'].min()),float(a['strain'].max())],requested_window=list(window),cells=rows,definition='Time mean of native local velocities in shared 8^3 bins; HemoCell node average and DPD particle average differ. Block scatter is descriptive, not a precise confidence interval.')

