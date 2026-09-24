"""Derived planar CFD port contract. Geometry/topology only; no FEM spaces."""
import numpy as np
from .cap_remesh import check_rim, project, cross2, coverage_check
from .mesh_qc import triangle_geometry


def polygon_metrics(rim, origin, basis):
    """Translation-stable contour formulas, longdouble intermediates."""
    d=np.asarray(rim,dtype=np.longdouble)-np.asarray(origin,dtype=np.longdouble)
    b=np.asarray(basis,dtype=np.longdouble)
    xy=d@b[:2].T
    c=cross2(xy,np.roll(xy,-1,axis=0)); area=c.sum()/2
    if area<=0 or not np.isfinite(area): raise ValueError('Reversed or degenerate projected polygon')
    centroid=((xy+np.roll(xy,-1,axis=0))*c[:,None]).sum(axis=0)/(6*area)
    vector=np.cross(d,np.roll(d,-1,axis=0)).sum(axis=0)/2
    return {'formal_projected_area_m2':float(area),'projected_centroid_xy_m':np.asarray(centroid,float).tolist(),
            'projected_centroid_m':np.asarray(np.asarray(origin,dtype=np.longdouble)+centroid@b[:2],float).tolist(),
            'vector_area_m2':np.asarray(vector,float).tolist()}


def verify_plane(port, origin, normal, basis):
    for key,a,b in [('origin',port['plane_origin_m'],origin),('normal',port['outward_normal'],normal),('basis',port['basis'],basis)]:
        if not np.array_equal(a,b): raise ValueError('Frozen plane '+key+' changed')


def verify_metrics(reference, measured, geometry):
    ar=reference['formal_projected_area_m2']; av=np.asarray(reference['vector_area_m2']); v=np.asarray(measured['vector_area_m2'])
    relative_area=abs(measured['formal_projected_area_m2']-ar)/ar
    relative_vector=float(np.linalg.norm(v-av)/np.linalg.norm(av))
    n=np.asarray(reference['outward_normal']); n=n/np.linalg.norm(n)
    dot=float(v@n/np.linalg.norm(v))
    shift=float(np.linalg.norm(np.asarray(measured['projected_centroid_m'])-reference['projected_centroid_m']))
    if not np.isfinite([relative_area,relative_vector,dot,shift]).all(): raise ValueError('Nonfinite port metric')
    if relative_area>geometry['projected_area_relative_error_max']: raise ValueError('Projected area changed')
    if relative_vector>geometry['vector_area_relative_error_max']: raise ValueError('Vector area changed or reversed')
    if dot<geometry['normal_dot_min']: raise ValueError('Vector normal reversed or changed')
    if shift>geometry['projected_centroid_displacement_m_max']: raise ValueError('Projected centroid changed')
    return {'relative_projected_area_error':relative_area,'relative_vector_area_error':relative_vector,'normal_dot_source':dot,'projected_centroid_displacement_m':shift}


def validate_port(source_points, points, triangles, port, interior_ids, policy, *, origin, normal, basis):
    verify_plane(port,origin,normal,basis)
    rim=check_rim(source_points,port['rim_edge_ids'],points,triangles)
    ids=np.asarray(port['rim_vertex_ids'],int)
    measured=polygon_metrics(points[ids],origin,basis)
    errors=verify_metrics(port,measured,policy['geometry'])
    # Independent implementation B: sum signed areas of ACTUAL projected triangles.
    xy,offset=project(points,origin,np.asarray(basis)); t=xy[triangles]
    weights=cross2(t[:,1]-t[:,0],t[:,2]-t[:,0])/2
    area_b=float(weights.sum()); error_b=abs(area_b-port['formal_projected_area_m2'])/port['formal_projected_area_m2']
    if error_b>policy['geometry']['projected_area_relative_error_max']: raise ValueError('Projected triangle area differs from polygon')
    centroid_b=np.sum(t.mean(axis=1)*weights[:,None],axis=0)/area_b
    cerror=float(np.linalg.norm(centroid_b-measured['projected_centroid_xy_m']))
    if cerror>policy['geometry']['projected_centroid_displacement_m_max']: raise ValueError('Projected triangle centroid mismatch')
    coverage=coverage_check(xy,triangles,ids,policy['geometry']['coverage_relative_roundoff_tolerance'])
    interior=np.asarray(interior_ids,dtype=int)
    deviation=float(np.abs(offset[interior]).max(initial=0))
    if deviation>policy['geometry']['interior_plane_deviation_m_max']: raise ValueError('Interior vertex off reference plane')
    a,c,v=triangle_geometry(points,triangles)
    if not np.isfinite(a).all() or (a<=0).any(): raise ValueError('Nonpositive scalar triangle area')
    scalar=float(a.sum()); legacy=port['legacy_scalar_area_m2']
    return {'status':'PASS',**measured,**errors,'rim':rim,'coverage':coverage,
            'A_vector_source':port['vector_area_m2'],'A_vector_candidate':measured['vector_area_m2'],
            'projected_triangle_sum_area_m2':area_b,'independent_projected_area_relative_error':error_b,
            'independent_projected_centroid_error_m':cerror,'new_interior_max_plane_deviation_m':deviation,
            'candidate_scalar_triangle_area_m2':scalar,'legacy_scalar_area_m2':legacy,
            'candidate_scalar_vs_legacy_relative_difference':(scalar-legacy)/legacy,
            'candidate_scalar_area_centroid_m':(np.sum(a[:,None]*c,axis=0)/scalar).tolist(),
            'legacy_scalar_area_is_gate':False}


def p2_proxy(tetra):
    """Count vertices and undirected edges, independent of element numbering."""
    t=np.asarray(tetra)
    if t.ndim!=2 or t.shape[1]!=4 or (np.diff(np.sort(t,axis=1),axis=1)==0).any(): raise ValueError('Invalid tetra topology')
    edges=np.concatenate([t[:,[i,j]] for i,j in ((0,1),(0,2),(0,3),(1,2),(1,3),(2,3))])
    nv=len(np.unique(t)); ne=len(np.unique(np.sort(edges,axis=1),axis=0))
    return {'N_vertex':nv,'N_edge':ne,'N_tetra':len(t),'N_P2_scalar_proxy':nv+ne,'N_P2_velocity_proxy':3*(nv+ne),'N_P1_pressure_proxy':nv}


def cap_gates(counts, q, policy):
    q=np.asarray(q); g=policy['cap_quality']; budgets=policy['density']
    checks={'finite_positive_quality':bool(np.isfinite(q).all() and (q>0).all()),
            'zero_q_below_0_1':bool(np.count_nonzero(q<.1)==0),
            'P5':bool(np.quantile(q,.05)>=g['P5_min']),'median':bool(np.median(q)>=g['median_min']),
            'per_port_density':all(counts[k]<=v for k,v in budgets['per_port_max'].items()),
            'total_density':bool(sum(counts.values())<=budgets['total_max'])}
    return {'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks}


def volume_gates(validity, quality, proxy, policy):
    g=policy['tetra_quality']; b=policy['cost']; q=quality['min_sicn']
    checks={'valid_tetra':all(v==0 for v in validity.values()),
            'cap_adjacent_low':quality['cap_adjacent_below_0_1']<=g['cap_adjacent_below_0_1_max'],
            'total_low':quality['total_below_0_1']<=g['total_below_0_1_max'],
            **{k:q[k]>=g[k+'_min'] for k in ('P1','P5','median')},
            'tetra_budget':proxy['N_tetra']<=b['tetra_max'],
            'P2_budget':proxy['N_P2_velocity_proxy']<=b['P2_velocity_ratio_max']*b['baseline']['N_P2_velocity_proxy']}
    checks={k:bool(v) for k,v in checks.items()}
    return {'status':'PASS' if all(checks.values()) else 'FAIL','checks':checks}


def select_candidate(records):
    eligible=[r for r in records if r['surface_status']=='PASS' and r.get('volume_status')=='PASS' and all(r['volume_gate_checks'].values())]
    key=lambda r:(r['proxy']['N_P2_velocity_proxy'],r['proxy']['N_tetra'],r['quality']['cap_adjacent_below_0_1'],r['quality']['total_below_0_1'],-r['quality']['min_sicn']['P1'])
    ordered=sorted(eligible,key=key)
    return {'selected_candidate':ordered[0]['candidate'] if ordered else None,'survivors':[r['candidate'] for r in ordered],
            'sort_keys':['N_P2_velocity_proxy','N_tetra','cap_adjacent_below_0_1','total_below_0_1','-P1'],
            'status':'CONDITIONAL PASS' if ordered else 'FAIL','records':records}
