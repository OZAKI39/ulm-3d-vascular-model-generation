"""Phase D1: affine field-only elevation of frozen V2, no mesher or CFD."""
from common import *
from p2_sampler import *
from representation import derived

def main():
    assert read(HERE/'P2_FIELD_SAMPLER_TESTS.json')['status']=='PASS'
    root=BASE/'production_workflow';case=root/'manufactured_pipe/R1_GMSH_V2/cfd_gpu'
    plan=read(root/'gmsh_reference/GMSH_REFERENCE_MESH_V2_PLAN.json');proto=read(root/'manufactured_pipe/PROTOCOL.json')
    R,L,U=plan['radius_m'],plan['length_m'],proto['mean_velocity_m_s']
    raw=np.load(case/'SV_MESH/mesh_arrays.npz');old=np.load(case/'analysis/FIELD_SAMPLES.npz')
    corners=raw['points_m'];tets=raw['tetra'].copy()
    # Exactly the native reader orientation, leaving corner coordinates and cells unchanged.
    det=np.linalg.det((corners[tets[:,:3]]-corners[tets[:,3,None]]).transpose(0,2,1))
    negative=det<0
    tets[negative,0],tets[negative,1]=tets[negative,1].copy(),tets[negative,0].copy()
    X,C,edges=elevate_affine(corners,tets)
    def exact(points):
        u=np.zeros_like(points);u[:,2]=2*U*(1-(points[:,0]**2+points[:,1]**2)/R**2)
        G=np.zeros((len(points),3,3));G[:,2,:2]=-4*U/R**2*points[:,:2]
        return u,G
    nodal,_=exact(X);p=-8*proto['mu_Pa_s']*U/R**2*(X[:,2]-L)
    sampler=P2Sampler(X,C,nodal,p)
    rng=np.random.default_rng(read(HERE/'PROTOCOL.json')['rng_seed'])
    random_ids=rng.integers(0,len(C),4096)
    bary=rng.dirichlet(np.ones(4),len(random_ids))
    random=np.einsum('ni,nij->nj',bary,corners[tets[random_ids]])
    query=np.vstack([old['points_m'],random])
    sampled=sampler.sample(query)
    assert sampled['inside_lumen'].all(),'strict containment failed; no probes moved or hidden'
    ue,Ge=exact(query)
    ids=sampled['element_id'];cellcorners=corners[tets[ids]]
    J=(cellcorners[:,:3]-cellcorners[:,3,None]).transpose(0,2,1)
    inv=np.linalg.inv(J);xi=np.einsum('nij,nj->ni',inv,query-cellcorners[:,3])
    weights=np.c_[xi,1-xi.sum(1)]
    up=np.einsum('ni,nij->nj',weights,nodal[tets[ids]])
    gp=np.einsum('nai,naj->nij',nodal[tets[ids]],np.einsum('ai,nij->naj',DL,inv))
    radial=np.linalg.norm(query[:,:2],axis=1)
    normal=-np.c_[query[:,:2]/radial[:,None],np.zeros(len(query))]
    def fields(u,G):
        gn=np.einsum('nij,nj->ni',G,normal)
        return dict(u=u,G=G,curl=np.stack([G[:,2,1]-G[:,1,2],G[:,0,2]-G[:,2,0],G[:,1,0]-G[:,0,1]],axis=1),
                    s=gn-(gn*normal).sum(1)[:,None]*normal)
    analytic=fields(ue,Ge);P1=fields(up,gp);P2=fields(sampled['u'],sampled['G'])
    groups=[];original=read(root/'FIELD_VALIDATION.json')['groups']
    for g,info in enumerate(original):
        mask=old['groups']==g
        row=dict(wall_distance_m=info['wall_distance_m'],role=info['role'],samples=int(mask.sum()),P1={},P2={})
        for name,data in [('P1',P1),('P2',P2)]:
            for k,expected in analytic.items():
                err=data[k][:960][mask]-expected[:960][mask]
                row[name][k]=dict(relative_L2=float(np.linalg.norm(err)/np.linalg.norm(expected[:960][mask])),max_abs=float(abs(err).max()))
        groups.append(row)
    # Physical scales declared from protocol, independent of vanishing near-wall u.
    scales={'u':2*U,'G':4*U/R,'s':4*U/R,'curl':4*U/R}
    error={k:float(abs(P2[k]-analytic[k]).max()/scales[k]) for k in analytic}
    tol=read(HERE/'PROTOCOL.json')['FP64_roundoff_tolerance']
    expected_p=-8*proto['mu_Pa_s']*U/R**2*(query[:,2]-L)
    assert max(error.values())<tol,error
    assert np.max(abs(sampled['pressure']-expected_p))<tol*max(abs(expected_p).max(),1.)
    # Deterministic global edge IDs, boundary midpoint inheritance, unchanged corners.
    boundary=raw['boundary_triangles'];tags=raw['facet_tags']
    lut={tuple(v):len(corners)+i for i,v in enumerate(edges)}
    bedges=np.stack([boundary[:,[0,1]],boundary[:,[1,2]],boundary[:,[2,0]]],axis=1)
    mids=np.array([lut[tuple(sorted(e))] for e in bedges.reshape(-1,2)]).reshape(-1,3)
    boundary6=np.column_stack([boundary,mids])
    assert np.array_equal(X[:len(corners)],corners)
    det2=np.linalg.det((X[C[:,:3]]-X[C[:,3,None]]).transpose(0,2,1))
    assert np.all(det2>0)
    rt=HERE/'runtime';rt.mkdir(exist_ok=True)
    np.savez_compressed(rt/'PIPE_AFFINE_ELEVATION.npz',points_m=X,tetra10=C,unique_edges=edges,boundary6=boundary6,facet_tags=tags)
    arrays=dict(points_m=query,cell=ids,groups=np.r_[old['groups'],np.full(4096,-1)],normal_inward=normal,inside=sampled['inside_lumen'])
    for label,data in [('exact',analytic),('P1',P1),('P2',P2)]:arrays.update({label+'_'+k:v for k,v in data.items()})
    np.savez_compressed(HERE/'P2_POISEUILLE_SAMPLES.npz',**arrays)
    inputs=read(HERE/'INPUT_SHA256.json')
    for path in [case/'SV_MESH/mesh_arrays.npz',case/'analysis/FIELD_SAMPLES.npz',root/'gmsh_reference/GMSH_REFERENCE_MESH_V2_PLAN.json',root/'manufactured_pipe/PROTOCOL.json']:
        inputs[str(path.relative_to(BASE))]=sha(path)
    put('INPUT_SHA256.json',inputs)
    result=dict(status='AFFINE_P2_REPRESENTATION_PASS',D1='PASS',D2='PENDING',new_CFD=0,
        R_m=R,L_m=L,Umean_m_s=U,rho=proto['rho_kg_m3'],mu=proto['mu_Pa_s'],frozen_probes=960,random_points=4096,
        invalid=0,groups=groups,physical_scale_max_errors=error,tolerance=tol,
        tolerance_definition='FP64 max component error / (2 Umean for u; 4 Umean/R for derivatives), predeclared 1e-11 threshold',
        mesh=dict(original_corners=len(corners),cells=len(C),unique_edges=len(edges),velocity_nodes=len(X),pressure_nodes=len(corners),
            independent_unknowns=3*len(X)+len(corners),stored_slots=4*len(X),native_orientation_swaps=int(negative.sum()),
            corners_unchanged=True,cell_order_unchanged=True,all_edge_midpoints_arithmetic=True,geometry_affine=True,
            min_positive_Jacobian=float(det2.min()),boundary_faces=len(boundary6),boundary_tags=np.unique(tags).tolist(),
            full_curved_mesh_generated=False,elevation_artifact='runtime/PIPE_AFFINE_ELEVATION.npz',elevation_sha256=sha(rt/'PIPE_AFFINE_ELEVATION.npz')),
        caveat='Exact nodal field representation only; not actual CFD, not wall-Dirichlet curved-geometry certification.')
    put('P2_POISEUILLE_REPRESENTATION.json',result)
    print('D1 PASS',error,'nodes',len(X),'cells',len(C))

if __name__=='__main__':main()
