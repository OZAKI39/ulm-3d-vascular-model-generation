"""Phase B: native basis polynomial reproduction and exact original WALL-edge topology."""
from common import *
import numpy as np
import pyvista as pv
import itertools

def derived(u,G,n=np.array([0.,0.,1.])):
    un=u@n
    Gn=G@n
    return dict(u=u,G=G,normal=un,tangent=u-un[:,None]*n,div=np.trace(G,axis1=1,axis2=2),
        curl=np.stack([G[:,2,1]-G[:,1,2],G[:,0,2]-G[:,2,0],G[:,1,0]-G[:,0,1]],axis=1),
        strain=(G+G.transpose(0,2,1))/2,s=Gn-(Gn@n)[:,None]*n)

def metrics(approx,exact):
    return {k:dict(max_abs_error=float(np.max(np.abs(approx[k]-v))),
                   scaled_max_error=float(np.max(np.abs(approx[k]-v))/max(1.,np.max(np.abs(v)))),
                   L2_error=float(np.linalg.norm(approx[k]-v))) for k,v in exact.items()}

def main():
    assert read(HERE/'NATIVE_P2P1_PATH_AUDIT.json')['status']=='NATIVE_P2P1_PATH_CONFIRMED'
    protocol=read(HERE/'PROTOCOL.json');tol=protocol['FP64_roundoff_tolerance']
    native=np.load(HERE/'NATIVE_BASIS_SAMPLES.npz')
    xi,N,D,P,DP,nodes,edges=(native[k] for k in ['xi','N','dN','P','dP','nodes','edges'])
    rng=np.random.default_rng(protocol['rng_seed'])
    out={};arrays={k:native[k] for k in native.files}
    # B1 uses a full quadratic with nonzero mixed and diagonal terms in all components.
    corners=np.array([[1.3,.1,.2],[.2,1.2,.3],[.1,.2,1.4],[-.2,.05,-.1]])
    J=(corners[:3]-corners[3]).T
    X=np.vstack([corners,corners[edges].mean(1)])
    points=P@corners
    c=rng.normal(size=3);L=rng.normal(size=(3,3));H=rng.normal(size=(3,3,3));H=(H+H.transpose(0,2,1))/2
    def exact(x):
        return derived(c+x@L.T+.5*np.einsum('ni,kij,nj->nk',x,H,x),L[None]+np.einsum('kij,nj->nki',H,x))
    U=exact(X)['u'];reference=exact(points)
    p2=derived(N@U,np.einsum('ai,naj->nij',U,D@np.linalg.inv(J)))
    p1=derived(P@U[:4],np.einsum('ai,naj->nij',U[:4],DP@np.linalg.inv(J)))
    m2,m1=metrics(p2,reference),metrics(p1,reference)
    out['B1']=dict(status='PASS' if max(v['scaled_max_error'] for v in m2.values())<tol else 'FAIL',
        points=len(points),random_internal_points=len(points)-10,units='dimensionless',P2=m2,P1=m1,
        geometry_det=float(np.linalg.det(J)),geometry_condition=float(np.linalg.cond(J)),
        coefficients=dict(constant=c.tolist(),linear=L.tolist(),symmetric_hessian=H.tolist()))
    arrays.update(B1_nodes=X,B1_points=points)
    for label,data in [('exact',reference),('P1',p1),('P2',p2)]:
        arrays.update({f'B1_{label}_{k}':v for k,v in data.items()})
    # B2 uses the native reference tetra: z=0 is the wall face [0,1,3].
    def near(x):
        a,b,z=x.T;u=np.column_stack([z*(1+a),z*(2+b),z*z]);g=np.zeros((len(x),3,3))
        g[:,0,0]=z;g[:,0,2]=1+a;g[:,1,1]=z;g[:,1,2]=2+b;g[:,2,2]=2*z
        return derived(u,g)
    U=near(nodes)['u']; reference=near(xi)
    p2=derived(N@U,np.einsum('ai,naj->nij',U,D))
    p1=derived(P@U[:4],np.einsum('ai,naj->nij',U[:4],DP))
    ranks1=np.linalg.matrix_rank(p1['G'][10:],tol=tol);ranks2=np.linalg.matrix_rank(p2['G'][10:],tol=tol)
    m2,m1=metrics(p2,reference),metrics(p1,reference)
    out['B2']=dict(status='PASS' if max(v['scaled_max_error'] for v in m2.values())<tol and np.all(ranks1==1) and np.all(ranks2==3) else 'FAIL',
        analytic_field=protocol['B2_analytic_field'],units='dimensionless',wall_normal_inward=[0,0,1],
        wall_local_nodes=np.flatnonzero(nodes[:,2]==0).tolist(),free_local_nodes=np.flatnonzero(nodes[:,2]!=0).tolist(),
        P1_gradient_ranks=np.unique(ranks1).tolist(),P2_gradient_ranks=np.unique(ranks2).tolist(),P2=m2,P1=m1,
        limitation='Representation capability only; this manufactured field is not divergence-free or a CFD solution.')
    arrays.update(B2_points=xi,B2_nodes=nodes,B2_node_u=U,B2_P1_rank=ranks1,B2_P2_rank=ranks2)
    for label,data in [('exact',reference),('P1',p1),('P2',p2)]:
        arrays.update({f'B2_{label}_{k}':v for k,v in data.items()})
    # B3: actual original input tetra numbering and GlobalNodeID-based WALL incidence.
    mp=BASE/'mesh/SV_MESH/mesh-complete.mesh.vtu';wp=BASE/'mesh/SV_MESH/mesh-surfaces/WALL.vtp'
    mesh=pv.read(mp);wall=pv.read(wp)
    assert np.all(mesh.celltypes==10)
    cells=mesh.cells.reshape(-1,5)[:,1:]
    f=wall.faces.reshape(-1,4);assert np.all(f[:,0]==3)
    gids=np.asarray(wall.point_data['GlobalNodeID'],dtype=np.int64)-1
    assert len(np.unique(gids))==len(gids)
    mismatch=float(np.max(np.abs(np.asarray(wall.points)-np.asarray(mesh.points)[gids])))
    assert mismatch<1e-14,mismatch
    tri=gids[f[:,1:]]
    for cell_id in [138584,388769]:
        ids=cells[cell_id];xyz=np.asarray(mesh.points,dtype=np.float64)[ids]
        rows=[]
        for local,(a,b) in enumerate(edges):
            pair=sorted([int(ids[a]),int(ids[b])])
            face_ids=np.flatnonzero(np.any(tri==pair[0],axis=1)&np.any(tri==pair[1],axis=1))
            rows.append(dict(local_midpoint_node=4+local,local_edge=[int(a),int(b)],
                mesh_point_ids_0based=pair,GlobalNodeIDs_1based=[v+1 for v in pair],
                midpoint_m=((xyz[a]+xyz[b])/2).tolist(),wall_triangle_ids_0based=face_ids.tolist(),
                wall_triangle_count=len(face_ids),fixed_velocity=bool(len(face_ids))))
        wall_vertices=np.isin(ids,gids)
        faces_on_wall=[]
        for face in itertools.combinations(range(4),3):
            face_ids=np.flatnonzero(np.all(np.isin(tri,ids[list(face)]),axis=1))
            if len(face_ids):faces_on_wall.append(dict(local_vertices=list(face),wall_triangle_ids_0based=face_ids.tolist()))
        free=4-int(wall_vertices.sum())+sum(not r['fixed_velocity'] for r in rows)
        info=dict(cell_id_original_0based=cell_id,point_ids_0based=ids.tolist(),GlobalNodeIDs_1based=(ids+1).tolist(),
            vertex_xyz_m=xyz.tolist(),centroid_m=xyz.mean(0).tolist(),wall_vertices=wall_vertices.tolist(),
            wall_faces=faces_on_wall,edges=rows,free_P1_velocity_nodes=int(4-wall_vertices.sum()),
            free_P2_velocity_nodes=free,free_scalar_velocity_DOF=3*free,
            P2_identically_zero_under_wall_BC=(free==0))
        out[f'real_{cell_id}']=info
        arrays[f'real_{cell_id}_nodes']=np.vstack([xyz,xyz[edges].mean(1)])
        arrays[f'real_{cell_id}_fixed_nodes']=np.r_[wall_vertices,[r['fixed_velocity'] for r in rows]]
        near_ids=np.flatnonzero(np.any(np.isin(tri,ids),axis=1))
        arrays[f'real_{cell_id}_wall_triangles_xyz_m']=np.asarray(mesh.points)[tri[near_ids]]
    real=out['real_138584']
    free_edges=[r['local_midpoint_node'] for r in real['edges'] if not r['fixed_velocity']]
    test_velocity=np.zeros((10,3))
    for node in free_edges:test_velocity[node]=[1.,2.,3.]
    xyz=np.asarray(real['vertex_xyz_m'])
    Gtest=np.einsum('ai,naj->nij',test_velocity,D@np.linalg.inv((xyz[:3]-xyz[3]).T))
    arrays.update(B3_unit_free_node_values=test_velocity,B3_unit_u=N@test_velocity,B3_unit_G=Gtest)
    real['unit_free_DOF_fixture']=dict(nodal_velocity=test_velocity.tolist(),
        max_speed=float(np.linalg.norm(N@test_velocity,axis=1).max()),
        nonzero_internal_samples=int(np.count_nonzero(np.linalg.norm((N@test_velocity)[10:],axis=1))),
        gradient_ranks=np.unique(np.linalg.matrix_rank(Gtest[10:],tol=1e-8)).tolist(),
        interpretation='Arbitrary unit coefficient illustrates admissible shape, NOT CFD or physical velocity. With only one free node, G remains rank <=1.')
    out['B3']=dict(status='TOPOLOGY_RESOLVED',mesh_sha256=sha(mp),wall_sha256=sha(wp),
        global_id_coordinate_max_mismatch_m=mismatch,real_target_free_P2_velocity_nodes=out['real_138584']['free_P2_velocity_nodes'],
        boundary_rule=protocol['B3_boundary_rule'],native_BC_reference='set_bc.cpp:1031-1144; face.gN includes all inherited TRI6 nodes',
        invented_interior_nodes=False,full_mesh_elevated=False)
    stop=out['B1']['status']!='PASS' or out['B2']['status']!='PASS' or out['real_138584']['free_P2_velocity_nodes']==0
    out.update(status='P2_REPRESENTATION_HYPOTHESIS_NOT_SUPPORTED' if stop else 'P2_REPRESENTATION_PASS',
        stop_gate='B' if stop else None,tolerance=tol,CFD_runs=0,
        distinction='A topology Stop B does not disprove quadratic approximation or the B2 rank improvement; it blocks this user-defined target route.')
    put('P2_REPRESENTATION_FIXTURES.json',out)
    np.savez_compressed(HERE/'P2_REPRESENTATION_FIXTURES.npz',**arrays)
    print(out['status'],'B1',out['B1']['status'],'B2',out['B2']['status'],'138584 free nodes',out['real_138584']['free_P2_velocity_nodes'])
    print('138584 edges',[(r['mesh_point_ids_0based'],r['wall_triangle_count']) for r in out['real_138584']['edges']])

if __name__=='__main__':main()
