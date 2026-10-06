"""D2: two native single-cell curved geometry probes, never a full curved mesh/CFD."""
from common import *
from p2_sampler import *
import os
import subprocess

def main():
    result=read(HERE/'P2_POISEUILLE_REPRESENTATION.json');assert result['D1']=='PASS'
    plan=read(BASE/'production_workflow/gmsh_reference/GMSH_REFERENCE_MESH_V2_PLAN.json')
    R,L=plan['radius_m'],plan['length_m']
    rt=HERE/'runtime'
    raw=np.load(BASE/'production_workflow/manufactured_pipe/R1_GMSH_V2/cfd_gpu/SV_MESH/mesh_arrays.npz')
    affine=np.load(rt/'PIPE_AFFINE_ELEVATION.npz')
    X,C=affine['points_m'],affine['tetra10'];tris=raw['boundary_triangles'];tags=raw['facet_tags']
    edge_tags={}
    for tag in np.unique(tags):
        f=tris[tags==tag]
        for pair in np.sort(f[:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1):
            edge_tags.setdefault(tuple(pair),set()).add(int(tag))
    # Protocol PhysicalGroups: WALL1, INLET2, OUTLET3. Verify actual cap planes.
    for tag,z in [(2,0.),(3,L)]:assert np.max(abs(X[tris[tags==tag],2]-z))<1e-12
    wall_ids=np.flatnonzero(tags==1);wtri=tris[wall_ids];cent=X[wtri].mean(1)
    central=wall_ids[np.argmin(np.linalg.norm(cent-np.array([R,0,L/2]),axis=1))]
    rim_candidates=wall_ids[np.sum(abs(X[wtri,2])<1e-12,axis=1)==2]
    assert len(rim_candidates)
    selected=[('central_wall',int(central)),('inlet_rim',int(rim_candidates[0]))]
    old=read(PREVIOUS/'PRODUCTION_SOURCE_IDENTITY.json')
    cmd=list(old['compile_command']);cmd[cmd.index('-c')+1]=str(HERE/'native_geometry.cpp');cmd[cmd.index('-o')+1]=str(rt/'native_geometry.o')
    cmd=[a.replace(old['source_checkout'],str(RETAINED)) for a in cmd]
    link=list(old['link_command']);link[next(i for i,a in enumerate(link) if a.endswith('/main.cpp.o'))]=str(rt/'native_geometry.o')
    link[link.index('-o')+1]=str(rt/'native_geometry')
    put('NATIVE_GEOMETRY_BUILD.json',dict(compile_command=cmd,link_command=link,cwd=old['link_cwd'],
        source_sha256=sha(HERE/'native_geometry.cpp'),production_objects_reused=True,production_modified=False))
    subprocess.run(cmd,cwd=old['link_cwd'],check=True);subprocess.run(link,cwd=old['link_cwd'],check=True)
    env=dict(os.environ);env['LD_LIBRARY_PATH']=':'.join([str(Path(old['PETSc']['path']).parent),
        '/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3l/external/gpu_mpi_fortran/lib',env.get('LD_LIBRARY_PATH','')])
    native=np.load(HERE/'NATIVE_BASIS_SAMPLES.npz');rows=[];arrays={}
    tol=read(HERE/'PROTOCOL.json')['FP64_roundoff_tolerance']
    for label,face_id in selected:
        tri=tris[face_id]
        incident=np.flatnonzero(np.isin(C[:,:4],tri).sum(1)==3);assert len(incident)==1
        cell=int(incident[0]);ids=C[cell];uncurved=X[ids].copy();curved=uncurved.copy()
        midpoint_tags=[]
        for e,(a,b) in enumerate(EDGES):
            tagset=edge_tags.get(tuple(sorted(ids[[a,b]])),set());midpoint_tags.append(sorted(tagset))
            if 1 in tagset:curved[4+e,:2]*=R/np.linalg.norm(curved[4+e,:2])
            if 2 in tagset:curved[4+e,2]=0.
            if 3 in tagset:curved[4+e,2]=L
        # Sufficient positivity certificate on entire reference tetra:
        # J0^-1*(J(xi)-J0) is affine; spectral norm's convexity bounds it by corner maxima.
        J0=(curved[:3]-curved[3]).T
        _,Dc,_=basis(native['nodes'][:4]);Jc=np.einsum('ai,naj->nij',curved-curved[3],Dc)
        perturb=np.linalg.inv(J0)@(Jc-J0)
        bound=float(np.linalg.norm(perturb,ord=2,axis=(1,2)).max())
        assert bound<1 and np.linalg.det(J0)>0
        # Same nodes/boundaries, affine and projected, to separate mapping defects from fixture setup.
        variants={}
        for variant,nodes in [('affine',uncurved),('curved',curved)]:
            prefix=f'{label}_{variant}';np.savetxt(rt/(prefix+'_nodes.txt'),nodes,fmt='%.17g')
            subprocess.run([str(rt/'native_geometry'),str(rt/(prefix+'_nodes.txt')),str(rt/(prefix+'_native.txt')),
                            str(rt/(prefix+'_native.json'))],env=env,check=True)
            d=np.loadtxt(rt/(prefix+'_native.txt'));summ=read(rt/(prefix+'_native.json'))
            opt,g,w,xi,jv,jq,used=d[:,0].astype(int),d[:,1].astype(int),d[:,2],d[:,3:6],d[:,6],d[:,7],d[:,8]
            N,D,LL=basis(xi);Y=nodes-nodes[3];J=np.einsum('ai,naj->nij',Y,D)
            world=N@Y+nodes[3];pressure_world=LL@nodes[:4]
            native_Dq=d[:,9:21].reshape(-1,4,3)
            common_Dq=np.einsum('ai,nij->naj',DL,np.linalg.inv(J))
            cont=opt==2
            predicted=3*np.einsum('n,ni->i',w[cont]*used[cont],LL[cont])
            common=3*np.einsum('n,ni->i',w[cont]*jv[cont],LL[cont])
            Rp=np.array(summ['native_construct_fluid_Rp'])
            match=float(np.max(abs(Rp-predicted))/np.max(abs(Rp)))
            assert match<tol,('native loop capture mismatch',label,variant,match)
            mismatch=float(np.linalg.norm(Rp-common)/np.linalg.norm(common))
            variants[variant]=dict(native_Rp=Rp.tolist(),native_weight_prediction=predicted.tolist(),
                common_geometry_same_quadrature_Rp=common.tolist(),native_prediction_relative_max_error=match,
                continuity_geometry_weight_relative_L2_mismatch=mismatch,
                max_pressure_coordinate_mismatch_m=float(np.linalg.norm(world-pressure_world,axis=1).max()),
                pressure_gradient_common_mapping_relative_L2=float(np.linalg.norm(native_Dq-common_Dq)/np.linalg.norm(common_Dq)),
                continuity_used_over_velocity_Jacobian=(used[cont]/jv[cont]).tolist(),
                momentum_used_over_velocity_Jacobian=(used[~cont]/jv[~cont]).tolist(),
                min_sampled_Jacobian=float(jv.min()),native_summary=summ)
            for k,v in dict(nodes=nodes,native_rows=d,world=world,pressure_world=pressure_world,native_Rp=Rp,
                            common_Rp=common,native_Dq=native_Dq,common_Dq=common_Dq).items():arrays[prefix+'_'+k]=v
        assert variants['affine']['continuity_geometry_weight_relative_L2_mismatch']<tol
        wallmid=[4+i for i,t in enumerate(midpoint_tags) if 1 in t]
        capmid=[4+i for i,t in enumerate(midpoint_tags) if 2 in t]
        rows.append(dict(label=label,original_cell_id=cell,boundary_face_id=face_id,original_corner_ids=ids[:4].tolist(),
            midpoint_tags=midpoint_tags,corners_unchanged=bool(np.array_equal(uncurved[:4],curved[:4])),
            all_reference_Jacobians_positive_certificate=dict(relative_perturbation_spectral_norm_upper_bound=bound,criterion='<1 with det(J0)>0',passed=True),
            wall_midpoint_radial_error_m=float(np.max(abs(np.linalg.norm(curved[wallmid,:2],axis=1)-R))),
            inlet_midpoint_plane_error_m=float(np.max(abs(curved[capmid,2]))) if capmid else None,
            maximum_midpoint_displacement_m=float(np.linalg.norm(curved-uncurved,axis=1).max()),variants=variants))
    blocked=any(r['variants']['curved']['continuity_geometry_weight_relative_L2_mismatch']>tol or
                r['variants']['curved']['pressure_gradient_common_mapping_relative_L2']>tol for r in rows)
    geometry=dict(status='P2_CURVED_GEOMETRY_NOT_SUPPORTED' if blocked else 'PASS',scope='Required common-geometry mixed semantics not validated; native velocity geometry is high order, not downgraded to linear.',
        fixtures=rows,full_curved_mesh_generated=False,CFD_steps=0,native_construct_fluid_executed=True,
        exact_native_control_flow='Momentum: pressure gnn then velocity gnn. Continuity: velocity gnn then pressure gnn; pressure lShpF=true so only g=0 overwrites Jacobian with affine value; g>0 uses curved velocity Jacobian.',
        interpretation='Mismatch is a deterministic element assembly inconsistency for curved geometry. It does not quantify error of an unrun P2 CFD solution or imply affine Taylor-Hood is unsupported.',
        no_solver_patch=True,harness_binary_sha256=sha(rt/'native_geometry'))
    np.savez_compressed(HERE/'P2_CURVED_GEOMETRY_SAMPLES.npz',**arrays)
    put('P2_CURVED_GEOMETRY_AUDIT.json',geometry)
    result['D2']=geometry['status'];result['curved_geometry_audit']='P2_CURVED_GEOMETRY_AUDIT.json'
    result['status']='AFFINE_REPRESENTATION_PASS_CURVED_IMPLEMENTATION_BLOCKED' if blocked else 'P2_POISEUILLE_REPRESENTATION_PASS'
    put('P2_POISEUILLE_REPRESENTATION.json',result)
    print(geometry['status'],[(r['label'],r['variants']['curved']['continuity_geometry_weight_relative_L2_mismatch']) for r in rows])

if __name__=='__main__':main()
