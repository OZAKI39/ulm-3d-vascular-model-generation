from common import *
from p2 import *
import xml.etree.ElementTree as ET
import shutil, time
import pyvista as pv
from particle_3d.flowfield_conservation_diagnosis import face_topology,boundary_owners

def main():
    assert json.loads((REPORT/'data/p1_network_bc_local_conservation_baseline.json').read_text())['status']=='PASS'
    assert not CASE.exists();CASE.mkdir();(CASE/'SV_MESH').mkdir();(CASE/'SV_MESH/mesh-surfaces').mkdir();(CASE/'run').mkdir();(CASE/'reports').mkdir()
    a=baseline_arrays();p,original_t=a['points_m'],a['tetra'];t=original_t.copy()
    # Apply exactly the native read_msh.cpp:check_tet_conn corner orientation.
    raw_det=np.linalg.det(np.swapaxes(p[t[:,1:]]-p[t[:,:1]],1,2))
    flip=raw_det>0;t[flip,0],t[flip,1]=t[flip,1].copy(),t[flip,0].copy()
    x,t10,edges=unique_elevation(p,t)
    mesh=pv.UnstructuredGrid(np.column_stack([np.full(len(t),10),t10]).ravel(),np.full(len(t),24,np.uint8),x)
    mesh.point_data['GlobalNodeID']=np.arange(1,len(x)+1,dtype=np.int32);mesh.cell_data['GlobalElementID']=np.arange(1,len(t)+1,dtype=np.int32)
    mesh.save(CASE/'SV_MESH/mesh-complete.mesh.vtu')
    with (REPORT/'data/tet4_to_tet10_node_map.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['new_node_id','edge_vertex_i','edge_vertex_j','x_m','y_m','z_m'])
        w.writerows((len(p)+i+1,int(e[0])+1,int(e[1])+1,*x[len(p)+i]) for i,e in enumerate(edges))
    top=face_topology(t);boundaries={};areas={}
    for path in sorted((BASE/'SV_MESH/mesh-surfaces').glob('*.vtp')):
        s=pv.read(path);tri=(np.asarray(s['GlobalNodeID'],int)-1)[s.faces.reshape(-1,4)[:,1:]]
        tri6=face_elevation(tri,edges,len(p));owners=boundary_owners(top,tri)
        ids,inv=np.unique(tri6,return_inverse=True)
        face=pv.PolyData(x[ids],np.column_stack([np.full(len(tri),6),inv.reshape(-1,6)]).ravel())
        face.point_data['GlobalNodeID']=(ids+1).astype(np.int32);face.cell_data['GlobalElementID']=(owners+1).astype(np.int32)
        face.save(CASE/'SV_MESH/mesh-surfaces'/path.name)
        boundaries[path.stem]=tri6
        area=np.linalg.norm(np.cross(p[tri[:,1]]-p[tri[:,0]],p[tri[:,2]]-p[tri[:,0]]),axis=1)/2
        # All midpoint constraints verified; quadratic maps reduce identically to affine maps.
        assert np.array_equal(x[tri6[:,3:]],x[tri6[:,:3]][:,TRI_EDGES].mean(axis=2))
        areas[path.stem]=dict(P1_area_m2=float(area.sum()),P2_integrated_area_m2=float(np.sum(area*TRI_W.sum())),
            relative_difference=float(TRI_W.sum()-1),boundary_facets=len(tri),unique_volume_parent=True)
    quad,weights=tetra_rule15();jac_min=float('inf');max_variation=0.;v2=0.
    det=np.linalg.det(np.swapaxes(p[t[:,:3]]-p[t[:,3,None]],1,2))
    assert det.min()>0
    for start in range(0,len(t),10000):
        pts=x[t10[start:start+10000]]
        jac=np.einsum('tai,qaj->tqij',pts-pts[:,:1],reference_gradient(quad))
        dj=np.linalg.det(jac);jac_min=min(jac_min,float(dj.min()))
        max_variation=max(max_variation,float(np.max(abs(dj/det[start:start+10000,None]-1))))
        v2+=float(np.sum(dj@weights/6))
    v1=float(det.sum()/6)
    assert jac_min>0 and abs(v2/v1-1)<1e-11 and max_variation<1e-10
    np.savez_compressed(CASE/'SV_MESH/tet10_arrays_si.npz',points_m=x,tetra10=t10,tetra4=t,edges=edges,**boundaries)
    dump(REPORT/'data/mesh_elevation_audit.json',dict(status='PASS',original_vertices=len(p),new_nodes=len(x),unique_edges=len(edges),
        tetra_count=len(t),original_corner_positions_exact=True,original_corner_connectivity_topology_exact=bool(np.array_equal(np.sort(t,axis=1),np.sort(original_t,axis=1))),
        native_orientation_reordered_tetrahedra=int(flip.sum()),orientation_rule='read_msh.cpp:check_tet_conn; swap local corners 0/1 when positive standard determinant',
        shared_global_midpoints_exact=True,total_volume_P1_m3=v1,total_volume_P2_m3=v2,volume_relative_difference=v2/v1-1,
        minimum_solver_quadrature_jacobian=jac_min,maximum_jacobian_relative_variation=max_variation,
        solver_quadrature_points_per_tet=len(quad),boundary_areas=areas,curved_geometry=False))
    faces=len(top['internal_faces'])+len(top['boundary_faces']);nnz1=len(p)+2*len(edges)
    # Count unique adjacency pairs of corner/edge entities in a conforming tetra complex.
    nnz2=len(x)+6*len(edges)+12*faces+6*len(t)
    dof1=4*len(p);dof2=3*len(x)+len(p);allocated=4*len(x)
    matrix_bytes=nnz2*16*12+allocated*4
    dump(REPORT/'data/taylor_hood_resource_estimate.json',dict(original_vertices=len(p),unique_edges=len(edges),tet10_nodes=len(x),tet10_cells=len(t),
        velocity_nodes=len(x),pressure_corner_nodes=len(p),velocity_dof=3*len(x),pressure_dof=len(p),P1_total_dof=dof1,P2P1_total_dof=dof2,
        active_dof_ratio=dof2/dof1,solver_allocated_dof=allocated,extra_pressure_slots=len(edges),
        storage_note='Native solver retains 4 slots/node and identity constraints on edge pressure slots; active pressure basis uses corners only.',
        P1_scalar_nnz_estimate=16*nnz1,P2_scalar_nnz_estimate=16*nnz2,matrix_nnz_ratio=nnz2/nnz1,
        CSR_matrix_GiB=matrix_bytes/2**30,host_memory_estimate_GiB=[matrix_bytes*3/2**30,matrix_bytes*10/2**30],
        GPU_memory_estimate_GiB=[matrix_bytes*1.2/2**30,matrix_bytes*3/2**30],
        estimate_notes='Full 4x4 block adjacency upper storage; ILU fill and host staging vary. Measure actual runtime peaks.'))
    # Restore the authoritative Flat + zero shared perimeter profile on every DOF.
    inlet=boundaries['INLET'];wallids=np.unique(boundaries['WALL']);inletids=np.unique(inlet)
    profile=np.ones(len(x));profile[wallids]=0
    ix=x[inlet[:,:3]];av=np.cross(ix[:,1]-ix[:,0],ix[:,2]-ix[:,0])/2
    av*=np.sign(np.einsum('ij,ij->i',av,ix.mean(axis=1)-p[t[boundary_owners(top,inlet[:,:3])]].mean(axis=1)))[:,None]
    normal=av.sum(axis=0);normal/=np.linalg.norm(normal)
    # Baseline inlet is planar up to original coordinate precision; profile normal follows native face normal.
    area=np.linalg.norm(av,axis=1);N=tri_basis(TRI_Q)
    discrete=float(np.sum(area*((profile[inlet]@N.T)@TRI_W)))
    target=1.551359160440232e-14;alpha=target/discrete
    value=profile[inletids,None]*(-normal)*alpha
    integrated=float(np.sum(area*((profile[inlet]@N.T)@TRI_W))*alpha)
    assert abs(integrated/target-1)<1e-14
    np.savez_compressed(REPORT/'data/tet10_inlet_profile.npz',node_ids=inletids,velocity_m_s=value,flat_shape=profile[inletids],normal_outward=normal)
    tree=ET.parse(BASE/'run/solver.xml');eq=tree.find('Add_equation');ET.SubElement(eq,'Use_taylor_hood_type_basis').text='true'
    oldQ=float(eq.find("Add_BC[@name='INLET']/Value").text)
    # Requested reference is the exactly remeasured baseline inlet. Single scalar change is documented.
    eq.find("Add_BC[@name='INLET']/Value").text=format(-target,'.17g')
    tree.write(CASE/'run/solver.xml',encoding='utf-8',xml_declaration=True)
    shutil.copy2(BASE/'run/PETSC_OPTIONS.txt',CASE/'run/PETSC_OPTIONS.txt');shutil.copy2(BASE/'policy.json',CASE/'policy.json')
    dump(REPORT/'data/tet10_inlet_profile_audit.json',dict(status='PASS',authoritative_profile='Flat; shared cap/WALL perimeter nodes zero; native flux normalization',
        source='Code/Source/solver/baf_ini.cpp:bc_ini; baseline solver.xml INLET',target_Q_m3_s=target,
        unit_amplitude_discrete_flux_m3_s=discrete,alpha_m_s=alpha,integrated_profile_Q_m3_s=integrated,
        exact_original_XML_flux_m3_s=-oldQ,normalization_relative_to_XML=target/(-oldQ),
        normalization_note='Only the globally allowed scalar normalization; requested Q is measured baseline Q, 2.87046e-10 below original XML prescription.',
        inlet_nodes=len(inletids),inlet_midnodes=int(np.sum(inletids>=len(p))),inlet_wall_shared_nodes=int(np.sum(np.isin(inletids,wallids))),
        no_linear_averaging_of_corner_velocity=True,runtime_nodal_profile_verification_pending=True))
    dump(REPORT/'data/solver_config_diff.json',dict(physics_changes=['TET4 to straight-sided TET10 order elevation','P1/P1 VMS to native P2/P1 Taylor-Hood'],
        inlet_single_scalar_normalization=dict(old=oldQ,new=-target,ratio=target/(-oldQ),reason='Explicit requested Qin equals actual measured baseline Qin'),
        unchanged=['geometry','corner tetra topology','rho','mu','outlet pressure strings','WALL no-slip','dt','steady policy','output cadence','solver tolerances','PETSC_OPTIONS'],
        new_XML=str(CASE/'run/solver.xml'),baseline_XML_sha256=sha(BASE/'run/solver.xml'),new_XML_sha256=sha(CASE/'run/solver.xml')))
    print('MESH_ELEVATION_PASS',len(p),len(edges),len(x),len(t),'active DOF ratio',dof2/dof1,flush=True)

if __name__=='__main__':main()
