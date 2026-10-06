"""Diagnostic constrained L2 RT0 projection using mature scikit-fem assembly.
Original WALL normal flux fixed; artificial ROI cuts explicitly free.
No pressure solve, no Navier-Stokes solve, no production particle field.
"""
from common import *
from skfem import MeshTet, Basis, ElementTetRT0, ElementTetP0, ElementTetP1, BilinearForm, LinearForm, asm
from skfem.helpers import dot, div
from scipy.sparse import bmat
from scipy.sparse.linalg import spsolve
import pyvista as pv

@BilinearForm
def mass(u,v,w):return dot(u,v)
@BilinearForm
def divergence(u,q,w):return div(u)*q
@LinearForm
def load(v,w):return dot(w.target,v)

def evaluate(basis,c,x,cell):
    """Library Piola basis on caller's verified original cell; no point relocation."""
    x=np.asarray(x).T;cell=np.asarray(cell,int)
    ref=basis.mapping.invF(x[:,:,None],tind=cell)
    value=np.zeros((3,len(cell),1))
    for k in range(basis.Nbfun):
        value+=basis.elem.gbasis(basis.mapping,ref,k,tind=cell)[0]*c[basis.element_dofs[k,cell]][None,:,None]
    return value[:,:,0].T

def face_geometry(m):
    tri=m.p[:,m.facets].transpose(2,1,0);center=tri.mean(1)
    area_vec=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])/2
    cell0=m.f2t[0];ct=m.p[:,m.t[:,cell0]].mean(1).T
    flip=np.einsum('ij,ij->i',area_vec,ct-center)>0;area_vec[flip]*=-1
    return center,area_vec

def project(m,u_nodes,fixed_faces=(),fixed_flux=None):
    r=Basis(m,ElementTetRT0(),intorder=2);p=Basis(m,ElementTetP0(),quadrature=r.quadrature)
    v=Basis(m,ElementTetP1(),quadrature=r.quadrature)
    uq=np.array([v.interpolate(u_nodes[:,j]) for j in range(3)])
    M=asm(mass,r);B=asm(divergence,r,p);f=asm(load,r,target=uq)
    center,av=face_geometry(m);fixed_faces=np.asarray(fixed_faces,int)
    # Do not assume library DOF equals integrated flux: actual library RT0
    # normalization is verified geometrically for every face.
    unit=evaluate(r,np.ones(r.N),center,m.f2t[0]);factor=np.einsum('ij,ij->i',unit,av)
    assert np.all(np.abs(factor)>.1)
    prescribed=np.zeros(r.N);fixed=r.facet_dofs[0,fixed_faces]
    flux=np.zeros(len(fixed_faces)) if fixed_flux is None else np.asarray(fixed_flux)
    prescribed[fixed]=flux/factor[fixed_faces]
    free=np.setdiff1d(np.arange(r.N),fixed);all_boundary_fixed=np.all(np.isin(m.boundary_facets(),fixed_faces))
    rows=np.arange(m.nelements)
    if all_boundary_fixed:
        if abs(flux.sum())>1e-11*max(float(np.abs(flux).sum()),1e-30):raise ValueError('INCOMPATIBLE_FIXED_BOUNDARY_FLUX')
        # Connected analytical test mesh: one redundant global divergence row.
        rows=rows[:-1]
    Bb=B[rows][:,free];rhs=np.r_[f[free]-M[free]@prescribed,-(B@prescribed)[rows]]
    A=bmat([[M[free][:,free],Bb.T],[Bb,None]],format='csc')
    sol=spsolve(A,rhs);assert np.isfinite(sol).all();c=prescribed.copy();c[free]=sol[:len(free)]
    residual=np.max(np.abs(A@sol-rhs))/max(np.max(np.abs(rhs)),1e-30)
    fluxes=np.einsum('ij,ij->i',evaluate(r,c,center,m.f2t[0]),av)
    shared=np.flatnonzero(m.f2t[1]>=0)
    jump=np.einsum('ij,ij->i',evaluate(r,c,center[shared],m.f2t[0,shared])-evaluate(r,c,center[shared],m.f2t[1,shared]),av[shared])
    net=B@c
    return dict(basis=r,c=c,B=B,M=M,center=center,area_vec=av,flux=fluxes,net=net,
      metrics=dict(linear_relative_residual=float(residual),max_cell_net_flux=float(np.abs(net).max()),max_shared_flux_jump=float(np.max(np.abs(jump),initial=0)),max_fixed_flux_error=float(np.max(np.abs(fluxes[fixed_faces]-flux),initial=0))))

def analytic_tests():
    m=MeshTet.init_tensor(np.linspace(0,1,3),np.linspace(0,1,3),np.linspace(0,1,3));x=m.p.T;ctr,av=face_geometry(m);bd=m.boundary_facets();results={}
    functions={'constant':lambda x:np.tile([.3,-.2,.1],(len(x),1)),
      'linear_divergence_free':lambda x:np.c_[x[:,0],-x[:,1],.2+0*x[:,2]],
      'nonconservative':lambda x:x}
    for name,fun in functions.items():
        flux=np.einsum('ij,ij->i',fun(ctr[bd]),av[bd]) if name!='nonconservative' else None
        fixed=bd if flux is not None else np.array([],int)
        result=project(m,fun(x),fixed,flux);r=result['basis'];value=evaluate(r,result['c'],m.p[:,m.t].mean(1).T,np.arange(m.nelements))
        error=np.max(np.abs(value-fun(m.p[:,m.t].mean(1).T)))
        assert result['metrics']['max_cell_net_flux']<1e-12
        assert result['metrics']['max_shared_flux_jump']<1e-12
        assert result['metrics']['max_fixed_flux_error']<1e-12
        if name=='constant':assert error<1e-12
        if name=='nonconservative':assert np.max(np.abs(value-m.p[:,m.t].mean(1).T))>.05
        results[name]={**result['metrics'],'centroid_max_component_difference':float(error),'assertion':'Conservation/trace preserved; exact reproduction required only for constant RT0 member'}
    try:project(m,x,bd,np.einsum('ij,ij->i',ctr[bd],av[bd]))
    except ValueError as e:assert str(e)=='INCOMPATIBLE_FIXED_BOUNDARY_FLUX';results['incompatible_boundary_rejected']=True
    else:raise AssertionError('Must reject incompatible all-fixed boundary')
    return results

def main():
    tests=analytic_tests();put('RT_ANALYTIC_TESTS.json',tests);print('RT_ANALYTIC_TESTS_PASS',flush=True)
    meshpath,flowpath=paths('M0');grid=pv.read(meshpath);flow=pv.read(flowpath)
    tet=grid.cells.reshape(-1,5)[:,1:];points=np.asarray(grid.points,float);u=np.asarray(flow.point_data['Velocity'],float)
    cells=np.load(PREVIOUS/'roi_1189.npz')['cells'];nodes,inverse=np.unique(tet[cells],return_inverse=True);localtet=inverse.reshape(-1,4)
    origin=points[nodes].mean(0);length=1e-4
    m=MeshTet(((points[nodes]-origin)/length).T,localtet.T,sort_t=False)
    assert np.array_equal(m.t.T,localtet)
    wall=pv.read(BASE/'mesh/SV_MESH/mesh-surfaces/WALL.vtp');wt=wall.point_data['GlobalNodeID'][wall.faces.reshape(-1,4)[:,1:]]-1;wallkeys={tuple(sorted(v)) for v in wt}
    facets=nodes[m.facets.T];bd=m.boundary_facets();wallfaces=np.array([f for f in bd if tuple(sorted(facets[f])) in wallkeys],int);cut=np.setdiff1d(bd,wallfaces)
    # ROI contains only physical WALL and artificial cuts, no inlet/outlet.
    for role in ['INLET','OUTLET_01','OUTLET_02','OUTLET_03']:
        z=pv.read(BASE/f'mesh/SV_MESH/mesh-surfaces/{role}.vtp');keys={tuple(sorted(v)) for v in z.point_data['GlobalNodeID'][z.faces.reshape(-1,4)[:,1:]]-1}
        assert not any(tuple(sorted(facets[f])) in keys for f in bd)
    rt=project(m,u[nodes],wallfaces);center=rt['center']*length+origin;av=rt['area_vec']*length**2
    originalflux=np.einsum('ij,ij->i',u[facets].mean(1),av);newflux=rt['flux']*length**2
    net=rt['net']*length**2
    # Reference with the same RT0 space and wall constraint, but without div=0,
    # separates space projection error from adding local conservation.
    rb=rt['basis'];p1=Basis(m,ElementTetP1(),quadrature=rb.quadrature);uq=np.array([p1.interpolate(u[nodes,j]) for j in range(3)])
    f=asm(load,rb,target=uq);fixed=rb.facet_dofs[0,wallfaces];free=np.setdiff1d(np.arange(rb.N),fixed)
    unconstrained=np.zeros(rb.N);unconstrained[free]=spsolve(rt['M'][free][:,free],f[free])
    ct=points[tet[cells]].mean(1);localids=np.arange(len(cells));ur=evaluate(rb,rt['c'],(ct-origin)/length,localids);uc=evaluate(rb,unconstrained,(ct-origin)/length,localids)
    mapping={int(c):i for i,c in enumerate(cells)};raw=np.load(CAUSE/'original_1189.npz');old=np.load(PREVIOUS/'M0_accepted_1189.npz');wallq=np.load(PREVIOUS/'wall_accepted_1189.npz')
    ids=np.array([mapping.get(int(c),-1) for c in old['cell']]);assert np.all(ids>=0)
    x=raw['x_m'];path_rt=evaluate(rb,rt['c'],(x-origin)/length,ids);path_un=evaluate(rb,unconstrained,(x-origin)/length,ids)
    oldpath=old['flow'][:,:3];normal=wallq['normal'];e=read(PREVIOUS/'EVENT_FIELD_PROTOCOL.json')['events']['1189'];reps={}
    for tag in ['A','B']:
        i=int(np.flatnonzero(np.all(x==e[tag]['x_m'],axis=1))[0]);n=np.array(e[tag]['normal']);reps[tag]={}
        for label,value in [('M0',oldpath[i]),('RT0_no_div_constraint',path_un[i]),('RT0_div0',path_rt[i])]:reps[tag][label]=dict(u_m_s=value.tolist(),speed_mm_s=float(np.linalg.norm(value)*1e3),inward_normal_mm_s=float(value@n*1e3),tangent_mm_s=float(np.linalg.norm(value-(value@n)*n)*1e3))
    cv=[];face_lookup={tuple(sorted(v)):i for i,v in enumerate(facets)}
    for layer in [0,1,2]:
        z=np.load(PREVIOUS/f'cv_1189_L{layer}.npz');fi=np.array([face_lookup[tuple(sorted(v))] for v in z['original_node_ids']]);sign=np.sign(np.einsum('ij,ij->i',av[fi],z['normal_out']));q=newflux[fi]*sign;qo=originalflux[fi]*sign
        assert np.isclose(q.sum(),net[[mapping[int(c)] for c in z['cells']]].sum(),rtol=0,atol=1e-20)
        cv.append(dict(layer=layer,cells=z['cells'].tolist(),original_Qnet_m3_s=float(qo.sum()),original_ratio=float(qo.sum()/(.5*np.abs(qo).sum())),RT_Qnet_m3_s=float(q.sum()),RT_throughput_m3_s=float(.5*np.abs(q).sum()),RT_ratio=float(q.sum()/(.5*np.abs(q).sum())),original_face_flux=qo.tolist(),RT_face_flux=q.tolist()))
    np.savez_compressed(HERE/'RT_REFERENCE.npz',cells=cells,nodes=nodes,local_tetra=localtet,points_m=points[nodes],original_velocity_m_s=u[nodes],facet_nodes=facets,face_center_m=center,area_vector_m2=av,original_flux_m3_s=originalflux,RT_flux_m3_s=newflux,wall_faces=wallfaces,cut_faces=cut,RT_net_m3_s=net,rt_coeff=rt['c'],RT0_no_div_coeff=unconstrained,
      original_x_m=x,time_s=raw['physical_time_s'],original_path_u=oldpath,RT_path_u=path_rt,RT0_no_div_path_u=path_un,inward_normals=normal,cell_centroid_m=ct,RT_centroid_u=ur,RT0_no_div_centroid_u=uc)
    selected=grid.extract_cells(cells);selected.cell_data['RT_div0_u_m_s']=ur;selected.cell_data['RT_no_div_u_m_s']=uc;selected.cell_data['RT_Qnet_m3_s']=net;selected.save(HERE/'RT_REFERENCE_DIAGNOSTIC.vtu')
    put('RT_REFERENCE.json',dict(status='DIAGNOSTIC_ONLY',RT_REFERENCE_NOT_PRODUCTION_FIELD=True,library='scikit-fem 12.0.2 ElementTetRT0 (alias ElementTetRT1)',
      method='Volume L2 closest RT0 velocity with per-tetra zero-divergence constraints, physical WALL normal flux zero; artificial cuts FREE. Additional same-space unconstrained-div control.',
      not_exact_Selzer_Cirpka='Their objective uses hydraulic gradients at centroids and Darcy conductivity; here velocity L2 norm, no Darcy conductivity.',
      coordinates='Original mesh SI; library algebra uses (x-origin)/1e-4; all delivered fluxes restored m3/s; velocities remain m/s.',
      cells=len(cells),wall_faces=len(wallfaces),artificial_cut_faces=len(cut),boundary_preservation='Only physical WALL normal flux is preserved. Tangential no-slip and all artificial cut fluxes are NOT constraints.',
      old_boundary_Qnet_m3_s=float(originalflux[bd].sum()),new_boundary_Qnet_m3_s=float(newflux[bd].sum()),cut_L1_change_m3_s=float(np.abs(newflux[cut]-originalflux[cut]).sum()),
      cannot_preserve_all_original_cut_fluxes=bool(abs(originalflux[bd].sum())>1e-15),
      solver_metrics_dimensionless_geometry=rt['metrics'],max_cell_Qnet_m3_s=float(np.abs(net).max()),max_wall_flux_m3_s=float(np.abs(newflux[wallfaces]).max()),
      representatives=reps,control_volumes=cv,
      limitations=['Zero-divergence RT0 is constant within each tet: interior G/curl/shear vanish, tangential jumps remain.','ROI artificial cuts changed; this is not a global production field nor a unique causal intervention.','Mass balance is enforced by construction, not an independent accuracy measurement.','No trajectory integrated.']))
    print('RT_COMPLETE',len(cells),reps['A'],flush=True)
if __name__=='__main__':main()
