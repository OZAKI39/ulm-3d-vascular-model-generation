"""Native FE checkpoint restoration and evaluation; no equilibrium solve here."""
import numpy as np
from scipy.spatial import cKDTree
from mpi4py import MPI
from petsc4py import PETSc
import basix,basix.ufl,ufl,dolfinx
from dolfinx import fem,mesh,geometry
from dolfinx.fem.petsc import assemble_vector
from .common import *
class Native:
    def __init__(self,label,degree=2,restore=True):
        self.label=label;self.m=dict(np.load(OLD/'geometry'/label/'mesh.npz'));self.X=self.m['X'];self.R=config()['physics']['R_m'];self.degree=degree
        self.domain=mesh.create_mesh(MPI.COMM_SELF,self.m['cells'],ufl.Mesh(basix.ufl.element('Lagrange','tetrahedron',1,shape=(3,))),self.X)
        d=self.domain;d.topology.create_connectivity(2,3);d.topology.create_connectivity(2,0)
        self.boundary=mesh.exterior_facet_indices(d.topology);self.wall=mesh.locate_entities_boundary(d,2,lambda x:np.abs(x[2])<1e-12*self.R)
        self.tags=np.where(np.isin(self.boundary,self.wall),2,3).astype(np.int32);self.mt=mesh.meshtags(d,2,self.boundary,self.tags)
        self.V=fem.functionspace(d,('Lagrange',degree,(3,)));self.u=fem.Function(self.V);self.coords=self.V.tabulate_dof_coordinates()
        self.fixed=np.zeros(len(self.coords),bool);self.fixed[fem.locate_dofs_topological(self.V,2,self.wall)]=True
        if restore:
            old=np.load(OLD/'reference'/label/'native_dofs.npz');self.restore_ids,self.mapping_error=bijective_map(old['X'],self.coords,self.R*1e-11)
            self.u.x.array[:]=old['u'][self.restore_ids].ravel();self.u.x.scatter_forward()
            assert np.array_equal(self.fixed,old['fixed'][self.restore_ids]);assert np.all(self.u.x.array.reshape(-1,3)[self.fixed]==0)
        self.F=ufl.Identity(3)+ufl.grad(self.u);self.gradspace=fem.functionspace(d,('DG',degree-1,(3,3)));self.gf=fem.Function(self.gradspace);self.update_gradient()
        self.tree=geometry.bb_tree(d,3);self.original_cells=d.topology.original_cell_index
        self.nc=d.topology.index_map(3).size_local;self.entities=np.arange(self.nc,dtype=np.int32)
        self.vertex_ref=np.array([[0.,0,0],[1.,0,0],[0.,1.,0],[0.,0,1.]])
        self.xn=fem.Expression(ufl.SpatialCoordinate(d),self.vertex_ref).eval(d,self.entities)
        dist,self.vertex_ids=cKDTree(self.X).query(self.xn);assert dist.max()<self.R*1e-11
        assert np.array_equal(np.sort(self.vertex_ids,axis=1),np.sort(self.m['cells'][self.original_cells],axis=1))
        self.cv=np.abs(np.linalg.det(self.xn[:,1:]-self.xn[:,0,None]))/6
    def update_gradient(self):self.gf.interpolate(fem.Expression(self.F,self.gradspace.element.interpolation_points))
    def locate(self,points,original=False):
        if not original:return self.locate_verified(points)
        points=np.ascontiguousarray(points,dtype=np.float64);a=geometry.compute_collisions_points(self.tree,points);b=geometry.compute_colliding_cells(self.domain,a,points)
        ids=np.full(len(points),-1,np.int32);multiplicity=np.diff(b.offsets).astype(np.int32)
        for i in np.flatnonzero(multiplicity):
            hits=b.links(int(i));ids[i]=hits[0] if original else hits[np.argmin(self.original_cells[hits])]
        canonical=np.full(len(points),-1,np.int32);valid=ids>=0;canonical[valid]=self.original_cells[ids[valid]]
        return ids,canonical,multiplicity
    def locate_verified(self,points):
        # Affine tetrahedra: certify bbox candidates with dimensionless barycentrics.
        # GJK collision lists alone can contain false positives at this SI scale.
        points=np.ascontiguousarray(points,dtype=np.float64);ids=np.full(len(points),-1,np.int32);multiplicity=np.zeros(len(points),np.int32)
        inv=np.linalg.inv(((self.xn[:,1:]-self.xn[:,0,None])/self.R).transpose(0,2,1));canonical=np.full(len(points),-1,np.int32)
        tolocal=np.empty(self.nc,np.int32);tolocal[self.original_cells]=np.arange(self.nc)
        for start in range(0,len(points),16384):
            p=points[start:start+16384];a=geometry.compute_collisions_points(self.tree,p);ci=a.array;qi=np.repeat(np.arange(len(p)),np.diff(a.offsets))
            b=np.einsum('kij,kj->ki',inv[ci],(p[qi]-self.xn[ci,0])/self.R);ok=(b.min(1)>=-1e-11)&(b.sum(1)<=1+1e-11)
            hit=np.full(len(p),self.nc,np.int32);np.minimum.at(hit,qi[ok],self.original_cells[ci[ok]]);count=np.bincount(qi[ok],minlength=len(p));good=hit<self.nc
            part=ids[start:start+len(p)];part[good]=tolocal[hit[good]];cp=canonical[start:start+len(p)];cp[good]=hit[good];multiplicity[start:start+len(p)]=count
        return ids,canonical,multiplicity
    def sample(self,points,ids,fun=None):
        f=self.u if fun is None else fun;shape=(3,) if f.function_space==self.V else (int(np.prod(f.ufl_shape)),)
        out=np.full((len(points),*shape),np.nan);valid=ids>=0
        for chunk in np.array_split(np.flatnonzero(valid),max(1,int(valid.sum()/32768)+1)):
            if len(chunk):out[chunk]=f.eval(np.ascontiguousarray(points[chunk]),ids[chunk])
        return out
    def forms(self,u=None,q=6):
        u=self.u if u is None else u;v=ufl.TestFunction(u.function_space);F=ufl.Identity(3)+ufl.grad(u);E=.5*(F.T*F-ufl.Identity(3));W=1000*ufl.inner(E,E)
        dx=ufl.Measure('dx',domain=self.domain,metadata={'quadrature_degree':q});ds=ufl.Measure('ds',domain=self.domain,subdomain_data=self.mt,metadata={'quadrature_degree':q})
        N=ufl.FacetNormal(self.domain);T=-3*N+3*(ufl.as_vector((1.,0.,0.))-N[0]*N)
        return W*dx,ufl.derivative(W*dx,u,v),ufl.dot(T,v)*ds(3)
    @staticmethod
    def vector(form):
        v=assemble_vector(fem.form(form));v.ghostUpdate(addv=PETSc.InsertMode.ADD,mode=PETSc.ScatterMode.REVERSE);out=v.array.copy();v.destroy();return out
    def p1(self,values,target=False):
        space=fem.functionspace(self.domain,('Lagrange',1,(3,)));p=fem.Function(space);ids,_=bijective_map(self.X,space.tabulate_dof_coordinates(),self.R*1e-11);p.x.array[:]=values[ids].ravel();p.x.scatter_forward()
        if not target:return p
        t=fem.Function(self.V);t.interpolate(p)
        if np.all(values[self.m['fixed']]==0):t.x.array.reshape(-1,3)[self.fixed]=0
        t.x.scatter_forward();return t
    def checkpoint(self,path):
        path=Path(path);path.mkdir(parents=True,exist_ok=True);imap=self.V.dofmap.index_map
        np.savez_compressed(path/'native_checkpoint.npz',mesh_X=self.X,mesh_cells=self.m['cells'],geometry_X=self.domain.geometry.x,geometry_dofmap=self.domain.geometry.dofmap,original_cell_ids=self.original_cells,facet_indices=self.boundary,facet_tags=self.tags,dofmap=self.V.dofmap.list,dof_coords=self.coords,u_coefficients=self.u.x.array,owned_global_dof_rows=imap.local_to_global(np.arange(imap.size_local,dtype=np.int32)),fixed=self.fixed,cell_vertex_ids=self.vertex_ids)
        dump(path/'CHECKPOINT_IDENTITY.json',dict(degree=self.degree,geometry_degree=1,element='continuous vector Lagrange blocked xyz',MPI_size=1,MPI_rank=0,owned_rows=imap.size_local,ghost_rows=imap.num_ghosts,block_size=self.V.dofmap.bs,dolfinx=dolfinx.__version__,basix=basix.__version__,PETSc=PETSc.Sys.getVersion(),dtype=str(self.u.x.array.dtype),configuration_sha256=digest(OUT/'AUDIT_CONFIGURATION.json'),source_sha256=digest(__file__)))
    def audit(self):
        path=OUT/'native'/self.label;path.mkdir(parents=True,exist_ok=True);self.checkpoint(path)
        trace=fem.Expression(self.F,self.vertex_ref).eval(self.domain,self.entities);F,A,B,disp,w=recover_traces(trace,self.vertex_ids,self.cv,len(self.X))
        ut=fem.Expression(self.u,self.vertex_ref).eval(self.domain,self.entities);un=np.zeros_like(self.X);np.add.at(un,self.vertex_ids.ravel(),(ut*self.cv[:,None,None]).reshape(-1,3));un/=w[:,None]
        old=np.load(OLD/'reference'/self.label/'state.npz');sc=config()['scales'];checks=dict(coefficient_bitwise=True,coordinate_error_over_R=self.mapping_error/self.R,cell_identity=True,volume_max_relative=float(np.max(abs(self.cv-self.m['cell_volume'][self.original_cells]))/self.cv.max()),lumped_volume_max_relative=float(np.max(abs(w/4-self.m['volume']))/self.m['volume'].max()),u_vertex_scaled=float(np.max(abs(un-old['u']))/sc['displacement_m']),F_vertex_scaled=float(np.max(abs(F-old['Fbar']))/sc['gradient']),Psi_A_scaled=float(np.max(abs(A-old['Psi_plus']))/sc['Psi_plus_Pa']))
        # Physical boundary facets: exact vertex sets and outward orientation from adjacent cell.
        topology_vertices=self.domain.topology.connectivity(2,0);dtop=self.domain.topology.connectivity(2,3);faces=[];normals=[]
        for fi in self.boundary:
            vids=topology_vertices.links(int(fi));gids=mesh.entities_to_geometry(self.domain,0,vids,False).ravel();coords=self.domain.geometry.x[gids];_,orig=cKDTree(self.X).query(coords);faces.append(orig)
            normal=np.cross(coords[1]-coords[0],coords[2]-coords[0]);normal/=np.linalg.norm(normal);center=coords.mean(0);cell=self.xn[dtop.links(int(fi))[0]].mean(0)
            if np.dot(normal,center-cell)<0:normal=-normal
            normals.append(normal)
        lookup={tuple(sorted(tri)):i for i,tri in enumerate(self.m['triangles'])};faceids=np.array([lookup[tuple(sorted(tri))] for tri in faces]);assert len(np.unique(faceids))==len(self.m['triangles'])
        checks['facet_identity']=True;checks['facet_normal_max_difference']=float(np.max(abs(np.asarray(normals)-self.m['normals'][faceids])));checks['facet_wall_identity']=bool(np.array_equal(self.tags==2,self.m['wall_faces'][faceids]));np.savez_compressed(path/'physical_mapping.npz',original_cell_ids=self.original_cells,original_facet_ids=faceids,native_to_saved_dof_rows=self.restore_ids,vertex_ids=self.vertex_ids)
        # DG1 exactness on affine cells, tested at fixed interior points via native Expression.
        points=np.array([[.1,.2,.3],[.25,.25,.25],[.6,.1,.1]])
        direct=fem.Expression(self.F,points).eval(self.domain,self.entities);dg=fem.Expression(self.gf,points).eval(self.domain,self.entities)
        checks['native_vs_DG_scaled']=float(np.max(abs(direct-dg))/sc['gradient'])
        energy,internal,external=self.forms();iv=self.vector(internal).reshape(-1,3);ev=self.vector(external).reshape(-1,3);ov=np.load(OLD/'reference'/self.label/'native_dofs.npz');r=iv-ev
        checks['external_scaled']=float(np.max(abs(ev-ov['external_N'][self.restore_ids]))/sc['force_N']);checks['reaction_scaled']=float(np.max(abs(r[self.fixed]-ov['reaction_N'][self.restore_ids][self.fixed]))/sc['force_N']);checks['free_residual_relative']=float(np.linalg.norm(r[~self.fixed])/np.linalg.norm(ev[~self.fixed]));checks['energy_J']=float(fem.assemble_scalar(fem.form(energy)))
        profiles=np.load(OLD/'reference'/self.label/'profiles.npz');profilechecks={}
        for key in ['axis','transverse']:
            xx=profiles[key+'_X'];ids,_,_=self.locate(xx,original=True);uu=self.sample(xx,ids);ff=self.sample(xx,ids,self.gf).reshape(-1,3,3)
            profilechecks[key]=dict(u_scaled=float(np.max(abs(uu-profiles[key+'_u']))/sc['displacement_m']),F_scaled=float(np.max(abs(ff-profiles[key+'_Fbar']))/sc['gradient']))
        checks['profiles']=profilechecks;checks['PASS']=bool(all(checks[k]<1e-10 for k in ['u_vertex_scaled','F_vertex_scaled','Psi_A_scaled','native_vs_DG_scaled','external_scaled','reaction_scaled']) and all(max(c.values())<1e-10 for c in profilechecks.values()) and checks['facet_wall_identity'] and checks['free_residual_relative']<1e-8)
        np.save(path/'native_F_cell_vertex_traces.npy',trace);np.savez_compressed(path/'recovery.npz',X=self.X,F_mean=F,A=A,B=B,trace_dispersion=disp,volume=self.m['volume'],ring_distance=ring_distance(self.X,self.R),fixed=self.m['fixed'],ROI=ring_distance(self.X,self.R)>=config()['ROI_width_m'],u=un)
        stats={}
        for k,values in [('A',A),('B',B)]:
            stats[k]={}
            for region,mask in [('full',np.ones(len(A),bool)),('ROI',ring_distance(self.X,self.R)>=config()['ROI_width_m'])]:
                s=weighted_stats(values,self.m['volume'],self.X,mask,save_order=path/f'{k}_{region}_ECDF.npz');tail=mask&(values>=s['P99']);s['tail_weight_fraction_on_wall']=float(self.m['volume'][tail&self.m['fixed']].sum()/self.m['volume'][tail].sum());s['tail_node_ids']=np.flatnonzero(tail);stats[k][region]=s
        dump(path/'FIDELITY.json',checks);dump(path/'RECOVERY_ORDER.json',stats);print(self.label,'NATIVE_FIDELITY',checks,flush=True)
        assert checks['PASS'];return checks
