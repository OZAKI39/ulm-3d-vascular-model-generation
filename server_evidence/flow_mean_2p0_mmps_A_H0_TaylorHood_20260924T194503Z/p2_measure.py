"""Lossless ten-node velocity reader and genuine quadratic integral monitor."""
from pathlib import Path
import math
import numpy as np
import pyvista as pv
from p2 import basis,tri_basis,tetra_rule15,TRI_Q,TRI_W

class P2Measurements:
    def __init__(self,case,target=1.551359160440232e-14):
        self.case=Path(case);mesh=pv.read(self.case/'SV_MESH/mesh-complete.mesh.vtu')
        self.points=np.asarray(mesh.points,float);self.tetra=mesh.cells.reshape(-1,11)[:,1:]
        assert np.all(mesh.celltypes==24)
        self.Q=target;self.ncorner=int(np.max(self.tetra[:,:4]))+1
        x=self.points[self.tetra[:,:4]];self.volumes=abs(np.linalg.det(x[:,1:]-x[:,:1]))/6
        centroids=x.mean(axis=1);self.boundary={}
        for path in (self.case/'SV_MESH/mesh-surfaces').glob('*.vtp'):
            f=pv.read(path);tri=(np.asarray(f['GlobalNodeID'],int)-1)[f.faces.reshape(-1,7)[:,1:]]
            owners=np.asarray(f['GlobalElementID'],int)-1;xyz=self.points[tri[:,:3]]
            av=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])/2
            av[np.einsum('ij,ij->i',av,xyz.mean(axis=1)-centroids[owners])<0]*=-1
            self.boundary[path.stem]=(tri,av)
        self.wall_nodes=np.unique(self.boundary['WALL'][0]);self.q,self.w=tetra_rule15();self.N=basis(self.q)

    def read(self,path):
        grid=pv.read(path)
        assert grid.n_points==len(self.points) and grid.n_cells==len(self.tetra)
        assert np.all(grid.celltypes==24), 'P2 connectivity must be retained'
        assert np.array_equal(grid.cells.reshape(-1,11)[:,1:],self.tetra)
        # Native writer quantizes only coordinates; authoritative input double points
        # remain the geometry used by the solver. Verify exact native rounding map.
        assert np.array_equal(np.asarray(grid.points),self.points.astype(np.float32))
        u=np.asarray(grid['Velocity'],float);p=np.asarray(grid['Pressure'],float).reshape(-1)
        assert u.shape==(len(self.points),3)
        return u,p

    def velocity_l2(self,u):
        total=0.
        for i in range(0,len(self.tetra),10000):
            v=np.einsum('qa,tai->tqi',self.N,u[self.tetra[i:i+10000]])
            total+=float(self.volumes[i:i+10000]@(np.sum(v*v,axis=-1)@self.w))
        return math.sqrt(total)

    def measure(self,u,p):
        flows={};N=tri_basis(TRI_Q)
        for name,(tri,av) in self.boundary.items():
            values=np.einsum('qa,tai->tqi',N,u[tri]);q=np.einsum('tqi,ti->tq',values,av)@TRI_W
            flows[name]=math.fsum(q)
        incoming=-flows['INLET'];out={name:flows[name] for name in ('OUTLET_01','OUTLET_02','OUTLET_03')}
        wall=float(np.linalg.norm(u[self.wall_nodes],axis=1).max())
        return dict(Q_target_m3_s=self.Q,Q_in_m3_s=incoming,outlet_flows_m3_s=out,Q_out_total_m3_s=math.fsum(out.values()),
            signed_outward_boundary_flows_m3_s=flows,outlet_fractions={k:v/incoming for k,v in out.items()} if incoming else None,
            epsilon_Q=abs(incoming-self.Q)/self.Q,epsilon_mass=abs(math.fsum(out.values())-incoming)/self.Q,
            velocity_L2=self.velocity_l2(u),velocity_max_m_s=float(np.linalg.norm(u,axis=1).max()),
            pressure_range_pa=[float(p[:self.ncorner].min()),float(p[:self.ncorner].max())],
            velocity_finite=bool(np.isfinite(u).all()),pressure_finite=bool(np.isfinite(p).all()),
            wall_velocity_max_m_s=wall,wall_noslip_pass=wall<=2e-13)
