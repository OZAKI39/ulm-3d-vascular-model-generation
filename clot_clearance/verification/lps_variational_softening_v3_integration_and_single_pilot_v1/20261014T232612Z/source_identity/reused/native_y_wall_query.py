"""Shared native sidewall, volume membership and open-port crossing queries."""
import time
import numpy as np
import vtk
from scipy.spatial import cKDTree
from .background_flow_sampler import BackgroundFlowSampler
from .y_mural_clot_geometry import locator,hit

class NativeYWallQuery:
    def __init__(self,sources=None):
        self.flow=BackgroundFlowSampler(sources);self.sources=self.flow.sources;s=self.sources;self.wall=s.faces['WALL'];self.wall_locator=locator(self.wall)
        self.tri=np.asarray(self.wall.points[self.wall.faces.reshape(-1,4)[:,1:]],float)
        # Orient using the opposite vertex of the actual adjacent lumen tetra.
        N=s.flow.n_points;faces=np.concatenate([np.sort(s.tet[:,[j for j in range(4) if j!=i]],axis=1) for i in range(4)])
        opposite=np.concatenate([s.tet[:,i] for i in range(4)]);keys=faces[:,0]*N*N+faces[:,1]*N+faces[:,2];order=np.argsort(keys)
        err,mapping=cKDTree(s.flow.points).query(self.wall.points);assert err.max()<1e-12
        ids=np.sort(mapping[self.wall.faces.reshape(-1,4)[:,1:]],axis=1);wk=ids[:,0]*N*N+ids[:,1]*N+ids[:,2];idx=np.searchsorted(keys[order],wk);assert np.array_equal(keys[order[idx]],wk)
        inward=s.flow.points[opposite[order[idx]]]-self.tri.mean(1)
        normals=np.cross(self.tri[:,1]-self.tri[:,0],self.tri[:,2]-self.tri[:,0]);normals/=np.linalg.norm(normals,axis=1)[:,None];normals[np.einsum('ij,ij->i',normals,inward)<0]*=-1;self.normals=normals
        self.ports={}
        for name,face in s.faces.items():
            if name=='WALL':continue
            centers=face.cell_centers().points;areas=face.compute_cell_sizes()['Area'];center=np.average(centers,axis=0,weights=areas)
            _,_,vt=np.linalg.svd(face.points-center,full_matrices=False);out=vt[-1]
            if np.dot(out,center-s.mesh.center)<0:out=-out
            self.ports[name]=dict(center=center,outward=out,locator=locator(face),surface=face)
        self.query_seconds=0.;self.query_points=0

    def query(self,positions,previous=None):
        start=time.perf_counter();X=np.atleast_2d(np.asarray(positions,float));P=np.empty_like(X);ids=np.empty(len(X),int)
        cp=[0.,0.,0.];cell=vtk.mutable(0);sub=vtk.mutable(0);d2=vtk.mutable(0.)
        for i,x in enumerate(X):
            self.wall_locator.FindClosestPoint(x,cp,cell,sub,d2);P[i]=cp;ids[i]=int(cell)
        distance=np.linalg.norm(X-P,axis=1);cells,_=self.flow.locate(X);inside=cells>=0
        clearance=np.where(inside,distance,-distance);normal=self.normals[ids]
        exit_ports=np.full(len(X),'',dtype='<U16')
        if previous is not None:
            prev=np.atleast_2d(previous)
            for i in np.flatnonzero(~inside):exit_ports[i]=self.port_crossing(prev[i],X[i]) or ''
        self.query_points+=len(X);self.query_seconds+=time.perf_counter()-start
        return dict(nearest_wall_point_m=P,inward_normal=normal,clearance_m=clearance,inside_lumen=inside,tetra_id=cells,wall_triangle_id=ids,exit_port=exit_ports,valid=np.isfinite(X).all(1))

    def port_crossing(self,previous,current):
        for name,p in self.ports.items():
            a=np.dot(previous-p['center'],p['outward']);b=np.dot(current-p['center'],p['outward'])
            if a<=1e-10 and b>0 and hit(p['locator'],previous,current)[0]:return name
        return None

    def swept_wall_crossing(self,a,b):
        return hit(self.wall_locator,a,b)

    def classification_crossing(self,a,b):
        # Native configured 10mm classification section, distinct from12mm outlet.
        distance=self.sources.config['geometry']['classification_distance_m'];out=[]
        for name,p in self.ports.items():
            if not name.startswith('OUTLET'):continue
            normal=p['outward'];center=normal*distance
            aa=np.dot(a-center,normal);bb=np.dot(b-center,normal)
            if aa<0<=bb:
                cross=a+(b-a)*(-aa/(bb-aa));v=cross-center
                if np.linalg.norm(v-np.dot(v,normal)*normal)<.0011:out.append(name.replace('OUTLET_','CLASSIFICATION_'))
        return out

    def minimal_checks(self):
        checks=[]
        for name,p in self.ports.items():
            a=p['center']-p['outward']*5e-6;b=p['center']+p['outward']*5e-6
            assert self.flow.locate([a])[0][0]>=0 and self.flow.locate([b])[0][0]<0
            assert self.port_crossing(a,b)==name
            checks.append(dict(port=name,inside_point_m=a.tolist(),outside_point_m=b.tolist(),legal_crossing=True))
        # The native sidewall is intentionally OPEN; containment uses tetrahedra.
        ids=np.arange(0,len(self.tri),max(1,len(self.tri)//100));center=self.tri[ids].mean(1)
        inside=self.flow.locate(center+self.normals[ids]*2e-6)[0]>=0;outside=self.flow.locate(center-self.normals[ids]*2e-6)[0]<0
        assert inside.all() and outside.all()
        return dict(status='PASS',ports=checks,orientation_triangles_checked=len(ids),native_wall_open_edges=int(self.wall.n_open_edges),containment='Native tetrahedra + strict barycentric weights; never closed-volume interpretation of open WALL',wall_normals='Opposite vertex of actual adjacent lumen tetra',classification_distance_m=self.sources.config['geometry']['classification_distance_m'],physical_outlet_distance_m=self.sources.config['geometry']['daughter_length_m'])
