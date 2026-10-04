"""One-way traction interface; all positions m, velocities m/s, stress Pa."""
from pathlib import Path
from typing import Protocol
import numpy as np
from scipy.spatial import Delaunay
from vendor.brava_wss_reference import p1_gradients


class StreamingProvider(Protocol):
    def traction(self, points, normals, time_s): ...


def fluid_traction(gradient, pressure, normals, viscosity):
    sigma = viscosity*(gradient+gradient.swapaxes(1,2))-pressure[:,None,None]*np.eye(3)
    return np.einsum('nij,nj->ni',sigma,normals)


def split_traction(traction,normals):
    normal=np.sum(traction*normals,axis=1)[:,None]*normals
    return normal,traction-normal


def pipe_field(points,c):
    R=c['radius_m'];L=c['length_m'];Q=c['flow_rate_m3_s'];mu=c['viscosity_Pa_s']
    y,z=points[:,1],points[:,2];umax=2*Q/(np.pi*R*R)
    u=np.zeros_like(points);u[:,0]=umax*(1-(y*y+z*z)/(R*R))
    grad=np.zeros((len(points),3,3));grad[:,0,1]=-2*umax*y/(R*R);grad[:,0,2]=-2*umax*z/(R*R)
    p=8*mu*Q/(np.pi*R**4)*(L/2-points[:,0])
    return u,p,grad


class AnalyticTestStreaming:
    def __init__(self,c,pipe,representative_frequency):
        self.c=c;self.pipe=pipe;self.frequency=representative_frequency

    def traction(self,points,normals,time_s):
        c=self.c
        scale=c.get('traction_scale_Pa')
        if scale is None:
            scale=self.pipe['viscosity_Pa_s']*c['streaming_velocity_scale_m_s']/c['streaming_length_scale_m']
        center=np.asarray(c['bubble_center_m'])
        envelope=np.exp(-.5*np.sum((points-center)**2,axis=1)/c['streaming_length_scale_m']**2)
        if c.get('uniform',False):envelope[:]=1
        tangent=np.tile([1.,0.,0.],(len(points),1))-normals[:,0,None]*normals
        cycle=c['mean_fraction']+c['oscillatory_amplitude']*np.sin(2*np.pi*self.frequency*time_s+c['phase_rad'])
        traction=envelope[:,None]*cycle*(scale*tangent+c['normal_traction_scale_Pa']*normals)
        if c['include_pipe_traction']:
            _,p,grad=pipe_field(points,self.pipe)
            traction+=fluid_traction(grad,p,normals,self.pipe['viscosity_Pa_s'])
        return traction


class FileStreaming:
    """Static SI snapshot, optionally multiplied by a configured slow envelope.

    VTU/legacy tetra grids use original-cell P1 interpolation and containing-cell
    rejection. CSV and point VTP use bounded Delaunay interpolation in their
    intrinsic affine dimension. CSV has no supplied solid/fluid domain topology.
    """
    def __init__(self,path,viscosity,arrays=None):
        self.path=Path(path);self.mu=viscosity;self.mesh=None;self.tetra=None
        arrays=arrays or {};self.gradient=None;self.pressure=None;self.direct=None
        if self.path.suffix=='.csv':
            d=np.genfromtxt(self.path,delimiter=',',names=True)
            names=d.dtype.names
            self.points=np.column_stack([d[k] for k in ['x','y','z']])
            if all(k in names for k in ['tx','ty','tz']):self.direct=np.column_stack([d[k] for k in ['tx','ty','tz']])
            else:
                self.velocity=np.column_stack([d[k] for k in ['ux','uy','uz']]);self.pressure=d['p']
        else:
            import pyvista as pv
            self.mesh=pv.read(path);self.points=np.asarray(self.mesh.points)
            tname=arrays.get('traction','traction_Pa');uname=arrays.get('velocity','Velocity');pname=arrays.get('pressure','Pressure')
            if tname in self.mesh.point_data:self.direct=np.asarray(self.mesh[tname])
            else:
                self.velocity=np.asarray(self.mesh.point_data[uname]);self.pressure=np.asarray(self.mesh.point_data[pname])
            if isinstance(self.mesh,pv.UnstructuredGrid) and np.all(self.mesh.celltypes==10):
                self.tetra=np.asarray(self.mesh.cells).reshape(-1,5)[:,1:]
                if self.direct is None:self.gradient=p1_gradients(self.points,self.tetra,self.velocity)
            elif self.direct is None:
                raise ValueError('Velocity/pressure VTK import requires linear tetrahedra; surface velocity alone cannot define a 3-D gradient')
        if not np.isfinite(self.points).all():raise ValueError('Nonfinite source points')
        for a in [self.direct,self.pressure,getattr(self,'velocity',None)]:
            if a is not None and not np.isfinite(a).all():raise ValueError('Nonfinite imported field')
        if self.tetra is None:
            self.origin=self.points.mean(axis=0);_,s,vh=np.linalg.svd(self.points-self.origin,full_matrices=False)
            self.rank=int(np.sum(s>s[0]*1e-10));self.basis=vh[:self.rank].T
            if self.rank<2 or (self.direct is None and self.rank!=3):
                raise ValueError('Traction needs a plane/volume; velocity gradients need volumetric samples')
            self.local=(self.points-self.origin)@self.basis;self.tri=Delaunay(self.local)
        self.mapping_description='original tetra P1' if self.tetra is not None else 'bounded intrinsic Delaunay; convex sample domain'

    def interpolate(self,points):
        if self.tetra is not None:
            cell=np.asarray(self.mesh.find_containing_cell(points))
            if np.any(cell<0):raise ValueError('Import query outside original CFD cells; extrapolation forbidden')
            ids=self.tetra[cell];xyz=self.points[ids]
            q=np.linalg.solve((xyz[:,1:]-xyz[:,:1]).swapaxes(1,2),(points-xyz[:,0])[...,None])[...,0]
            weights=np.column_stack((1-q.sum(axis=1),q))
            grad=self.gradient[cell] if self.direct is None else None
        else:
            local=(points-self.origin)@self.basis
            residual=points-self.origin-local@self.basis.T
            if np.max(np.linalg.norm(residual,axis=1))>1e-9:raise ValueError('Query outside sampled plane')
            cell=self.tri.find_simplex(local,tol=1e-10)
            if np.any(cell<0):raise ValueError('Import query outside sampled convex hull; extrapolation forbidden')
            transform=self.tri.transform[cell];q=np.einsum('nij,nj->ni',transform[:,:self.rank],local-transform[:,self.rank])
            weights=np.column_stack((q,1-q.sum(axis=1)));ids=self.tri.simplices[cell]
            grad=None
            if self.direct is None:
                xyz=self.points[ids];du=self.velocity[ids[:,1:]]-self.velocity[ids[:,:1]]
                grad=np.linalg.solve(xyz[:,1:]-xyz[:,:1],du).swapaxes(1,2)
        if self.direct is not None:return np.einsum('ni,nij->nj',weights,self.direct[ids]),None
        pressure=np.sum(weights*self.pressure[ids],axis=1)
        return pressure,grad

    def traction(self,points,normals,time_s):
        value,gradient=self.interpolate(points)
        if self.direct is not None:return value
        return fluid_traction(gradient,value,normals,self.mu)


def make_provider(c,root=None):
    s=c['streaming']
    if s['provider']=='analytic_test':return AnalyticTestStreaming(s,c['pipe'],c['simulation']['representative_frequency_Hz'])
    if s['provider']=='file':
        path=Path(s['path']);path=path if path.is_absolute() else Path(root or '.')/path
        return FileStreaming(path,c['pipe']['viscosity_Pa_s'],s.get('arrays'))
    raise ValueError('Unknown streaming provider')


def map_forces(cloud,traction):
    force=np.zeros_like(cloud.X)
    np.add.at(force,cloud.face_particle,traction*cloud.face_area[:,None])
    return force
