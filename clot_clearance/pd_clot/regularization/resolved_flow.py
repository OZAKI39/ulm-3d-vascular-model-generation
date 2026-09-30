"""Explicit-unit, bounded static velocity/pressure import for the same PD API.

Tetrahedral VTU/VTK uses original-cell P1 interpolation. CSV and VTP use bounded
intrinsic Delaunay interpolation, which requires a declared convex-hull domain.
A planar field cannot supply a normal velocity derivative without explicit data.
"""
import hashlib,json
from pathlib import Path
import numpy as np
import pyvista as pv
from scipy.spatial import Delaunay
from ..streaming import fluid_traction
from .flow import FlowSample


UNITS={'length_unit':{'m':1.,'mm':1e-3,'um':1e-6},
       'velocity_unit':{'m/s':1.,'mm/s':1e-3,'um/s':1e-6},
       'pressure_unit':{'Pa':1.,'kPa':1e3},'time_unit':{'s':1.,'ms':1e-3,'us':1e-6}}


class ResolvedStreamingFieldProvider:
    def __init__(self,path,metadata,viscosity,surface_points,report_path=None):
        self.path=Path(path);self.mu=float(viscosity)
        self.metadata=json.loads(Path(metadata).read_text()) if isinstance(metadata,(str,Path)) else dict(metadata)
        self.report={'source':str(self.path.resolve()),'status':'VALIDATING','field_time_model':'Static snapshot',
                     'normal_convention_check':'Declared outward solid normal; this is not a geometric orientation proof'}
        try:
            self._initialize(surface_points)
            self.report['status']='ACCEPTED_VERIFICATION_CHECKS'
        except Exception as error:
            self.report.update(status='REJECTED',error=repr(error))
            raise
        finally:
            if report_path is not None:
                output=Path(report_path)
                if output.exists():raise FileExistsError(output)
                output.parent.mkdir(parents=True,exist_ok=True)
                output.write_text(json.dumps(self.report,indent=2,allow_nan=False)+'\n')

    def _initialize(self,surface_points):
        m=self.metadata;scales={}
        for key,known in UNITS.items():
            if key not in m or m[key] not in known:raise ValueError(f'Explicit supported {key} required; no unit inference')
            scales[key]=known[m[key]]
        if m.get('normal_convention')!='outward_solid':raise ValueError('PD traction requires declared outward_solid normal convention')
        if self.mu<=0:raise ValueError('Positive viscosity required')
        self.report.update(metadata=m,source_sha256=hashlib.sha256(self.path.read_bytes()).hexdigest(),conversion_to_SI=scales)
        arrays=m.get('arrays',{});self.mesh=None;self.tetra=None;self.supplied_gradient=None;self.gradient_location=None
        if self.path.suffix.lower()=='.csv':
            d=np.atleast_1d(np.genfromtxt(self.path,delimiter=',',names=True))
            self.points=np.column_stack([d[k] for k in ['x','y','z']])
            self.velocity_data=np.column_stack([d[k] for k in ['ux','uy','uz']]);self.pressure_data=np.asarray(d['p'])
            names=[f'd{v}_d{x}' for v in ['u','v','w'] for x in ['x','y','z']]
            present=[n in d.dtype.names for n in names]
            if any(present) and not all(present):raise ValueError('Partial velocity gradient is not sufficient')
            if all(present):self.supplied_gradient=np.column_stack([d[k] for k in names]).reshape(-1,3,3);self.gradient_location='point'
        elif self.path.suffix.lower() in ['.vtu','.vtp','.vtk']:
            self.mesh=pv.read(self.path)
            self.points=np.asarray(self.mesh.points,dtype=float).copy()
            self.velocity_data=np.asarray(self.mesh.point_data[arrays.get('velocity','velocity')],dtype=float)
            self.pressure_data=np.asarray(self.mesh.point_data[arrays.get('pressure','pressure')],dtype=float).reshape(-1)
            key=arrays.get('velocity_gradient','velocity_gradient')
            if key in self.mesh.point_data:
                self.supplied_gradient=np.asarray(self.mesh.point_data[key],dtype=float).reshape(-1,3,3);self.gradient_location='point'
            elif key in self.mesh.cell_data:
                self.supplied_gradient=np.asarray(self.mesh.cell_data[key],dtype=float).reshape(-1,3,3);self.gradient_location='cell'
            if isinstance(self.mesh,pv.UnstructuredGrid):
                if not np.all(self.mesh.celltypes==10):raise ValueError('Volume import currently requires linear tetrahedral cells')
                self.tetra=np.asarray(self.mesh.cells).reshape(-1,5)[:,1:]
        else:raise ValueError('Supported fields: VTU, VTP, legacy VTK, CSV')
        self.points=self.points*scales['length_unit'];self.velocity_data=self.velocity_data*scales['velocity_unit']
        self.pressure_data=self.pressure_data*scales['pressure_unit']
        if self.supplied_gradient is not None:
            if m.get('gradient_unit','velocity_unit/length_unit')!='velocity_unit/length_unit':
                raise ValueError('Gradient unit must be explicitly velocity_unit/length_unit')
            self.supplied_gradient=self.supplied_gradient*(scales['velocity_unit']/scales['length_unit'])
        n=len(self.points)
        if self.points.shape!=(n,3) or self.velocity_data.shape!=(n,3) or self.pressure_data.shape!=(n,):raise ValueError('Invalid field array shapes')
        for value in [self.points,self.velocity_data,self.pressure_data,self.supplied_gradient]:
            if value is not None and not np.isfinite(value).all():raise ValueError('NaN/Inf in imported field')
        if self.tetra is not None:
            self.mesh.points=self.points
            self.rank=3;cells=self.tetra
            self.mapping='Original tetrahedral P1; containing-cell rejection, no extrapolation'
        else:
            if m.get('sample_domain')!='convex_hull':raise ValueError('CSV/VTP requires explicit sample_domain=convex_hull; holes are not inferred')
            self.origin=self.points.mean(axis=0)
            _,singular,basis=np.linalg.svd(self.points-self.origin,full_matrices=False)
            self.rank=int(np.sum(singular>singular[0]*1e-10));self.basis=basis[:self.rank].T
            if self.rank<2:raise ValueError('At least a planar field is required')
            self.local=(self.points-self.origin)@self.basis
            self.tri=Delaunay(self.local);cells=self.tri.simplices
            self.mapping='Bounded intrinsic Delaunay over explicitly declared convex hull'
            if self.gradient_location=='cell':raise ValueError('Scattered surface/point import requires a point velocity gradient')
        if self.rank==2 and self.supplied_gradient is None:
            raise ValueError('A surface velocity field cannot determine the 3-D gradient; supply the full gradient or a volume field')
        self.cells=cells
        self.p1_gradient=None
        if self.rank==3:
            xyz=self.points[cells];du=self.velocity_data[cells[:,1:]]-self.velocity_data[cells[:,:1]]
            self.p1_gradient=np.linalg.solve(xyz[:,1:]-xyz[:,:1],du).swapaxes(1,2)
        if self.supplied_gradient is None:
            gradients=self.p1_gradient;self.gradient_location='derived_cell'
        elif self.gradient_location=='point':gradients=self.supplied_gradient[cells].mean(axis=1)
        else:gradients=self.supplied_gradient
        if self.gradient_location=='cell' and len(gradients)!=len(cells):raise ValueError('Gradient cell count mismatch')
        if self.p1_gradient is not None and self.supplied_gradient is not None:
            error=np.linalg.norm(gradients-self.p1_gradient,axis=(1,2))/np.maximum(np.maximum(np.linalg.norm(gradients,axis=(1,2)),np.linalg.norm(self.p1_gradient,axis=(1,2))),1e-30)
            self.report['maximum_relative_gradient_inconsistency']=float(error.max())
            if error.max()>m.get('validation',{}).get('maximum_relative_gradient_inconsistency',.25):raise ValueError('Supplied gradient is inconsistent with volume velocity interpolation')
        checked_gradients=np.concatenate((gradients,self.supplied_gradient),axis=0) if self.gradient_location=='point' else gradients
        divergence=np.trace(checked_gradients,axis1=1,axis2=2)
        relative=np.abs(divergence)/np.maximum(np.linalg.norm(checked_gradients,axis=(1,2)),1e-30)
        tol=m.get('validation',{}).get('maximum_relative_divergence',.05)
        if np.max(relative)>tol:raise ValueError(f'Imported relative divergence exceeds {tol}')
        limits=m.get('validation',{})
        if 'pressure_range_Pa' in limits:
            lo,hi=limits['pressure_range_Pa']
            if self.pressure_data.min()<lo or self.pressure_data.max()>hi:raise ValueError('Pressure outside declared validation range')
        if np.linalg.norm(self.velocity_data,axis=1).max()>limits.get('maximum_speed_m_s',float('inf')):raise ValueError('Velocity outside declared validation range')
        if np.linalg.norm(checked_gradients,axis=(1,2)).max()>limits.get('maximum_gradient_s_inv',float('inf')):raise ValueError('Gradient outside declared validation range')
        self.report.update(mapping=self.mapping,intrinsic_rank=self.rank,point_count=n,cell_count=len(cells),
            bounds_m=[self.points.min(axis=0).tolist(),self.points.max(axis=0).tolist()],
            pressure_range_Pa=[float(self.pressure_data.min()),float(self.pressure_data.max())],
            speed_range_m_s=[float(np.linalg.norm(self.velocity_data,axis=1).min()),float(np.linalg.norm(self.velocity_data,axis=1).max())],
            gradient_frobenius_range_s_inv=[float(np.linalg.norm(gradients,axis=(1,2)).min()),float(np.linalg.norm(gradients,axis=(1,2)).max())],
            divergence_range_s_inv=[float(divergence.min()),float(divergence.max())],maximum_relative_divergence=float(relative.max()),
            snapshot_time_s=m.get('snapshot_time',0.)*scales['time_unit'],gradient_source=self.gradient_location,
            normal_derivative_validation='SUPPLIED_ONLY_NOT_RECOVERABLE_FROM_PLANAR_DATA' if self.rank==2 else 'VOLUME_FIELD')
        surface_points=np.asarray(surface_points,float)
        if surface_points.ndim!=2 or surface_points.shape[1]!=3 or len(surface_points)==0:raise ValueError('Nonempty SI surface coverage points are required before accepting a field')
        value=self.sample(surface_points,0.)
        self.report.update(surface_query_count=len(surface_points),surface_coverage_fraction=1.,all_validation_arrays_finite=True)

    def _coordinates(self,points):
        points=np.asarray(points,float)
        if points.ndim!=2 or points.shape[1]!=3 or not np.isfinite(points).all():raise ValueError('Invalid query points')
        if self.tetra is not None:
            cell=np.asarray(self.mesh.find_containing_cell(points))
            if np.any(cell<0):raise ValueError('Query outside original CFD cells; extrapolation forbidden')
            ids=self.tetra[cell];xyz=self.points[ids]
            q=np.linalg.solve((xyz[:,1:]-xyz[:,:1]).swapaxes(1,2),(points-xyz[:,0])[...,None])[...,0]
            weights=np.column_stack((1-q.sum(axis=1),q))
        else:
            local=(points-self.origin)@self.basis
            residual=points-self.origin-local@self.basis.T
            tolerance=max(float(np.ptp(self.points,axis=0).max())*1e-9,1e-14)
            if np.any(np.linalg.norm(residual,axis=1)>tolerance):raise ValueError('Query outside source affine plane')
            cell=self.tri.find_simplex(local,tol=1e-12)
            if np.any(cell<0):raise ValueError('Query outside declared convex sample hull')
            tr=self.tri.transform[cell];q=np.einsum('nij,nj->ni',tr[:,:self.rank],local-tr[:,self.rank])
            weights=np.column_stack((q,1-q.sum(axis=1)));ids=self.tri.simplices[cell]
        return cell,ids,weights

    def sample(self,points,time_s):
        cell,ids,w=self._coordinates(points)
        velocity=np.einsum('ni,nij->nj',w,self.velocity_data[ids]);pressure=np.sum(w*self.pressure_data[ids],axis=1)
        if self.supplied_gradient is None:gradient=self.p1_gradient[cell]
        elif self.gradient_location=='point':gradient=np.einsum('ni,nijk->njk',w,self.supplied_gradient[ids])
        else:gradient=self.supplied_gradient[cell]
        return FlowSample(velocity,pressure,gradient,np.zeros_like(velocity))

    def traction(self,points,normals,time_s):
        value=self.sample(points,time_s)
        normals=np.asarray(normals,float)
        if not np.allclose(np.linalg.norm(normals,axis=1),1,rtol=0,atol=1e-8):raise ValueError('Unit outward-solid normals required')
        return fluid_traction(value.gradient,value.pressure,normals,self.mu)
