"""Read-only SI geometry and a vectorized equivalent of the frozen native sampler."""
from pathlib import Path
import numpy as np
import h5py, vtk
from vtk.util.numpy_support import numpy_to_vtk,numpy_to_vtkIdTypeArray,vtk_to_numpy

class FrozenSampler:
    VALID, OUTSIDE, SOLID, MISSING, NONFINITE_POSITION = range(5)
    def __init__(self,path):
        with h5py.File(path,'r') as f:
            self.dims=f['dims'][:].astype(np.int64);self.origin=f['origin_m'][:];self.dx=float(f['dx_m'][()])
            self.indices=f['linear_index'][:];self.fluid=f['fluid_linear_index'][:];self.velocity=f['Velocity_m_s'][:]
            for name,want in [('coordinate_units','m'),('velocity_units','m/s')]:
                value=f.attrs[name];value=value.decode() if isinstance(value,bytes) else value
                assert value==want
    def query(self,positions):
        positions=np.asarray(positions,dtype=float).reshape(-1,3)
        q=(positions-self.origin)/self.dx;finite=np.isfinite(q).all(axis=1)
        valid=finite & (q>=0).all(axis=1) & (q<self.dims-1).all(axis=1)
        status=np.where(finite,self.OUTSIDE,self.NONFINITE_POSITION).astype(np.int8);status[valid]=self.VALID
        base=np.floor(np.where(np.isfinite(q),q,0)).astype(np.int64);frac=q-base;out=np.zeros_like(positions)
        for z in [0,1]:
            for y in [0,1]:
                for x in [0,1]:
                    active=np.flatnonzero(status==self.VALID)
                    if not len(active):continue
                    b=base[active]+[x,y,z];ids=b[:,0]+self.dims[0]*(b[:,1]+self.dims[1]*b[:,2]);fi=np.searchsorted(self.fluid,ids);vi=np.searchsorted(self.indices,ids)
                    fluid=(fi<len(self.fluid)) & (self.fluid[np.minimum(fi,len(self.fluid)-1)]==ids)
                    present=(vi<len(self.indices)) & (self.indices[np.minimum(vi,len(self.indices)-1)]==ids)
                    status[active[~fluid]]=self.SOLID;status[active[fluid & ~present]]=self.MISSING
                    ok=fluid & present;idx=active[ok];weights=np.prod(np.where(np.array([x,y,z]),frac[idx],1-frac[idx]),axis=1)
                    out[idx]+=weights[:,None]*self.velocity[vi[ok]]
        out[status!=self.VALID]=0 # Exactly matches native status+zero return; never valid physical zero.
        return status,out

def poly(points,faces):
    p=vtk.vtkPolyData();v=vtk.vtkPoints();v.SetData(numpy_to_vtk(points,deep=True));p.SetPoints(v);c=vtk.vtkCellArray();c.SetCells(len(faces),numpy_to_vtkIdTypeArray(np.c_[np.full(len(faces),3),faces].astype(np.int64).ravel(),deep=True));p.SetPolys(c);return p

class Geometry:
    def __init__(self,path):
        self.g=dict(np.load(path));self.closed=poly(self.g['points'],self.g['faces']);self.walls=poly(self.g['points'],self.g['faces'][self.g['classes']<=2])
        self.locator=vtk.vtkStaticCellLocator();self.locator.SetDataSet(self.walls);self.locator.BuildLocator()
    def distance(self,p):
        q=[0.,0.,0.];cell=vtk.reference(0);sub=vtk.reference(0);d2=vtk.reference(0.)
        self.locator.FindClosestPoint(p,q,cell,sub,d2);return np.sqrt(float(d2)),q,int(cell)
    def section(self,center,normal):
        plane=vtk.vtkPlane();plane.SetOrigin(center);plane.SetNormal(normal)
        cutter=vtk.vtkCutter();cutter.SetCutFunction(plane);cutter.SetInputData(self.closed);cutter.SetOutputPointsPrecision(vtk.vtkAlgorithm.DOUBLE_PRECISION);cutter.Update()
        strip=vtk.vtkStripper();strip.SetInputConnection(cutter.GetOutputPort());strip.JoinContiguousSegmentsOn();strip.Update()
        m=strip.GetOutput();p=vtk_to_numpy(m.GetPoints().GetData());cells=vtk_to_numpy(m.GetLines().GetData());j=0;loops=[]
        while j<len(cells):
            n=cells[j];ids=cells[j+1:j+n+1];j+=n+1
            if n>=4 and np.linalg.norm(p[ids[0]]-p[ids[-1]])<1e-10:loops.append(p[ids[:-1]])
        assert loops,'STOP_INJECTION_SECTION_NOT_CLOSED'
        loop=min(loops,key=lambda x:np.linalg.norm(x.mean(axis=0)-center))
        e1=np.cross(normal,np.eye(3)[np.argmin(abs(normal))]);e1/=np.linalg.norm(e1);e2=np.cross(normal,e1)
        xy=(loop-center)@np.c_[e1,e2]
        p=vtk.vtkPoints();p.SetData(numpy_to_vtk(loop,deep=True));lines=vtk.vtkCellArray()
        for i in range(len(loop)):
            lines.InsertNextCell(2);lines.InsertCellPoint(i);lines.InsertCellPoint((i+1)%len(loop))
        contour=vtk.vtkPolyData();contour.SetPoints(p);contour.SetLines(lines)
        triang=vtk.vtkContourTriangulator();triang.SetInputData(contour);triang.Update();m=triang.GetOutput()
        faces=vtk_to_numpy(m.GetPolys().GetData()).reshape(-1,4)[:,1:];tri=vtk_to_numpy(m.GetPoints().GetData())[faces]
        area=.5*np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1).sum()
        polygon_area=.5*abs(np.sum(xy[:,0]*np.roll(xy[:,1],-1)-xy[:,1]*np.roll(xy[:,0],-1)))
        assert abs(area-polygon_area)<1e-8*area,(area,polygon_area)
        return {'points':loop,'triangles':tri,'center':np.array(center),'normal':np.array(normal),'axis1':e1,'axis2':e2,'area':area}

def refine(tri):
    a,b,c=tri[:,0],tri[:,1],tri[:,2];ab=(a+b)/2;bc=(b+c)/2;ca=(c+a)/2
    return np.concatenate([np.stack(x,axis=1) for x in [(a,ab,ca),(ab,b,bc),(ca,bc,c),(ab,bc,ca)]])

def quadrature(tri):
    weights=.5*np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/3
    bary=np.array([[2/3,1/6,1/6],[1/6,2/3,1/6],[1/6,1/6,2/3]])
    return np.einsum('qa,tad->tqd',bary,tri).reshape(-1,3),np.repeat(weights,3)
