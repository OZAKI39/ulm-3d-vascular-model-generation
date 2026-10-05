"""Optional static-mesh query implementations; the original VTK path stays available."""
from collections import OrderedDict
from time import perf_counter

import numpy as np


class PersistentVTKWallBackend:
    """One vtkStaticCellLocator, reusable query objects, signed mesh distance.

    Uses VTK triangle closest points and angle-weighted pseudonormals at shared
    edges/vertices, following the documented vtkImplicitPolyDataDistance semantics.
    This is an experimental query backend until the separate trajectory gate passes.
    No tree, triangle-distance solver or physical projection tolerance is replaced.
    Instances are single-worker objects, not shared between threads.
    """
    def __init__(self, surface):
        import vtk
        self.vtk = vtk
        begin=perf_counter()
        clean=vtk.vtkCleanPolyData()
        clean.SetInputData(surface)
        triangles=vtk.vtkTriangleFilter()
        triangles.PassVertsOff()
        triangles.PassLinesOff()
        triangles.SetInputConnection(clean.GetOutputPort())
        triangles.Update()
        self.mesh=vtk.vtkPolyData()
        self.mesh.DeepCopy(triangles.GetOutput())
        self.mesh.BuildLinks()
        self.locator=vtk.vtkStaticCellLocator()
        self.locator.SetDataSet(self.mesh)
        self.locator.BuildLocator()
        self.cell=vtk.vtkGenericCell()
        self.adjacent=vtk.vtkGenericCell()
        self.ids=vtk.vtkIdList()
        self.cell_id=vtk.reference(0)
        self.sub_id=vtk.reference(0)
        self.dist2=vtk.reference(0.)
        self.evaluated_dist2=vtk.reference(0.)
        self.raw_point=[0.,0.,0.]
        self.closest=[0.,0.,0.]
        self.pcoords=[0.,0.,0.]
        self.weights=[0.,0.,0.]
        self.normals=np.empty((self.mesh.GetNumberOfCells(),3))
        existing=self.mesh.GetCellData().GetNormals()
        normal=[0.,0.,0.]
        for index in range(len(self.normals)):
            if existing is not None:
                existing.GetTuple(index,normal)
            else:
                self.mesh.GetCell(index,self.adjacent)
                vtk.vtkPolygon.ComputeNormal(self.adjacent.GetPoints(),normal)
            self.normals[index]=normal
        self.pseudonormals={}
        self.build_seconds=perf_counter()-begin
        self.calls=self.points_queried=0
        self.query_seconds=0.

    def _normal(self, face):
        zeros=[i for i,w in enumerate(self.weights) if abs(w)<1e-12]
        if not zeros:
            return self.normals[face]
        ids=[self.cell.GetPointId(i) for i in range(3)]
        if len(zeros)==1:
            edge=zeros[0]
            a,b=ids[(edge+1)%3],ids[(edge+2)%3]
            key=("edge",min(a,b),max(a,b))
            if key not in self.pseudonormals:
                self.mesh.GetCellEdgeNeighbors(-1,a,b,self.ids)
                normal=np.zeros(3)
                for i in range(self.ids.GetNumberOfIds()):
                    normal+=self.normals[self.ids.GetId(i)]
                self.vtk.vtkMath.Normalize(normal)
                self.pseudonormals[key]=normal
            return self.pseudonormals[key]
        if len(zeros)==2:
            vertex=next(ids[i] for i,w in enumerate(self.weights) if abs(w)>1e-12)
            key=("vertex",vertex)
            if key not in self.pseudonormals:
                self.mesh.GetPointCells(vertex,self.ids)
                normal=np.zeros(3)
                point=np.asarray(self.mesh.GetPoint(vertex))
                for i in range(self.ids.GetNumberOfIds()):
                    adjacent_id=self.ids.GetId(i)
                    self.mesh.GetCell(adjacent_id,self.adjacent)
                    b,c=self.adjacent.GetPointId(0),self.adjacent.GetPointId(1)
                    if vertex==b:
                        b=self.adjacent.GetPointId(2)
                    elif vertex==c:
                        c=self.adjacent.GetPointId(2)
                    pb=np.asarray(self.mesh.GetPoint(b))-point
                    pc=np.asarray(self.mesh.GetPoint(c))-point
                    self.vtk.vtkMath.Normalize(pb)
                    self.vtk.vtkMath.Normalize(pc)
                    alpha=np.arccos(self.vtk.vtkMath.Dot(pb,pc))
                    normal+=alpha*self.normals[adjacent_id]
                self.vtk.vtkMath.Normalize(normal)
                self.pseudonormals[key]=normal
            return self.pseudonormals[key]
        raise RuntimeError("Degenerate closest triangle: no signed distance accepted")

    def query(self, points):
        begin=perf_counter()
        x=np.atleast_2d(np.asarray(points,float))
        distances=np.empty(len(x))
        closest=np.empty_like(x)
        gradients=np.empty_like(x)
        cell_ids=np.empty(len(x),np.int64)
        for i,point in enumerate(x):
            self.locator.FindClosestPoint(point,self.raw_point,self.cell,
                self.cell_id,self.sub_id,self.dist2)
            face=int(self.cell_id)
            if face<0:
                raise RuntimeError("Static wall closest point unavailable")
            # The locator may use the generic cell as scratch for a final rejected
            # candidate. Bind the returned winning ID before barycentric evaluation.
            self.mesh.GetCell(face,self.cell)
            distance=np.sqrt(float(self.dist2))
            gradient=(np.asarray(self.raw_point)-point)/(distance if distance else 1.)
            self.cell.EvaluatePosition(self.raw_point,self.closest,self.sub_id,self.pcoords,
                self.evaluated_dist2,self.weights)
            normal=self._normal(face)
            if distance==0.:
                gradient=normal.copy()
            distance*=1. if self.vtk.vtkMath.Dot(gradient,normal)<0. else -1.
            if distance>0.:
                gradient=-gradient
            distances[i],closest[i],gradients[i],cell_ids[i]=distance,self.closest,gradient,face
        self.calls+=1
        self.points_queried+=len(x)
        self.query_seconds+=perf_counter()-begin
        return distances,closest,gradients,cell_ids


class ExactCoordinateWallCache:
    """Bounded cache of original double-precision query outputs; no coordinate bins.

    Exact bit-pattern keys safely survive ragged projection subsets. Cache entries
    are invalidated if the bound VTK distance object changes. This supplements the
    original previous-batch cache, whose shape/values and certificate remain intact.
    """
    def __init__(self, geometry, capacity=16384, query_backend=None):
        self.geometry=geometry
        self.capacity=capacity
        self.query_backend=query_backend
        self.cache=OrderedDict()
        self.mesh=None
        self.hits=self.misses=0

    def level_and_gradient(self, positions):
        g=self.geometry
        x=np.atleast_2d(np.asarray(positions,float))
        if g._mesh_distance is not self.mesh:
            self.cache.clear()
            self.mesh=g._mesh_distance
        previous=g._distance_cache
        if previous is not None and previous[0].shape==x.shape:
            changed=np.any(x!=previous[0],axis=1)
            values,gradients=previous[1].copy(),previous[2].copy()
            g.distance_cache_hits+=int((~changed).sum())
        else:
            changed=np.ones(len(x),bool)
            values,gradients=np.empty(len(x)),np.empty_like(x)
        g.distance_cache_misses+=int(changed.sum())
        missing=[]
        for row in np.flatnonzero(changed):
            key=x[row].tobytes()
            hit=self.cache.get(key)
            if hit is not None:
                values[row],gradients[row]=hit
                self.hits+=1
            else:
                missing.append((row,key))
        if missing:
            if self.query_backend is not None:
                distances,_,normals,_=self.query_backend.query(x[[row for row,_ in missing]])
            for index,(row,key) in enumerate(missing):
                if self.query_backend is None:
                    value=g._mesh_distance.EvaluateFunction(x[row])
                    normal=[0.,0.,0.]
                    g._mesh_distance.EvaluateGradient(x[row],normal)
                else:
                    value,normal=distances[index],normals[index]
                values[row],gradients[row]=value,normal
                self.cache[key]=(float(value),np.array(normal,copy=True))
                self.misses+=1
                if len(self.cache)>self.capacity:
                    self.cache.popitem(last=False)
        g._distance_cache=(x.copy(),values.copy(),gradients.copy())
        return values,gradients
