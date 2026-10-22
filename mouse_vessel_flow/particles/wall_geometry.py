"""Stationary zero-thickness WALL triangles and BVH candidate enumeration only."""
from pathlib import Path
import json
import numpy as np
import pyvista as pv
import vtk
from .audit import check_hash
from .particle_shapes import roundoff_length


class WallGeometry:
    def __init__(self,triangles_m,*,provenance=None):
        triangles=np.asarray(triangles_m,dtype=float)
        if triangles.ndim!=3 or triangles.shape[1:]!=(3,3) or not len(triangles) or not np.isfinite(triangles).all():
            raise ValueError('nonempty finite wall triangles required')
        cross=np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]);sizes=np.linalg.norm(cross,axis=1)
        if np.any(sizes==0):raise ValueError('degenerate wall triangle')
        self.triangles=np.frombuffer(triangles.tobytes(),dtype=float).reshape(triangles.shape)
        self.normal_out=cross/sizes[:,None];self.normal_in=-self.normal_out;self.areas_m2=sizes/2
        self.triangle_centers=triangles.mean(axis=1)
        self.triangle_bounding_radii=np.max(np.linalg.norm(triangles-self.triangle_centers[:,None,:],axis=2),axis=1)
        self.provenance=provenance or {'role':'SYNTHETIC_VALIDATION_GEOMETRY'}
        self.roundoff_m=roundoff_length(triangles)
        faces=np.column_stack((np.full(len(triangles),3),np.arange(3*len(triangles)).reshape(-1,3)))
        self._surface=pv.PolyData(triangles.reshape(-1,3),faces.ravel())
        self._locator=vtk.vtkStaticCellLocator();self._locator.SetDataSet(self._surface);self._locator.BuildLocator()

    @classmethod
    def from_frozen(cls,root):
        root=Path(root).resolve()
        mesh=json.loads((root/'frozen_reference/mesh_manifest.json').read_text())
        boundaries=json.loads((root/mesh['boundary_manifest']).read_text())
        spec=boundaries['boundaries']['WALL'];path=(root/spec['path']).resolve()
        if not path.is_relative_to(root) or spec['sv_face_id']!=1 or path.name!='WALL.vtp':
            raise ValueError('Only manifest WALL.vtp is solid')
        check_hash(path,spec['sha256'])
        wall=pv.read(path);faces=wall.faces.reshape(-1,4)
        if not np.all(faces[:,0]==3) or wall.n_cells!=spec['facets']:raise ValueError('WALL triangles changed')
        obj=cls(np.asarray(wall.points)[faces[:,1:]],provenance=dict(path=str(path),sha256=spec['sha256'],triangle_count=wall.n_cells,
                boundary_ids=mesh['face_ids'],solid_boundaries=['WALL'],wall_rigid=True,wall_stationary=True,physical_coating_thickness_m=0.))
        obj.global_node_ids=np.asarray(wall.point_data['GlobalNodeID'])[faces[:,1:]]-1
        from .open_boundary_rim import build_rim_topology
        caps={};hashes={'WALL':spec['sha256']}
        for role,cap in boundaries['boundaries'].items():
            if role!='INLET' and not role.startswith('OUTLET_'):continue
            cap_path=(root/cap['path']).resolve()
            if not cap_path.is_relative_to(root):raise ValueError('Cap outside frozen root')
            check_hash(cap_path,cap['sha256']);surface=pv.read(cap_path)
            cf=surface.faces.reshape(-1,4)
            if not np.all(cf[:,0]==3):raise ValueError('Triangular open cap required')
            caps[role]=np.asarray(surface.point_data['GlobalNodeID'])[cf[:,1:]]-1
            hashes[role]=cap['sha256']
        obj.open_boundary_topology=build_rim_topology(obj.global_node_ids,caps,hashes)
        return obj

    def candidates(self,center_m,radius_m):
        center=np.asarray(center_m,dtype=float);pad=roundoff_length(center,radius_m,self.triangles)
        bounds=np.column_stack((center-radius_m-pad,center+radius_m+pad)).ravel()
        cells=vtk.vtkIdList();self._locator.FindCellsWithinBounds(bounds,cells)
        return np.array(sorted(cells.GetId(i) for i in range(cells.GetNumberOfIds())),dtype=np.int64)

    def nearest_center_triangle(self,center_m):
        closest=[0.,0.,0.];cell=vtk.mutable(0);sub=vtk.mutable(0);distance2=vtk.mutable(0.)
        self._locator.FindClosestPoint(center_m,closest,cell,sub,distance2)
        # VTK chooses the candidate only; final finite-triangle distance is ours.
        from .convex_triangle import triangle_closest_many
        index=int(cell)
        point,_=triangle_closest_many(np.asarray(center_m),self.triangles[index:index+1])
        return index,float(np.linalg.norm(np.asarray(center_m)-point[0]))

    def swept_candidates(self,start,end,radius):
        center=.5*(np.asarray(start)+np.asarray(end))
        return self.candidates(center,radius+.5*np.linalg.norm(np.asarray(end)-start))
