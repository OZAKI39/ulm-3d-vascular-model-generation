"""STL precision cleanup for newly designed support, never source anatomy."""
import numpy as np,io
import trimesh
import manifold3d
from .support_custom_qc import canonical_stl

def exact_stl(mesh):
    result=trimesh.load_mesh(io.BytesIO(mesh.export(file_type='stl')),file_type='stl')
    # A positive-area skinny seam is not a degenerate triangle. Deleting all
    # faces below trimesh's default height threshold opens otherwise closed
    # Boolean seams. Remove only genuinely zero-area triangles.
    nonzero=result.area_faces>0.
    if not nonzero.all():result.update_faces(nonzero);result.remove_unreferenced_vertices()
    if not result.is_volume:raise ValueError('EXPORTED_SUPPORT_STL_NOT_CLOSED_VOLUME')
    return result

def support_stl(mesh):
    try:
        return exact_stl(mesh)
    except ValueError:
        solid=manifold3d.Manifold(manifold3d.Mesh(np.asarray(mesh.vertices,dtype=np.float32),np.asarray(mesh.faces,dtype=np.uint32)))
        if solid.is_empty():
            raise ValueError('SUPPORT_UNION_NOT_MANIFOLD')
        for tolerance in [0.00002,0.00005,0.0001,0.0005,0.001]:
            m=solid.simplify(tolerance).to_mesh()
            candidate=trimesh.Trimesh(m.vert_properties[:,:3],m.tri_verts,process=False)
            try:
                result=exact_stl(candidate)
                result.metadata['support_only_numerical_cleanup_tolerance_mm']=tolerance
                return result
            except ValueError:
                continue
        raise ValueError('SUPPORT_STL_PRECISION_REPAIR_FAILED')
