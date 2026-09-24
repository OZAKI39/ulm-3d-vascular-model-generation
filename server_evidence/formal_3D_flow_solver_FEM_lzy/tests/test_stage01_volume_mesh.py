import json
import numpy as np
from fem3d.mesh_qc import tetra_volumes, boundary_partition


def test_tetrahedra_are_finite_positive_and_connected(stage01_mesh):
    q=stage01_mesh["qc"]
    assert q["topology"]["connected_fluid_components"]==1
    assert q["validity"]=={"zero_volume":0,"negative_volume":0,"nonfinite_volume":0}
    assert q["quality"]["gmsh_min_sicn"]["minimum"]>0
    assert q["volume_closure"]["relative_error"]<1e-10
    assert q["cell_tags"]=={"100":q["tetrahedron_count"]}
    metadata=json.loads((stage01_mesh["base"]/"metadata/meshing.json").read_text())
    assert metadata["actual_tetra_count"]==q["tetrahedron_count"]
    assert metadata["actual_vertex_count"]==q["vertex_count"]
    assert metadata["exit_status"]==0 and metadata["gpu_used"] is False
    assert metadata["threads"]==metadata["mpi_ranks"]==1


def test_signed_volume_detects_inversion_and_collapse():
    points=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],dtype=float)
    volume=tetra_volumes(points,np.array([[0,1,2,3],[0,2,1,3],[0,1,2,2]]))
    assert volume[0]>0 and volume[1]<0 and volume[2]==0


def test_disconnected_components_are_not_mistaken_for_one():
    points=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]],dtype=float)
    points=np.vstack([points,points+3])
    faces=np.array([[1,2,3],[0,2,3],[0,1,3],[0,1,2]])
    _,_,q=boundary_partition(points,np.array([[0,1,2,3],[4,5,6,7]]),np.vstack([faces,faces+4]),np.ones(8,dtype=int))
    assert q["connected_fluid_components"]==2
