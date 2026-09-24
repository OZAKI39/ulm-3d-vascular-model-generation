import numpy as np
from particle_3d.wall_geometry import WallGeometry
from particle_3d.particle82a_geometry import cylinder_triangles
from particle_3d.particle82a_entry_audit import swept_entry


def test_open_cap_crossing_swept_sphere_is_clear():
    wall=WallGeometry(cylinder_triangles(radius=5e-6)[0])
    result=swept_entry(np.array([[0,0,0],[0,0,1e-6],[0,0,2e-6]]),1e-6,wall)
    assert result['handoff_clear'] and result['full_sweep_evaluated']


def test_legal_endpoint_does_not_erase_illegal_anchor():
    wall=WallGeometry(cylinder_triangles(radius=5e-6)[0])
    result=swept_entry(np.array([[4.5e-6,0,0],[0,0,2e-6]]),1e-6,wall)
    assert not result['wall_geometrically_clear'] and not result['handoff_clear']
    assert result['first_failure']=='WALL_INTERSECTION_AT_ANCHOR'
