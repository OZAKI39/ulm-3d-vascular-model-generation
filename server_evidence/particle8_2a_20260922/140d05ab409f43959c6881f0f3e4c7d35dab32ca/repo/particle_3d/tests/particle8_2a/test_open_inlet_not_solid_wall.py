from types import SimpleNamespace
import numpy as np
from particle_3d.injection_admission import FiniteSizeAdmission
from particle_3d.wall_geometry import WallGeometry
from particle_3d.particle82a_geometry import cylinder_triangles


def test_current_checker_accepts_sphere_straddling_open_cap():
    wall, _ = cylinder_triangles()
    checker = FiniteSizeAdmission(None, None, wall=WallGeometry(wall))
    event = dict(species='MB', radius_m=1e-6, particle_id=1, q=[1, 0, 0, 0])
    p, status, _ = checker.check(event, np.zeros(3), {})
    assert status == 'ACCEPTED' and p is not None
    # The sphere extends upstream; the actual current checker allows this.
    assert p.position[2]-event['radius_m'] < 0


def test_positive_control_solid_cap_rejects_same_sphere():
    wall, cap = cylinder_triangles()
    checker = FiniteSizeAdmission(None, None, wall=WallGeometry(np.concatenate([wall, cap])))
    event = dict(species='MB', radius_m=1e-6, particle_id=1, q=[1, 0, 0, 0])
    assert checker.check(event, np.zeros(3), {})[1] == 'WALL_REJECTED'
