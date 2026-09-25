import numpy as np
from particle_3d.injection_admission import FiniteSizeAdmission
from particle_3d.wall_geometry import WallGeometry
from particle_3d.particle82a_geometry import cylinder_triangles

def test_sphere_straddles_open_cap_accepts():
 wall,_=cylinder_triangles();check=FiniteSizeAdmission(None,None,wall=WallGeometry(wall))
 p,status,_=check.check(dict(species='MB',radius_m=1e-6,particle_id=1,q=[1,0,0,0]),np.zeros(3),{})
 assert p is not None and status=='ACCEPTED' and p.position[2]-p.radius<0
