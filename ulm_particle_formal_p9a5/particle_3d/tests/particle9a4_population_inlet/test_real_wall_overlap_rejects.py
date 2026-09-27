import numpy as np
from particle_3d.injection_admission import FiniteSizeAdmission
from particle_3d.wall_geometry import WallGeometry
from particle_3d.particle82a_geometry import cylinder_triangles

def test_rim_wall_overlap_rejects():
 wall,_=cylinder_triangles(radius=5e-6);check=FiniteSizeAdmission(None,None,wall=WallGeometry(wall))
 p,status,_=check.check(dict(species='MB',radius_m=1e-6,particle_id=1,q=[1,0,0,0]),[4.8e-6,0,0],{})
 assert p is None and status=='WALL_REJECTED'
