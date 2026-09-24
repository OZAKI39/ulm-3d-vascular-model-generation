import numpy as np
from particle_3d.particle65_cases import synthetic,MU
from particle_3d.particle65_motion import assemble_v1
from particle_3d.wall_gap import wall_gap
from particle_3d.nearfield_regularization import contract

def test_original_wall_and_geometry(policy):
 shapes,provider,wall,_=synthetic('wall');before=wall.triangles.copy();x=shapes[17].center_m.copy()
 free={17:np.r_[provider(17,shapes[17],0)[0],np.zeros(3)]};r=assemble_v1(shapes,free,MU,wall,policy=policy)
 assert np.array_equal(wall.triangles,before) and np.array_equal(shapes[17].center_m,x)
 assert r.blocks[0]['h_geom_m']==wall_gap(shapes[17],wall).gap_m
 assert contract()['GLYCOCALYX_WALL_MODEL']=='NOT_FROZEN' and contract()['geometry_offset_m']==0
