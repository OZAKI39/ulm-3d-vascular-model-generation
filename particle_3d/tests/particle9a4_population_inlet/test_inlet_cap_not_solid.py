import numpy as np
from particle_3d.injection_admission import FiniteSizeAdmission
from particle_3d.wall_geometry import WallGeometry
from particle_3d.particle82a_geometry import cylinder_triangles

def test_same_sphere_with_solid_cap_rejects():
 wall,cap=cylinder_triangles();check=FiniteSizeAdmission(None,None,wall=WallGeometry(np.concatenate([wall,cap])))
 assert check.check(dict(species='MB',radius_m=1e-6,particle_id=1,q=[1,0,0,0]),[0,0,0],{})[1]=='WALL_REJECTED'

def test_inlet_cap_not_in_wall_bvh(real_env):
 wall={tuple(sorted(x)) for x in real_env.wall.global_node_ids}
 cap=real_env.boundaries['INLET'];ids=np.asarray(cap.point_data['GlobalNodeID'],int)-1
 assert not wall&{tuple(sorted(x)) for x in ids[cap.faces.reshape(-1,4)[:,1:]]}
 assert real_env.wall.provenance['solid_boundaries']==['WALL']

def test_center_on_inlet_field_sample_is_valid(real_env):
 for tri in real_env.sampler.triangles:assert real_env.field.sample(tri.mean(0)).inside_lumen

def test_edge_and_vertex_boundary_ownership_stable(real_env):
 tri=real_env.sampler.triangles
 pts=np.vstack([tri[:,0],(tri[:,0]+tri[:,1])/2])
 first=[real_env.field.locate(p)[0] for p in pts]
 second=[real_env.field.locate(p)[0] for p in pts[::-1]][::-1]
 assert first==second and min(first)>=0
