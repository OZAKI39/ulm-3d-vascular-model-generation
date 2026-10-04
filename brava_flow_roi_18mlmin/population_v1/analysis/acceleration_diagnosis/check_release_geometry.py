from pathlib import Path
import json,sys,numpy as np,csv
out=Path(__file__).resolve().parent;root=out.parents[1];base=root/'analysis/realization_0000_contactfix'
sys.path.insert(0,'/home/lzy/projects/ulm_particle_formal_p9a5/particle_3d/src')
from particle_3d.population_v1.environment import load_environment
from particle_3d.particle_shapes import Sphere
from particle_3d.wall_gap import wall_gap
from particle_3d.convex_triangle import triangle_closest_many
from particle_3d.nearfield_handoff import handoff_constraints
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
env=load_environment(json.loads((root/'data/environment.local.json').read_text()));policy=NearFieldRegularizationV1();records=[]
for r in json.loads((out/'CPU_STATE_REPLAY.json').read_text()):
 p=np.array(r['position_m']);f=env.field.sample(p);n=np.array(r['normal']);G=f.velocity_gradient_s_inv
 ids=env.field.tetra[f.tetra_id];record=dict(particle_id=r['particle_id'],time_s=r['time_s'],tetra=f.tetra_id,normal_fluid_speed_mm_s=float(f.velocity_m_s@n*1000),tangent_fluid_speed_mm_s=float(np.linalg.norm(f.velocity_m_s-(f.velocity_m_s@n)*n)*1000),
  divergence_s_inv=float(np.trace(G)),gradient_norm_s_inv=float(np.linalg.norm(G)),divergence_over_gradient_norm=float(np.trace(G)/np.linalg.norm(G)),tetra_vertex_speeds_m_s=np.linalg.norm(env.field.velocity_nodes_m_s[ids],axis=1).tolist(),
  tetra_edge_lengths_um=[float(np.linalg.norm(env.field.points_m[ids[i]]-env.field.points_m[ids[j]])*1e6) for i in range(4) for j in range(i)],
  nearest_feature=wall_gap(Sphere(p,r['radius_m']),env.wall).wall_feature)
 records.append(record)
release=[]
for pid in [636,1189,393]:
 r=next(r for r in json.loads((out/'CPU_STATE_REPLAY.json').read_text()) if r['particle_id']==pid);radius=r['radius_m'];z=np.load(base/'tracks'/f'mb_{pid:06d}.npz');t=z['physical_time_s'];x=z['x_m'];gap=z['wall_gap_m'];vel=z['velocity_m_s']
 mid=next(float(row['physical_time_s']) for row in csv.DictReader((out/'visible_speedup_candidates.csv').open()) if int(row['particle_id'])==pid)
 possible=np.flatnonzero((t>mid-.1)&(t<mid+.1)&(gap>2e-9+10*env.wall.roundoff_m))
 k=int(possible[0]);ix=np.unique(np.clip([k-3,k-2,k-1,k,k+1,k+3,k+10],0,len(t)-2))
 for i in ix:
  g=wall_gap(Sphere(x[i],radius),env.wall);f=env.field.sample(x[i]);tri=env.wall.triangles[g.wall_triangle_id];_,bary=triangle_closest_many(x[i],tri[None]);contacts,_=handoff_constraints({pid:Sphere(x[i],radius)},env.wall,policy)
  release.append(dict(particle_id=pid,time_s=float(t[i]),video_s=float(t[i]/.5526601319117461),wall_triangle=g.wall_triangle_id,wall_feature=g.wall_feature,
   barycentric=bary[0].tolist(),gap_nm=float(g.gap_m*1e9),tetra=f.tetra_id,contact_count=len(contacts),normal=g.normal_inward.tolist(),
   bubble_speed_mm_s=float(np.linalg.norm(vel[i+1])*1e3),fluid_speed_mm_s=float(np.linalg.norm(f.velocity_m_s)*1e3),
   normal_bubble_speed_mm_s=float(vel[i+1]@g.normal_inward*1e3),normal_fluid_speed_mm_s=float(f.velocity_m_s@g.normal_inward*1e3)))
from itertools import combinations
cell_id=388769;nodes=env.field.tetra[cell_id];pts=env.field.points_m[nodes];uu=env.field.velocity_nodes_m_s[nodes]
direct_G=np.linalg.solve(pts[1:]-pts[0],uu[1:]-uu[0]).T
volume=abs(np.linalg.det(pts[1:]-pts[0]))/6;flux=0.
for face in combinations(range(4),3):
 ix=np.array(face);p=pts[ix];area_normal=np.cross(p[1]-p[0],p[2]-p[0])/2
 if area_normal@(p.mean(0)-pts.mean(0))<0:area_normal=-area_normal
 flux+=float(uu[ix].mean(0)@area_normal)
independent=dict(tetra_id=cell_id,direct_nodal_linear_solve_gradient_max_error_s_inv=float(np.max(abs(direct_G-env.field.gradients_s_inv[cell_id]))),
 divergence_from_direct_solve_s_inv=float(np.trace(direct_G)),divergence_from_four_face_fluxes_s_inv=float(flux/volume),
 tetra_volume_m3=float(volume),net_outward_flux_m3_s=flux,net_flux_over_inlet_Q=float(flux/3e-7),
 interpretation='This is the actual frozen P1 field sampled by microbubbles; not a recomputation of the original CFD basis or a global boundary mass-balance claim')
result=dict(state_near_wall_FEM=records,release_closeup=release,independent_local_divergence_check=independent)
(out/'WALL_RELEASE_AND_FEM.json').write_text(json.dumps(result,indent=2)+'\n')
for r in records:print('FEM',json.dumps({k:v for k,v in r.items() if k!='tetra_edge_lengths_um'}),flush=True)
for r in release:print('RELEASE',json.dumps(r),flush=True)
