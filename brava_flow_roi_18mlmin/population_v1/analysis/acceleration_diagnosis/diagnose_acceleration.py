"""Frozen-field and original CPU physics replay at measured speedup states; no integration."""
from pathlib import Path
import numpy as np,json,csv,sys
out=Path(__file__).resolve().parent;root=out.parents[1];base=root/'analysis/realization_0000_contactfix'
sys.path.insert(0,'/home/lzy/projects/ulm_particle_formal_p9a5/particle_3d/src')
from particle_3d.population_v1.environment import load_environment
from particle_3d.population_v1.components import assemble_components,solve_components
from particle_3d.population_v1.pairs import candidates,sphere_gap
from particle_3d.particle6_stepper import bind_query_dependency
from particle_3d.nearfield_handoff import handoff_constraints
from particle_3d.nearfield_regularization import NearFieldRegularizationV1
from particle_3d.particle_shapes import Sphere
from particle_3d.wall_gap import wall_gap
from particle_3d.planar_wall_hydrodynamics import wall_weight
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
policy=NearFieldRegularizationV1();env=load_environment(json.loads((root/'data/environment.local.json').read_text()));print('ENV_LOADED',flush=True)
events=list(csv.DictReader((out/'visible_speedup_candidates.csv').open()));catalog={int(r['particle_id']):r for r in csv.DictReader((base/'trajectory_catalog.csv').open())};pairs=list(csv.DictReader((base/'pair_interactions.csv').open()))
pair_by_id={}
for r in pairs:
 for k in ['particle_i','particle_j']:pair_by_id.setdefault(int(r[k]),[]).append(r)
geom=json.loads((root.parent/'inputs/geometry.json').read_text());direction=np.array([-.72,.69,.11]);direction/=np.linalg.norm(direction)
summary=[];replays=[];selected=[636,1189,393,1329];curves={}
for n,row in enumerate(events):
 pid=int(row['particle_id']);mid=float(row['physical_time_s']);win=float(row['window_s']);radius=float(catalog[pid]['radius_m'])
 with np.load(base/'tracks'/f'mb_{pid:06d}.npz') as z:
  t=z['physical_time_s'];x=z['x_m'];v=z['velocity_m_s'];gap=z['wall_gap_m'];tid=z['tetra_id']
  # Saved velocity at row i+1 was evaluated at position/time of row i.
  tc=t[:-1];vs=v[1:];nodes=env.field.tetra[tid[:-1]];G=env.field.gradients_s_inv[tid[:-1]]
  uf=env.field.velocity_nodes_m_s[nodes[:,0]]+np.einsum('nij,nj->ni',G,x[:-1]-env.field.points_m[nodes[:,0]])
  speed=np.linalg.norm(vs,axis=1);fluid=np.linalg.norm(uf,axis=1)
  masks=[(tc>=mid-win)&(tc<mid),(tc>=mid)&(tc<mid+win)]
  record=dict(row)
  for tag,mask in zip(['before','after'],masks):
   record.update({tag+'_gap_median_m':float(np.median(gap[:-1][mask])),tag+'_fluid_speed_median_m_s':float(np.median(fluid[mask])),tag+'_raw_bubble_speed_median_m_s':float(np.median(speed[mask])),tag+'_wall_weight_median':float(np.median([wall_weight(h/radius) for h in gap[:-1][mask]]))})
  record['radius_m']=radius;record['status']=catalog[pid]['status'];record['outlet']=catalog[pid]['outlet'];record['pair_rows_in_window']=sum(mid-win<=float(r['physical_time_s'])<=mid+win for r in pair_by_id.get(pid,[]));record['pair_ever']=pid in pair_by_id
  record['slow_position_m']=x[np.searchsorted(t,mid-win)].tolist();record['fast_position_m']=x[np.searchsorted(t,mid+win)].tolist()
  summary.append(record)
  if pid not in selected:continue
  mask=(tc>=mid-.16)&(tc<=mid+.16);ix=np.flatnonzero(mask)
  table=[]
  for k in ix:
   r=dict(time_s=float(tc[k]),video_time_s=float(tc[k]/(.5526601319117461)),bubble_speed_mm_s=float(speed[k]*1e3),fluid_speed_mm_s=float(fluid[k]*1e3),wall_gap_um=float(gap[k]*1e6),tetra_id=int(tid[k]))
   for pose in ['0','15']:
    R=np.array(geom['pose_transforms'][pose]['source_mm_to_candidate_mm'])[:3,:3];u=R@vs[k]
    r['pose'+pose+'_projected_speed_mm_s']=float(np.linalg.norm(u-direction*np.dot(u,direction))*1e3)
   table.append(r)
  with (out/f'case_{pid}_timeseries.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
  curves[pid]=(mid,table,radius)
  for target in [mid-win*.8,mid,mid+win*.8]:
   k=int(np.clip(np.searchsorted(tc,target),0,len(tc)-1));p=x[k];sample=env.field.sample(p);s={pid:Sphere(p,radius)};free={pid:np.r_[sample.velocity_m_s,.5*sample.vorticity_s_inv]}
   pair_now=[r for r in pair_by_id.get(pid,[]) if abs(float(r['physical_time_s'])-t[k+1])<1e-10]
   if pair_now:raise RuntimeError(('Selected replay needs coupled partner',pid,tc[k]))
   system=assemble_components(s,free,env.mu,env.wall,policy=policy,gradient_at=lambda p:env.field.sample(p).velocity_gradient_s_inv)
   constraints,_=bind_query_dependency(handoff_constraints,dict(all_pairs=candidates,pair_gap=sphere_gap))(s,env.wall,policy)
   sol=solve_components(system,particles=s,constraints=constraints);wg=wall_gap(s[pid],env.wall)
   replay=dict(particle_id=pid,time_s=float(tc[k]),position_m=p.tolist(),tetra=int(tid[k]),triangle=int(wg.wall_triangle_id),
    radius_m=radius,gap_m=float(wg.gap_m),wall_weight=system.systems[0].planar_diagnostics[0]['wall_weight'],
    normal=wg.normal_inward.tolist(),gradient_norm_s_inv=float(np.linalg.norm(sample.velocity_gradient_s_inv)),
    saved_speed_mm_s=float(speed[k]*1e3),fluid_speed_mm_s=float(fluid[k]*1e3),unconstrained_speed_mm_s=float(np.linalg.norm(sol.unconstrained[:3])*1e3),
    original_CPU_speed_mm_s=float(np.linalg.norm(sol.velocity[:3])*1e3),CPU_recorded_error_m_s=float(np.max(abs(sol.velocity[:3]-vs[k]))),
    contact_count=len(constraints),contact_ids=[c.canonical_id for c in constraints],pair_interaction=False,solver=sol.record)
   replays.append(replay);print('REPLAY',pid,round(tc[k],6),'bubble/fluid',round(speed[k]*1e3,6),round(fluid[k]*1e3,6),'gap_um',round(wg.gap_m*1e6,8),'contacts',len(constraints),'error',replay['CPU_recorded_error_m_s'],flush=True)
 if n%50==0:print('EVENTS',n,flush=True)
with (out/'speedup_diagnosis_all_candidates.csv').open('w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
(out/'CPU_STATE_REPLAY.json').write_text(json.dumps(replays,indent=2)+'\n')
fig,ax=plt.subplots(4,2,figsize=(13,12),layout='constrained')
for j,pid in enumerate(selected):
 mid,rows,radius=curves[pid];tt=np.array([r['video_time_s'] for r in rows]);left=ax[j,0];right=ax[j,1]
 left.plot(tt,[r['bubble_speed_mm_s'] for r in rows],label='Microbubble speed',color='#cc5d27');left.plot(tt,[r['fluid_speed_mm_s'] for r in rows],label='FEM fluid speed at center',color='#167ca0',ls='--');left.set_yscale('log');left.set_ylabel('Speed (mm/s)');left.set_title(f'ID {pid}: same recorded trajectory, CPU replay checked');left.legend(fontsize=8)
 right.plot(tt,[r['wall_gap_um'] for r in rows],color='#167ca0');right.axhline(radius*1e6,color='k',ls='--',label='Gap = bubble radius');right.axhline(.002,color='gray',ls=':',label='2 nm handoff');right.set_yscale('log');right.set_ylabel('Surface-to-wall gap (um)');right.legend(fontsize=8)
 for a in [left,right]:a.set_xlabel('Video playback time (s)');a.axvline(mid/.5526601319117461,color='gray',ls=':',alpha=.5)
fig.savefig(out/'speedup_examples.png',dpi=300);fig.savefig(out/'speedup_examples.pdf');plt.close(fig)
agg=dict(candidates=len(summary),after_video_3s=sum(float(r['video_time_s'])>=3 for r in summary),
 no_pair_ever=sum(not r['pair_ever'] for r in summary),no_pair_during_speedup=sum(r['pair_rows_in_window']==0 for r in summary),
 slow_wall_gap_below_radius=sum(r['before_gap_median_m']<r['radius_m'] for r in summary),
 fast_gap_above_radius=sum(r['after_gap_median_m']>r['radius_m'] for r in summary),
 maximum_CPU_replay_velocity_error_m_s=max(r['CPU_recorded_error_m_s'] for r in replays),selected=selected,
 analysis_scope='Recorded-state replay and same-flow sampling only, no time integration; selected events are an explicit threshold subset, not all perceived speedups')
(out/'DIAGNOSIS_SUMMARY.json').write_text(json.dumps(agg,indent=2)+'\n');print(json.dumps(agg),flush=True)
