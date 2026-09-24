#!/usr/bin/env python3
"""Existing RK23 tracer resolution check on unresolved saved seeds; no new sampler."""
from particle_3d.routing_stationary_audit import *
from particle_3d.particle8_replay import REPO
from particle_3d.particle82_tracers import trace_positions

def main():
 R=REPO/'particle_3d/reports/particle9a1_routing_stationary_audit';env=environment();out=R/'diagnostic_outputs/point_resolution';out.mkdir(parents=True,exist_ok=True);result=[];start=time.time()
 rows=read(R/'data/point_completed.json')['results']
 # Unresolved seeds plus up to five already resolved seeds per outlet, selected by stable ID.
 selected=[r for r in rows if r['point_outlet']=='NO_EXIT']
 for role in ROLES[:3]:selected += [r for r in rows if r['point_outlet']==role][:5]
 for r in selected:
  pid=r['particle_id'];trials=[]
  for level,config in enumerate([dict(step_m=.1e-6,error=1e-12,horizon_m=2e-3),dict(step_m=.05e-6,error=1e-13,horizon_m=2e-3)]):
   t=trace_positions([r['anchor_m']],**config)[0];path=t.pop('path');f=out/f'point_{pid:06d}_level{level+1}.npz';np.savez_compressed(f,path=path)
   field=env.field.sample(path[-1,1:]);tri,dist=env.wall.nearest_center_triangle(path[-1,1:])
   trials.append(dict(config=config,point_outlet=t['outlet'] or 'NO_EXIT',end_reason=t['end_reason'],terminal_time_s=float(path[-1,0]),path_sha256=sha(f),terminal_position_m=path[-1,1:].tolist(),terminal_speed_m_s=float(np.linalg.norm(field.velocity_m_s)),wall_distance_m=dist,wall_triangle_id=tri,sample_count=len(path)))
  result.append(dict(particle_id=pid,accepted=r['accepted'],baseline_outlet=r['point_outlet'],baseline_reason=r['point_end_reason'],trials=trials));dump(R/'data/point_resolution_progress.json',result);print(pid,[(x['point_outlet'],x['end_reason']) for x in trials],flush=True)
 dump(R/'data/point_resolution.json',dict(results=result,wall_seconds=time.time()-start,role='AUDIT_ONLY_ORIGINAL_FEM_AND_TRACER; BASELINE_ROWS_PRESERVED',native_source_sha256=sha(REPO/'particle_3d/src/particle_3d/particle82_point_native.cpp')))
if __name__=='__main__':main()
