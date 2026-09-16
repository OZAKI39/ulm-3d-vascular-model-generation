from pathlib import Path
import numpy as np,json,sys,hashlib,time
P=Path(__file__).resolve().parents[1];sys.path.insert(0,str(P/'src'))
from reference_passive_transport import ReferenceField,WallReference,velocity
ct=json.loads((P/'contracts/PASSIVE_MICROBUBBLE_TRANSPORT_V0_CONTRACT.json').read_text());pop=json.loads((P/'provenance/POPULATION_IDENTITY.json').read_text());f=ReferenceField(P/'fields/FROZEN_FLOW_FIELD_V0.h5');wall=WallReference(P/'provenance/closed_geometry_m.stl');dx=f.dx
ids=f.interior_lower_ids();points=f.xyz(ids)+.5*dx;vel,status=f.sample(points);assert np.all(status==0);distances=np.array([wall.gap(x) for x in points]);speeds=np.linalg.norm(vel,axis=1);order=np.lexsort((ids,-speeds))
np.savez_compressed(P/'validation/GEOMETRY_CANDIDATES.npz',ids=ids,xyz_m=points,center_wall_m=distances,speed_m_s=speeds)
d=pop['case2']['diameters_um'][0]*1e-6;eligible=order[distances[order]>d/2+5*dx];eligible=eligible[wall.enclosed(points[eligible])];dt=.125*ct['field_speed_statistics']['real']['dt_base_s'];n=ct['case_duration_policy']['real_coarse_steps']*4;attempts=[];chosen=None
print('CANDIDATES',len(ids),'SAFE_FOR_CASE2',len(eligible),'diameter_um',d/1e-6,flush=True)
for index in eligible[:64]:
 x0=points[index].copy();x=x0.copy();reason='MAX_TIME';maxdisp=0.;steps=0
 for j in range(n):
  try:
   u=velocity(f,x);mid=x+.5*dt*u
   if wall.gap(mid)<d/2+3*dx:reason='WALL_SAFETY';break
   k2=velocity(f,mid);new=x+dt*k2
   if not wall.segment_safe(x,new,d/2+3*dx):reason='WALL_SAFETY';break
   velocity(f,new)
  except ValueError:reason='INVALID_QUERY';break
  x=new;steps=j+1;maxdisp=float(np.linalg.norm(x-x0))
  if maxdisp>=10*dx:reason='TEN_DX_REACHED';break
 attempts.append({'lower_cell_id':int(ids[index]),'candidate_index':int(index),'steps':steps,'displacement_m':maxdisp,'reason':reason})
 print('PATH',attempts[-1],flush=True)
 if reason=='TEN_DX_REACHED':chosen=int(index);break
safe={'status':'PASS' if chosen is not None else 'BLOCKED_SHORT_SAFE_PATH','diameter_um':d/1e-6,'source_bubble_id':pop['case2']['selected_source_ids'][0],'attempts':attempts,'candidate_count':len(ids),'safe_start_count':len(eligible),'search_did_not_change_diameter_seed_margins':True}
if chosen is not None:safe.update(spawn_m=points[chosen].tolist(),center_wall_distance_m=float(distances[chosen]),initial_velocity_m_s=vel[chosen].tolist(),lower_cell_id=int(ids[chosen]))
(P/'contracts/CASE2_GEOMETRY_CONTRACT.json').write_text(json.dumps(safe,indent=2)+'\n')
# Fixed original samples, largest-first deterministic placement. No new diameter draws.
ds=np.array(pop['case5']['diameters_um'])*1e-6;placed=[];skipped=[]
for particle in np.lexsort((np.arange(len(ds)),-ds)):
 options=order[distances[order]>=ds[particle]/2+5*dx];options=options[wall.enclosed(points[options])];selected=None
 for idx in options:
  if all(np.linalg.norm(points[idx]-np.array(q['position_m']))>=(ds[particle]+q['diameter_m'])/2 for q in placed):selected=int(idx);break
 if selected is None:skipped.append(int(particle));continue
 placed.append({'source_bubble_id':int(particle),'diameter_m':float(ds[particle]),'position_m':points[selected].tolist(),'center_wall_distance_m':float(distances[selected]),'lower_cell_id':int(ids[selected])})
placed.sort(key=lambda q:q['source_bubble_id']);p5={'status':'PASS' if len(placed)>=3 else 'BLOCKED_INSUFFICIENT_SAFE_STARTS','N_requested':10,'N':len(placed),'placements':placed,'unplaced_source_ids':skipped,'method':ct['geometry_search']['case5_packing'],'frozen_before_case5_run':True}
(P/'contracts/CASE5_GEOMETRY_CONTRACT.json').write_text(json.dumps(p5,indent=2)+'\n')
# Independent distances to validate new C++ wall query without using it for truth.
rng=np.random.default_rng(2026092800);pick=rng.choice(len(ids),2000,replace=False);xyz=points[pick]+rng.uniform(-.35,.35,(len(pick),3))*dx;dd=np.array([wall.gap(x) for x in xyz]);np.savetxt(P/'validation/WALL_QUERY_REFERENCE.csv',np.column_stack((np.arange(len(xyz)),xyz,dd)),delimiter=',',header='id,x,y,z,wall_distance_m',comments='',fmt='%.17g')
print(json.dumps({'CASE2':safe,'CASE5':p5},indent=2))
