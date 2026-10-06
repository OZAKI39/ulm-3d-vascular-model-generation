from pathlib import Path
import hashlib,json
import h5py,numpy as np,pyvista as pv
from mevo_mb.bubbles.metrics import aggregation_metrics
from mevo_mb.coupling.bifurcation import BifurcationGeometry
root=Path.cwd()
delivery=root/'results/stage4/lubrication_upgrade/visualization/B2_B3_clean_enlarged_v8_20261006'
out=root/'results/stage4/lubrication_upgrade/visualization_audit/20261006_v8_geometry'
out.mkdir(parents=True,exist_ok=False)
manifest=json.loads((delivery/'RENDER_MANIFEST.json').read_text())
report={'scope':'Read-only check of every stored primary-only state, not unsaved accepted steps. Positive wall gap denotes an interior sphere. Display transform exactly matches the delivered 60-fold connected-component scaling.', 'magnification':60.0,'cases':{}}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
for case in ['B2','B3']:
 source=manifest['videos'][case+'_primary']['source'];run=Path(source['run'])
 identity=json.loads((run/'RUN_IDENTITY.json').read_text());cfg=identity['config']
 mesh_path=root/str(cfg['geometry']['mesh_surface']).removeprefix('/workspace/mevo_acoustic_microbubbles/')
 assert sha(mesh_path)==identity['inputs_sha256'][cfg['geometry']['mesh_surface']]
 assert sha(source['trajectory'])==source['trajectory_sha256']
 with h5py.File(source['trajectory'],'r') as h:
  t,x,ids,r,active=[h[k][:] for k in ['time_s','positions_m','particle_id','radii_m','active']]
 geom=BifurcationGeometry(cfg['geometry']);geom.bind_mesh_surface(pv.read(mesh_path))
 shown=np.full_like(x,np.nan);centers=np.full_like(x,np.nan);membership=[]
 for k in range(len(t)):
  ix=np.flatnonzero(active[k]);_,labels,_,_=aggregation_metrics(x[k,ix],r[ix],cfg['metrics'])
  members={}
  for label in np.unique(labels):
   group=ix[labels==label];center=x[k,group].mean(axis=0)
   shown[k,group]=center+60*(x[k,group]-center);centers[k,group]=center
   tup=tuple(int(v) for v in ids[group])
   for i in group:members[int(i)]=tup
  membership.append(members)
 radius=np.broadcast_to(r,active.shape)
 real_gap,_,valid=geom.wall_gap_and_normal(x[active],radius[active]);assert valid.all()
 display_gap,_,valid=geom.wall_gap_and_normal(shown[active],60*radius[active]);assert valid.all()
 real_center_distance=real_gap+radius[active]
 display_center_distance=display_gap+60*radius[active]
 real_grid=np.full(active.shape,np.nan);real_grid[active]=real_gap
 display_grid=np.full(active.shape,np.nan);display_grid[active]=display_center_distance
 common=active[1:] & active[:-1]
 physical=np.linalg.norm(np.diff(x,axis=0),axis=2)
 apparent=np.linalg.norm(np.diff(shown,axis=0),axis=2)
 changed=np.zeros_like(common)
 for k,i in zip(*np.where(common)):
  changed[k,i]=membership[k][int(i)]!=membership[k+1][int(i)]
 # The additive jump in display offset is separate from the actual physical displacement.
 offset=np.linalg.norm(np.diff(shown-x,axis=0),axis=2)
 candidates=np.argwhere(common & changed & (apparent>0.0002) & (physical<0.00001))
 candidates=sorted(candidates,key=lambda z:apparent[tuple(z)],reverse=True)
 examples=[]
 for k,i in candidates[:8]:
  examples.append({'particle_id':int(ids[i]),'time_start_s':float(t[k]),'time_end_s':float(t[k+1]),'video_time_s':float(t[k]/0.3),'physical_displacement_um':float(physical[k,i]*1e6),'display_displacement_um':float(apparent[k,i]*1e6),'display_offset_change_um':float(offset[k,i]*1e6),'group_size_before':len(membership[k][int(i)]),'group_size_after':len(membership[k+1][int(i)]),'group_members_before':membership[k][int(i)],'group_members_after':membership[k+1][int(i)],'real_positions_m':[x[k,i].tolist(),x[k+1,i].tolist()],'display_positions_m':[shown[k,i].tolist(),shown[k+1,i].tolist()],'real_surface_gaps_m':[float(real_grid[k,i]),float(real_grid[k+1,i])],'display_center_wall_distance_m':[float(display_grid[k,i]),float(display_grid[k+1,i])]})
 outside=np.argwhere(active & (display_grid<0))
 worst=outside[np.argmin(display_grid[tuple(outside.T)])] if len(outside) else None
 stats={'saved_states':len(t),'active_particle_state_samples':int(active.sum()),'raw_center_outside_count':int((real_center_distance<0).sum()),'raw_min_center_wall_distance_um':float(real_center_distance.min()*1e6),'raw_min_surface_gap_m':float(real_gap.min()),'raw_gap_below_minus_1pm_count':int((real_gap<-1e-12).sum()),'display_center_outside_count':int((display_center_distance<0).sum()),'frames_with_display_centers_outside':int(np.any(display_grid<0,axis=1).sum()),'display_max_center_outside_depth_um':float(max(0,-display_center_distance.min())*1e6),'display_sphere_crosses_wall_count':int((display_gap<0).sum()),'max_display_sphere_wall_protrusion_um':float(max(0,-display_gap.min())*1e6),'raw_inactive_to_active_transitions':int(((~active[:-1])&active[1:]).sum()),'membership_changed_particle_intervals':int(changed.sum()),'max_physical_step_um':float(physical[common].max()*1e6),'max_display_step_um':float(apparent[common].max()*1e6),'max_display_step_with_changed_membership_um':float(apparent[changed].max()*1e6),'small_real_motion_large_display_jump_events':len(candidates),'jump_examples':examples,'source':source,'mesh_path':str(mesh_path),'mesh_sha256':sha(mesh_path)}
 if worst is not None:
  k,i=worst;stats['worst_display_outside']={'time_s':float(t[k]),'video_time_s':float(t[k]/0.3),'particle_id':int(ids[i]),'group_size':len(membership[k][int(i)]),'real_surface_gap_um':float(real_grid[k,i]*1e6),'display_center_outside_depth_um':float(-display_grid[k,i]*1e6),'real_position_m':x[k,i].tolist(),'display_position_m':shown[k,i].tolist()}
 np.savez_compressed(out/f'{case}_display_audit.npz',time_s=t,particle_id=ids,active=active,physical_positions_m=x,display_positions_m=shown,group_centers_m=centers,physical_surface_gap_m=real_grid,display_center_wall_distance_m=display_grid,membership_changed=changed)
 report['cases'][case]=stats
 print(json.dumps({k:v for k,v in stats.items() if k not in ['jump_examples','source']},indent=2),flush=True)
 print('Jump example',json.dumps(examples[0] if examples else {}),flush=True)
report['render_code_evidence']={'group_transform':'scripts/stage4/render_single_overview.py: center + 60*(real-center), radius*60','regrouping':'Connected components are recomputed independently at each 20 ms saved state. A changed group changes the local magnification origin.','temporal_blend':'Primary frames blend two saved raster images at alpha 0, .25, .5, .75, without reconstructed intermediate positions.'}
(out/'AUDIT.json').write_text(json.dumps(report,indent=2)+'\n')
import shutil
shutil.copy2(__file__,out/'audit_reproduce.py')
print('REPORT',out/'AUDIT.json')
