from common import *
from graph_ops import topology,boundary,compactness,proposals,internal_components
import time
def main():
 start=time.perf_counter();ref=reference();n=int(ref['clot_count']);b=ref['bonds'];bt=ref['bond_types'];snapshots=[];candidates=[];screened=[];all_actual=[];best_rejected=[]
 for path,s in keyframes():
  step=int(s['step']);t=float(s['time']);labels,sizes,att=topology(ref,s['d']);free=sizes[~att];assert np.array_equal(~att[labels],s['detached']);assert sizes.sum()==n
  comps=[]
  for j in range(len(sizes)):
   nodes=np.flatnonzero(labels==j);out,inside=boundary(ref,nodes)
   comps.append(dict(label=j,status='ATTACHED' if att[j] else 'ACTUAL DETACHED',atom_ids=nodes+1,nodes=len(nodes),mass_kg=len(nodes)*RHO*V,volume_m3=len(nodes)*V,equivalent_diameter_m=(6*len(nodes)*V/np.pi)**(1/3),internal_live_bulk_bonds=int(np.sum((bt[inside]==1)&(s['d'][inside]<1))),external_live_bulk_bonds=int(np.sum((bt[out]==1)&(s['d'][out]<1))),external_live_interface_bonds=int(np.sum((bt[out]==2)&(s['d'][out]<1)))))
  all_actual.append(dict(step=step,time_s=t,components=comps));snapshots.append(dict(step=step,time_s=t,attached_nodes=int(sizes[att].sum()),detached_nodes=int(free.sum()),detached_components=len(free),largest_detached_nodes=int(free.max(initial=0)),N_ge200=int(np.sum(free>=34)),N_ge300=int(np.sum(free>=114)),total_material_mass_kg=n*RHO*V))
  proposed,thresholds=proposals(ref,s);screened.append(dict(step=step,time_s=t,threshold_screens=thresholds,raw_proposals=len(proposed)))
  accepted=0;rejected=[]
  for ids,methods in proposed.items():
   nodes=np.asarray(ids);out,inside=boundary(ref,nodes);inside=inside[bt[inside]==1];live=s['d']<1;external=out[live[out]];initial_external=len(out);internal_survival=float(np.mean(live[inside])) if len(inside) else 0.;internal_damage=float(np.mean(s['d'][inside])) if len(inside) else 1.;cut_loss=float(np.mean(s['d'][out])) if len(out) else 1.;qcut=float(np.sum(1-s['d'][out]));qinside=float(np.sum(1-s['d'][inside]));conductance=qcut/max(2*qinside+qcut,1e-30);shape=compactness(ref['X'][nodes],DX)
   nc=internal_components(ref,nodes,s['d']);connected=nc==1
   qualifies=connected and len(inside)>0 and internal_survival>=.8 and cut_loss>=.5 and cut_loss>internal_damage and conductance<=.15 and shape['voxel_bbox_fill']>=.2 and np.any(att[labels[nodes]])
   metrics=dict(nodes=len(nodes),internal_survival=internal_survival,mean_internal_d=internal_damage,boundary_integrity_loss=cut_loss,weighted_conductance=conductance,reference_fill=shape['voxel_bbox_fill'])
   if not qualifies:rejected.append(dict(**metrics,search_methods=methods));continue
   accepted+=1;stretch=np.linalg.norm(s['x'][b[out,1]]-s['x'][b[out,0]],axis=1)/ref['r0'][out]-1
   candidates.append(dict(candidate_id='step_'+str(step)+'_candidate_'+str(accepted),status='TETHERED_CANDIDATE',step=step,time_s=t,search_methods=methods,atom_ids=nodes+1,reference_volume_m3=len(nodes)*V,equivalent_diameter_m=(6*len(nodes)*V/np.pi)**(1/3),reference_geometry=shape,current_geometry=compactness(s['x'][nodes],DX),internal_reference_bulk_bonds=len(inside),internal_alive_bulk_bonds=int(live[inside].sum()),external_reference_bonds=initial_external,external_alive_bonds=len(external),external_alive_to_bulk=int(np.sum(bt[external]==1)),external_alive_to_anchors=int(np.sum(bt[external]==2)),external_bond_ids=out+1,external_bond_types=bt[out],external_kappa=s['kappa'][out],external_d=s['d'][out],external_actual_stretch=stretch,external_alive=live[out],free_surface_nodes=int(np.sum(ref['areas'][nodes]>0)),**metrics))
  best_rejected.append(dict(step=step,time_s=t,top_diagnostic_proposals=sorted(rejected,key=lambda x:(-x['boundary_integrity_loss'],x['weighted_conductance']))[:5]))
  print('KEY',step,'detached',int(free.sum()),'proposals',len(proposed),'qualified',accepted,flush=True)
 dump(O/'ACTUAL_KEYFRAME_TOPOLOGY.json',dict(status='PASS',contact_edges_excluded=True,anchors_excluded_from_mass=True,frames=all_actual,summary=snapshots))
 dump(O/'V3_CHUNK_CANDIDATES.json',dict(status='TETHERED_CANDIDATE' if candidates else 'NO_COHERENT_CHUNK_CANDIDATE',scope='Only diagnostic graph searches on original accepted keyframes; absence is not an exhaustive mathematical proof',criteria=dict(min_nodes=3,max_nodes=n//2,minimum_internal_survival=.8,minimum_boundary_integrity_loss=.5,max_weighted_conductance=.15,minimum_reference_bbox_fill=.2,weights='1-d is a graph heuristic, NOT certified force capacity or fracture energy'),candidate_count=len(candidates),candidate_ge200=sum(c['nodes']>=34 for c in candidates),candidate_ge300=sum(c['nodes']>=114 for c in candidates),candidates=candidates,threshold_screens=screened,best_rejected_diagnostic_proposals=best_rejected,wall_s=time.perf_counter()-start))
 print('DONE candidates',len(candidates),'ge200',sum(c['nodes']>=34 for c in candidates),'wall',time.perf_counter()-start,flush=True)
if __name__=='__main__':main()
