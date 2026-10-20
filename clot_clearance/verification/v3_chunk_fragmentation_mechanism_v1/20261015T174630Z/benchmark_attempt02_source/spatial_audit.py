from common import *
from graph_ops import topology,boundary,internal_components
from scipy.spatial import cKDTree
import gzip,time,csv

def main():
 start=time.perf_counter();ref=reference();X=ref['X'];b=ref['bonds'];bt=ref['bond_types'];n=int(ref['clot_count']);nb=len(b)
 hist=np.genfromtxt(B/'V3_NATIVE_PILOT_HISTORY.csv',delimiter=',',names=True)
 soft=np.full(nb,np.inf);fail=np.full(nb,np.inf);failstep=np.full(nb,2**30,int);events=[]
 for line in (B/'full/events.jsonl').open():
  e=json.loads(line)
  if e['kind']!='ACCEPTED':continue
  events.append(e);si=np.asarray(e['new_soft_ids'],int)-1;fi=np.asarray(e['new_failed_ids'],int)-1
  soft[si]=e['time_s'];fail[fi]=e['time_s'];failstep[fi]=e['step']
 assert len(events)==len(hist)
 dtype=np.dtype([('step','<i4'),('bond_id','<i4'),('kappa','<f8'),('d','<f8')])
 with gzip.open(B/'full/damage_updates.bin.gz','rb') as f:updates=np.frombuffer(f.read(),dtype=dtype)
 def damage(step):
  a=updates[updates['step']<=step];d=np.zeros(nb);d[a['bond_id']-1]=a['d'];return d
 keys=list(keyframes())
 for p,s in keys:assert np.array_equal(damage(int(s['step'])),s['d']),p
 maxstep=int(hist['step'][np.argmax(hist['largest_detached_nodes'])]);dmax=damage(maxstep);lab,sizes,att=topology(ref,dmax);ci=np.flatnonzero(~att)[np.argmax(sizes[~att])];nodes=np.flatnonzero(lab==ci);out,inside=boundary(ref,nodes);inside=inside[bt[inside]==1]
 track=[]
 for step in sorted(set([int(s['step']) for _,s in keys]+[maxstep]+failstep[inside].tolist())):
  if step>10502:continue
  d=damage(step);ll,ss,aa=topology(ref,d);u,c=np.unique(ll[nodes],return_counts=True)
  track.append(dict(step=step,time_s=float(hist['time_s'][step]),internal_live=int(np.sum(d[inside]<1)),internal_mean_d=float(np.mean(d[inside])),pieces=c,actual_component_sizes=ss[u],piece_attached=aa[u],external_live=int(np.sum(d[out]<1))))
 maxgroup=dict(status='ACTUAL DETACHED',step=maxstep,time_s=float(hist['time_s'][maxstep]),atom_ids=nodes+1,nodes=len(nodes),reference_internal_bonds=len(inside),live_internal_bonds_at_formation=int(np.sum(dmax[inside]<1)),mean_internal_d_at_formation=float(np.mean(dmax[inside])),external_live_at_formation=int(np.sum(dmax[out]<1)),internal_bond_ids=inside+1,internal_fail_time_s=[None if not np.isfinite(v) else float(v) for v in fail[inside]],trajectory=track)
 # Exact topology reconstructed at every accepted full-failure change; reuse otherwise.
 rows=[];live_d=np.zeros(nb);ll,ss,aa=topology(ref,live_d)
 for e in events:
  if e['new_failed_ids']:
   live_d[np.asarray(e['new_failed_ids'])-1]=1;ll,ss,aa=topology(ref,live_d)
  st=e['step'];free=ss[~aa]
  assert free.sum()==hist['detached_nodes'][st] and free.max(initial=0)==hist['largest_detached_nodes'][st]
  rows.append([st,e['time_s'],int(free.sum()),len(free),int(free.max(initial=0))])
 np.savetxt(O/'derived/all_accepted_topology.csv',rows,delimiter=',',header='step,time_s,detached_nodes,detached_components,largest_detached_nodes',comments='')
 cand=json.loads((O/'V3_CHUNK_CANDIDATES.json').read_text());cand['large_chunk_status']='TETHERED_CANDIDATE' if cand['candidate_ge200'] else 'NO_COHERENT_CHUNK_CANDIDATE'
 for c in cand['candidates']:
  cn=np.asarray(c['atom_ids'])-1;co,ci=boundary(ref,cn);c['evolution']=[]
  for _,s in keys:
   ll,ss,aa=topology(ref,s['d']);ext=co[s['d'][co]<1]
   c['evolution'].append(dict(step=int(s['step']),time_s=float(s['time']),internal_components=internal_components(ref,cn,s['d']),internal_alive=int(np.sum(s['d'][ci]<1)),external_alive_bond_ids=ext+1,attached_atom_count=int(np.sum(aa[ll[cn]]))))
  c['remaining_bridge_failure_times']=[dict(bond_id=int(q+1),time_s=None if not np.isfinite(fail[q]) else float(fail[q])) for q in co[np.array(c['external_alive'])]]
 dump(O/'V3_CHUNK_CANDIDATES.json',cand)
 # Original free surface, not the subset of load-bearing points.
 import pyvista as pv
 surface=R/'results/y_clot_realistic_morphology_wall_contact_and_flow_transport_v1/paraview/exposed_surface.vtp'
 mesh=pv.read(surface);surf_distance=np.abs(pv.PolyData(X).compute_implicit_distance(mesh)['implicit_distance'])
 offsets=np.stack(np.meshgrid(*[np.arange(-3,4)*DX]*3,indexing='ij'),-1).reshape(-1,3);rr=np.linalg.norm(offsets,axis=1);rr=rr[(rr>0)&(rr<=HORIZON)]
 hh=np.minimum(1,(HORIZON+DX/2-rr)/DX);ideal_m=float(np.sum(hh*rr*V));degree=np.bincount(b.ravel(),minlength=len(X));bulkdeg=np.bincount(b[bt==1].ravel(),minlength=len(X));interdeg=np.bincount(b[bt==2].ravel(),minlength=len(X))
 load0=ref['areas']>0;early=(bt==1)&(fail<=.0003225);earlynode=np.unique(b[early]);ever=np.unique(b[(bt==1)&np.isfinite(fail)]);mid=.5*(X[b[:,0]]+X[b[:,1]])
 def stats(ids):
  return dict(nodes=len(ids),neighbor_count_mean=float(degree[ids].mean()),missing_fraction_mean=float((1-degree[ids]/len(rr)).mean()),m_ratio_mean=float((ref['m'][ids]/ideal_m).mean()),surface_distance_um_quantiles=np.quantile(surf_distance[ids]*1e6,[0,.25,.5,.75,1]),within_one_horizon_fraction=float(np.mean(surf_distance[ids]<=HORIZON)),load_area_node_fraction=float(load0[ids].mean()),interface_neighbor_fraction=float(np.mean(interdeg[ids]>0)))
 bins=[]
 for lo,hi in [(0,.5),(.5,1),(1,2),(2,100)]:
  ids=np.flatnonzero((surf_distance[:n]/HORIZON>=lo)&(surf_distance[:n]/HORIZON<hi));bins.append(dict(distance_delta_bin=[lo,hi],nodes=len(ids),early_failure_node_fraction=float(np.isin(ids,earlynode).mean()) if len(ids) else None,load_area_node_fraction=float(load0[ids].mean()) if len(ids) else None))
 np.savez_compressed(O/'derived/spatial_fields.npz',soft_time_s=soft,fail_time_s=fail,fail_step=failstep,free_surface_distance_m=surf_distance,degree=degree,bulk_degree=bulkdeg,interface_degree=interdeg,m_ratio=ref['m']/ideal_m,early_failure_nodes=earlynode,max_fragment_nodes=nodes)
 surfaceaudit=dict(source=str(surface.relative_to(R)),sha256=sha(surface),method='absolute closest triangle distance to archived exposed_surface.vtp; not load-node distance',ideal_neighbors=len(rr),ideal_m=ideal_m,all_material=stats(np.arange(n)),early_bulk_failed_nodes=stats(earlynode),ever_bulk_failed_nodes=stats(ever),distance_bins=bins,interpretation='Association only: free surface, incomplete horizon, applied surface load and geometry are confounded; no horizon correction justified by this comparison alone.')
 first=[]
 for name,times,typ in [('bulk_soft',soft,1),('bulk_failure',fail,1),('interface_soft',soft,2),('interface_failure',fail,2)]:
  ids=np.flatnonzero(bt==typ);q=ids[np.argmin(times[ids])];first.append(dict(event=name,time_s=times[q],bond_id=int(q+1),atom_ids=b[q]+1,reference_midpoint_m=mid[q],endpoint_free_surface_distance_m=surf_distance[b[q]],endpoint_original_neighbors=degree[b[q]],endpoint_m_ratio=ref['m'][b[q]]/ideal_m))
 dump(O/'COHESIVE_ADHESIVE_COMPETITION.json',dict(status='PASS',first_events=first,first_detachment_time_s=float(hist['time_s'][np.flatnonzero(hist['detached_nodes']>0)[0]]),surface_horizon_audit=surfaceaudit,historical_maximum_fragment=maxgroup,all_accepted_topology_checked=len(events),final_bulk_failure_fraction=float(np.mean(np.isfinite(fail[bt==1]))),final_interface_failure_fraction=float(np.mean(np.isfinite(fail[bt==2]))),interpretation='Bulk softening/failure precede interface failure and detachment. Earliest bulk failure is near exposed loaded surface; spatial evidence does not establish a resolved interior crack front. Interface network remains largely intact. The 13-node free component retains some internal bonds but subsequently splits; see exact identities/times.',wall_s=time.perf_counter()-start))
 # Read accepted states recovered losslessly, never trial states.
 states=sorted([load_npz(p) for p in (O/'derived/recovered_accepted').glob('*.npz') if int(load_npz(p)['step'])<2000],key=lambda s:int(s['step']))
 pair=np.array([573,1428]);near=np.flatnonzero((np.linalg.norm(X-midpoint,axis=1)<=HORIZON).any(0)) if False else np.flatnonzero(np.min(np.linalg.norm(X[:,None,:]-X[pair][None,:,:],axis=2),axis=1)<=HORIZON)
 localbulk=np.flatnonzero((bt==1)&(np.isin(b[:,0],near)|np.isin(b[:,1],near)));contact_c=15*(18*K/(np.pi*HORIZON**4))*V**2/HORIZON;cr=[]
 for s in states:
  x=s['x'];vec=x[pair[1]]-x[pair[0]];dist=np.linalg.norm(vec);force=contact_c*(dist-DX)*vec/dist;tree=cKDTree(x[:n]);local=[]
  for q in pair:
   ids=tree.query_ball_point(x[q],HORIZON);rids=cKDTree(X[:n]).query_ball_point(X[q],HORIZON)
   local.append(dict(atom_id=int(q+1),reference_position_m=X[q],displacement_m=x[q]-X[q],velocity_m_s=s['velocity'][q],total_force_N=s['force'][q],material_force_N=s['material_force'][q],load_N=s['load'][q],current_material_points_within_delta=len(ids),reference_material_points_within_delta=len(rids),local_count_ratio=len(ids)/len(rids),attached=not bool(s['detached'][q])))
  cr.append(dict(step=int(s['step']),time_s=float(s['time']),distance_m=dist,overlap_dx=(DX-dist)/DX,pair_force_on_first_N=force,pair_force_magnitude_N=np.linalg.norm(force),pair_contact_energy_J=.5*contact_c*(DX-dist)**2,local_nodes=local,neighbor_bulk_failed=int(np.sum(s['d'][localbulk]==1)),neighbor_bulk_mean_d=float(np.mean(s['d'][localbulk]))))
 windows=[]
 for lo,hi in [(0,.2305e-3),(.2305e-3,.37575e-3),(.37575e-3,.37625e-3),(.37625e-3,.45e-3),(.45e-3,2.331e-3)]:windows.append(dict(interval_s=[lo,hi],local_bulk_new_failures=int(np.sum((fail[localbulk]>=lo)&(fail[localbulk]<hi)))))
 dump(O/'CONTACT_0376MS_AUDIT.json',dict(status='PASS_ACCEPTED_STATES_RECOVERED',atom_ids=pair+1,original_direct_bond_exists=bool(np.any(np.all(b==pair,axis=1))),contact_stiffness_N_m=contact_c,states=cr,reference_neighborhood_atom_ids=near+1,local_bulk_failure_windows=windows,interpretation='Exact local contact force and accepted node data. Local damage predates the peak; time ordering alone does not identify contact as the primary fragmentation cause. No contact stiffness modified.',local_density_limit='Counts within a current Euclidean delta-ball; diagnostic compression proxy, not a continuum density/Jacobian measurement.'))
 print('DONE spatial audit',time.perf_counter()-start,flush=True)
if __name__=='__main__':main()
