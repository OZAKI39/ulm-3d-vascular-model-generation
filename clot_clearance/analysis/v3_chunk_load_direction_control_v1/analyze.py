"""Accepted-state analysis only; immutable original identities and graph criteria."""
from common import *
from graph_ops import topology,boundary,compactness,internal_components
import gzip,time,csv
DTYPE=np.dtype([('step','<i4'),('bond_id','<i4'),('kappa','<f8'),('d','<f8')])
class Archive:
 def __init__(self,root):
  self.root=root;self.ref=load_npz(root/'provenance/reference.npz');self.hist=np.genfromtxt(root/'V3_NATIVE_PILOT_HISTORY.csv',delimiter=',',names=True)
  with gzip.open(root/'full/damage_updates.bin.gz','rb') as f:self.updates=np.frombuffer(f.read(),DTYPE)
  self.keys=[load_npz(p) for p in sorted((root/'full/keyframes').glob('step_*.npz'))];self.index=np.genfromtxt(root/'full/native_frames.csv',delimiter=',',names=True)
  self.frames=np.memmap(root/'full/native_frames.bin',dtype='<f8',mode='r',shape=(len(self.index),4,len(self.ref['X']),3))
  self.events=[json.loads(l) for l in (root/'full/events.jsonl').open() if '"kind":"ACCEPTED"' in l]
  assert len(self.events)==len(self.hist)
  for s in self.keys:assert np.array_equal(self.damage(int(s['step'])),s['d']),('Accepted damage stream/keyframe mismatch',root,int(s['step']))
 def damage(self,step):
  a=self.updates[self.updates['step']<=step];d=np.zeros(len(self.ref['bonds']));d[a['bond_id']-1]=a['d'];return d
 def state(self,t):
  ki=int(np.argmin([abs(float(s['time'])-t) for s in self.keys]));fi=int(np.argmin(abs(self.index['time_s']-t)));s=self.keys[ki]
  if abs(float(s['time'])-t)<=abs(self.index['time_s'][fi]-t):return s
  raw=self.frames[fi];step=int(self.index['step'][fi]);return dict(x=raw[0],velocity=raw[1],force=raw[2],material_force=raw[3],d=self.damage(step),step=step,time=float(self.index['time_s'][fi]))
def group_metrics(ref,d,nodes,s=None):
 b=ref['bonds'];bt=ref['bond_types'];outside,inside=boundary(ref,nodes);inside=inside[bt[inside]==1];alive=d<1
 z=dict(atom_ids=nodes+1,nodes=len(nodes),reference_volume_m3=len(nodes)*V,mass_kg=len(nodes)*MASS,equivalent_diameter_m=(6*len(nodes)*V/np.pi)**(1/3),internal_reference_bulk_bonds=len(inside),internal_alive_bulk_bonds=int(alive[inside].sum()),internal_survival=float(alive[inside].mean()) if len(inside) else None,mean_internal_d=float(d[inside].mean()) if len(inside) else None,external_alive_bond_ids=outside[alive[outside]]+1,external_reference_bond_ids=outside+1,external_alive_to_bulk=int(np.sum(alive[outside]&(bt[outside]==1))),external_alive_to_anchors=int(np.sum(alive[outside]&(bt[outside]==2))),reference_shape=compactness(ref['X'][nodes],DX))
 if s is not None:
  vel=s['velocity'][nodes];com=vel.mean(0);z.update(current_shape=compactness(s['x'][nodes],DX),COM_position_m=s['x'][nodes].mean(0),COM_velocity_m_s=com,K_internal_J=float(.5*MASS*np.sum((vel-com)**2)),K_translation_J=float(.5*len(nodes)*MASS*np.dot(com,com)),nodal_state_actual_time_s=float(s['time']))
 return z
def topology_history(ar,label):
 ref=ar.ref;bt=ref['bond_types'];b=ref['bonds'];nb=len(b);d=np.zeros(nb);soft=np.full(nb,np.inf);fail=np.full(nb,np.inf);rows=[];groups={};active=set();at=0;u=ar.updates;lab,sz,att=topology(ref,d);observations={}
 stream=ar.root/'full/component_observations.jsonl'
 if stream.exists():
  for line in stream.open():
   q=json.loads(line)
   for g in q['components']:observations[(q['step'],tuple(g['atom_ids']))]=g
 for idx,e in enumerate(ar.events):
  step=e['step'];t=e['time_s'];end=np.searchsorted(u['step'],step,side='right');a=u[at:end];d[a['bond_id']-1]=a['d'];at=end
  soft[np.array(e['new_soft_ids'],int)-1]=t;fail[np.array(e['new_failed_ids'],int)-1]=t
  if e['new_failed_ids'] or idx==0:lab,sz,att=topology(ref,d)
  free=sz[~att];row=dict(step=step,time_s=t,detached_nodes=int(free.sum()),detached_components=len(free),largest_detached_nodes=int(free.max(initial=0)),tiny_nodes=int(free[free<=2].sum()),tiny_components=int(np.sum(free<=2)),detached_mass_kg=float(free.sum()*MASS),N_ge14=int(np.sum(free>=14)),N_ge34=int(np.sum(free>=34)),N_ge114=int(np.sum(free>=114)))
  for k in ['step','detached_nodes','detached_components','largest_detached_nodes','tiny_nodes']:assert ar.hist[k][idx]==row[k],(label,k,step)
  rows.append(row)
  if e['new_failed_ids']:
   present={tuple(np.flatnonzero(lab==j)+1) for j in np.flatnonzero(~att)}
   for ids in present-active:
    nodes=np.array(ids)-1;g=group_metrics(ref,d,nodes);g.update(first_separation_step=step,first_separation_time_s=t,first_split_time_s=None)
    obs=observations.get((step,ids));g['formation_motion']=obs if obs else {'status':'NOT_MEASURED_AT_EXACT_EVENT','reason':'Historical full nodal stream is sampled; topology event remains exact.'};groups[ids]=g
   for ids in active-present:
    groups[ids]['first_split_time_s']=t;groups[ids]['first_split_step']=step;groups[ids]['first_daughter_sizes']=np.unique(lab[np.array(ids)-1],return_counts=True)[1]
   active=present
 final=ar.keys[-1];assert int(final['step'])==int(ar.hist['step'][-1]);assert np.array_equal(d,final['d'])
 for ids,g in groups.items():
  nodes=np.array(ids)-1;g['final_same_identity_set']=group_metrics(ref,d,nodes,final);g['final_daughter_sizes']=np.unique(lab[nodes],return_counts=True)[1];g['retains_single_component_until_end']=g['first_split_time_s'] is None
 out=O/'derived';out.mkdir(exist_ok=True)
 with (out/(label+'_all_accepted_topology.csv')).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 np.savez_compressed(out/(label+'_event_fields.npz'),soft_time_s=soft,fail_time_s=fail,final_damage=d)
 largest=max(groups.values(),key=lambda x:(x['nodes'],-x['first_separation_time_s']),default=None)
 maxima={k:max(v[k] for v in rows) for k in ['detached_nodes','detached_components','largest_detached_nodes','tiny_nodes','N_ge14','N_ge34','N_ge114']}
 summary=dict(label=label,all_accepted_states_checked=len(rows),final=rows[-1],history_maxima=maxima,largest_historical_group=largest,unique_free_component_identity_sets=len(groups),all_free_components=list(groups.values()),attached_final_components=[group_metrics(ref,d,np.flatnonzero(lab==j),final) for j in np.flatnonzero(att)],threshold_note='14/34/114 nodes are discrete volume thresholds, not certified resolution or clinical size limits; 34/114 correspond to >=200/300 micrometres.',contact_is_not_material_connectivity=True)
 return summary,soft,fail
def failure_space(ref,ar,soft,fail,label):
 import pyvista as pv
 b=ref['bonds'];bt=ref['bond_types'];X=ref['X'];n=int(ref['clot_count']);mid=X[b].mean(1);top=pv.read(SURFACE/'exposed_surface.vtp');attach=pv.read(SURFACE/'attachment_surface.vtp')
 ds=np.abs(pv.PolyData(mid).compute_implicit_distance(top)['implicit_distance']);da=np.abs(pv.PolyData(mid).compute_implicit_distance(attach)['implicit_distance'])
 # Fixed geometric categories shared across the two histories, independent of outcomes.
 # Surface cells use proximity to the ORIGINAL loaded nodes, never current damage.
 from scipy.spatial import cKDTree
 distload=cKDTree(X[ref['areas']>0]).query(mid)[0];regions=dict(loaded_free_surface=(bt==1)&(ds<=HORIZON)&(distload<=HORIZON),other_free_surface=(bt==1)&(ds<=HORIZON)&(distload>HORIZON),interior=(bt==1)&(ds>HORIZON),interface=bt==2)
 req=json.loads((ar.root/'provenance/request.json').read_text());radial=np.linalg.norm(mid-np.asarray(req['load_centers_m'][0]),axis=1);surface=(bt==1)&(ds<=HORIZON)
 radial_regions=dict(surface_within_one_gaussian_sigma=surface&(radial<=req['gaussian_sigma_m']),surface_outside_one_gaussian_sigma=surface&(radial>req['gaussian_sigma_m']),interior=(bt==1)&(ds>HORIZON),interface=bt==2)
 counts=[];radial_counts=[]
 for t in [0.00015,.0003,.0005,.0008,.00125,float(ar.hist['time_s'][-1])]:
  if t>ar.hist['time_s'][-1]+1e-15:continue
  counts.append(dict(time_s=t,regions={k:dict(reference_bonds=int(mask.sum()),softened=int(np.sum(mask&(soft<=t+1e-15))),failed=int(np.sum(mask&(fail<=t+1e-15)))) for k,mask in regions.items()}))
  radial_counts.append(dict(time_s=t,regions={k:dict(reference_bonds=int(mask.sum()),softened=int(np.sum(mask&(soft<=t+1e-15))),failed=int(np.sum(mask&(fail<=t+1e-15)))) for k,mask in radial_regions.items()}))
 first=json.loads((ar.root/'FIRST_EVENTS.json').read_text());positions=[]
 for kind,q in first.items():
  if q is None or 'bond_ids' not in q:continue
  ids=np.array(q['bond_ids'])-1;positions.append(dict(event=kind,time_s=q['time_s'],bond_ids=ids+1,atom_ids=b[ids]+1,midpoints_m=mid[ids],distance_exposed_m=ds[ids],distance_attachment_m=da[ids],region={k:mask[ids] for k,mask in regions.items()}))
 np.savez_compressed(O/'derived'/f'{label}_spatial_fields.npz',midpoints_m=mid,distance_exposed_m=ds,distance_attachment_m=da,**regions)
 # Damage clusters are graph-connected damaged original bulk endpoints, not crack surfaces.
 from scipy.sparse import coo_matrix,csgraph
 damaged=(bt==1)&np.isfinite(fail);e=b[damaged];nc,ll=csgraph.connected_components(coo_matrix((np.ones(2*len(e)),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(n,n)),directed=False);used=np.unique(e);u,cc=np.unique(ll[used],return_counts=True)
 return dict(first_events=first,first_event_geometry=positions,region_histories=counts,region_definition='Mutually exclusive bulk categories: reference midpoint within delta of exposed mesh and within delta of an original loaded point; other exposed within delta; interior farther than delta. Interface is typed separately. Fixed common reference geometry.',gaussian_radial_region_histories=radial_counts,gaussian_radial_definition='Secondary spatial diagnostic only: exposed-layer bond reference midpoint lies within or outside one ORIGINAL Gaussian sigma of the ORIGINAL center. This compares the central loading patch with peripheral/side surface material. Original loaded IDs cover the exposed mesh, so the primary other_free_surface category is empty. This spatial partition does not change loading, fracture or candidate thresholds.',failed_bulk_endpoint_cluster_sizes=sorted(cc.tolist(),reverse=True),interpretation='Clusters of failed original-bond endpoints are a damage-organization diagnostic, not reconstructed continuous crack surfaces.')
def main():
 start=time.perf_counter();status=json.loads((O/'FINAL_STATUS.json').read_text());assert status['full_protocol_completed'] or status['last_accepted_time_s'] is not None
 old=Archive(B);new=Archive(O);ref=new.ref
 # Execute the exact previous screening implementation with ONLY archive bindings changed.
 audit=module_from_file('unchanged_candidate_screen',PREV_A/'topology_audit.py');audit.reference=lambda:ref;audit.keyframes=lambda:((None,s) for s in new.keys);audit.O=O;audit.main()
 candidates=json.loads((O/'V3_CHUNK_CANDIDATES.json').read_text());assert candidates['criteria']==json.loads((PREV/'V3_CHUNK_CANDIDATES.json').read_text())['criteria']
 alltop={};spaces={};energy={}
 for label,ar in [('historical_tangent',old),('reference_normal',new)]:
  topo,soft,fail=topology_history(ar,label);alltop[label]=topo;spaces[label]=failure_space(ref,ar,soft,fail,label);h=ar.hist
  fields=['W_load_J','W_drag_J','W_numeric_J','W_wall_elastic_J','W_wall_damping_J','W_J','U_J','K_J','D_J','Uc_J','Uw_J','Jc_J','R_J','relative_R','energy_scale_J','drag_background_supply_J','drag_relative_loss_J']
  maxima=['max_speed_m_s','max_dd','max_contact_penetration_dx','max_wall_overlap_m','relative_R','load_sum_norm_N','load_max_N']
  energy[label]=dict(final_time_s=float(h['time_s'][-1]),final={k:float(h[k][-1]) for k in fields},maxima={k:dict(value=float(h[k].max()),time_s=float(h['time_s'][np.argmax(h[k])])) for k in maxima})
  if label=='reference_normal':
   for c in candidates['candidates']:
    nodes=np.array(c['atom_ids'])-1;outside,inside=boundary(ref,nodes);c['remaining_bridge_failure_times']=[dict(bond_id=int(q+1),time_s=None if not np.isfinite(fail[q]) else float(fail[q])) for q in outside[np.asarray(c['external_alive'])]];c['evolution']=[]
    for s in new.keys:
     lab,sz,att=topology(ref,s['d']);c['evolution'].append(dict(step=int(s['step']),time_s=float(s['time']),internal_components=internal_components(ref,nodes,s['d']),attached_nodes=int(att[lab[nodes]].sum()),**group_metrics(ref,s['d'],nodes,s)))
 candidates['unique_identity_sets']={str(k):len({tuple(c['atom_ids']) for c in candidates['candidates'] if c['nodes']>=k}) for k in [3,14,34,114]}
 candidates['count_semantics']='candidate_count and candidate_ge200/ge300 count qualifying saved-state instances; unique_identity_sets deduplicates original IDs. Unequal keyframe sampling is not an exhaustive count of all possible regions.'
 dump(O/'V3_CHUNK_CANDIDATES.json',candidates);dump(O/'V3_DIRECTION_CONTROL_TOPOLOGY.json',alltop);dump(O/'BULK_INTERFACE_FAILURE_AUDIT.json',spaces)
 dump(O/'CHUNK_COHESION_AUDIT.json',dict(normal_actual=alltop['reference_normal']['all_free_components'],historical_largest=alltop['historical_tangent']['largest_historical_group'],candidates=candidates,formation_motion_note='New-case motion is stored at exact accepted topology changes. Historical formation motion is NOT MEASURED unless it coincides with an original full nodal record; no interpolated event velocities.'))
 ratio=energy['reference_normal']['final']['W_load_J']/energy['historical_tangent']['final']['W_load_J']
 dump(O/'DYNAMIC_ENERGY_LEDGER.json',dict(formula='W_load+W_drag+W_numeric+W_wall_damping=deltaU+deltaK+D+deltaUc+deltaUw-Jc+R',runs=energy,normal_over_tangent_load_work=ratio,ENERGY_COMPARABILITY_STATUS='SAME_SCALAR_RULE_NOT_EQUAL_INPUT_WORK',same_input_energy_efficiency='NOT_TESTED',no_energy_matching_second_case=True))
 matches=[]
 targets=[0.,.0005,.0008,.00125,min(float(new.hist['time_s'][-1]),END)]
 firsttimes=[v['time_s'] for v in [spaces['historical_tangent']['first_events']['bulk_failure'],spaces['reference_normal']['first_events']['bulk_failure']] if v is not None]
 if firsttimes:targets.append(min(firsttimes))
 for t in sorted(set(targets)):
  if t>new.hist['time_s'][-1]+1e-15:continue
  row=dict(requested_time_s=t)
  for label,ar in [('historical_tangent',old),('reference_normal',new)]:
   s=ar.state(t);p=O/'derived'/f'{label}_t{t*1e6:010.3f}us.npz';save_state(p,{k:s[k] for k in ['x','velocity','force','material_force','d','step','time']});row[label]=dict(step=int(s['step']),actual_accepted_time_s=float(s['time']),offset_s=float(s['time'])-t,state_path=str(p.relative_to(O)))
  matches.append(row)
 dump(O/'VISUAL_STATE_MATCHING.json',dict(method='Nearest actual accepted full nodal/keyframe state; no interpolation. Damage exactly reconstructed from original accepted updates when using a regular node frame.',matches=matches))
 dump(O/'HISTORICAL_DIRECTION_COMPARISON.json',dict(physical_cases_new=1,old_case_repeated=False,normal_final=alltop['reference_normal']['final'],historical_final=alltop['historical_tangent']['final'],normal_maxima=alltop['reference_normal']['history_maxima'],historical_maxima=alltop['historical_tangent']['history_maxima'],normal_candidate_ge14=sum(c['nodes']>=14 for c in candidates['candidates']),normal_candidate_ge34=candidates['candidate_ge200'],normal_candidate_ge114=candidates['candidate_ge300'],normal_over_tangent_W_load=ratio,first_events={k:v['first_events'] for k,v in spaces.items()},equal_energy_inference='NOT_SUPPORTED_BY_THIS_SINGLE_DIRECTION_CONTROL'))
 dump(O/'ANALYSIS_STATUS.json',dict(status='PASS',wall_s=time.perf_counter()-start,candidate_classifier_sha256=sha(PREV_A/'topology_audit.py'),candidate_graph_sha256=sha(PREV_A/'graph_ops.py'),actual_states_only=True));print('ANALYSIS_DONE',time.perf_counter()-start,flush=True)
if __name__=='__main__':main()
