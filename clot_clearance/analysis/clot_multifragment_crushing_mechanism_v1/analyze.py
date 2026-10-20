from common import *
from indenter import *
import csv,gzip
from scipy.sparse import coo_matrix,csgraph
from scipy.spatial import cKDTree
graph=module('analysis_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
def read(p):return json.loads(Path(p).read_text())
def histories(case):
 rows=[]
 for r in csv.DictReader((case/'history.csv').open()):
  rows.append({k:(v if k=='phase' else v=='True' if k=='instant_target' else float(v) if v else None) for k,v in r.items()})
 return rows
def deformation(ref,s,ids):
 """Retained-neighbor local affine fit, diagnostic only; no stress calibration."""
 n=int(ref['clot_count']);b=ref['bonds'];keep=(ref['bond_types']==1)&~ref['initial_cut']&(s['d']<1);edge=b[keep];out=[]
 model=energy_model(ref);state=model.evaluate(s['x'],s['d']);dy=s['x'][b[:,1]]-s['x'][b[:,0]];length=np.linalg.norm(dy,axis=1);tensor=.5/V*np.einsum('bi,bj,b->bij',dy,dy,state['t']/length);virial=np.zeros((len(ref['X']),3,3));np.add.at(virial,b[:,0],tensor);np.add.at(virial,b[:,1],tensor)
 for i in sorted(set(ids)):
  if i>=n:continue
  nb=np.r_[edge[edge[:,0]==i,1],edge[edge[:,1]==i,0]];xi=ref['X'][nb]-ref['X'][i];y=s['x'][nb]-s['x'][i];rank=np.linalg.matrix_rank(xi) if len(xi) else 0
  if rank<3:out.append(dict(atom_id=int(i+1),retained_neighbors=len(nb),rank=rank,resolved=False));continue
  Ft=np.linalg.lstsq(xi,y,rcond=None)[0];E=(Ft@Ft.T-np.eye(3))/2;eig=np.linalg.eigvalsh(E);res=np.linalg.norm(xi@Ft-y)/max(np.linalg.norm(y),1e-30)
  sigma=virial[i];mean=np.trace(sigma)/3;dev=sigma-mean*np.eye(3)
  out.append(dict(atom_id=int(i+1),retained_neighbors=len(nb),rank=rank,resolved=True,green_principal_strains=eig,max_green_shear=float((eig[-1]-eig[0])/2),volume_ratio=float(np.linalg.det(Ft)),affine_fit_relative_residual=float(res),material_virial_per_reference_volume_Pa=sigma,material_virial_principal_Pa=np.linalg.eigvalsh(sigma),material_virial_mean_tension_positive_Pa=float(mean),material_virial_von_mises_Pa=float(np.sqrt(1.5*np.sum(dev*dev))),stress_scope='Nonlocal material-edge virial per original volume; contact excluded; not a calibrated local Cauchy stress'))
 return out
def crack_audit(case,ref):
 n=int(ref['clot_count']);b=ref['bonds'];bt=ref['bond_types'];parent=np.arange(n);active=np.zeros(n,bool);allfailed=set();cluster_ids={};births=[];merges=[];events=[];largest_clusters=0
 def root(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return int(i)
 def join(a,z):
  a=root(a);z=root(z)
  if a!=z:parent[max(a,z)]=min(a,z)
 with gzip.open(case/'damage_events.jsonl.gz','rt') as stream:
  for line in stream:
   z=json.loads(line);ids=np.asarray(z['original_bond_ids'],int)-1;d=np.asarray(z['d']);new=ids[(d==1)&(bt[ids]==1)]
   if not len(new):continue
   previous_roots={root(i):set(v) for i,v in cluster_ids.items()};newnodes=set(b[new].ravel().tolist());existing={root(i) for i in newnodes if active[i]}
   for a,c in b[new]:join(int(a),int(c));active[a]=active[c]=True
   combined={}
   for old,labels in previous_roots.items():combined.setdefault(root(old),set()).update(labels)
   born_roots={root(i) for i in newnodes if root(i) not in combined}
   for rt in sorted(born_roots):
    cid=len(births)+1;combined[rt]={cid};involved=[int(q+1) for q in new if root(int(b[q,0]))==rt];mid=ref['X'][b[np.array(involved)-1]].mean((0,1));births.append(dict(cluster_id=cid,time_s=z['time_s'],step=z['step'],first_bond_ids=involved,reference_centroid_m=mid))
   for rt,labels in combined.items():
    parents=[old for old in previous_roots if root(old)==rt]
    if len(parents)>1:merges.append(dict(time_s=z['time_s'],step=z['step'],n_merged_components=len(parents),original_nucleation_ids=sorted(labels),new_bridge_bond_ids=[int(q+1) for q in new if root(int(b[q,0]))==rt]))
   cluster_ids=combined;allfailed.update((new+1).tolist());largest_clusters=max(largest_clusters,len(combined));events.append(dict(time_s=z['time_s'],step=z['step'],separate_failed_bulk_networks=len(combined),new_full_failure_bonds=len(new)))
 return dict(definition='Naturally fully failed bulk bonds connected through shared original endpoints; precuts excluded. Discrete damage-network births/merges, not proof of continuous crack surfaces.',births=births,merges=merges,history=events,max_simultaneous_networks=largest_clusters,n_births=len(births),n_merge_events=len(merges),new_failed_bulk_bonds=len(allfailed))
def analyze_case(name,ref):
 case=O/'cases'/name;summary=read(case/'SUMMARY.json');rows=histories(case);registry={int(k):v for k,v in read(case/'COMPONENT_REGISTRY.json').items()};final=load_npz(case/'keys'/('final.npz' if summary['status']=='COMPLETED' else 'last_accepted.npz'));n=int(ref['clot_count']);lab,sizes,att=graph.topology(ref,final['d']);free=sizes[~att];snap={};first={};last={};trace={};midcohort=None
 with gzip.open(case/'CONNECTED_COMPONENT_HISTORY.jsonl.gz','rt') as stream:
  for line in stream:
   z=json.loads(line);t=z['time_s']
   if abs(t-CLEAR)<1e-13:midcohort=z
   for q in z['components']:
    cid=q['component_id'];first.setdefault(cid,dict(time_s=t,**q));last[cid]=dict(time_s=t,**q)
    if registry[cid]['nodes']>=2:trace.setdefault(cid,[]).append(dict(time_s=t,internal_survival=q['internal_survival'],mean_internal_damage=q['mean_internal_damage'],internal_relative_kinetic_J=q['internal_relative_kinetic_J'],centroid_velocity_m_s=q['centroid_velocity_m_s'],internal_initial_bonds=q['internal_initial_bonds'],internal_alive_bonds=q['internal_alive_bonds'],remaining_live_external_connections=q['remaining_live_external_connections']))
   snap['final']=z
 contact_rows=[q for q in rows if q['contact_nodes']>0];lastcontact=contact_rows[-1]['time_s'] if contact_rows else None;contacts_after=[q for q in rows if q['time_s']>=CLEAR-1e-14 and q['contact_nodes']>0]
 stability=[]
 for cid,c in registry.items():
  ids=np.array(c['atom_ids'])-1;connected=len(set(lab[ids].tolist()))==1 and not att[lab[ids[0]]];same=connected and sizes[lab[ids[0]]]==len(ids);lasttime=last.get(cid,{}).get('time_s',c['birth_time_s'])
  stability.append(dict(component_id=cid,nodes=c['nodes'],birth_s=c['birth_time_s'],last_observed_s=lasttime,observed_lifetime_s=lasttime-c['birth_time_s'],same_original_members_connected_at_final=bool(same),status='CONNECTED_THROUGH_OBSERVED_WINDOW' if same else 'SPLIT_WITHIN_OBSERVED_WINDOW',internal_survival_at_birth=first.get(cid,{}).get('internal_survival'),internal_survival_last=last.get(cid,{}).get('internal_survival'),post_withdrawal_observed_s=max(0,lasttime-max(c['birth_time_s'],CLEAR))))
 stableids={q['component_id'] for q in stability if q['same_original_members_connected_at_final']};cohort_ids=[q['component_id'] for q in midcohort['components'] if 20<=registry[q['component_id']]['nodes']<=200] if midcohort else [];clear_row=next((q for q in rows if abs(q['time_s']-CLEAR)<1e-13),None)
 target=bool(clear_row and clear_row['instant_target'] and sum(cid in stableids for cid in cohort_ids)>=5 and summary['status']=='COMPLETED')
 classification=[];end=rows[-1]
 if not end['softened']:classification.append('NO_DAMAGE')
 elif not end['detached_nodes']:classification.append('DAMAGE_WITHOUT_DETACHMENT')
 if end['largest_fraction']>.3:classification.append('SINGLE_LARGE_DETACHMENT')
 if end['N20_200']>=2:classification.append('MULTIPLE_COHESIVE_FRAGMENTS')
 if end['detached_nodes'] and end['debris12_detached_fraction']>.2:
  classification.append('DEBRIS_DOMINATED_DETACHMENT')
  classification.append('EXCESSIVE_PULVERIZATION' if end['detached_nodes']/n>=.1 else 'SURFACE_EROSION')
 elif end['detached_nodes'] and not classification:classification.append('SURFACE_EROSION')
 end['debris12_initial_fraction']=end['debris12_nodes']/n
 if summary['secondary_events']:classification.append('SECONDARY_FRAGMENTATION')
 if summary['status']!='COMPLETED':classification.append('NUMERICALLY_UNRESOLVED')
 finalcs=[]
 for q in snap['final']['components']:
  cid=q['component_id'];finalcs.append(dict(**registry[cid],**q))
 failed=ref['bonds'][(final['d']==1)&~ref['initial_cut']&(ref['bond_types']==1)];cracks=crack_audit(case,ref);diagnostics={}
 for ev in ['first_softening','first_failure']:
  if ev in summary['events']:
   e=summary['events'][ev];s=load_npz(case/e['evidence']);diagnostics[ev]=dict(event=e,local_affine=deformation(ref,s,np.asarray(e['atom_ids']).ravel()-1),note='Local retained-neighbor Green strain and nonaffinity; not calibrated Cauchy stress or a new compression/shear failure rule')
 result=dict(name=name,status=summary['status'],classification=classification,exploratory_target='EXPLORATORY_MULTIFRAGMENT_TARGET_MET' if target else 'EXPLORATORY_MULTIFRAGMENT_TARGET_NOT_MET',initial_material_nodes=n,first_softening_s=summary['events'].get('first_softening',{}).get('time_s'),first_failure_s=summary['events'].get('first_failure',{}).get('time_s'),first_multinode_s=summary['events'].get('first_multinode',{}).get('time_s'),history_max_N20_200=summary['max_medium_components'],history_max_free_nodes=summary['history_max_free_nodes'],final=end,final_components=finalcs,first_contact_s=contact_rows[0]['time_s'] if contact_rows else None,last_actual_contact_s=lastcontact,contact_after_scheduled_withdrawal=bool(contacts_after),actual_contact_free_tail_s=END-max(CLEAR,lastcontact) if lastcontact is not None and summary['status']=='COMPLETED' else None,peak_reaction_N=max(q['reaction_norm_N'] for q in rows),max_overlap_m=max(q['max_overlap_m'] for q in rows),max_overlap_dx=max(q['max_overlap_m'] for q in rows)/DX,max_relative_energy_residual=summary['max_relative_energy_residual'],crack_networks=dict(n_births=cracks['n_births'],n_merges=cracks['n_merge_events'],max_simultaneous=cracks['max_simultaneous_networks']),withdrawal_medium_cohort=cohort_ids,withdrawal_medium_cohort_surviving=sum(q in stableids for q in cohort_ids),secondary_after_withdrawal=sum(q['time_s']>=CLEAR for q in summary['secondary_events']),evidence_path=str(case.relative_to(R)),scope='Same prescribed kinematic trajectory, not equal work or validated physical hammer/bubble')
 result['history_max_size_threshold_counts']={key:int(max(row[key] for row in rows)) for key in ['N_ge34','N_ge114','N_ge200um','N_ge300um']}
 result['classification_policy']='Debris-dominated means >20% of detached mass is 1–2 nodes; call it excessive pulverization only if >=10% of initial material is detached, otherwise surface erosion. This disclosed descriptive convention does not change the frozen exploratory target.'
 dump(case/'CRACK_NETWORK_AUDIT.json',cracks);dump(case/'LOCAL_DEFORMATION_AUDIT.json',diagnostics);dump(case/'POST_UNLOAD_FRAGMENT_STABILITY.json',dict(fragments=stability,withdrawal_medium_cohort=cohort_ids,secondary_events=summary['secondary_events'],actual_contact_free_tail_s=result['actual_contact_free_tail_s'],contact_after_withdrawal=bool(contacts_after)));dump(case/'FRAGMENT_TRACE.json',trace);dump(case/'ANALYSIS.json',result)
 return result,summary,rows
def main():
 ref=load_npz(O/'provenance/reference.npz');cases={};summaries={};rows={}
 for name in ['uniform','heterogeneous']:cases[name],summaries[name],rows[name]=analyze_case(name,ref)
 dump(O/'UNIFORM_CRUSHING_RESULTS.json',cases['uniform']);dump(O/'HETEROGENEOUS_CRUSHING_RESULTS.json',cases['heterogeneous']);dump(O/'FRACTURE_EVENTS.json',{k:v['events'] for k,v in summaries.items()})
 dump(O/'CONNECTED_COMPONENT_HISTORY.json',{k:dict(history=f'cases/{k}/CONNECTED_COMPONENT_HISTORY.jsonl.gz',registry=f'cases/{k}/COMPONENT_REGISTRY.json',every_accepted_state=True) for k in cases});dump(O/'POST_UNLOAD_FRAGMENT_STABILITY.json',{k:read(O/'cases'/k/'POST_UNLOAD_FRAGMENT_STABILITY.json') for k in cases})
 dump(O/'DYNAMIC_CONTACT_AND_ENERGY_LEDGER.json',dict(definition='W_indenter = Delta U_LPS + Delta K + D + Delta U_pair - J_pair_activation + Delta U_indenter + R; numerical damping work = 0',external_work='Integral of partial U_indent / partial center dot center increment on accepted piecewise-linear paths; Gauss quadrature split at exact contact entry/exit. Not force times final displacement; not counted twice.',cases={k:dict(history=f'cases/{k}/history.csv',final=v['final'],peak_reaction_N=v['peak_reaction_N'],max_overlap_m=v['max_overlap_m']) for k,v in cases.items()}))
 distributions=[]
 for k,v in cases.items():
  for q in v['final_components']:distributions.append(dict(case=k,component_id=q['component_id'],nodes=q['nodes'],reference_volume_m3=q['reference_volume_m3'],equivalent_diameter_m=q['equivalent_diameter_m'],mass_kg=q['nodes']*RHO*V,internal_survival=q['internal_survival'],mean_internal_damage=q['mean_internal_damage'],centroid_speed_m_s=float(np.linalg.norm(q['centroid_velocity_m_s'])),internal_relative_kinetic_J=q['internal_relative_kinetic_J']))
 columns=['case','component_id','nodes','reference_volume_m3','equivalent_diameter_m','mass_kg','internal_survival','mean_internal_damage','centroid_speed_m_s','internal_relative_kinetic_J']
 with (O/'MULTIFRAGMENT_SIZE_DISTRIBUTION.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=columns);w.writeheader();w.writerows(distributions)
 resources=dict(cases={k:{q:v[q] for q in ['wall_s','cpu_s','peak_RSS_bytes','accepted_steps','predictive_substeps','predictive_rejections']} for k,v in summaries.items()},native_wall_s=sum(v['wall_s'] for v in summaries.values()),native_cpu_s=sum(v['cpu_s'] for v in summaries.values()),peak_RSS_bytes=max(v['peak_RSS_bytes'] for v in summaries.values()),case_count=len(summaries),rank_per_case=1,threads_per_case=1,build_wall_s=read(O/'BUILD_AUDIT.json')['wall_s'],key_check_wall_s=read(O/'KEY_CONSISTENCY_CHECKS.json')['wall_s'],new_dependencies=[],full_V3_runs=0,CFD_runs=0)
 oldpath=PREVIOUS/'cases/local_impingement_1000pa_r2/SUMMARY.json';old=read(oldpath)
 dump(O/'HISTORICAL_LOADING_CONTEXT.json',dict(source=str(oldpath),source_sha256=sha(oldpath),old_protocol='1000 Pa fixed reference inward dead traction, 2.5 ms; original loaded nodes remain loaded after detachment',old_W_J=old['last_row']['W_J'],old_K_J=old['last_row']['K_J'],old_D_J=old['last_row']['D_J'],old_history_largest=old['history_max_free_nodes'],current={k:dict(W_J=v['final']['W_J'],K_J=v['final']['K_J'],D_J=v['final']['D_J'],W_ratio_to_old=v['final']['W_J']/old['last_row']['W_J']) for k,v in cases.items()},comparison='Context only; different duration, force/area path, work and contact rule, not a single-factor test of H1'))
 resources['measurement_scope']='Each formal driver wall and process CPU, including native stepping, online consistency/energy audits and evidence output; build and short checks listed separately. Postprocessing is not included in native totals.'
 timestep={}
 for k,hist in rows.items():
  tiny=[]
  for i,q in enumerate(hist):
   if q['dt_s']<MIN_H*(1-1e-10):
    previous_time=hist[i-1]['time_s'] if i else 0.
    tiny.append(dict(step=int(q['step']),dt_s=q['dt_s'],time_s=q['time_s'],previous_time_s=previous_time,max_dd=q['max_dd'],at_prescribed_boundary=any(abs(q['time_s']-t)<1e-15 for t in [APPROACH,PEAK,RETRACT,CLEAR,END]),below_event_floor_reason='Prescribed-boundary residual clipping; MIN_H limits predictive damage subdivision, not endpoint alignment'))
  timestep[k]=dict(accepted_steps=len(hist),base_dt_s=BASE,predictive_subdivision_floor_s=MIN_H,min_actual_accepted_dt_s=min(q['dt_s'] for q in hist),min_non_roundoff_dt_s=min(q['dt_s'] for q in hist if q['dt_s']>1e-12),boundary_residual_steps=tiny)
 dump(O/'TIMESTEP_BOUNDARY_AUDIT.json',dict(cases=timestep,interpretation='Actual accepted histories are preserved. Endpoint floating-point residual steps are explicitly distinguished from damage-driven substeps; no rerun or rewritten clock.'))
 dump(O/'RUN_RESOURCES.json',resources);print('ANALYZED',[(k,v['classification'],v['final']['N20_200']) for k,v in cases.items()]);return cases
if __name__=='__main__':main()
