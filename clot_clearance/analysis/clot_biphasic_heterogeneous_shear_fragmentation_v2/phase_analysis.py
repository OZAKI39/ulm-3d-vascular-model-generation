"""Original-ID phase evidence; failed-edge regions never substitute for free fragments."""
from common import *
import gzip,csv
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
EARLY_WINDOW=.0001
SITE_SEPARATION=2*HORIZON
def analyze_spatial(ref,results,topology):
 field=load_npz(O/'BIPHASIC_FIELD.npz');b=ref['bonds'];nb=len(b);n=int(ref['clot_count']);cut=ref['initial_cut'];bulk=(ref['bond_types']==1)&~cut;phase=field['material_phase'];domain=field['strong_domain_id'];mid=ref['X'][b].mean(1);out={};eventrows=[];domainrows=[]
 for mode in CASES:
  case=O/'cases'/mode;summary=results[mode]['native'];final=load_npz(case/'keys'/('final.npz' if summary['status']=='COMPLETED' else 'last_accepted.npz'));d=final['d'];soft=np.full(nb,np.nan);fail=soft.copy()
  with gzip.open(case/'damage_events.jsonl.gz','rt') as f:
   for line in f:
    row=json.loads(line);ids=np.array(row['original_bond_ids'])-1;v=np.array(row['d']);new=ids[np.isnan(soft[ids])&(v>0)];soft[new]=row['time_s'];new=ids[np.isnan(fail[ids])&(v==1)];fail[new]=row['time_s']
  np.savez_compressed(case/'bond_event_times.npz',first_softening_s=soft,first_failure_s=fail,original_bond_ids=np.arange(1,nb+1),final_damage=d)
  kinds=np.full(nb,'INTERFACE',dtype='<U16');kinds[ref['bond_types']==3]='ANCHOR_SCAFFOLD';kinds[(ref['bond_types']==1)&(phase[b]==0).all(1)]='WEAK_WEAK';kinds[(ref['bond_types']==1)&(phase[b]==1).all(1)]='STRONG_STRONG';kinds[field['cross_phase']]='WEAK_STRONG'
  groups=[]
  for kind in ['WEAK_WEAK','WEAK_STRONG','STRONG_STRONG','INTERFACE']:
   sel=(kinds==kind)&~cut;groups.append(dict(kind=kind,eligible_bonds=int(sel.sum()),softened=int(np.isfinite(soft[sel]).sum()),failed=int((d[sel]==1).sum()),failed_fraction=float(np.mean(d[sel]==1)),mean_damage=float(np.mean(d[sel])),live=int((d[sel]<1).sum())))
  selected=np.flatnonzero(np.isfinite(soft)|np.isfinite(fail))
  for i in selected:eventrows.append(dict(case=mode,original_bond_id=int(i+1),atom_i=int(b[i,0]+1),atom_j=int(b[i,1]+1),phase_pair=kinds[i],domain_i=int(domain[b[i,0]]),domain_j=int(domain[b[i,1]]),first_softening_s=None if np.isnan(soft[i]) else soft[i],first_failure_s=None if np.isnan(fail[i]) else fail[i],final_damage=d[i]))
  first={}
  for key,times in [('softening',soft),('failure',fail)]:
   good=np.flatnonzero(np.isfinite(times));first_t=float(times[good].min()) if len(good) else None;ids=good[times[good]==first_t] if len(good) else good
   first[key]=dict(time_s=first_t,original_bond_ids=ids+1,phase_pairs=kinds[ids],reference_midpoints_m=mid[ids],domain_endpoint_ids=domain[b[ids]])
  # Identical descriptive early-site rule for both cases; not a continuum-crack count.
  good=np.flatnonzero(bulk&np.isfinite(soft));sites=[]
  first_bulk_t=float(soft[good].min()) if len(good) else None
  first_bulk_ids=good[soft[good]==first_bulk_t] if len(good) else good
  first['bulk_softening']=dict(time_s=first_bulk_t,original_bond_ids=first_bulk_ids+1,phase_pairs=kinds[first_bulk_ids],reference_midpoints_m=mid[first_bulk_ids],domain_endpoint_ids=domain[b[first_bulk_ids]])
  if len(good):
   early=good[soft[good]<=soft[good].min()+EARLY_WINDOW]
   for i in early[np.lexsort((early,soft[early]))]:
    if not sites or min(np.linalg.norm(mid[i]-np.array(q['reference_midpoint_m'])) for q in sites)>SITE_SEPARATION:
     sites.append(dict(original_bond_id=int(i+1),time_s=float(soft[i]),reference_midpoint_m=mid[i],phase_pair=kinds[i]))
  domains=[]
  member_records={r['component_id']:r for r in topology[mode]['fragment_lifetimes']}
  containing_free_size=np.zeros(n,dtype=int)
  for c in topology[mode]['final_components']:
   r=member_records[c['component_id']];containing_free_size[np.asarray(r['atom_ids'])-1]=r['nodes']
  for c in range(1,11):
   ids=np.flatnonzero(domain==c);inside=np.isin(b,ids).all(1);live=inside&(d<1);remap=np.full(len(ref['X']),-1);remap[ids]=np.arange(len(ids));bb=remap[b[live]];g=csr_matrix((np.ones(len(bb)*2),(np.r_[bb[:,0],bb[:,1]],np.r_[bb[:,1],bb[:,0]])),shape=(len(ids),len(ids)));nc,lab=connected_components(g,directed=False)
   cross=field['cross_phase']&np.isin(b,ids).any(1)&~cut;otherstrong=bulk&(phase[b]==1).all(1)&(domain[b[:,0]]!=domain[b[:,1]])&np.isin(b,ids).any(1)
   destination=containing_free_size[ids]
   row=dict(case=mode,strong_domain_id=c,original_nodes=len(ids),final_internal_components=nc,largest_internal_connected_nodes=int(np.bincount(lab).max()),original_internal_live_bonds=int((inside&~cut).sum()),final_internal_live_bonds=int(live.sum()),remaining_cross_phase_bonds=int((cross&(d<1)).sum()),initial_other_strong_domain_bonds=int(otherstrong.sum()),remaining_other_strong_domain_bonds=int((otherstrong&(d<1)).sum()),attached_nodes=int((destination==0).sum()),free_nodes=int((destination>0).sum()),nodes_in_target_components=int(((destination>=20)&(destination<=200)).sum()),nodes_in_gt200_components=int((destination>200).sum()),nodes_in_1_2_components=int(((destination>0)&(destination<=2)).sum()));domains.append(row);domainrows.append(row)
  stretch=np.linalg.norm(final['x'][b[:,1]]-final['x'][b[:,0]],axis=1)/ref['r0']-1
  on=field['bond_s0'] if mode=='biphasic' else S0[ref['bond_types']];bridge_state=[]
  for name,mask in [('weak_weak',bulk&(phase[b]==0).all(1)),('weak_strong',field['cross_phase']&~cut),('strong_strong',bulk&(phase[b]==1).all(1)),('direct_interdomain_strong_strong',bulk&(phase[b]==1).all(1)&(domain[b[:,0]]!=domain[b[:,1]]))]:
   live=mask&(d<1);count=int(live.sum());ratio=final['kappa'][live]/on[live]
   bridge_state.append(dict(group=name,live_bonds=count,currently_nonextended=int((live&(stretch<=1e-12)).sum()),never_exceeded_onset=int((live&(final['kappa']<=on)).sum()),mean_kappa_over_s0=float(ratio.mean()) if count else None,maximum_kappa_over_s0=float(ratio.max()) if count else None,interpretation='Geometric extension and tensile-history evidence; not a per-bond force measurement or proof that a different damage law is calibrated.'))
  comp=[]
  for r in topology[mode]['fragment_lifetimes']:
   ids=np.array(r['atom_ids'])-1;comp.append(dict(component_id=r['component_id'],nodes=len(ids),strong_nodes=int((phase[ids]==1).sum()),weak_nodes=int((phase[ids]==0).sum()),strong_domain_ids=np.unique(domain[ids][domain[ids]>0]),birth_s=r['birth_s'],lifetime_s=r['lifetime_s'],right_censored=r['right_censored']))
  out[mode]=dict(phase_label_scope='H uses actual biphasic properties; U uses the same H reference-region labels only as spatial covariates, with uniform bulk thresholds .15/.18 everywhere.',surviving_bridge_geometric_state=bridge_state,phase_damage=groups,first_events=first,early_spatially_separated_softening_sites=sites,early_site_window_s=EARLY_WINDOW,early_site_separation_m=SITE_SEPARATION,multisite_spatial_evidence=len(sites)>=2,domain_retention=domains,fragment_phase_composition=comp,remaining_live_cross_phase_bonds=int((field['cross_phase']&~cut&(d<1)).sum()),remaining_direct_strong_to_other_strong_domain_bonds=int((bulk&(phase[b]==1).all(1)&(domain[b[:,0]]!=domain[b[:,1]])&(d<1)).sum()),true_graph_separated_free_components=len(results[mode]['final_free_sizes']),target_components_ge5=results[mode]['objective']['n_fragments_20_200']>=5,continuum_3D_crack_surface_certification='NOT_ESTABLISHED_FROM_POINT_GRAPH',interpretation='Early separated softened bonds support multiple spatial initiation locations, not independent continuum crack surfaces. Strong regions are reference-grid labels, not free fragments. Actual free components certify separation in the complete surviving material graph. One spatial realization only.')
 u=results['uniform']['native']['events'];h=results['biphasic']['native']['events'];out['difference']={key+'_same_original_bonds':u.get(key,{}).get('original_bond_ids')==h.get(key,{}).get('original_bond_ids') for key in ['first_softening','first_failure']}
 # Follow the first released >200-node body by its original member IDs. This
 # distinguishes erosion of a real cohesive body from a collection of fines.
 large_lineage={}
 for mode in CASES:
  first_large=results[mode]['first_large_release']
  if first_large is None:
   large_lineage[mode]=dict(status='NOT_OBSERVED');continue
  parent=set(first_large['atom_ids']);life={r['component_id']:r for r in topology[mode]['fragment_lifetimes']};children=[];represented=set()
  for c in topology[mode]['final_components']:
   r=life[c['component_id']];members=set(r['atom_ids'])
   if not members.intersection(parent):continue
   assert members<=parent,'Disconnected material components cannot merge'
   represented.update(members);ids=np.array(sorted(members))-1
   children.append(dict(component_id=c['component_id'],nodes=len(ids),internal_alive_bonds=c['internal_alive_bonds'],target_20_200=20<=len(ids)<=200,strong_nodes=int((phase[ids]==1).sum()),weak_nodes=int((phase[ids]==0).sum()),strong_domain_ids=np.unique(domain[ids][domain[ids]>0]),original_atom_ids=ids+1))
  assert represented==parent,'All previously free original IDs must remain free'
  large_lineage[mode]=dict(status='PASS_ORIGINAL_ID_DESCENDANT_CLOSURE',parent_component_id=first_large['component_id'],parent_nodes=len(parent),parent_birth_s=first_large['birth_s'],parent_original_atom_ids=sorted(parent),last_accepted_time_s=results[mode]['native']['time_s'],final_descendants=sorted(children,key=lambda r:-r['nodes']),largest_descendant_nodes=max(r['nodes'] for r in children),largest_descendant_fraction_of_parent=max(r['nodes'] for r in children)/len(parent),target_descendant_count=sum(r['target_20_200'] for r in children),interpretation='Exact original-member ancestry. A >200-node body is a real cohesive free component but lies outside the specified target band. U phase labels are reference covariates only.')
 dump(O/'FIRST_LARGE_COMPONENT_DESCENDANTS.json',large_lineage)
 for name,rows in [('NATIVE_EVENT_HISTORY.csv',eventrows),('STRONG_DOMAIN_RETENTION.csv',domainrows)]:
  with (O/name).open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['status']);w.writeheader();w.writerows(rows)
 dump(O/'SPATIAL_DAMAGE_ANALYSIS.json',out);return out
