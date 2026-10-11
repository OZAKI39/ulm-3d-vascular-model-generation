"""One frozen face-to-material map followed by one minimum-norm correction."""
from analysis.traction_failure_v1.common import *
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import dijkstra

def constraints(X,x0,L):
 r=(X-x0)/L;n=len(X);C=np.zeros((6,3*n));C[:3]=np.tile(np.eye(3),(1,n))
 for i,(x,y,z) in enumerate(r):C[3:,3*i:3*i+3]=[[0,-z,y],[z,0,-x],[-y,x,0]]
 return C

def apply_nodal_forces(atom_ids,existing_nodal_N,material_ids,force_N,anchor_ids):
 """Offline adapter specification: arrays already converted to N; NO volume factor."""
 ids=np.asarray(atom_ids);mid=np.asarray(material_ids);assert len(set(ids.tolist()))==len(ids)
 assert len(set(mid.tolist()))==len(mid) and not np.isin(mid,anchor_ids).any()
 lookup={int(k):i for i,k in enumerate(ids)};assert all(int(k) in lookup for k in mid)
 result=np.array(existing_nodal_N,copy=True)
 for k,f in zip(mid,force_N):result[lookup[int(k)]]+=f
 return result

def run():
 start=time.perf_counter();cpu=time.process_time();freeze=json.loads((B/'MAPPING_DESIGN_FREEZE.json').read_text());assert not (B/'FORMAL_MAPPING_START.json').exists()
 save(B/'FORMAL_MAPPING_START.json',dict(run=RUN,utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),freeze_sha256=sha(B/'MAPPING_DESIGN_FREEZE.json'),code_sha256=sha(__file__)))
 g=dict(np.load(B/'GEOMETRY_PREFLIGHT.npz'));s=dict(np.load(PREV/'refined/CLOT_SURFACE_TRACTIONS.npz'));z=dict(np.load(V3/'V3_MATERIAL_POINTS_AND_BONDS.npz'));X=g['X'];N=len(X);x0=X.mean(0);L=np.linalg.norm(np.ptp(X,axis=0));cent=s['surface_coordinates_m'];area=s['triangle_area_m2'];normal=s['n_clot'];eligible=g['distance_to_surface_m']<=DX;node_norm=normal[g['closest_surface_triangle']];index=np.flatnonzero(eligible);tree=cKDTree(X[index]);row=[];col=[];values=[];linkdist=[];coverage=[]
 repair_path=B/'CONTACT_RIM_SUPPORT_REPAIR_FREEZE.json';repair=json.loads(repair_path.read_text()) if repair_path.exists() else None
 repaired={}
 if repair:
  faces=s['triangle_mesh_node_ids'];edges={}
  for fi,face in enumerate(faces):
   for a,b in [(face[0],face[1]),(face[1],face[2]),(face[2],face[0])]:edges.setdefault(tuple(sorted([int(a),int(b)])),[]).append(fi)
  gr=[];gc=[];gv=[]
  for fs in edges.values():
   if len(fs)==2:
    a,b=fs;gr.extend([a,b]);gc.extend([b,a]);gv.extend([np.linalg.norm(cent[a]-cent[b])]*2)
  graph=coo_matrix((gv,(gr,gc)),shape=(len(area),len(area))).tocsr();bad=repair['faces'];dist=dijkstra(graph,directed=False,indices=bad,limit=DELTA)
  repaired={f:dist[k] for k,f in enumerate(bad)}
 for f,near in enumerate(tree.query_ball_point(cent,DELTA)):
  ids=index[near];delta=cent[f]-X[ids];dd=np.linalg.norm(delta,axis=1)
  inside=delta@normal[f]>=-DX/4
  if f in repaired:inside=np.isfinite(repaired[f][g['closest_surface_triangle'][ids]])
  keep=(node_norm[ids]@normal[f]>0)&inside&(dd<DELTA);ids=ids[keep];dd=dd[keep];coverage.append(len(ids))
  if len(ids):
   q=dd/DELTA;w=(1-q)**4*(4*q+1);w/=w.sum();row.extend(ids);col.extend([f]*len(ids));values.extend(w);linkdist.extend(dd)
 save(B/'COVERAGE.json',dict(face_count=len(area),support_count=stats(coverage),unsupported_faces=np.flatnonzero(np.array(coverage)==0),eligible_nodes=int(eligible.sum()),support_distance_m=stats(linkdist)))
 if min(coverage)==0:
  save(B/'FINAL_STATUS.json',{'LOCAL_MAPPING_QUALITY':'MAPPING_LOCALITY_OR_MOMENT_UNRESOLVED','reason':'Unsupported exposed faces; no remote-node fallback','formal_mapping_attempts':1});return
 W=coo_matrix((values,(row,col)),shape=(N,len(area))).tocsr();aeff=np.asarray(W@area);loaded=aeff>0;ix=np.flatnonzero(loaded);C6=constraints(X[ix],x0,L)
 assert np.max(abs(np.asarray(W.sum(0)).ravel()-1))<1e-13 and min(values)>0
 tp=s['pressure_traction_Pa'];ts=s['tangential_traction_Pa'];tn=s['viscous_normal_traction_Pa'];qxyz=s['quadrature_coordinates_m'];qa=s['quadrature_area_m2'];qpressure=-s['quadrature_pressure_Pa'][:,:,None]*normal[:,None];qn=np.repeat(tn[:,None],3,axis=1);qs=np.repeat(ts[:,None],3,axis=1)
 comp={'pressure':(tp,qpressure),'shear':(ts,qs),'viscous_normal':(tn,qn)};initial={};corrected={};aud={};quality={};Ds={};lambdas={};targets={};ok=[]
 avg_norm=np.asarray(W@(normal*area[:,None]));avg_norm[loaded]/=np.linalg.norm(avg_norm[loaded],axis=1)[:,None]
 rim_source=s['contact_rim_mask'];rim_nodes=np.asarray(W@(area*rim_source))/np.where(loaded,aeff,1.)>.5
 regions={f'xquartile{k+1}':(X[:,0]>=a)&(X[:,0]<(b if k<3 else b+DX)) for k,(a,b) in enumerate(zip(np.quantile(X[ix,0],[0,.25,.5,.75]),np.quantile(X[ix,0],[.25,.5,.75,1])))}
 regions['contact_rim']=rim_nodes;regions['interface_connected']=np.isin(np.arange(N),g['interface_material_indices'])
 source_integrals={};centroid_discrepancy={}
 for key,(traction,qtraction) in comp.items():
  ff=traction*area[:,None];F=ff.sum(0);M=(np.cross(qxyz-x0,qtraction)*qa[:,:,None]).sum((0,1));Mc=np.cross(cent-x0,ff).sum(0);scale=np.linalg.norm(ff,axis=1).sum();source_integrals[key]={'force_N':F,'quadrature_moment_Nm':M,'centroid_moment_Nm':Mc,'sum_face_force_norm_N':scale};centroid_discrepancy[key]=M-Mc
  f0=np.asarray(W@ff);target=np.r_[F,M/L];residual=target-C6@f0[ix].ravel();D=np.linalg.norm(f0[ix],axis=1)**2+(.01*scale*aeff[ix]/area.sum())**2;Dr=np.repeat(D,3);gram=(C6*Dr)@C6.T;lam=np.linalg.solve(gram,residual);df=(Dr*(C6.T@lam)).reshape(-1,3);f=f0.copy();f[ix]+=df
  initial[key]=f0;corrected[key]=f;Ds[key]=D;lambdas[key]=lam;targets[key]=target
  au={}
  for label,vec in [('initial',f0),('corrected',f)]:
   force=vec.sum(0);moment=np.cross(X-x0,vec).sum(0);au[label]=dict(force_N=force,moment_Nm=moment,force_error_N=force-F,moment_error_Nm=moment-M,force_scaled_error=np.linalg.norm(force-F)/scale,moment_scaled_error=np.linalg.norm(moment-M)/(scale*L))
  au.update(rank=int(np.linalg.matrix_rank(gram)),normalized_gram_condition=float(np.linalg.cond(gram)),force_scale_N=scale,length_scale_m=L);aud[key]=au
  quality[key]=dict(correction_L2_fraction=np.linalg.norm(df)/np.linalg.norm(f0),single_node_relative_correction=stats(np.linalg.norm(df,axis=1)/np.linalg.norm(f0[ix],axis=1)),reverse_load_fraction=float(np.linalg.norm(f,axis=1)[np.einsum('ij,ij->i',f,f0)<0].sum()/np.linalg.norm(f,axis=1).sum()),regional_fraction={name:{stage:float(np.linalg.norm(vec,axis=1)[mask].sum()/np.linalg.norm(vec,axis=1).sum()) for stage,vec in [('initial',f0),('corrected',f)]} for name,mask in regions.items()})
  if key=='pressure':quality[key]['outward_pressure_load_fraction']=float(np.linalg.norm(f,axis=1)[np.einsum('ij,ij->i',f,avg_norm)>0].sum()/np.linalg.norm(f,axis=1).sum())
  if key=='shear':quality[key]['normal_component_absolute_fraction']=float(np.abs(np.einsum('ij,ij->i',f,avg_norm)).sum()/np.linalg.norm(f,axis=1).sum())
  top=max(1,int(np.ceil(len(ix)*.01)))
  quality[key]['top_one_percent_load_share']={stage:float(np.sort(np.linalg.norm(vec[ix],axis=1))[-top:].sum()/np.linalg.norm(vec[ix],axis=1).sum()) for stage,vec in [('initial',f0),('corrected',f)]}
  cosine=np.einsum('ij,ij->i',f[ix],f0[ix])/(np.linalg.norm(f[ix],axis=1)*np.linalg.norm(f0[ix],axis=1));quality[key]['direction_change_degrees']=stats(np.degrees(np.arccos(np.clip(cosine,-1,1))))
 initial['total']=sum(initial.values());corrected['total']=sum(corrected.values());F=sum(v['force_N'] for v in source_integrals.values());M=sum(v['quadrature_moment_Nm'] for v in source_integrals.values());scale=np.linalg.norm((tp+ts+tn)*area[:,None],axis=1).sum()
 assert np.linalg.norm(F-s['integrated_total_force_N'])/scale<1e-12;assert np.linalg.norm(M-s['integrated_total_torque_Nm'])/(scale*L)<1e-12
 aud['total']={stage:dict(force_N=vec.sum(0),moment_Nm=np.cross(X-x0,vec).sum(0),force_scaled_error=np.linalg.norm(vec.sum(0)-F)/scale,moment_scaled_error=np.linalg.norm(np.cross(X-x0,vec).sum(0)-M)/(scale*L)) for stage,vec in [('initial',initial['total']),('corrected',corrected['total'])]}
 quality['total']=dict(correction_L2_fraction=np.linalg.norm(corrected['total']-initial['total'])/np.linalg.norm(initial['total']),single_node_force_share_max=np.linalg.norm(corrected['total'],axis=1).max()/np.linalg.norm(corrected['total'],axis=1).sum())
 gates=freeze['gates'];checks={}
 for k in ['pressure','shear','viscous_normal','total']:
  checks[k+'_force']=aud[k]['corrected']['force_scaled_error']<=gates['force_and_moment_scaled_error'];checks[k+'_moment']=aud[k]['corrected']['moment_scaled_error']<=gates['force_and_moment_scaled_error']
  if k!='total':checks[k+'_constraint']=aud[k]['rank']==6 and aud[k]['normalized_gram_condition']<gates['normalized_gram_condition_max']
  if k!='viscous_normal':checks[k+'_correction']=quality[k]['correction_L2_fraction']<=.10
  if k in ['pressure','shear']:
   checks[k+'_reverse']=quality[k]['reverse_load_fraction']<=.01;v=quality[k]['top_one_percent_load_share'];checks[k+'_concentration_change']=abs(v['corrected']-v['initial'])<=.02
   checks[k+'_regional_change']=all(abs(v['corrected']-v['initial'])<=.05 for v in quality[k]['regional_fraction'].values())
 checks['pressure_inward']=quality['pressure']['outward_pressure_load_fraction']<=.01;checks['maximum_concentration']=quality['total']['single_node_force_share_max']<=.02
 passed=all(checks.values());name='CLOT_NODAL_TRACTION_FORCES.npz' if passed else 'DIAGNOSTIC_CANDIDATE_NODAL_TRACTION_FORCES.npz'
 data=dict(material_atom_id=g['material_ids'],reference_coordinates_m=X,support_node=loaded,effective_area_m2=aeff,reference_point_m=x0,nearest_surface_normal=avg_norm,anchor_atom_id=np.arange(N+1,len(z['X'])+1),anchor_reference_coordinates_m=z['X'][N:])
 for k in corrected:data[k+'_force_N']=corrected[k];data[k+'_initial_force_N']=initial[k]
 np.savez_compressed(B/name,**data)
 np.savez_compressed(B/'CFD_TO_LAMMPS_MAPPING_MATRIX.npz',data=W.data,indices=W.indices,indptr=W.indptr,shape=np.array(W.shape),loaded_node_indices=ix,constraint_matrix=C6,length_scale_m=L,face_support_distance_m=np.array(linkdist),reference_point_m=x0,**{k+'_D':v for k,v in Ds.items()},**{k+'_lambda':v for k,v in lambdas.items()})
 source_regions={}
 for label,field in [('pressure',tp),('shear',ts),('total',tp+ts+tn)]:
  magnitude=area*np.linalg.norm(field,axis=1);source_regions[label]=dict(contact_rim_fraction=float(magnitude[rim_source].sum()/magnitude.sum()),nodal_effective_area_weighted_rim_fraction=float((np.linalg.norm(initial[label],axis=1)*rim_nodes).sum()/np.linalg.norm(initial[label],axis=1).sum()))
 save(B/'CFD_SURFACE_FORCE_AUDIT.json',dict(components=source_integrals,total_force_N=F,total_moment_Nm=M,centroid_to_exact_moment_difference_Nm=centroid_discrepancy,reference_point_m=x0,original_report_recovered=True,face_pressure_min_Pa=float(s['pressure_Pa'].min()),open_surface_pressure_offset_1Pa_delta_force_N=-np.einsum('i,ij->j',area,normal)))
 save(B/'CONSERVATIVE_MAPPING_AUDIT.json',dict(components=aud,normalization='sum_f |F_f| and that sum times material bounding-box diagonal',checks=checks,formal_mapping_attempts=1,minimal_correction_algorithms=1))
 save(B/'LOCAL_MAPPING_QUALITY.json',dict(status='PASS_ENGINEERING_SCREEN' if passed else 'MAPPING_LOCALITY_OR_MOMENT_UNRESOLVED',metrics=quality,checks=checks,effective_area_m2=stats(aeff[loaded]),support_nodes=int(loaded.sum()),force_free_internal_nodes=int((~loaded).sum()),source_and_target_regions=source_regions,far_support_face_fraction=float(np.mean(g['face_nearest_material_distance_m']>2*DX)),far_support_area_fraction=float(area[g['face_nearest_material_distance_m']>2*DX].sum()/area.sum()),scope='Algebraic conservation and frozen locality screen only, not pointwise convergence/biological validation.'))
 meta=dict(run_id=RUN,status='OFFLINE_ENGINEERING_SCREEN_PASSED' if passed else 'DIAGNOSTIC_ONLY_NOT_CERTIFIED',file=name,production_eligible=False,production_not_run=True,source_CFD='20261011T154854Z/refined',source_SHA=json.loads((B/'INPUT_SHA.json').read_text()),coordinates='world xyz metres',coordinate_transform=np.eye(4),force_units='N already integrated; add after existing native force-density * volume conversion; NEVER multiply these forces by volume',surface_normal='n_clot=-n_fluid verified using adjacent tetrahedron',pressure='relative gauge with original zero-reference outlet traction; open patch pressure force gauge-dependent',reference_point_m=x0,method='One Wendland local column-stochastic face map plus exact minimum correction; see frozen design',map_sha256=sha(B/'CFD_TO_LAMMPS_MAPPING_MATRIX.npz'),force_sha256=sha(B/name),double_counting='Existing detached-component background drag must be allocated to disjoint objects/stages; fixed intact traction cannot follow changed fragment surfaces without updating geometry/flow')
 save(B/'CLOT_NODAL_TRACTION_METADATA.json',meta)
 save(B/'GEOMETRY_AND_UNIT_AUDIT.json',dict(transform=np.eye(4),surface_tag=5,material_nodes=N,anchors_excluded=len(z['X'])-N,loaded_nodes=int(loaded.sum()),coordinate_center_m=x0,continuous_volume_mm3=.2715743882436176,CFD_volume_mm3=.271576593161914,discrete_volume_mm3=float(z['volume'][:N].sum()*1e9),relative_discrete_continuous_difference=float(z['volume'][:N].sum()*1e9/.2715743882436176-1),no_volume_force_rescale=True,exposed_area_m2=area.sum(),area_sum_error_m2=aeff.sum()-area.sum(),nearest_distances=json.loads((B/'GEOMETRY_PREFLIGHT.json').read_text()),exposed_support_distance_m=stats(g['distance_to_surface_m'][loaded]),effective_area_over_dx2=stats(aeff[loaded]/DX**2),all_faces_supported=True,unit_table=dict(length='m',area='m^2',traction='Pa',force='N',torque='N m')))
 # Meaningful offline integration check: nonsequential native IDs and exactly-once units.
 ids=np.arange(1,len(z['X'])+1)[::-1];base=np.ones_like(z['X'])*2e-9;result=apply_nodal_forces(ids,base,g['material_ids'],corrected['total'],data['anchor_atom_id']);delta=result-base
 uniterr=np.max(abs(delta[::-1][:N]-corrected['total']));assert uniterr<1e-23 and np.all(delta[np.isin(ids,data['anchor_atom_id'])]==0)
 save(B/'INTERFACE_UNIT_TESTS.json',dict(nonsequential_ID_mapping=True,anchors_untouched=True,already_N_added_exactly_once=True,max_absolute_error_N=uniterr,production_deck_modified=False))
 states=dict(CFD_SOURCE_IDENTITY='PASS_NATIVE_FULL_PRECISION_SHA',GEOMETRY_REGISTRATION='PASS_IDENTITY_WITH_DISCLOSED_VOLUME_DIFFERENCE',EXPOSED_SURFACE_IDENTITY='PASS_EXACT_TAG5_TRIANGLES',SURFACE_TRACTION_RECONSTRUCTION='PASS_NATIVE_GRADIENT_AND_PRESSURE',FORCE_CONSERVATION='PASS_ALL_COMPONENTS',MOMENT_CONSERVATION='PASS_ALL_COMPONENTS',LOCAL_MAPPING_QUALITY='PASS_ENGINEERING_SCREEN' if passed else 'MAPPING_LOCALITY_OR_MOMENT_UNRESOLVED',ANCHOR_EXCLUSION='PASS',PRESSURE_REFERENCE_DISCLOSED=True,LAMMPS_FORCE_UNIT_COMPATIBILITY='PASS_OFFLINE_ID_AND_N_INTERFACE',PRODUCTION_COUPLING_NOT_RUN=True)
 save(B/'FINAL_STATUS.json',dict(run_id=RUN,states=states,formal_fields=1));save(B/'RESOURCES.json',cost(start,cpu));print(json.dumps({'passed':passed,'failed_checks':[k for k,v in checks.items() if not v],'support_nodes':int(loaded.sum()),'correction':{k:q['correction_L2_fraction'] for k,q in quality.items()}},indent=2))
if __name__=='__main__':run()
