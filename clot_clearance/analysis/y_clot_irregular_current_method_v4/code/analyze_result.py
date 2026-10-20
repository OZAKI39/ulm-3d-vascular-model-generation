"""Compare saved native trajectories; graph membership, never spatial clustering."""
from common import *
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components,dijkstra
import csv,gzip,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASE=O/'cases/H_Y_IRREGULAR_LPS_V4';OLD=Y_OLD/'cases/biphasic'
BINS=[(1,2),(3,13),(14,19),(20,33),(34,113),(114,200),(201,1000000)]
BIN_LABELS=['1–2','3–13','14–19','20–33','34–113','114–200','>200']
def table(path):
 with Path(path).open() as f:return list(csv.DictReader(f))
def topology(ref,d):
 b=ref['bonds'][d<1];N=len(ref['X']);nc,lab=connected_components(csr_matrix((np.ones(2*len(b)),(np.r_[b[:,0],b[:,1]],np.r_[b[:,1],b[:,0]])),shape=(N,N)),directed=False);attached=np.zeros(nc,bool);attached[lab[int(ref['clot_count']):]]=True
 return lab,attached
def state_paths(case):
 paths=[case/'keys/initial.npz']+list((case/'keys').glob('frame_*.npz'))
 events=['pilot','final','first_softening','first_failure','first_free_component','first_multinode','first_target_observed','max_target_count']
 for p in [case/'keys'/(name+'.npz') for name in events]+[case/'last_accepted.npz']:
  if p.exists():paths.append(p)
 found={}
 for p in paths:
  with np.load(p) as q:found[float(q['time'])]=p
 return sorted(found.items())
def choose(paths,t):
 valid=[q for q in paths if q[0]<=t+1e-14];return valid[-1] if valid else paths[0]
def fragments(ref,state):
 lab,attached=topology(ref,state['d']);order=np.argsort(lab,kind='stable');starts=np.r_[0,np.cumsum(np.bincount(lab))];b=ref['bonds'];alive=state['d']<1;nb=np.bincount(lab[b[alive,0]],minlength=len(attached));rows=[]
 boundary_mask=(lab[b[:,0]]!=lab[b[:,1]])&(ref['bond_types']==1)&~ref['initial_cut'];boundary=b[boundary_mask];directions=(ref['X'][boundary[:,1]]-ref['X'][boundary[:,0]])/ref['r0'][boundary_mask,None];angular=np.sum(directions**4,axis=1);bl=lab[boundary].ravel();boundary_count=np.bincount(bl,minlength=len(attached));boundary_sum=np.bincount(bl,weights=np.repeat(angular,2),minlength=len(attached))
 for k in np.flatnonzero(~attached):
  ids=np.sort(order[starts[k]:starts[k+1]]);vol=ref['V'][ids];w=vol/vol.sum();x=state['x'][ids];center=w@x;dx=x-center;cov=(dx*w[:,None]).T@dx;eig,axes=np.linalg.eigh(cov);eig=eig[::-1];axes=axes[:,::-1];r=np.sqrt(np.maximum(eig,0));ratio=r/max(r[0],1e-300)
  rows.append(dict(atom_ids=(ids+1).tolist(),nodes=len(ids),reference_volume_m3=float(vol.sum()),mass_kg=float(RHO*vol.sum()),equivalent_nominal_nodes=float(vol.sum()/V),equivalent_sphere_diameter_m=float((6*vol.sum()/np.pi)**(1/3)),centroid_m=center,principal_axes=axes,principal_radii_m=r,principal_inertia_kg_m2=RHO*vol.sum()*(eig.sum()-eig),aspect_ratios=ratio,radius_of_gyration_m=float(np.sqrt(eig.sum())),internal_alive_bonds=int(nb[k]),local_graph_density=float(2*nb[k]/(len(ids)*(len(ids)-1))) if len(ids)>1 else None,boundary_failed_original_bonds=int(boundary_count[k]),boundary_reference_axis_fourth_moment=float(boundary_sum[k]/boundary_count[k]) if boundary_count[k] else None,boundary_metric_scope='Angular statistic of severed incident original bulk connections; not a continuum surface roughness',morphology='SMALL_DEBRIS' if len(ids)<20 else 'FILAMENT' if ratio[1]<.2 else 'THIN_SHEET' if ratio[2]<.2 else 'COHESIVE_BLOCK'))
 return rows
def crack_metrics(ref,state):
 b=ref['bonds'];mask=(state['d']==1)&~ref['initial_cut']&(ref['bond_types']==1);broken=b[mask];r=ref['r0'][mask]
 if not len(broken):return dict(bulk_failed=0,direction_tensor=None,axis_alignment_fourth_moment=None,broken_edge_graph_path_ratio=None)
 directions=(ref['X'][broken[:,1]]-ref['X'][broken[:,0]])/r[:,None];tensor=directions.T@directions/len(directions)
 N=len(ref['X']);graph=csr_matrix((np.r_[r,r],(np.r_[broken[:,0],broken[:,1]],np.r_[broken[:,1],broken[:,0]])),shape=(N,N));_,lab=connected_components(graph,directed=False);counts=np.bincount(lab);members=np.flatnonzero(lab==counts.argmax());source=members[np.argmin(ref['X'][members,0])];dist=dijkstra(graph,directed=False,indices=source);ids=members[np.isfinite(dist[members])];end=ids[np.argmax(dist[ids])];chord=float(np.linalg.norm(ref['X'][end]-ref['X'][source]));ratio=float(dist[end]/chord) if chord>0 else None
 return dict(bulk_failed=len(broken),direction_tensor=tensor,direction_eigenvalues=np.linalg.eigvalsh(tensor),axis_alignment_fourth_moment=float(np.mean(np.sum(directions**4,axis=1))),broken_edge_graph_path_ratio=ratio,path_endpoint_original_ids=[int(source+1),int(end+1)],path_graph_distance_m=float(dist[end]),path_chord_m=chord,path_interpretation='Shortest path on broken original edges in the largest failed-edge network, between an X-extreme seed and its farthest reachable node. A discrete network diagnostic, NOT an identified continuum crack path or surface.',midpoint_covariance_m2=np.cov(ref['X'][broken].mean(1).T) if len(broken)>1 else np.zeros((3,3)))
def lineage_data(name):
 if name=='V4':return json.loads((CASE/'LINEAGE.json').read_text())
 reg=json.loads((OLD/'COMPONENT_REGISTRY.json').read_text());splits=json.loads((OLD/'SUMMARY.json').read_text())['secondary_events'];return dict(registry=reg,splits=splits)
def survival(data,final_time):
 death={int(s['parent_component_id']):float(s['time_s']) for s in data['splits']};target=[]
 for cid,q in data['registry'].items():
  if 20<=q['nodes']<=200 and q['birth_time_s']<=final_time+1e-14:
   actual=death.get(int(cid));stop=actual if actual is not None else float('inf');censored=stop>final_time;target.append(dict(component_id=int(cid),nodes=q['nodes'],birth_time_s=q['birth_time_s'],split_time_s=actual if not censored else None,observed_duration_s=min(stop,final_time)-q['birth_time_s'],right_censored=censored,at_drive_off=final_time>=OFF and q['birth_time_s']<=OFF and stop>OFF,survives_to_end=censored))
 cohort=[q for q in target if q['at_drive_off']];return dict(target_member_sets=target,drive_off_reached=final_time>=OFF,unload_cohort_count=len(cohort) if final_time>=OFF else None,unload_cohort_survived=sum(q['survives_to_end'] for q in cohort) if final_time>=OFF else None,unload_cohort_split=sum(not q['survives_to_end'] for q in cohort) if final_time>=OFF else None,final_time_s=final_time)
def csv_write(path,rows):
 if not rows:path.write_text('status\nNO_OBSERVED_FRAGMENTS\n');return
 keys=list(dict.fromkeys(k for row in rows for k in row));f=path.open('w');w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
 for r in rows:w.writerow({k:json.dumps(v,default=cv) if isinstance(v,(list,dict,np.ndarray)) else v for k,v in r.items()})
 f.close()

def evolution(ref,paths,end,phase):
 b=ref['bonds'];bulk=(ref['bond_types']==1)&~ref['initial_cut'];cross=bulk&(phase[b[:,0]]!=phase[b[:,1]]);directions=(ref['X'][b[:,1]]-ref['X'][b[:,0]])/ref['r0'][:,None];fourth=np.sum(directions**4,axis=1);rows=[]
 for t,path in paths:
  if t>end+1e-14:continue
  with np.load(path) as q:d=q['d']
  failed=(d==1)&bulk;lab,attached=topology(ref,d);sizes=np.bincount(lab)[~attached];volumes=np.bincount(lab,weights=ref['V'])[~attached]
  spectra=[dict(node_bin=label,count=int(((sizes>=lo)&(sizes<=hi)).sum()),reference_volume_m3=float(volumes[(sizes>=lo)&(sizes<=hi)].sum()),mass_kg=float(RHO*volumes[(sizes>=lo)&(sizes<=hi)].sum()),equivalent_nominal_nodes=float(volumes[(sizes>=lo)&(sizes<=hi)].sum()/V)) for (lo,hi),label in zip(BINS,BIN_LABELS)]
  rows.append(dict(time_s=t,source=str(path.relative_to(R)),bulk_failed=int(failed.sum()),bulk_softened=int(((d>0)&bulk).sum()),cross_phase_reference_bonds=int(cross.sum()),cross_phase_surviving_bonds=int((cross&(d<1)).sum()),cross_phase_surviving_fraction=float((cross&(d<1)).sum()/cross.sum()),reference_axis_fourth_moment=float(fourth[bulk].mean()),failed_axis_fourth_moment=float(fourth[failed].mean()) if failed.any() else None,failed_minus_reference_axis_moment=float(fourth[failed].mean()-fourth[bulk].mean()) if failed.any() else None,free_component_spectra=spectra))
 return rows

def main():
 result=json.loads((O/'NATIVE_V4_RESULT.json').read_text())
 assert result['status'] in ['COMPLETED','NUMERICAL_FAILURE']
 tf=result['accepted_time_s']; (F/'figures').mkdir(parents=True,exist_ok=True)
 refs={'Y regular':load_npz(Y_OLD/'provenance/reference.npz'),'Y V4':load_npz(O/'provenance/reference.npz')}
 paths={'Y regular':state_paths(OLD),'Y V4':state_paths(CASE)}
 states={k:load_npz(choose(paths[k],tf)[1]) for k in refs}
 groups={k:fragments(refs[k],states[k]) for k in refs}
 history={'Y regular':table(OLD/'history.csv'),'Y V4':table(CASE/'history.csv')}
 phase=load_npz(O/'BIPHASIC_FIELD.npz')['material_phase']
 ev={k:evolution(refs[k],paths[k],tf,phase) for k in refs}
 dump(O/'CRACK_EVOLUTION.json',ev)
 dump(O/'CRACK_METRICS.json',{k:crack_metrics(refs[k],states[k]) for k in refs})
 csv_write(O/'FRAGMENT_MORPHOLOGY.csv',[dict(case=k,time_s=float(states[k]['time']),**q) for k in groups for q in groups[k]])
 csv_write(O/'FRAGMENT_SPECTRA_HISTORY.csv',[dict(case=k,time_s=q['time_s'],**b) for k in ev for q in ev[k] for b in q['free_component_spectra']])
 csv_write(O/'ENERGY_LEDGER.csv',[{k:v for k,v in q.items() if k in ['step','time_s','dt_s','W_J','W_drive_J','W_drag_J','W_wall_J','U_J','K_J','D_J','Uc_J','Jc_J','R_J','relative_R','wall_potential_J','min_wall_clearance_m']} for q in history['Y V4']])
 surv={'Y regular':survival(lineage_data('V3'),float(states['Y regular']['time'])),'Y V4':survival(lineage_data('V4'),tf)}
 dump(O/'FRAGMENT_SURVIVAL.json',surv)
 stats={};spectra={}
 for k,rows in groups.items():
  n=int(refs[k]['clot_count']);mass=RHO*refs[k]['V'][:n].sum()
  stats[k]=dict(time_s=float(states[k]['time']),free_nodes=sum(q['nodes'] for q in rows),free_components=len(rows),largest_nodes=max((q['nodes'] for q in rows),default=0),target_20_200=sum(20<=q['nodes']<=200 for q in rows),peak_target_count=max(int(float(q['N20_200'])) for q in history[k] if float(q['time_s'])<=tf+1e-14),peak_largest_nodes=max(int(float(q['largest_nodes'])) for q in history[k] if float(q['time_s'])<=tf+1e-14),free_mass_kg=sum(q['mass_kg'] for q in rows),target_mass_kg=sum(q['mass_kg'] for q in rows if 20<=q['nodes']<=200),debris_1_2_mass_fraction=sum(q['mass_kg'] for q in rows if q['nodes']<=2)/mass,total_material_mass_kg=float(mass))
  spectra[k]=dict(node_bins=BIN_LABELS,counts=[sum(lo<=q['nodes']<=hi for q in rows) for lo,hi in BINS],mass_kg=[sum(q['mass_kg'] for q in rows if lo<=q['nodes']<=hi) for lo,hi in BINS])
 dump(O/'FRAGMENT_SIZE_SPECTRA.json',spectra)
 h=history['Y V4'];energy=max(h,key=lambda q:float(q['relative_R']))
 outcome=dict(status=result['status'],accepted_time_s=tf,formal_trajectories=1,old_case_reused_without_rerun=True,statistics=stats,maximum_normalized_energy_residual=float(energy['relative_R']),energy_peak_time_s=float(energy['time_s']),minimum_wall_clearance_m=min(float(q['min_wall_clearance_m']) for q in h),events=result['events'],wall_s=result['wall_s'],peak_RSS_bytes=result['peak_RSS_bytes'],accepted_steps=result['accepted_steps'],unload_survival=surv,scope='Single exploratory conditioned irregular discretization; engineered original Y geometry. Mechanical surface-shear surrogate, no new CFD/acoustics/FSI. Native agreement is numerical verification, not material validation.')
 dump(O/'SCIENTIFIC_OUTCOME.json',outcome)
 dump(O/'COMPUTE_COST.json',{k:result[k] for k in ['wall_s','cpu_s','peak_RSS_bytes','accepted_steps','rejected_predictions','substeps','accepted_time_s']})
 plt.rcParams.update({'font.size':14,'axes.grid':True,'grid.alpha':.45,'grid.color':'#8e9aa4','savefig.facecolor':'white'})
 fig,ax=plt.subplots(2,2,figsize=(15,10),layout='constrained')
 for j,k in enumerate(history):
  rows=[q for q in history[k] if float(q['time_s'])<=tf+1e-14];t=np.array([float(q['time_s']) for q in rows])*1e3;c=['#778c9b','#a31f3a'][j]
  for a,key,title in [(ax[0,0],'N20_200','20–200-node free components'),(ax[0,1],'largest_nodes','Largest free component (nodes)'),(ax[1,0],'actual_peak_slip_shear_Pa','Peak slip traction (Pa)'),(ax[1,1],'effective_loaded_area_m2','Effective loaded area (mm²)')]:
   if key in rows[0]:a.plot(t,np.array([float(q[key]) for q in rows])*(1e6 if key.endswith('m2') else 1),color=c,label=k)
   a.set(xlabel='Physical time (ms)',ylabel=title);a.legend()
 fig.savefig(F/'figures/FIG_04_LOAD_AND_FRAGMENT_HISTORY.png',dpi=180);plt.close(fig)
 fig,ax=plt.subplots(2,1,figsize=(15,9),layout='constrained');t=np.array([float(q['time_s']) for q in h])*1e3
 for key in ['W_drive_J','W_drag_J','W_wall_J','U_J','K_J','D_J','Uc_J','Jc_J']:
  ax[0].plot(t,np.array([float(q[key]) for q in h])*1e6,label=key.replace('_J',''))
 ax[0].set(ylabel='Energy (µJ)',title='Wdrive + Wdrag + Wwall = ΔU + ΔK + D + ΔUc − Jc + R');ax[0].legend(ncol=4)
 ax[1].plot(t,[100*float(q['relative_R']) for q in h],c='#a31f3a');ax[1].axhline(1,c='#333',ls='--',label='1% numerical gate');ax[1].set(xlabel='Physical time (ms)',ylabel='Normalized energy residual (%)');ax[1].legend()
 fig.savefig(F/'figures/FIG_05_ENERGY_LEDGER.png',dpi=180);plt.close(fig)
 fig,ax=plt.subplots(1,2,figsize=(16,6),layout='constrained')
 for j,k in enumerate(spectra):
  for z,key,scale in [(0,'counts',1),(1,'mass_kg',1e6)]:
   vals=np.array(spectra[k][key])*scale;bars=ax[z].bar(np.arange(len(BINS))+(j-.5)*.37,vals,width=.37,color=['#778c9b','#a31f3a'][j],label=k);ax[z].bar_label(bars,labels=[f'{v:.3g}' for v in vals],fontsize=11,padding=4,rotation=90 if z else 0)
 for a in ax:a.set_xticks(np.arange(len(BINS)),BIN_LABELS);a.set_xlabel('Original nodes per free component');a.legend();a.margins(y=.25)
 ax[0].set_yscale('symlog',linthresh=1);ax[0].set_ylabel('Free component count (symlog)');ax[1].set_ylabel('Reference mass (mg)')
 fig.savefig(F/'figures/FIG_06_FRAGMENT_SIZE_DISTRIBUTION.png',dpi=180);plt.close(fig)
 print('ANALYSIS_COMPLETE',json.dumps(stats,default=cv),flush=True)

if __name__=='__main__':main()

