"""Compare saved native trajectories; graph membership, never spatial clustering."""
from common import *
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components,dijkstra
import csv,gzip,shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

CASE=O/'cases/H_IRREGULAR_LPS_V4';OLD=V3_O/'cases/biphasic'
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

def plot_spectra(spectra,tf):
 fig,axes=plt.subplots(1,2,figsize=(16,6),layout='constrained')
 for j,k in enumerate(['V3','V4']):
  counts=spectra[k]['counts'];mass=np.array(spectra[k]['mass_kg'])*1e6
  bars0=axes[0].bar(np.arange(7)+(j-.5)*.36,counts,width=.36,color=['#718c9e','#a41c3b'][j],label=k)
  bars1=axes[1].bar(np.arange(7)+(j-.5)*.36,mass,width=.36,color=['#718c9e','#a41c3b'][j],label=k)
  axes[0].bar_label(bars0,labels=[str(n) for n in counts],padding=4,fontsize=9)
  axes[1].bar_label(bars1,labels=[f'{m:.3g}' for m in mass],padding=5,fontsize=10,rotation=90)
 for ax in axes:ax.set_xticks(np.arange(7),BIN_LABELS);ax.set_xlabel('Original material nodes / surviving component');ax.legend()
 axes[0].set_yscale('symlog',linthresh=1);axes[0].set_ylabel('Free component count (symlog scale)');axes[1].set_ylabel('Reference mass (mg)');axes[0].margins(y=.2);axes[1].margins(y=.2)
 fig.suptitle(f'Surviving original material graph; matched accepted time ≈ {tf*1e3:.4f} ms');fig.savefig(F/'figures/FIG_07_FRAGMENT_SIZE_DISTRIBUTION.png',dpi=180);plt.close(fig)

def main():
 result=json.loads((O/'NATIVE_V4_RESULT.json').read_text());assert result['status'] in ['COMPLETED','NUMERICAL_FAILURE'];tf=result['accepted_time_s'];complete=result['status']=='COMPLETED'
 paths={'V3':state_paths(OLD),'V4':state_paths(CASE)};refs={'V3':load_npz(V3_O/'provenance/reference.npz'),'V4':load_npz(O/'provenance/reference.npz')};states={k:load_npz(choose(p,tf)[1]) for k,p in paths.items()};groups={k:fragments(refs[k],s) for k,s in states.items()};cracks={k:crack_metrics(refs[k],s) for k,s in states.items()}
 survival_results={k:survival(lineage_data(k),float(states[k]['time'])) for k in refs};history={'V3':table(OLD/'history.csv'),'V4':table(CASE/'history.csv')};out=F/'figures';out.mkdir(parents=True,exist_ok=True)
 dump(O/'CRACK_METRICS.json',dict(comparison_at_accepted_end=True,times_s={k:float(s['time']) for k,s in states.items()},metrics=cracks));dump(O/'FRAGMENT_SURVIVAL.json',survival_results);dump(O/'FRAGMENT_LINEAGE.json',lineage_data('V4'))
 phase=load_npz(O/'BIPHASIC_FIELD.npz')['material_phase'];evolutions={k:evolution(refs[k],paths[k],tf,phase) for k in refs};dump(O/'CRACK_EVOLUTION.json',evolutions)
 csv_write(O/'FRAGMENT_SPECTRA_HISTORY.csv',[dict(case=k,time_s=row['time_s'],**spectrum) for k,rows in evolutions.items() for row in rows for spectrum in row['free_component_spectra']])
 first_events=dict(V4=result['events'],V3=json.loads((OLD/'SUMMARY.json').read_text())['events']);dump(O/'FIRST_FRACTURE_EVENTS.json',{k:{name:{key:value for key,value in event.items() if key in ['time_s','step','evidence']} for name,event in events.items() if name in ['first_softening','first_failure','first_multinode','first_free_component']} for k,events in first_events.items()})
 morphology=[]
 for k,rows in groups.items():morphology.extend(dict(case=k,time_s=float(states[k]['time']),**row) for row in rows)
 csv_write(O/'FRAGMENT_MORPHOLOGY.csv',morphology);csv_write(O/'FRAGMENT_HISTORY.csv',history['V4']);csv_write(O/'ENERGY_LEDGER.csv',[{k:v for k,v in r.items() if k in ['step','time_s','dt_s','W_J','W_drive_J','W_drag_J','U_J','K_J','D_J','Uc_J','Jc_J','R_J','relative_R','numerical_damping_work_J']} for r in history['V4']])
 energy_keys=['step','time_s','dt_s','W_J','U_J','K_J','D_J','Uc_J','Jc_J','R_J','relative_R']
 energy_extrema={name:{k:float(row[k]) for k in energy_keys} for name,row in [('max_relative',max(history['V4'],key=lambda r:float(r['relative_R']))),('max_absolute',max(history['V4'],key=lambda r:abs(float(r['R_J']))))]}
 energy_extrema['normalization']='Absolute R divided by cumulative max of 1e-24 J, abs(W), and abs(delta U)+K+D+Uc; exactly the frozen runner definition.'
 dump(O/'ENERGY_EXTREMA.json',energy_extrema)
 plt.rcParams.update({'font.size':12,'axes.grid':True,'grid.alpha':.4,'grid.color':'#9ba5ae','savefig.facecolor':'white'})
 fig=plt.figure(figsize=(14,13),layout='constrained')
 for j,k in enumerate(['V3','V4']):
  ref=refs[k];s=states[k];b=ref['bonds'];mask=(s['d']==1)&~ref['initial_cut']&(ref['bond_types']==1);ax=fig.add_subplot(3,2,j+1,projection='3d');X=ref['X'][:int(ref['clot_count'])]*1e3;ax.scatter(*X.T,s=1,c='#b6bdc5',alpha=.12,linewidths=0)
  mid=ref['X'][b[mask]].mean(1)*1e3
  if len(mid):ax.scatter(*mid.T,s=1,c='#9c1536',alpha=.3,linewidths=0)
  ax.set(xlim=(-.65,.65),ylim=(-.65,.65),zlim=(-.35,.35),xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)',title=f'{k}: failed bulk-bond midpoints, {float(s["time"])*1e3:.4f} ms');ax.view_init(15,-58);ax.set_box_aspect((1.2,1.2,.6))
  ax=fig.add_subplot(3,2,j+3)
  if mask.any():
   d=(ref['X'][b[mask,1]]-ref['X'][b[mask,0]])/ref['r0'][mask,None];ax.hist(abs(d).ravel(),bins=np.linspace(0,1,31),density=True,color=['#718c9e','#a41c3b'][j],alpha=.8);ax.set_title(f'Axis fourth moment = {cracks[k]["axis_alignment_fourth_moment"]:.3f}')
  else:ax.text(.5,.5,'No bulk failures in accepted window',ha='center',transform=ax.transAxes)
  ax.set(xlabel='|failed-bond direction · axis|',ylabel='Probability density')
 evolution_axes=[fig.add_subplot(3,2,5),fig.add_subplot(3,2,6)]
 for j,k in enumerate(['V3','V4']):
  rows=evolutions[k];t=np.array([r['time_s'] for r in rows])*1e3;color=['#718c9e','#a41c3b'][j]
  evolution_axes[0].plot(t,[r['bulk_failed'] for r in rows],label=k,color=color)
  evolution_axes[1].plot(t,[r['failed_minus_reference_axis_moment'] if r['failed_minus_reference_axis_moment'] is not None else np.nan for r in rows],label=k,color=color)
 evolution_axes[0].set(xlabel='Physical time (ms)',ylabel='New fully failed bulk bonds',title='Accepted-state damage development')
 evolution_axes[1].set(xlabel='Physical time (ms)',ylabel='Failed A4 − original reference A4',title='Axis alignment relative to each reference graph')
 for ax in evolution_axes:ax.legend()
 fig.savefig(out/'FIG_06_CRACK_PATH_COMPARISON.png',dpi=180);plt.close(fig)
 spectra={}
 for j,k in enumerate(['V3','V4']):
  sizes=np.array([q['nodes'] for q in groups[k]]);masses=np.array([q['mass_kg'] for q in groups[k]]);counts=[int(((sizes>=lo)&(sizes<=hi)).sum()) for lo,hi in BINS];mass=[float(masses[(sizes>=lo)&(sizes<=hi)].sum()) for lo,hi in BINS];spectra[k]=dict(node_bins=BIN_LABELS,counts=counts,mass_kg=mass)
 plot_spectra(spectra,tf);dump(O/'FRAGMENT_SIZE_SPECTRA.json',spectra)
 fig,axes=plt.subplots(1,2,figsize=(14,5),layout='constrained')
 for j,k in enumerate(['V3','V4']):
  q=survival_results[k];ax=axes[j];items=q['target_member_sets']
  for i,r in enumerate(items):ax.plot([r['birth_time_s']*1e3,(r['birth_time_s']+r['observed_duration_s'])*1e3],[i,i],color='#a41c3b',lw=2);ax.scatter((r['birth_time_s']+r['observed_duration_s'])*1e3,i,marker='>' if r['right_censored'] else 'x',color='#a41c3b')
  if not items:ax.text(.5,.5,'No observed 20–200-node free member sets',ha='center',transform=ax.transAxes)
  ax.axvline(OFF*1e3,color='#333',ls='--');ax.set(xlim=(0,END*1e3),xlabel='Physical time (ms)',ylabel='Exact target member set',title=k+(' — unload reached' if q['drive_off_reached'] else ' — unload not reached'))
 fig.suptitle('Exact member-set survival; × split, > right-censored at accepted end');fig.savefig(out/'FIG_08_FRAGMENT_SURVIVAL.png',dpi=180);plt.close(fig)
 fig,axes=plt.subplots(2,1,figsize=(13,8),layout='constrained');rows=history['V4'];t=np.array([float(r['time_s']) for r in rows])*1e3
 for key in ['W_J','U_J','K_J','D_J','Uc_J','Jc_J']:
  axes[0].plot(t,np.array([float(r[key]) for r in rows])*1e6,label=key.replace('_J',''),lw=1.5)
 axes[0].set(ylabel='Energy (µJ)',title='V4: W = ΔU + ΔK + D + ΔUc − Jc + R');axes[0].legend(ncol=6)
 axes[1].plot(t,100*np.array([float(r['relative_R']) for r in rows]),color='#a41c3b');axes[1].axhline(1,color='#333',ls='--',label='Frozen 1% gate');axes[1].set(xlabel='Physical time (ms)',ylabel='Absolute normalized residual (%)');axes[1].legend();fig.savefig(out/'FIG_09_COMPLETE_ENERGY_LEDGER.png',dpi=180);plt.close(fig)
 target_counts={k:sum(20<=q['nodes']<=200 for q in groups[k]) for k in refs};target_masses={k:sum(q['mass_kg'] for q in groups[k] if 20<=q['nodes']<=200) for k in refs};geo=json.loads((O/'IRREGULAR_REFERENCE_GEOMETRY.json').read_text());elastic=json.loads((O/'ELASTIC_PATCH_CHECK.json').read_text());performance=json.loads((O/'PERFORMANCE_OPTIMIZATION.json').read_text())
 mass_statistics={k:dict(free_mass_kg=sum(q['mass_kg'] for q in rows),debris12_mass_kg=sum(q['mass_kg'] for q in rows if q['nodes']<=2),free_volume_m3=sum(q['reference_volume_m3'] for q in rows),largest_nodes=max((q['nodes'] for q in rows),default=0),target_mass_fraction=target_masses[k]/refs[k]['V'][:int(refs[k]['clot_count'])].sum()/RHO) for k,rows in groups.items()}
 outcome=dict(status='COMPLETED_SINGLE_EXPLORATORY_SAMPLE' if complete else 'FORMAL_TRAJECTORY_STOPPED_AT_NUMERICAL_GATE',numerical_status=result['status'],full_3p6ms_completed=complete,accepted_time_s=tf,final_target_count=target_counts,target_mass_kg=target_masses,mass_statistics=mass_statistics,history_max_target_count=result['max_N20_200'],unload_survival=survival_results,performance_speedup_damaged_window=performance['comparisons']['damaged']['speedup'],elastic_relative_changes=elastic['relative_modulus_changes'],more_stable_than_V3='NOT_ASSESSABLE' if not complete else 'SEE_EXACT_COHORT_COUNTS_SINGLE_SAMPLE',physical_conclusion='No calibrated clot microstructure or acoustic streaming conclusion; positions, neighbor topology and quadrature all changed.',no_added_cases=True)
 outcome['matched_window_max_target_count']={k:max((int(float(row['N20_200'])) for row in rows if float(row['time_s'])<=tf+1e-14),default=0) for k,rows in history.items()}
 outcome['final_axis_alignment_excess_over_reference']={k:rows[-1]['failed_minus_reference_axis_moment'] for k,rows in evolutions.items()}
 outcome['endpoint_sample_times_s']={k:float(s['time']) for k,s in states.items()}
 dump(O/'SCIENTIFIC_OUTCOME.json',outcome);dump(O/'COMPUTE_COST.json',{k:result[k] for k in ['wall_s','cpu_s','peak_RSS_bytes','accepted_steps','rejected_predictions','substeps','accepted_time_s']})
 fig,ax=plt.subplots(figsize=(13,7));ax.axis('off');lines=[('Numerical trajectory',f'{result["status"]} — {tf*1e3:.6f} / 3.6 ms'),('Elastic discretization',f'Bulk unchanged; shear changes {100*min(elastic["relative_modulus_changes"][mode] for mode in ["xy","xz","yz"]):.2f}% to {100*max(elastic["relative_modulus_changes"][mode] for mode in ["xy","xz","yz"]):.2f}%'),('Reference topology',f'{geo["original_bonds"]:,} → {geo["new_reference_bonds"]:,} original bonds'),('20–200 node free blocks',f'Matched-time V3: {target_counts["V3"]}; V4: {target_counts["V4"]}'),('V4 trajectory cost',f'{result["wall_s"]/3600:.3f} wall h; {result["peak_RSS_bytes"]/2**30:.2f} GiB peak RSS'),('Exact V3 optimization',f'{performance["comparisons"]["damaged"]["speedup"]:.3f}× damaged short window; no full-run speedup claim'),('Scientific scope','One conditioned seed; exploratory, uncalibrated material'),('Load scope','Prescribed mechanical shear surrogate; no solved acoustic streaming')]
 ax.table(cellText=lines,colLabels=['Question / evidence','Measured result / limitation'],colWidths=[.28,.72],cellLoc='left',loc='center').auto_set_font_size(False)
 for tab in ax.tables:
  tab.set_fontsize(12);tab.scale(1,2.2)
 fig.savefig(out/'FIG_10_SCIENTIFIC_OUTCOME.png',dpi=180,bbox_inches='tight');plt.close(fig)
 print('POSTPROCESS_ANALYSIS_COMPLETE',json.dumps({k:v for k,v in outcome.items() if k!='unload_survival'},default=cv),flush=True)
if __name__=='__main__':
 import argparse
 parser=argparse.ArgumentParser();parser.add_argument('--spectrum-only',action='store_true');args=parser.parse_args()
 if args.spectrum_only:
  plt.rcParams.update({'font.size':12,'axes.grid':True,'grid.alpha':.4,'grid.color':'#9ba5ae','savefig.facecolor':'white'})
  plot_spectra(json.loads((O/'FRAGMENT_SIZE_SPECTRA.json').read_text()),json.loads((O/'NATIVE_V4_RESULT.json').read_text())['accepted_time_s'])
 else:main()
