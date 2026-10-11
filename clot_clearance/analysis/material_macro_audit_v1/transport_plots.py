"""Actual-step / saved-state distinction, physical-time weighted comparisons."""
import sys,json
from pathlib import Path
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.sparse.csgraph import connected_components
sys.path.insert(0,str(Path(__file__).parent))
from postprocess import save,js,figsave
R=1.78e-6

def measure(path):
 with h5py.File(path/'trajectory.h5') as h:
  t=h['time_s'][:];x=h['positions_m'][:,0];ids=h['particle_id'][:];dt=h['trajectory/physical_dt_s'][:];quality=h['quality/raw_F1_diagonal'][:,0];r=h['radii_m'][:]
 n=len(ids);delta=x[:,:,None,:]-x[:,None,:,:];distance=np.linalg.norm(delta,axis=-1);distance[:,np.arange(n),np.arange(n)]=np.inf;gap=distance-2*R;near=(gap<=.02*R).any(2).mean(1);aggregation=[];largest=[]
 for g in gap:
  count,lab=connected_components(g<=.1*R,directed=False);size=np.bincount(lab);largest.append(size.max()/n);aggregation.append(np.mean(size[lab]>=2))
 largest=np.array(largest);aggregation=np.array(aggregation);weights=np.diff(t)
 weighted=lambda a:float(np.sum(.5*(a[:-1]+a[1:])*weights)/(t[-1]-t[0]))
 cum=np.r_[0,np.cumsum(dt)];ix=np.searchsorted(cum,t);ix=np.clip(ix,0,len(cum)-1);clockerr=np.minimum(abs(cum[ix]-t),abs(cum[np.maximum(ix-1,0)]-t)).max()
 d=dict(saved_frames=len(t),accepted_steps=len(dt),dt_quantiles_s=np.quantile(dt,[0,.5,.9,.99,1]),saved_interval_quantiles_s=np.quantile(np.diff(t),[0,.5,1]),max_saved_clock_vs_accumulated_accepted_dt_error_s=clockerr,total_accepted_time_s=dt.sum(),initial_zero_dt_rows=int((dt==0).sum()),timeweighted_largest_fraction=weighted(largest),timeweighted_aggregation_fraction=weighted(aggregation),timeweighted_nearcontact_fraction=weighted(near),weighting='Trapezoidal integration on identical actual 0.5ms saved states; no interpolation; substep-dense-region bias avoided; nearcontact inter-save behavior is not fully sampled in baseline',minimum_saved_gap_over_R=float(gap.min()/R),summary=js(path/'SUMMARY.json'))
 return d,dict(t=t,x=x,dt=dt,quality=quality,largest=largest,aggregation=aggregation,near=near,ids=ids,r=r)

def run(base,figs):
 mb=base/'microbubble';baseline=Path(js(mb/'test_plan.json')['baseline']);b,bdata=measure(baseline);h,hdata=measure(mb/'dt_half');stats=np.load(mb/'dt_half/accepted_step_audit.npz')['values']
 assert np.array_equal(bdata['ids'],hdata['ids']) and np.array_equal(bdata['r'],hdata['r']) and np.array_equal(bdata['x'][0],hdata['x'][0]);assert np.max(abs(bdata['t']-hdata['t']))<1e-13
 c0=js(baseline/'CONFIG.json');c1=js(mb/'dt_half/CONFIG.json');changes={k:[c0[k],c1[k]] for k in c0 if c0[k]!=c1[k]};assert list(changes)==['dt_s']
 diffs={k:abs(h[k]-b[k]) for k in ['timeweighted_largest_fraction','timeweighted_aggregation_fraction','timeweighted_nearcontact_fraction']}
 inst=js(mb/'INSTRUMENTATION.json');maxpos=float(np.linalg.norm(hdata['x']-bdata['x'],axis=-1).max());comdiff=float(np.linalg.norm(hdata['x'].mean(1)-bdata['x'].mean(1),axis=-1).max())
 result=dict(baseline=b,half_step=h,actual_config_changes=changes,initial_ids_radii_positions_identical=True,core_normalized_differences=diffs,cluster_metric_screen_pass=all(v<=.05 for v in diffs.values()),maximum_same_ID_saved_position_difference_over_R=maxpos/R,maximum_saved_COM_difference_over_R=comdiff/R,half_step_instrumentation=inst,half_timeweighted_all_accepted_step_nearcontact_fraction=float(np.sum(stats[:,7]*hdata['dt'])/sum(hdata['dt'])),half_actual_displacement_over_R_quantiles=np.quantile(stats[:,0],[0,.5,.9,.99,1]),half_start_prediction_over_R_quantiles=np.quantile(stats[:,1],[0,.5,.9,.99,1]),half_endpoint_projection_over_R_quantiles=np.quantile(stats[:,2],[0,.5,.9,.99,1]),half_raw_positive_gap_closing_s_quantiles=np.quantile(stats[:,5][stats[:,5]<1e8],[0,.5,.9,.99,1]),positive_gap_interpretation='No denominator floor; near machine-zero positive gaps produce min closing time3.4e-14s and very large dt/closing ratios. Projection enforces a geometric budget; this does not establish temporally resolved lubrication/contact physics.',baseline_actual_displacement_per_step='UNAVAILABLE; baseline diagnostic field max_physical_substep_displacement is actually pre-projection prediction, verified source kernels.py predict(); never relabel it accepted displacement',RK_stage_scope='NOT_APPLICABLE: this actual F6/F7 branch is Euler prediction +12 contact passes, not old RK tracker. Per-pass and endpoint-net projection explicitly distinct.',stage_model_scope='Mature initial cluster persists across 0–1s; no de-novo cluster-formation or whole6s claim',physics_validity='UNCALIBRATED_NEAR_FIELD_EXPLORATORY; independent of timestep screen',not_run='No full baseline rerun, no max_displacement_per_radius sweep: this actual backend has no such cap',files_missing_or_schema_replaced='F7 uses CONFIG/IDENTITY/SUMMARY plus trajectory.h5; old particles/events.jsonl/PHYSICS_BANNER names not present. Equivalent actual metadata audited, no invented files.')
 # Remove only rigid pose for a diagnostic of internal rearrangement; never alter saved states.
 aligned=[];pairdiff=[];edge_disagreement=[]
 for xb,xh in zip(bdata['x'],hdata['x']):
  xb=xb-xb.mean(0);xh=xh-xh.mean(0);u,_,vt=np.linalg.svd(xh.T@xb);rotation=u@vt
  if np.linalg.det(rotation)<0:u[:,-1]*=-1;rotation=u@vt
  aligned.append(float(np.sqrt(np.mean(np.sum((xh@rotation-xb)**2,axis=1)))/R))
  db=np.linalg.norm(xb[:,None]-xb[None,:],axis=-1);dh=np.linalg.norm(xh[:,None]-xh[None,:],axis=-1);tri=np.triu_indices(len(xb),1)
  pairdiff.append(float(np.sqrt(np.mean((db[tri]-dh[tri])**2))/R));eb=db[tri]<=2.1*R;eh=dh[tri]<=2.1*R;edge_disagreement.append(float(np.count_nonzero(eb!=eh)/max(np.count_nonzero(eb|eh),1)))
 result['max_rigid_pose_aligned_same_ID_RMS_over_R']=max(aligned)
 result['max_pair_distance_RMS_difference_over_R']=max(pairdiff)
 result['max_cluster_edge_Jaccard_disagreement']=max(edge_disagreement)
 result['internal_trajectory_timestep_status']='NOT_ESTABLISHED: multi-radius same-ID divergence despite stable cluster size/COM; no certified cluster AI labels'
 result['registered_geometric_gate_pass']=all(s['summary']['quality']['error_bits']==0 and s['summary']['quality']['max_accepted_penetration_m']/R<=.01 and s['summary']['quality']['max_contact_pass_correction_m']/R<=.1 for s in [b,h])
 result['allowed_AI_labels']='NO_CERTIFIED_CLUSTER_AI_LABELS; retain audit observables as exploratory records only'
 save(mb/'MICROBUBBLE_AUDIT.json',result)
 np.savez_compressed(mb/'COMPARISON_SERIES.npz',time_s=bdata['t'],largest_baseline=bdata['largest'],largest_half=hdata['largest'],near_baseline=bdata['near'],near_half=hdata['near'],baseline_quality=bdata['quality'],half_quality=hdata['quality'],half_step_time_s=np.cumsum(hdata['dt']),half_step_audit=stats)
 fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained');ax=axs[0,0]
 for d,label in [(bdata,'Baseline'),(hdata,'dt/2')]:ax.plot(np.cumsum(d['dt']),d['dt']*1e6,label=label,lw=1)
 ax.set(xlabel='Physical time (s)',ylabel='Actual accepted dt (µs)',title='10 → 5 µs; saved interval = 500 µs');ax.legend()
 ax=axs[0,1]
 for d,label in [(bdata,'Baseline'),(hdata,'dt/2')]:ax.plot(d['t'],d['quality'][:,:,6].max(1)/R,label=label+' cumulative max/pass')
 display=np.unique(np.r_[np.arange(0,len(stats),100),len(stats)-1,np.argmax(stats[:,:3],axis=0),np.argmax(stats[:,6])]);ax.plot(np.cumsum(hdata['dt'])[display],stats[display,2],alpha=.6,lw=.7,label='dt/2 endpoint net (stride100 + maxima)');ax.set(xlabel='Physical time (s)',ylabel='Projection / R',title='Per-pass maximum ≠ endpoint net');ax.legend(fontsize=8)
 ax=axs[1,0]
 for d,label in [(bdata,'Baseline'),(hdata,'dt/2')]:ax.plot(d['t'],d['largest'],label=label+' largest',lw=1.2);ax.plot(d['t'],d['near'],'--',label=label+' near-contact',lw=1)
 ax.set(xlabel='Physical time (s)',ylabel='Particle fraction',ylim=(-.03,1.08),title='Time weighted, matched actual saved times');ax.legend(fontsize=8)
 ax=axs[1,1];t=np.cumsum(hdata['dt']);ax.plot(t[display],stats[display,0],label='Actual accepted displacement');ax.plot(t[display],stats[display,1],alpha=.6,label='Step-start prediction');ax.set(xlabel='Physical time (s)',ylabel='Maximum displacement / R',title='dt/2 audit; display stride100 + maxima');ax.legend(fontsize=8)
 fig.suptitle('F7 replay M0, 32 bubbles; fixed model / controls / initial IDs; near-field EXPLORATORY');figsave(fig,figs,'06_microbubble_timestep')
 fig,axs=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
 same=np.linalg.norm(hdata['x']-bdata['x'],axis=-1).max(1)/R
 axs[0].plot(bdata['t'],same,label='Maximum same-ID separation / R');axs[0].plot(bdata['t'],aligned,label='RMS after rigid-pose alignment / R');axs[0].plot(bdata['t'],pairdiff,label='Pair-distance RMS difference / R');axs[0].set(xlabel='Physical time (s)',ylabel='Difference / R',title='Internal trajectories are sensitive');axs[0].legend(fontsize=8)
 axs[1].plot(bdata['t'],edge_disagreement,color='#b33132');axs[1].set(xlabel='Physical time (s)',ylabel='Contact-neighbor edge Jaccard disagreement',ylim=(0,1),title='Cluster size alone hides rearrangement')
 tm=np.cumsum(hdata['dt']);axs[2].plot(tm[display],stats[display,6],lw=.7,color='#8b4693',label='Stride100 + true maximum');ix=int(np.argmax(stats[:,6]));axs[2].scatter([tm[ix]],[stats[ix,6]],s=22,c='crimson',label=f'Max = {stats[ix,6]:.3g}');axs[2].legend(fontsize=8);axs[2].set_yscale('log');axs[2].set(xlabel='Physical time (s)',ylabel='dt / minimum positive gap-closing time',title='Raw positive gaps; no denominator floor')
 fig.suptitle('No certified cluster AI labels: finite geometric tolerance does not establish internal time convergence')
 figsave(fig,figs,'09_microbubble_internal_sensitivity')

 y=js(base/'production_audit/Y_AUDIT.json');native=np.genfromtxt(base/'production_audit/Y_NATIVE_SPEED_TIMELINE.csv',delimiter=',',names=True,dtype=None,encoding='utf8');tr=np.load(base/'production_audit/Y_TRANSPORT_SPEED_TIMELINE.npz')
 fig,axs=plt.subplots(1,3,figsize=(15,4),layout='constrained')
 for source,label in [('saved_frame_force_audit','All nodes: saved-record maxima'),('detached_force_audit','Detached: each native step')]:
  sel=native['source']==source;axs[0].plot(native['time_s'][sel]*1e3,native['max_node_speed_m_s'][sel],label=label,lw=.9)
 axs[0].set(xlabel='Native time (ms)',ylabel='Full-step speed (m/s)',title='Native 0–2.33 ms');axs[0].legend(fontsize=8)
 axs[1].scatter(tr['relative_time_s'][::100],tr['speed_m_s'][::100],s=2,alpha=.2,color='#286f96');axs[1].scatter([y['transport']['max_accepted_endpoint_COM_speed']['start_relative_time_s']],[y['transport']['max_accepted_endpoint_COM_speed']['m_s']],s=35,color='crimson',marker='x',label='Maximum including handoff');axs[1].legend(fontsize=8);axs[1].set(xlabel='Time after handoff (s)',ylabel='Accepted COM speed (m/s)',title='232,276 actual steps; display stride100')
 axs[2].hist(tr['dt_s']*1e3,bins=60,color='#328692');axs[2].set(xlabel='Transport accepted dt (ms)',ylabel='Step count',title='Maximum 5 ms; displacement ≤20 µm');axs[2].set_yscale('log');figsave(fig,figs,'07_Y_native_transport_audit')
 save(base/'TIMESTEP_AUDIT.json',dict(material=js(base/'MATERIAL_COMPARISONS.json') if (base/'MATERIAL_COMPARISONS.json').exists() else 'PENDING',microbubble=result,production=y,production_crosschecks=js(base/'production_audit/Y_CROSSCHECKS.json')))
 print(json.dumps({k:v for k,v in result.items() if k in ['core_normalized_differences','cluster_metric_screen_pass','registered_geometric_gate_pass','maximum_same_ID_saved_position_difference_over_R','maximum_saved_COM_difference_over_R']},indent=2))
if __name__=='__main__':
 b=Path(sys.argv[1]);f=Path(sys.argv[2]);f.mkdir(parents=True,exist_ok=True);run(b,f)
