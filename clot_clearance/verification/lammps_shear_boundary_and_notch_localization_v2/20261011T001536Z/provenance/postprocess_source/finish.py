"""Analyze the two completed cases; no solver invocation or outcome-selected refits."""
import sys,json,csv,hashlib,os,time
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from audit import shear,notch,save,js,fracture_area,regions
sys.path.insert(0,str(Path(__file__).parents[1]/'material_macro_audit_v1'))
from postprocess import figsave
BLUE='#236c9d';RED='#b43c37';TEAL='#16818b';GOLD='#d7a52f'

def main(old,base,figs):
 figs.mkdir(parents=True,exist_ok=True)
 assert all(js(base/n/'EXECUTION.json')['status']=='COMPLETED' for n in ['M2P','M5B'])
 original_s=shear(old/'M2',base/'baseline/M2');original_n=notch(old/'M5',base/'baseline/M5')
 new_s=shear(base/'M2P',base/'analysis/M2P');new_n=notch(base/'M5B',base/'analysis/M5B')
 previous=js(old/'ELASTIC_MAPPING.json')['measured'];infer=previous['E_Pa']/(2*(1+previous['nu']));moment=js(base/'DISCRETE_SHEAR_MOMENT.json')
 h={n:np.genfromtxt(p/'history.csv',delimiter=',',names=True) for n,p in [('M2',old/'M2'),('M2P',base/'M2P'),('M5',old/'M5'),('M5B',base/'M5B')]}
 native=np.genfromtxt(base/'M2P/boundary_audit.csv',names=True,delimiter=',');a=h['M2P'];sel=(native['time_s']>=new_s['fit_time_s'][0])&(native['time_s']<=new_s['fit_time_s'][1]);bal=np.max(abs(native['top_Rx_N'][sel]+native['bottom_Rx_N'][sel])/abs(native['top_Rx_N'][sel]))
 new_s['actual_native_fit_window_grip_force_balance_relative_max']=float(bal)
 sa=dict(baseline=original_s,periodic=new_s,M1_inferred_G_Pa=infer,M1_inferred_vs_M2_difference=original_s['G_ROI_Pa']/infer-1,M1_inferred_vs_M2P_difference=new_s['G_ROI_Pa']/infer-1,discrete_moment=moment,M2P_vs_discrete_prediction=new_s['G_ROI_Pa']/moment['predicted_full_neighborhood_small_strain_G_Pa']-1,interpretation='Free-surface/grip-induced heterogeneous field and global-traction/local-strain mismatch are major causes. Original force conversion, area and symmetric shear-strain formula reproduced. Periodic residual to940 assessed against frozen cubic-stencil fourth moment, not eliminated by parameter tuning.',internal_consistency='M1/M2/M3P values cannot be assembled as one isotropic calibrated material tensor; finite boundaries, measurement fields and cubic quadrature differ.')
 save(base/'SHEAR_BOUNDARY_AUDIT.json',sa)
 # Late imposed collar damage is a protocol limitation, independent of the outcome.
 limit=js(base/'KINEMATIC_VALIDITY_LIMIT.json');p=base/'M5B';g=dict(np.load(p/'geometry.npz'));X=g['X'];b=g['bonds'];ev=np.load(p/'break_events.npz')['events'];mask=g['prescribed_mask'][:,2];collar=mask[b].all(1);times=ev[:,0]*.5e-6;cid=ev[:,1]
 early_candidates=[]
 for statepath in p.glob('frame_*.npz'):
  z=np.load(statepath)
  if float(z['time_s'])<limit['threshold_time_s']:early_candidates.append((float(z['time_s']),statepath))
 t,ep=max(early_candidates);z=np.load(ep);dead=~z['alive'].astype(bool);reg,_,_,band=regions(X,b);cross=((X[b[:,0],2]<-1e-12)!=(X[b[:,1],2]<-1e-12));area=fracture_area(p,g,z,ev)
 matrix=coo_matrix((np.ones(2*np.sum(~dead)),(np.r_[b[~dead,0],b[~dead,1]],np.r_[b[~dead,1],b[~dead,0]])),shape=(len(X),len(X)));nc,comp=connected_components(matrix,directed=False);sizes=np.sort(np.bincount(comp))[::-1]
 early=dict(actual_state_file=str(ep),time_s=t,selection='Latest predeclared1ms accepted snapshot strictly before analytically forced collar-break threshold; no interpolated state or new stopping time',broken=int(dead.sum()),band_fraction=float(band[dead].mean()),grip_fraction=float((reg[dead]==1).mean()),off_projected_plane_fraction=float((~cross[dead]).mean()),area_attempt=area,connected_component_sizes=sizes,component_count=nc)
 save(base/'analysis/M5B/EARLY_AREA_ATTEMPT.json',early)
 n=dict(baseline=original_n,improved=new_n,kinematic_protocol_limit=limit,early_before_forced_collar=early,collar_bonds_initial=int(collar.sum()),collar_bonds_broken_end=int(collar[cid].sum()),collar_bonds_broken_before_limit=int(np.sum(collar[cid]&(times<limit['threshold_time_s']))),collar_bonds_broken_after_limit=int(np.sum(collar[cid]&(times>=limit['threshold_time_s']))),fracture_work_status='GC_UNAVAILABLE',reason='Initiation can improve while finite-band fragmentation, uncertain new surface and imposed late collar stretching still prevent reliable fracture work per crack area. Signed topology energy is retained without clamping negative increments.')
 save(base/'NOTCH_FAILURE_DIAGNOSTIC.json',n)
 rows=[];checks={}
 for name,p in [('M2',old/'M2'),('M2P',base/'M2P'),('M5',old/'M5'),('M5B',base/'M5B')]:
  a=h[name];d=js(p/'DIAGNOSTIC_CHECKS.json');peak=int(np.argmax(a['stress_Pa']));row=dict(case=name,source=str(p),G_ROI_Pa=sa['baseline' if name=='M2' else 'periodic']['G_ROI_Pa'] if name in ['M2','M2P'] else None,G_nominal_Pa=sa['baseline' if name=='M2' else 'periodic']['G_nominal_Pa'] if name in ['M2','M2P'] else None,peak_nominal_stress_Pa=float(a['stress_Pa'][peak]),peak_time_s=float(a['time_s'][peak]),peak_opening_m=float(a['displacement_m'][peak]),peak_observed=bool(peak<len(a)-1 and a['stress_Pa'][-1]<.9*a['stress_Pa'][peak]) if name.startswith('M5') else None,first_break_s=float(a['time_s'][np.flatnonzero(a['broken_bulk']>0)[0]]) if a['broken_bulk'][-1] else None,first_region=(n['baseline' if name=='M5' else 'improved']['first_event'][0]['region'] if name.startswith('M5') else None),broken_end=int(a['broken_bulk'][-1]),net_topology_energy_J=float(a['D_break_independent_J'][-1]),positive_topology_energy_J=float(np.maximum(np.diff(a['D_break_independent_J']),0).sum()),negative_topology_energy_J=float(np.minimum(np.diff(a['D_break_independent_J']),0).sum()),negative_topology_events=len(d['negative_topology_releases']),energy_residual_work_fraction=float(max(abs(a['R_energy_J']))/max(abs(a['W_ext_J']))),wall_seconds=d['wall_seconds'],Gc_J_m2=None)
  if name.startswith('M5'):
   audit=n['baseline' if name=='M5' else 'improved'];row.update(peak_target_band_fraction=audit['snapshots']['peak']['target_band_fraction'],end_off_plane_fraction=audit['snapshots']['end']['off_projected_plane_fraction'],KE_U_immediately_before_first=audit['first_break_KE_U_before'],peak_KE_U=float(a['KE_over_U'][peak]));assert int(np.load(p/'peak.npz')['step'])==peak
  rows.append(row)
  if name in ['M2P','M5B']:
   assert not list(p.glob('*ERROR*'));assert len(a)==js(p/'test_plan.json')['steps']+1;assert np.isfinite(np.column_stack([a[k] for k in a.dtype.names])).all();assert np.all(np.diff(a['broken_bulk'])>=0)
   fc=d['force_checks'];force_ok=len(fc)==3 and max(v['native_force_relative_error'] for v in fc)<1e-9;fd_ok=len(fc)==3 and max(v['finite_difference_relative_error'] for v in fc)<1e-5
   assert force_ok and fd_ok
   checks[name]=dict(finite=True,actual_time_index=bool(np.max(abs(a['time_s']-a['step']*.5e-6))<1e-14),native_force_gradient_match=force_ok,finite_difference_match=fd_ok,energy_residual_below5percent=row['energy_residual_work_fraction']<=.05,energy_residual=row['energy_residual_work_fraction'],topology_exact_at_saved_states=True,unintended_breaks_absent=bool(a['broken_bulk'][-1]==0) if name=='M2P' else 'Fracture intended; late collar rupture separately identified')
   idx=np.unique(np.r_[np.arange(0,len(a),20),len(a)-1,peak,np.flatnonzero(np.diff(a['broken_bulk'])>0)+1]);np.savez_compressed(p/'review_curves.npz',records=a[idx],selection='Every20 steps plus every topology event, exact maximum and endpoint; full accepted history.csv remains server-side')
 fields=list(dict.fromkeys(k for row in rows for k in row))
 with (base/'CASE_COMPARISON.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)
 save(base/'IMPLEMENTATION_CHECKS.json',checks)
 plot_shear(base,figs,h,sa);plot_notch(old,base,figs,h,n);write_report(old,base,figs,sa,n,rows,checks)

def plot_shear(base,figs,h,s):
 fig,axs=plt.subplots(2,2,figsize=(12,9),layout='constrained')
 for name,color,sub in [('M2',BLUE,'baseline/M2'),('M2P',TEAL,'analysis/M2P')]:
  a=h[name];axs[0,0].plot(a['gamma_xz']*100,a['stress_Pa'],color=color,label=f'{name}: nominal stress / ROI shear');axs[0,1].plot(a['boundary_strain']*100,a['stress_Pa'],color=color,label=f'{name}: nominal strain')
  z=np.load(base/sub/'peak_fields.npz');pr=z['profile'];axs[1,0].plot(pr[:,1]*1e6,pr[:,0]*1e3,color=color,label=name+' mean ux')
  axs[1,1].plot(pr[:,4]*100,pr[:,0]*1e3,color=color,label=name+' dux/dZ');axs[1,1].plot(pr[:,5]*100,pr[:,0]*1e3,'--',color=color,label=name+' duz/dX')
 x=np.array([0,max(h['M2P']['gamma_xz'].max(),h['M2']['gamma_xz'].max())]);axs[0,0].plot(x*100,940*x,'k--',label='Continuum input G=940 Pa');axs[0,0].plot(x*100,s['discrete_moment']['predicted_full_neighborhood_small_strain_G_Pa']*x,':',color=GOLD,label='Discrete full-stencil prediction996.03 Pa')
 axs[0,0].set(xlabel='Central fitted symmetric shear (%)',ylabel='Nominal shear traction (Pa)',title=f'M2 {s["baseline"]["G_ROI_Pa"]:.2f} Pa; M2P {s["periodic"]["G_ROI_Pa"]:.2f} Pa');axs[0,0].legend(fontsize=8)
 axs[0,1].set(xlabel='Grip displacement / 850µm (%)',ylabel='Nominal shear traction (Pa)',title='Different strain definition; not local G');axs[0,1].legend()
 axs[1,0].set(xlabel='Mean accepted displacement ux (µm)',ylabel='Reference height Z (mm)',title='Accepted peaks: M2 9.9855 ms; M2P 9.9760 ms');axs[1,0].legend()
 axs[1,1].set(xlabel='Local displacement gradient (%)',ylabel='Reference height Z (mm)',title='Shear includes both gradient components');axs[1,1].legend(fontsize=8)
 for ax in axs[1]:ax.axhspan(-.1995,.1995,color='purple',alpha=.07,label='ROI height extent')
 fig.suptitle('A1 | Fixed 600×600×1000µm LPS; only lateral boundary differs; unchanged20–70% loading fit')
 figsave(fig,figs,'A1_shear_definitions_and_profiles')
 fig=plt.figure(figsize=(14,5.4),layout='constrained')
 for j,(name,sub) in enumerate([('M2','baseline/M2'),('M2P','analysis/M2P')]):
  z=np.load(base/sub/'peak_fields.npz');X=z['X']*1e3;roi=z['ROI'];boundary=z['prescribed_mask'].any(1);ax=fig.add_subplot(1,3,j+1,projection='3d',computed_zorder=False);sc=ax.scatter(*X.T,c=z['gamma']*100,cmap='viridis',s=6,alpha=.7,vmin=0,vmax=1.2,depthshade=False);ax.scatter(*X[roi].T,s=9,facecolors='none',edgecolors='purple',linewidths=.35);ax.view_init(17,-55);ax.set(xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)',title=f'{name} actual shear, t={float(z["time_s"])*1e3:.4f} ms',xticks=[-.3,0,.3],yticks=[-.3,0,.3]);ax.set_box_aspect([.6,.6,1]);ax.tick_params(labelsize=8)
  bounds=np.array([.1495,.1495,.1995]);corners=np.array([[xx,yy,zz] for xx in [-1,1] for yy in [-1,1] for zz in [-1,1]])*bounds
  edges=[corners[[ii,jj]] for ii in range(8) for jj in range(ii+1,8) if np.sum(corners[ii]!=corners[jj])==1];ax.add_collection3d(Line3DCollection(edges,colors='purple',linewidths=1.1,zorder=12))
  for zz in [-.35,.35]:ax.plot([-.3,.3,.3,-.3,-.3],[-.3,-.3,.3,.3,-.3],[zz]*5,c='#236c9d',ls='--',lw=1,zorder=11)
 ax=fig.add_subplot(1,3,3)
 for name,color in [('M2',BLUE),('M2P',TEAL)]:a=h[name];ax.plot(a['time_s']*1e3,a['KE_over_U']*100,color=color,label=name)
 ax.set(xlabel='Physical time (ms)',ylabel='Free-DOF kinetic / bulk energy (%)',title='Full startup/unloading disclosed');ax.set_yscale('symlog',linthresh=.1);ax.axhline(5,c='black',ls='--');ax.legend();fig.colorbar(sc,ax=fig.axes[:2],shrink=.65,label='Local symmetric shear (%)');fig.suptitle('A2 | Real node positions in reference space; purple box=measurement ROI; blue planes=grip transitions (x/z prescribed)')
 figsave(fig,figs,'A2_shear_ROI_and_inertia',False)

def plot_notch(old,base,figs,h,n):
 fig=plt.figure(figsize=(13,6));fig.subplots_adjust(top=.78,bottom=.12,wspace=.2)
 for j,(name,sub) in enumerate([('M5','baseline/M5'),('M5B','analysis/M5B')]):
  z=np.load(base/sub/'first_break_risk.npz');X=z['X'];b=z['bonds'];alive=z['alive'];util=z['utilization'];eligible=alive.copy().astype(bool);eligible[z['first_bond_indices']]=True;node=np.zeros(len(X));np.maximum.at(node,b[eligible].ravel(),np.repeat(util[eligible],2));ax=fig.add_subplot(1,2,j+1,projection='3d');sc=ax.scatter(*(X*1e3).T,c=node,vmin=0,vmax=1.05,cmap='magma',s=9,depthshade=False,alpha=.9)
  cut=(abs(X[b[:,0],1])<55e-6)&(abs(X[b[:,1],1])<55e-6)&alive.astype(bool);ax.add_collection3d(Line3DCollection(X[b[cut]]*1e3,linewidths=.25,colors='#b8b4b4',alpha=.25))
  mids=X[b[z['first_bond_indices']]].mean(1)*1e3;ax.scatter(*mids.T,s=80,c='cyan',marker='*',edgecolors='black',linewidths=.5,depthshade=False)
  ax.plot([-.1,-.1],[-.3,.3],[0,0],c='limegreen',lw=2,label='Original notch front')
  for zz in [-.35,.35]:ax.plot([-.3,.3],[0,0],[zz,zz],ls='--',c='#237ac0',lw=2)
  ax.set(xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)',title=f'{name} first-break accepted state\nt={float(z["time_s"])*1e3:.4f} ms',xticks=[-.3,0,.3],yticks=[-.3,0,.3]);ax.set_box_aspect([.6,.6,1]);ax.view_init(17,-55);ax.tick_params(labelsize=9)
 fig.colorbar(sc,ax=fig.axes,shrink=.7,pad=.06,label='Max incident live/new-break stretch / s00')
 fig.suptitle('B1 | Actual point lattice and surviving central-slab bonds; reference coordinates\nCyan stars=first native breaks; green=notch front; blue=grip transitions',fontsize=13,y=.98);figsave(fig,figs,'B1_notch_vs_grip_first_break',False)
 fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
 for col,(name,sub) in enumerate([('M5','baseline/M5'),('M5B','analysis/M5B')]):
  z=np.load(base/sub/'EVENT_LOCALIZATION.npz');v=z['cumulative'];ax=axs[0,col]
  for k,label,color in [(4,'Notch-tip region',TEAL),(5,'Grip transition',BLUE),(6,'Other surface',GOLD),(7,'Free interior',RED)]:ax.plot(v[:,0]*1e3,v[:,k],label=label,c=color)
  ax.set(xlabel='Physical time (ms)',ylabel='Cumulative native broken bonds',title=name+' fixed spatial classification');ax.legend(fontsize=8)
  ids=z['events'][:,1];mid=z['midpoint_reference_m'][ids];ax=axs[1,col];step=z['events'][:,0]*.5e-6
  ax.scatter(step*1e3,mid[:,2]*1e3,s=2,c=np.where(z['target_band'][ids],0,1),cmap='coolwarm',alpha=.25);ax.axhspan(-.07525,.07525,color=TEAL,alpha=.08);ax.set(xlabel='Actual native break time (ms)',ylabel='Reference bond midpoint Z (mm)',title='All native events; shaded target-plane band')
  if name=='M5B':
   for row in axs:row[col].axvline(n['kinematic_protocol_limit']['threshold_time_s']*1e3,c='black',ls='--',lw=1)
 fig.suptitle('B2 | Late dashed threshold: affine collar vertical bonds necessarily exceed s00; not material Gc')
 figsave(fig,figs,'B2_native_event_localization')
 fig,axs=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
 for name,color in [('M5',BLUE),('M5B',TEAL)]:
  a=h[name];axs[0].plot(a['displacement_m']*1e6,a['reaction_N']*1e6,label=name+' resultant top reaction',c=color)
  for key,line in [('notch_tip','-'),('grip_transition','--')]:
   if name=='M5':
    q=n['baseline']['snapshots']['first_break'];axs[1].scatter([q['time_s']*1e3],[q['region_max_utilization'][key]],c=color,marker='x' if key=='notch_tip' else 's',label='M5 first-state '+key)
 reg=np.genfromtxt(base/'M5B/regional_utilization.csv',delimiter=',',names=True)
 for key,color in [('notch',TEAL),('grip',BLUE)]:axs[1].plot(reg['time_s']*1e3,reg[key+'_including_newbreak_max'],label='M5B '+key,c=color,lw=.9)
 axs[1].axhline(1,c='black',ls=':');axs[1].set(xlim=(0,7.2),xlabel='Physical time (ms)',ylabel='Max current eligible stretch / s00',title='No synthetic baseline pre-break field');axs[1].legend(fontsize=7)
 q=np.genfromtxt(base/'M5B/boundary_audit.csv',delimiter=',',names=True);disp=h['M5B']['displacement_m'][q['step'].astype(int)];axs[0].plot(disp*1e6,q['generalized_external_Q_N']*1e6,c=GOLD,ls='--',lw=1,label='M5B generalized load Q');axs[0].set(xlabel='Actual mean grip opening (µm)',ylabel='Force (µN)',title='Strained collars: Q differs from net reaction');axs[0].legend(fontsize=8)
 vals=[n['baseline']['snapshots']['peak']['target_band_fraction'],n['improved']['snapshots']['peak']['target_band_fraction'],n['early_before_forced_collar']['band_fraction']];axs[2].bar(['M5 peak','M5B peak','M5B 6ms'],np.array(vals)*100,color=[BLUE,TEAL,GOLD]);axs[2].axhline(70,c='black',ls='--');axs[2].set(ylabel='Broken bonds in frozen target band (%)',title='Localization criterion is not crack area',ylim=(0,105))
 figsave(fig,figs,'C1_loading_and_localization_comparison')
 # Actual first, peak and final accepted states; same limits and physical aspect.
 fig=plt.figure(figsize=(15,11));fig.subplots_adjust(left=.04,right=.96,bottom=.06,top=.90,hspace=.25,wspace=.08);states={}
 for name,p in [('M5',old/'M5'),('M5B',base/'M5B')]:
  for label in ['first_break','peak','end']:states[name,label]=dict(np.load(p/(label+'.npz')))
 xx=np.concatenate([z['x_m'] for z in states.values()])*1e3;lo=xx.min(0);hi=xx.max(0);pad=.04*(hi-lo);lo-=pad;hi+=pad
 for row,(name,p) in enumerate([('M5',old/'M5'),('M5B',base/'M5B')]):
  g=np.load(p/'geometry.npz');X=g['X'];b=g['bonds']
  for col,label in enumerate(['first_break','peak','end']):
   z=states[name,label];x=z['x_m']*1e3;alive=z['alive'].astype(bool);ax=fig.add_subplot(2,3,row*3+col+1,projection='3d');ax.scatter(*x.T,c=RED,s=5,depthshade=False,alpha=.75);cut=alive&(abs(X[b[:,0],1])<55e-6)&(abs(X[b[:,1],1])<55e-6);ax.add_collection3d(Line3DCollection(x[b[cut]],colors='#aaa19c',linewidths=.25,alpha=.3));ax.view_init(16,-55);ax.set(xlabel='x (mm)',ylabel='y (mm)',zlabel='z (mm)',title=f'{name} {label.replace("_"," ")}\nt={float(z["time_s"])*1e3:.4f} ms',xlim=(lo[0],hi[0]),ylim=(lo[1],hi[1]),zlim=(lo[2],hi[2]),xticks=[-.3,0,.3],yticks=[-.3,0,.3]);ax.set_box_aspect(hi-lo);ax.tick_params(labelsize=8)
 fig.suptitle('C2 | Actual accepted states, unchanged node geometry; surviving bonds; common physical axes\nM5B final collar damage includes externally prescribed stretching and cannot certify a fracture surface',y=.98,fontsize=13);figsave(fig,figs,'C2_actual_fracture_states',False)
 fig,axs=plt.subplots(1,2,figsize=(13,5),layout='constrained')
 for ax,name in zip(axs,['M5','M5B']):
  a=h[name]
  for key,label in [('W_ext_J','Actual boundary work'),('U_bulk_J','Recoverable bulk U'),('U_contact_J','Contact U'),('K_free_J','Free-DOF K'),('D_num_J','Numerical damping'),('D_break_independent_J','Signed topology transfer'),('R_energy_J','Independent residual')]:ax.plot(a['time_s']*1e3,a[key]*1e9,label=label,lw=1)
  if name=='M5B':ax.axvline(n['kinematic_protocol_limit']['threshold_time_s']*1e3,c='black',ls='--',lw=1)
  ax.set(xlabel='Physical time (ms)',ylabel='Energy (nJ)',title=name+' | Gc unavailable');ax.legend(fontsize=8)
 fig.suptitle('D | Unchanged source-consistent LPS energy; negative topology increments retained')
 figsave(fig,figs,'D_independent_energy_and_limits')

def write_report(old,base,figs,s,n,rows,checks):
 # Report generation is kept separate from the scientific computations.
 save(base/'REPORT_INPUTS.json',dict(shear=s,notch=n,rows=rows,checks=checks))
 print(json.dumps(dict(shear_G=[s['baseline']['G_ROI_Pa'],s['periodic']['G_ROI_Pa']],M2P_Gnom=s['periodic']['G_nominal_Pa'],M2P_vs_stencil=s['M2P_vs_discrete_prediction'],first_M5B=n['improved']['first_event'],early=n['early_before_forced_collar'],rows=rows,checks=checks),default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x.item() if isinstance(x,np.generic) else x,indent=2))

if __name__=='__main__':main(Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]))
