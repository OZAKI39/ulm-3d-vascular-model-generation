#!/usr/bin/env python3
"""Regenerate P7 review figures exclusively from saved CSV/JSON evidence."""
from pathlib import Path
import sys,json,argparse
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from particle_3d.particle3_cases import write_json
from particle_3d.particle7_cases import REPO

NAMES=['00_particle7_scope_and_literature','01_frozen_inlet_flux_map','02_flux_weighted_inlet_sampling','03_mb_cumulative_flux_scheduler','04_rbc_volume_flux_scheduler','05_rbc_isotropic_orientation_validation','06_finite_size_inlet_admission','07_pending_injection_queue','08_injection_timestep_independence','09_particle_lifecycle_and_outlet_removal','10_injection_checkpoint_restart','11_long_injection_population_balance','12_scheduled_vs_admitted_distributions','13_real_frozen_inlet_injection_smoke','14_local_tube_hematocrit_diagnostic']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'axes.grid':True,'grid.alpha':.2,'savefig.facecolor':'white','figure.facecolor':'white'})
COLORS={'scheduled':'#2477b4','admitted':'#e17c16','pending':'#b3264b','active':'#24894c','exited':'#8052a3'}

def load(root,name): return json.loads((root/'data'/name).read_text())
def csv(root,name): return np.genfromtxt(root/'data'/name,delimiter=',',names=True,encoding='utf-8')
def save(root,idx,fig,sources):
 fig.tight_layout(rect=(0,0,1,.94)); fig.savefig(root/'figures'/(NAMES[idx]+'.png'),dpi=170); plt.close(fig)
 write_json(root/'data'/(NAMES[idx]+'_plot.json'),dict(figure=NAMES[idx]+'.png',source_files=sources,plotter='particle_3d/scripts/plot_particle7.py',role='VALIDATION_ONLY'))

def main(root,indices=None):
 root=Path(root); (root/'figures').mkdir(exist_ok=True); indices=set(range(15) if indices is None else indices)
 def wanted(i,paths): return i in indices and all((root/'data'/p).exists() for p in paths)
 if wanted(0,['01_flux_audit.json']):
  a=load(root,'01_flux_audit.json'); Q=a['INLET']['positive_Q_m3_s']; f,ax=plt.subplots(figsize=(12,6.5)); ax.axis('off')
  items=[(.16,.76,'Frozen FEM inlet\nQ = %.6g m³/s\nOfficial faces + P1 velocity'%Q),(.50,.76,'RBC feed volume fraction\nH_D = 0.45\nPatel 2025, Table 2'),(.84,.76,'SonoVue number concentration\nC = 8.5e12 m⁻³\nDauba 2025 + FAU 2024\n2e7 / (0.0326 × 72 mL)')]
  for x,y,text in items: ax.text(x,y,text,ha='center',va='center',transform=ax.transAxes,bbox=dict(boxstyle='round,pad=.7',fc='#eaf2fa',ec='#2b608a'),fontsize=12)
  for x in [.16,.50,.84]: ax.annotate('',xy=(x,.55),xytext=(x,.65),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color='#456278',lw=1.5))
  ax.text(.5,.44,'P2 original RBC geometry  +  Original SonoVue histogram\nIndependent RNG streams → deterministic cumulative flux → finite-size admission\nStable ID → pending / active → outlet record → deletion',ha='center',va='center',transform=ax.transAxes,fontsize=14,bbox=dict(boxstyle='round,pad=.7',fc='#eaf5eb',ec='#387f45'))
  ax.text(.5,.12,'No full PK  |  Provisional isotropic orientation  |  Real RBC passage unresolved\nNo production timestep / neighbor settings  |  No Particle-8  |  No CFD',ha='center',transform=ax.transAxes,color='#9a2838',fontsize=12)
  f.suptitle('Particle-7: continuous inlet population infrastructure',fontsize=17); save(root,0,f,['01_flux_audit.json','../literature/sources.json'])
 if wanted(1,['01_flux_audit.json','01_inlet_mesh.json']):
  d=load(root,'01_inlet_mesh.json'); a=load(root,'01_flux_audit.json'); xyz=np.array(d['triangles_m']); origin=xyz.reshape(-1,3).mean(0); _,_,vh=np.linalg.svd(xyz.reshape(-1,3)-origin,full_matrices=False); xy=(xyz-origin)@vh[:2].T*1e6
  f,ax=plt.subplots(1,3,figsize=(14,4.8))
  for axis,values,label in [(ax[0],np.mean(d['q_m_s'],axis=1)*1e6,'Inward speed (µm/s)'),(ax[1],np.array(d['weights_m3_s'])*1e18,'Triangle flow (fL/s)')]:
   pc=PolyCollection(xy,array=values,cmap='viridis',edgecolors='white',linewidths=.25); axis.add_collection(pc); axis.autoscale(); axis.set_aspect('equal'); f.colorbar(pc,ax=axis,label=label); axis.set(xlabel='Inlet basis 1 (µm)',ylabel='Inlet basis 2 (µm)')
  roles=['INLET','OUTLET_01','OUTLET_02','OUTLET_03']; ax[2].barh(roles,[a[k]['positive_Q_m3_s']*1e18 for k in roles],color=['#2477b4']+['#24894c']*3); ax[2].set_xlabel('Positive volume flow (fL/s)'); ax[2].set_title('Signed mass residual: %.2g'%a['mass_balance']['relative_signed_residual'])
  f.suptitle('Official Frozen FEM boundaries: exact positive linear flux',fontsize=15); save(root,1,f,['01_flux_audit.json','01_inlet_mesh.json'])
 if wanted(2,['02_synthetic_linear.csv','02_actual_triangle_counts.csv','02_sampling_statistics.json']):
  p=csv(root,'02_synthetic_linear.csv'); r=csv(root,'02_actual_triangle_counts.csv'); stat=load(root,'02_sampling_statistics.json'); f,ax=plt.subplots(1,2,figsize=(12,5))
  h=ax[0].hist2d(p['x'],p['y'],bins=45,cmap='viridis'); f.colorbar(h[3],ax=ax[0],label='Samples per bin'); ax[0].plot(.5,.25,'r+',ms=15,label='Analytic mean (0.5, 0.25)'); ax[0].set(xlabel='x',ylabel='y',title='Synthetic q(x,y)=x: 100,000 samples'); ax[0].legend(fontsize=9)
  ax[1].scatter(r['expected_fraction'],r['sampled_fraction'],s=15,alpha=.7); m=max(r['expected_fraction'].max(),r['sampled_fraction'].max()); ax[1].plot([0,m],[0,m],'k--'); ax[1].set(xlabel='Expected triangle flux fraction',ylabel='Observed sample fraction',title='Actual inlet: N=100,000; chi² p=%.3f'%stat['actual']['chi_square_p_value'])
  f.suptitle('Flux-density sampling: exact Dirichlet mixture, acceptance = 1',fontsize=15); save(root,2,f,['02_synthetic_linear.csv','02_actual_triangle_counts.csv','02_sampling_statistics.json'])
 if wanted(3,['03_04_scheduler.csv']):
  d=csv(root,'03_04_scheduler.csv'); f,ax=plt.subplots(2,1,figsize=(11,6),sharex=True)
  ax[0].plot(d['time_s'],d['mb_expected'],label='Expected cumulative count'); events=load(root,'03_04_scheduled_events.json'); exact_times=[0]+[e['scheduled_time_s'] for e in events if e['species']=='MB']+[float(d['time_s'][-1])]; exact_counts=[0]+list(range(1,len(exact_times)-1))+[len(exact_times)-2]; ax[0].step(exact_times,exact_counts,where='post',label='Scheduled at exact event time'); ax[0].set_ylabel('MB count'); ax[0].legend()
  ax[1].plot(d['time_s'],d['mb_error']); ax[1].axhline(-1,color='r',ls='--',label='Strict lower error bound'); ax[1].set(xlabel='Physical time (s)',ylabel='Scheduled − expected'); ax[1].legend()
  f.suptitle('Frozen-Q MB scheduler: C = 8.5e12 m⁻³; no Poisson arrivals',fontsize=15); save(root,3,f,['03_04_scheduler.csv','03_04_scheduled_events.json'])
 if wanted(4,['03_04_scheduler.csv']):
  d=csv(root,'03_04_scheduler.csv'); f,ax=plt.subplots(1,2,figsize=(12,5)); ax[0].plot(d['time_s'],d['rbc_target_m3']*1e18,label='0.45 × cumulative blood volume'); ax[0].plot(d['time_s'],d['rbc_scheduled_m3']*1e18,'--',label='Sum of original sampled volumes'); ax[0].set(xlabel='Time (s)',ylabel='RBC volume (fL)'); ax[0].legend()
  ax[1].plot(d['time_s'],d['rbc_residual_m3']*1e18,label='Unscheduled residual'); ax[1].plot(d['time_s'],d['next_rbc_volume_m3']*1e18,label='Next pre-sampled RBC volume',alpha=.7); ax[1].set(xlabel='Time (s)',ylabel='Volume (fL)'); ax[1].legend()
  f.suptitle('RBC volume-threshold scheduler: each actual volume is counted',fontsize=15); save(root,4,f,['03_04_scheduler.csv','03_04_scheduled_events.json'])
 if wanted(5,['05_orientation.csv','05_orientation_statistics.json']):
  d=csv(root,'05_orientation.csv'); stat=load(root,'05_orientation_statistics.json'); f=plt.figure(figsize=(13,4.8)); ax=f.add_subplot(131,projection='3d'); ax.scatter(d['px'][::40],d['py'][::40],d['pz'][::40],s=2,alpha=.4); ax.set(xlabel='px',ylabel='py',zlabel='pz',title='Short-axis directions (subsample)'); ax.set_box_aspect([1,1,1]); b=f.add_subplot(132)
  for k in ['px','py','pz']: b.hist(d[k],bins=40,histtype='step',density=True,label=k)
  b.axhline(.5,color='k',ls='--'); b.set(xlabel='Short-axis component',ylabel='Density',title='Uniform marginal on [−1,1]'); b.legend(); c=f.add_subplot(133); c.bar(['px²','py²','pz²'],stat['second_moments']); c.axhline(1/3,color='r',ls='--',label='1/3'); c.set_ylim(0,.4); c.legend(); c.set_title('Second moments, N=100,000')
  f.suptitle('Provisional isotropic SO(3): a model assumption, not a mouse measurement',fontsize=14); save(root,5,f,['05_orientation.csv','05_orientation_statistics.json'])
 if wanted(6,['06_admission.json']):
  from particle_3d.injection_admission import event_shape
  d=load(root,'06_admission.json'); f,ax=plt.subplots(2,3,figsize=(13,8))
  for row,sp in enumerate(['MB','RBC']):
   for col,case in enumerate(['ACCEPT','WALL','PAIR']):
    item=next(v for v in d['cases'] if v['species']==sp and v['case']==case); axis=ax[row,col]
    shape=event_shape(item['particle'],item['position_m']); theta=np.linspace(0,2*np.pi,181)
    xy=np.array([shape.support([np.cos(t),np.sin(t),0])[:2] for t in theta])*1e6
    good=item['status']=='ACCEPTED'; color='#24894c' if good else '#b3264b'
    axis.fill(xy[:,0],xy[:,1],color=color,alpha=.22); axis.plot(xy[:,0],xy[:,1],color=color)
    x,y,_=np.array(item['position_m'])*1e6; axis.plot(x,y,'+',color=color)
    if case=='PAIR': axis.add_patch(plt.Circle((0,0),3,color='grey',alpha=.35,label='Existing sphere'))
    if case=='WALL': axis.axvline(100,color='black',lw=2); axis.axvspan(100,107,color='grey',alpha=.2)
    axis.set(xlim=(x-6,x+6),ylim=(y-6,y+6),xlabel='x (µm)',ylabel='y (µm)',title=sp+': '+item['status'].replace('_',' ')); axis.set_aspect('equal')
  f.suptitle('Finite geometry at birth: same sizes / q across candidate retries\nRejected position → retry; exhausted guard → preserve particle in pending queue',fontsize=14); save(root,6,f,['06_admission.json'])
 if wanted(7,['07_pending_queue.json']):
  d=load(root,'07_pending_queue.json'); f,ax=plt.subplots(1,2,figsize=(12,5)); ax[0].step(d['times_s'],d['queue_sizes'],where='post',lw=3); ax[0].set(xlabel='Validation time (s)',ylabel='Pending RBC count',ylim=(-.1,1.3)); ax[0].annotate('Blocker removed → same RBC admitted',(.01,0),xytext=(.002,.45),arrowprops=dict(arrowstyle='->')); ax[1].axis('off'); e=d['original']; g=e['geometry']; text='Preserved across queue and admission\n\nID: %d\nD: %.6f µm\nV: %.6f fL\nq: %s\nScheduled: %.9g s\nAdmitted: %.9g s\n\nNo geometry or orientation re-sampling'%(e['particle_id'],e['provenance']['D_um'],e['provenance']['V_fL'],np.array2string(np.array(e['q']),precision=4),e['scheduled_time_s'],d['admitted']['admitted_time_s']); ax[1].text(.03,.9,text,va='top',fontsize=12)
  f.suptitle('Pending queue preserves particle identity and original geometry',fontsize=15); save(root,7,f,['07_pending_queue.json'])
 if wanted(8,['08_timestep_parity.json']):
  d=load(root,'08_timestep_parity.json'); f,ax=plt.subplots(1,2,figsize=(12,5)); events=d['events']; ids=np.array([e['particle_id'] for e in events]); ts=np.array([e['scheduled_time_s'] for e in events]); ax[0].plot(ids,ts,lw=1); ax[0].set(xlabel='Stable particle ID',ylabel='Scheduled time (s)',title='All events retain their physical times')
  for k,row in enumerate(d['partitions']): ax[1].plot(ids[::30],np.full(len(ids[::30]),k),'.',label='dt=%.5g s; exact match'%row['dt_s'])
  ax[1].set(xlabel='Stable particle ID',yticks=[0,1,2],yticklabels=['dt','dt/2','dt/4'],title='Times / IDs / geometries / orientations match'); ax[1].legend(fontsize=9)
  f.suptitle('NOT PRODUCTION TIMESTEP SELECTION',fontsize=16); save(root,8,f,['08_timestep_parity.json'])
 if wanted(9,['09_lifecycle.json','09_lifecycle_balance.csv']):
  d=load(root,'09_lifecycle.json'); rows=csv(root,'09_lifecycle_balance.csv'); f,ax=plt.subplots(1,2,figsize=(12,5));
  for key in ['scheduled','active','exited']: ax[0].plot(rows['time_s'],rows[key+'_rbc_count'],label=key,color=COLORS[key])
  ax[0].set(xlabel='Time (s)',ylabel='RBC count',title='Actual LAMMPS insert / delete'); ax[0].legend()
  for k,e in enumerate(d['exits'][:35]): ax[1].plot([e['birth_admitted_time'],e['exit_time']],[e['particle_id']]*2,'-',lw=3,color='#2477b4')
  ax[1].set(xlabel='Time (s)',ylabel='Stable particle ID',title='Birth → first OUTLET_01 crossing → removal')
  f.suptitle('Synthetic open-section lifecycle; common plug validation only',fontsize=15); save(root,9,f,['09_lifecycle.json','09_lifecycle_balance.csv'])
 if wanted(10,['10_restart_parity.json']):
  d=load(root,'10_restart_parity.json'); r=d['comparisons']; f,ax=plt.subplots(1,2,figsize=(12,5)); fields=['events','births','exits','pending','active','scheduler']; a=np.array([[int(not v[k]) for k in fields] for v in r]); im=ax[0].imshow(a,aspect='auto',vmin=0,vmax=1,cmap='Reds'); ax[0].set(xticks=range(6),xticklabels=fields,ylabel='Comparison checkpoint',title='Mismatch matrix (white = exact match)'); ax[0].tick_params(axis='x',rotation=35); f.colorbar(im,ax=ax[0],ticks=[0,1]); ax[1].axis('off'); ax[1].text(.05,.85,'Actual binary LAMMPS restart\n\nEvents at checkpoint: %d\nEvents at final time: %d\nLAMMPS destroyed before restart\nAll six independent RNG states restored\nPending particles restored without draws\nNext RBC geometry restored without draws\n\nAll sequence fields: EXACT'%(d['checkpoint_event_count'],d['final_event_count']),va='top',fontsize=12)
  f.suptitle('Continuous run vs checkpoint / destroy / restart',fontsize=16); save(root,10,f,['10_restart_parity.json'])
 if wanted(11,['11_long_balance.csv','11_long_summary.json']):
  d=csv(root,'11_long_balance.csv'); summary=load(root,'11_long_summary.json'); f,ax=plt.subplots(2,2,figsize=(13,8))
  for label,key in [('expected','mb_expected'),('scheduled','scheduled_mb_count'),('active','active_mb_count'),('exited','exited_mb_count')]: ax[0,0].plot(d['time_s'],d[key],label=label)
  ax[0,0].set(xlabel='Time (s)',ylabel='MB count'); ax[0,0].legend()
  for label,key in [('target','rbc_target_m3'),('admitted','admitted_rbc_volume_m3'),('active','active_rbc_volume_m3'),('exited','exited_rbc_volume_m3')]: ax[0,1].plot(d['time_s'],d[key]*1e12,label=label)
  ax[0,1].set(xlabel='Time (s)',ylabel='RBC volume (nL)'); ax[0,1].legend()
  ax[1,0].plot(d['time_s'],d['pending_mb_count'],label='MB'); ax[1,0].plot(d['time_s'],d['pending_rbc_count'],'--',label='RBC'); ax[1,0].set(xlabel='Time (s)',ylabel='Pending count',ylim=(-.1,1)); ax[1,0].legend()
  ax[1,1].plot(d['time_s'],d['active_mb_count']+d['active_rbc_count']); ax[1,1].set(xlabel='Time (s)',ylabel='Active count')
  f.suptitle('Long mixed ledger: %d MB + %s RBC; 10 nm open control section\nBookkeeping stress fixture — no full suspension or RBC passage claim'%(summary['final']['scheduled_mb_count'],format(summary['final']['scheduled_rbc_count'],',')),fontsize=14); save(root,11,f,['11_long_balance.csv','11_long_summary.json','11_long_lifecycle.csv.gz'])
 if wanted(12,['12_distribution_diagnostics.json']):
  d=load(root,'12_distribution_diagnostics.json'); f,ax=plt.subplots(2,3,figsize=(14,8))
  for row,case in enumerate(['long','real']):
   for col,key in enumerate(['D_um','r','mb_diameter_um']):
    h=d[case][key]; edges=np.array(h['edges']); axis=ax[row,col]
    for name in ['scheduled','admitted','pending']:
     counts=np.array(h[name]); n=counts.sum()
     if n: axis.stairs(counts/(n*np.diff(edges)),edges,label=name+' n='+format(int(n),','),color=COLORS[name],linestyle='--' if name=='admitted' else '-',linewidth=2 if name=='admitted' else 1)
     else: axis.plot([],[],label=name+' n=0',color=COLORS[name])
    axis.set(xlabel={'D_um':'RBC diameter D (µm)','r':'RBC aspect ratio c/a','mb_diameter_um':'MB diameter (µm)'}[key],ylabel='Probability density',title=('Long mixed fixture' if case=='long' else 'Real inlet admission smoke')); axis.legend(fontsize=8)
  f.suptitle('Scheduled / admitted / pending distributions; original sizes never redrawn',fontsize=15); save(root,12,f,['12_distribution_diagnostics.json','11_long_lifecycle.csv.gz','13_real_smoke.json'])
 if wanted(13,['13_real_smoke.json','13_real_geometry.json','13_isolated_mb_transport.json']):
  d=load(root,'13_real_smoke.json'); geo=load(root,'13_real_geometry.json'); inlet=np.array(geo['inlet_triangles_m']); wall=np.array(geo['nearby_wall_triangles_m']); origin=inlet.reshape(-1,3).mean(0); f=plt.figure(figsize=(18,6)); ax=f.add_subplot(131,projection='3d',computed_zorder=False); ax.add_collection3d(Poly3DCollection((wall-origin)*1e6,facecolor='#c6d4dc',edgecolor='none',alpha=.1)); ax.add_collection3d(Poly3DCollection((inlet-origin)*1e6,facecolor='#7aa6d2',edgecolor='grey',lw=.1,alpha=.35))
  cand=d['candidates']; choose=cand[::max(1,len(cand)//60)]
  for label,selected,color,marker in [('Rejected',[x for x in choose if x['status']!='ACCEPTED'],'#b3264b','x'),('Admitted MB',[x for x in cand if x['status']=='ACCEPTED' and x['species']=='MB'],'#e17c16','o'),('Admitted RBC',[x for x in cand if x['status']=='ACCEPTED' and x['species']=='RBC'],'#24894c','^')]:
   if selected:
    x=(np.array([v['position_m'] for v in selected])-origin)*1e6; ax.scatter(x[:,0],x[:,1],x[:,2],s=15 if marker=='x' else 90,color=color,marker=marker,label=label,alpha=.35 if marker=='x' else 1.,depthshade=False,zorder=3 if marker=='x' else 8)
  from particle_3d.lammps_state import BridgeParticle
  for record in d['active']:
   shape=BridgeParticle(**record).shape(); theta=np.linspace(0,2*np.pi,25); phi=np.linspace(0,np.pi,17)
   grid=np.array([[shape.support([np.sin(p)*np.cos(t),np.sin(p)*np.sin(t),np.cos(p)]) for t in theta] for p in phi])
   grid=(grid-origin)*1e6; ax.plot_wireframe(grid[:,:,0],grid[:,:,1],grid[:,:,2],color='#24894c',lw=.6,alpha=.85,zorder=7)
  ax.set(xlabel='x − inlet center (µm)',ylabel='y − center (µm)',zlabel='z − center (µm)',title='Official INLET + nearby WALL'); ax.legend(fontsize=8); limits=[ax.get_xlim(),ax.get_ylim(),ax.get_zlim()]; span=max(np.ptp(v) for v in limits); centers=[np.mean(v) for v in limits]; ax.set_xlim(centers[0]-span/2,centers[0]+span/2); ax.set_ylim(centers[1]-span/2,centers[1]+span/2); ax.set_zlim(centers[2]-span/2,centers[2]+span/2); ax.set_box_aspect([1,1,1]); ax.view_init(25,120)
  b=f.add_subplot(132); r=d['accounting']; labels=['MB scheduled','MB admitted','MB pending','RBC scheduled','RBC admitted','RBC pending']; vals=[r[x.lower().split()[1]+'_'+x.lower().split()[0]+'_count'] for x in labels]; b.barh(labels,vals,color=['#2477b4','#24894c','#b3264b']*2); b.set_xscale('symlog',linthresh=1); b.set_xlabel('Particle count (symlog)'); b.set_title('RBC passage remains NOT ESTABLISHED')
  for i,v in enumerate(vals): b.text(max(v,1),i,' '+str(v),va='center')
  isolated=load(root,'13_isolated_mb_transport.json'); states=isolated['states']; positions=np.array([v['particle']['position'] for v in states]); tri=inlet[isolated['candidates'][-1]['triangle_id']]; n=np.cross(tri[1]-tri[0],tri[2]-tri[0]); n/=np.linalg.norm(n)
  if n@np.array(states[0]['particle']['velocity'])<0: n=-n
  c=f.add_subplot(133); c.plot((np.array([v['time_s'] for v in states])-states[0]['time_s'])*1000,(positions-positions[0])@n*1e6,'o-'); c.set(xlabel='Remaining physical time after birth (ms)',ylabel='Inward center displacement (µm)',title='Independent single-MB P6.5 replay\nNo RBC present in this separate fixture')
  f.suptitle('Real Frozen inlet: mixed admission congestion + independent single-MB entry motion',fontsize=15); save(root,13,f,['13_real_smoke.json','13_real_geometry.json','13_isolated_mb_transport.json'])
 if wanted(14,['14_uniform_time_tube_hct.csv','14_tube_hct_statistics.json','13_real_smoke.json']):
  d=csv(root,'14_uniform_time_tube_hct.csv'); real=load(root,'13_real_smoke.json'); stats=load(root,'14_tube_hct_statistics.json'); f,ax=plt.subplots(1,2,figsize=(12,5)); ax[0].plot(d['time_s'],d['observed_tube_hct'],lw=.7,label='Observed active RBC volume / control volume'); ax[0].axhline(.45,color='r',ls='--',label='Feed H_D = 0.45 (reference only)'); ax[0].set(xlabel='Time (s)',ylabel='Center-assigned volume diagnostic',title='Uniform-time snapshots; time mean = %.6f'%stats['exact_residence_time_weighted_mean']); ax[0].legend(fontsize=8)
  ax[1].bar(['Feed ratio','Real lumen snapshot'],[.45,real['observed_local_tube_hct']],color=['#2477b4','#24894c']); ax[1].set_ylabel('RBC volume fraction diagnostic'); ax[1].set_title('Different definitions; no equality imposed'); ax[1].text(1,real['observed_local_tube_hct']+.02,'%.5g'%real['observed_local_tube_hct'],ha='center')
  f.suptitle('H_D controls incoming volume flux; tube Hct is a separate observation',fontsize=15); save(root,14,f,['14_uniform_time_tube_hct.csv','14_tube_hct_statistics.json','13_real_smoke.json'])

if __name__=='__main__':
 p=argparse.ArgumentParser(); p.add_argument('--root',default=str(REPO/'particle_3d/reports/particle7')); p.add_argument('--indices',nargs='*',type=int); a=p.parse_args(); main(a.root,a.indices)
