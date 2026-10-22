"""Publication figures from computed fields/trajectories, with explicit model scope."""
import time
import numpy as np,h5py,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .common import *
from .admission import time_windows
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelcolor':'#202a35','text.color':'#202a35','xtick.color':'#202a35','ytick.color':'#202a35','axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':220,'pdf.fonttype':42})
COL=['#356fa4','#b75a51','#5a9474','#9270a7'];OUT=BASE/'figures'
def save(fig,i,name,evidence):
 fig.suptitle(f'F10AG | {name}',fontsize=15);fig.text(.5,.015,'EXPLORATORY · Homogeneous-water proxy · Zero background flow · Uncalibrated physics',ha='center',fontsize=8,color='#4d5965');fig.tight_layout(rect=[0,.04,1,.94]);stem=f'Fig_F10AG_{i:02d}_{name.lower().replace(" ","_").replace("/","_")}'
 for ext in ['png','pdf']:fig.savefig(OUT/(stem+'.'+ext),facecolor='white')
 plt.close(fig);return dict(index=i,name=name,png=stem+'.png',pdf=stem+'.pdf',evidence=evidence)
def series(kind,id):return np.load(BASE/'series'/(kind+'_'+id+'.npz'))
def xyz(ax,x,title):
 ax.scatter(*(x.T*1e6),s=8,color=COL[0],depthshade=True);ax.set(xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)',title=title);ax.set_xlim(-550,550);ax.set_ylim(-450,450);ax.set_zlim(-80,80);ax.set_zticks([-60,60]);ax.zaxis.set_tick_params(labelsize=8);ax.set_box_aspect([1100,900,160]);ax.view_init(24,-60)
def main():
 begin=time.perf_counter();OUT.mkdir(exist_ok=True);s=read(BASE/'SCIENTIFIC_RESULTS.json');rows=s['rows'];hf=[r for r in rows if r['kind']=='HF'];good=[r for r in hf if r['status']=='PASS'];st=read(BASE/'STATIC_RESULTS.json');sing=read(BASE/'SINGLE_RESULTS.json');rot=read(BASE/'ROTATION_STATIC.json');manifest=[]
 rootrow=next(v for v in st['rows'] if v['frequency_Hz']==2500000 and v['phase_deg']==0);roots=np.array([r['position_m'][0] for r in rootrow['roots']]);reference=next(r for r in good if r['id']=='baseline_long50_X_PAIR')
 f,ax=plt.subplots(1,2,figsize=(12,5));old=read(AF/'SCIENTIFIC_RESULTS.json')['rows'];old=[r for r in old if r['id'].startswith('long50_') and not r['id'].endswith(('DT05','FINE'))];ax[0].barh([r['id'].replace('long50_','') for r in old],[r['max_fixed_ID_displacement_m']*1e6 for r in old],color=COL);ax[0].set(xlabel='Maximum displacement at 50 ms (µm)',title='Frozen F10AF reference')
 ax[1].axis('off');ax[1].text(.03,.9,'Historical replay: bitwise equal\n\nX-restoring spacing: 304.487 µm\n\nNew questions:\n• Full 3D stability and actual capture\n• Time to enter and remain in a band\n• Full F1/F2/contact vs guarded RK45\n• Tangential force vs actual rotation',va='top',fontsize=12);manifest.append(save(f,1,'Historical baseline and scope',['BASELINE_RESULTS.json','f10af/SCIENTIFIC_RESULTS.json']))
 f=plt.figure(figsize=(12,5));a=f.add_subplot(121,projection='3d');b=f.add_subplot(122)
 with h5py.File(BASE/'FIELD_3D_2500000.h5') as h:q=h['grid/positions_m'][:][::15];v=h['grid/drift_m_s'][:][::15];a.quiver(*(q.T*1e6),*(v.T*1e6),length=.4,normalize=False,color=COL[0],linewidth=.6);a.set(xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)',title='Instantaneous F1 / Stokes drift vectors');a.set_zticks([-60,60]);a.zaxis.set_tick_params(labelsize=8);a.scatter(roots*1e6,np.zeros(len(roots)),np.zeros(len(roots)),s=25,color=COL[1],marker='*',label='Stable points');a.legend(fontsize=8,loc='upper left');a.set_box_aspect([1100,900,140]);q=h['yz_center/positions_m'][:];v=h['yz_center/drift_m_s'][:];b.quiver(q[::10,1]*1e6,q[::10,2]*1e6,v[::10,1]*1e6,v[::10,2]*1e6,color=COL[0]);b.set(xlabel='Y (µm)',ylabel='Z (µm)',title='Y/Z drift at central X-restoring location');b.set_aspect('equal')
 manifest.append(save(f,2,'Three-dimensional acoustic drift',['FIELD_3D_2500000.h5']))
 f,ax=plt.subplots(1,2,figsize=(12,5));J=np.array(rootrow['roots'][1]['jacobians_s1'][-1]);im=ax[0].imshow(J,cmap='RdBu_r',vmin=-3,vmax=3);ax[0].set(xticks=range(3),xticklabels=['X','Y','Z'],yticks=range(3),yticklabels=['vX','vY','vZ'],title='Central drift Jacobian (s⁻¹)');f.colorbar(im,ax=ax[0]);
 for i in range(3):
  for j in range(3):ax[0].text(j,i,f'{J[i,j]:.4g}',ha='center',va='center',color='white' if abs(J[i,j])>1.5 else '#222222')
 for j,r in enumerate(rootrow['roots']):ax[1].plot([1,2,3],np.sort(-np.array(r['eigenvalues_real_s1'])),'o-',label=f'X = {r["position_m"][0]*1e6:.1f} µm')
 ax[1].set(yscale='log',xticks=[1,2,3],xticklabels=['Slow','Slow','Fast'],ylabel='Restoring rate −Re(λ) (s⁻¹)',title='Stable points; unequal restoring rates');ax[1].legend();manifest.append(save(f,3,'Local stability and Jacobian',['STATIC_RESULTS.json']))
 with h5py.File(BASE/'single/atlas.h5') as h:at=h['time_s'][:];xx=h['positions_m'][:]
 f=plt.figure(figsize=(14,5));a=f.add_subplot(131,projection='3d');b=f.add_subplot(132,projection='3d');c=f.add_subplot(133)
 outcomes=sing['rows'][0]['capture']['20'];colors=[COL[r['final_band']] if not r['right_censored'] else '#9ca5ad' for r in outcomes]
 for axx,q,title in [(a,xx[0],'Initial locations; final capture outcome'),(b,xx[-1],'F1-only positions at 250 ms')]:
  axx.scatter(*(q.T*1e6),s=7,c=colors);axx.set(xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)',title=title,xlim=(-550,550),ylim=(-450,450),zlim=(-80,80));axx.set_zticks([-60,60]);axx.zaxis.set_tick_params(labelsize=8);axx.set_box_aspect([1100,900,160]);axx.view_init(24,-60)
 with h5py.File(BASE/'single/local_3d_point.h5') as h:
  qt=h['time_s'][:];qx=h['positions_m'][:];distance=np.sqrt((qx[:,:,0,None]-roots)**2+np.sum(qx[:,:,1:]**2,axis=2)[:,:,None]).min(2)*1e6
 c.plot(qt*1000,distance.mean(1),color=COL[0],label='Mean Euclidean distance');c.fill_between(qt*1000,distance.min(1),distance.max(1),color=COL[0],alpha=.2);c.axhline(20,color=COL[1],ls='--',label='20 µm point region');c.set(xlabel='Physical time (ms)',ylabel='Distance to 3D stable point (µm)',title='24 prescribed local 3D starts');c.legend(fontsize=8)
 f.text(.05,.045,'3D panels: color = held X-band assignment; grey = right-censored. Point capture is evaluated separately.',fontsize=8)
 manifest.append(save(f,4,'Three-dimensional capture atlas',['single/atlas.h5','SINGLE_RESULTS.json','POINT_CAPTURE_RESULTS.json']))
 f,ax=plt.subplots(1,2,figsize=(12,5));atlas=sing['rows'][0]
 for j,th in enumerate(['10','20','30']):
  rr=atlas['capture'][th];d=np.array([r['initial_distance_m']*1e6 for r in rr]);tt=np.array([r['first_held_entry_s'] if r['first_held_entry_s'] is not None else r['completed_s'] for r in rr])*1000;ok=np.array([not r['right_censored'] for r in rr]);ax[0].scatter(d[ok]+j,tt[ok],s=12,color=COL[j],label=f'{th} µm + 10 ms hold');ax[0].scatter(d[~ok]+j,tt[~ok],s=20,marker='v',facecolors='none',edgecolors=COL[j]);ax[1].plot([15,45,100],[sum(not r['right_censored'] for r in rr if abs(r['initial_distance_m']*1e6-d0)<.01)/sum(abs(r['initial_distance_m']*1e6-d0)<.01 for r in rr) for d0 in [15,45,100]],'o-',color=COL[j],label=th+' µm')
 ax[0].set(xlabel='Initial band distance (µm; display offsets ≤2 µm)',ylabel='Held-entry time / censor time (ms)',title='Open triangles: right-censored at 250 ms');ax[0].legend(fontsize=8);ax[1].set(xlabel='Initial band distance (µm)',ylabel='Held-capture fraction',ylim=(-.05,1.05),title='Single-bubble F1 diagnostic');ax[1].legend();manifest.append(save(f,5,'Time to capture and censoring',['SINGLE_RESULTS.json']))
 f,ax=plt.subplots(1,2,figsize=(12,5))
 for j,freq in enumerate([2400000,2500000,2600000]):
  rr=[r for r in st['rows'] if r['frequency_Hz']==freq]
  for r in rr:ax[0].scatter([r['phase_deg']]*len(r['roots']),[z['position_m'][0]*1e6 for z in r['roots']],color=COL[j],s=20,label=f'{freq/1e6:.1f} MHz' if r['phase_deg']==0 else None)
  rr=[r for r in sing['rows'] if r['id'].startswith(f'freq_{freq}_')];ax[1].plot([r['phase_deg'] for r in rr],[np.mean([q['final_distance_m'] for q in r['capture']['20']])*1e6 for r in rr],'o-',color=COL[j],label=f'{freq/1e6:.1f} MHz')
 ax[0].set(xlabel='X+ phase relative to X− (deg)',ylabel='X-restoring location (µm)',title='Finite-disc force roots');ax[0].legend();ax[1].set(xlabel='Relative phase (deg)',ylabel='Mean final band distance (µm)',title='Same absolute starts; F1-only at 250 ms');ax[1].legend();manifest.append(save(f,6,'Frequency and phase dependence',['STATIC_RESULTS.json','SINGLE_RESULTS.json']))
 selected=[]
 for seed in range(3):
  candidates=[r for r in good if r['id'] in [f'long150_random_{seed}',f'long250_random_{seed}']];selected.append(max(candidates,key=lambda r:r['completed_s']))
 a0=series('HF',selected[0]['id']);f=plt.figure(figsize=(13,9));a=f.add_subplot(221,projection='3d');b=f.add_subplot(222,projection='3d');c=f.add_subplot(223);d=f.add_subplot(224);xyz(a,a0['positions_m'][0],'Frozen random cloud A');xyz(b,a0['positions_m'][-1],f'Full dynamics: {selected[0]["completed_s"]*1e3:.0f} ms')
 for r in selected:
  z=series('HF',r['id']);label=r['state'].replace('AG32_','');c.plot(z['time_s']*1000,z['S'],label=label);dist=abs(z['positions_m'][:,:,0,None]-z['roots_m']);physical_occupied=(dist.min(axis=1)<=20e-6).sum(axis=1);d.step(z['time_s']*1000,physical_occupied,where='post',label=label)
 c.axhline(12,color=COL[1],ls='--',label='S threshold');c.set(xlabel='Physical time (ms)',ylabel='Peak structure factor S',title='Actual evolution; no inferred pattern');c.legend(fontsize=8);d.set(xlabel='Physical time (ms)',ylabel='Occupied physical bands (<20 µm)',title='Actual force roots; proximity is not a stripe');d.legend(fontsize=8);manifest.append(save(f,7,'Multibubble stripe dynamics',['SCIENTIFIC_RESULTS.json','series/HF_long*_random_*.npz']))
 f,ax=plt.subplots(2,2,figsize=(13,9));pairs=[r for r in s['paired'] if 'max_fixed_ID_m' in r];dur=[r['paired_duration_s']*1000 for r in pairs]
 for key,label,color in [('max_fixed_ID_m','Fixed ID',COL[0]),('COM_error_m','COM',COL[1]),('Rg_error_m','Final Rg',COL[2])]:ax[0,0].scatter(dur,[r[key]*1e9 for r in pairs],s=20,color=color,label=label)
 ax[0,0].axhline(50,ls='--',color='black',label='Fixed-ID limit');ax[0,0].set(xlabel='Paired duration (ms)',ylabel='Difference (nm)',title='Full HF vs guarded LF');ax[0,0].legend(fontsize=8)
 conv=s['convergence'];ax[0,1].barh([r['id'].replace('long150_','') for r in conv],[r['error_m']*1e9 for r in conv],color=COL[2]);ax[0,1].axvline(50,ls='--',color=COL[1]);ax[0,1].set(xlabel='Maximum fixed-ID difference (nm)',title='Independent half-step / fine-cache checks');ax[0,1].tick_params(axis='y',labelsize=8)
 ax[1,0].scatter(dur,[r['structure_peak_difference'] for r in pairs],color=COL[0]);ax[1,0].set(xlabel='Paired duration (ms)',ylabel='Final peak-S difference',title='Periodic structure comparison')
 categories=['Band assignment','Held capture','Neighbor graph'];values=[sum(r[k] for r in pairs) for k in ['band_disagreements','capture_outcome_disagreements','neighbor_disagreements']];ax[1,1].bar(categories,values,color=COL[:3]);ax[1,1].set(ylabel='Disagreements over compared samples',title='Fixed-ID classification and contacts',ylim=(0,max(1,max(values)*1.15)))
 for j,v in enumerate(values):ax[1,1].text(j,v+.02,f'{v}',ha='center')
 manifest.append(save(f,8,'High and low fidelity validation',['SCIENTIFIC_RESULTS.json']))
 f,ax=plt.subplots(1,2,figsize=(12,5))
 for kind,color in [('HF',COL[0]),('LF',COL[2])]:
  rr=[r for r in rows if r['kind']==kind];ax[0].scatter([r['completed_s']*1000 for r in rr],[r['wall_s'] for r in rr],color=color,label=kind,alpha=.7)
 rr=[r for r in rows if r['kind']=='LF'];ax[1].scatter([r['completed_s']*1000 for r in rr],[r['HF_fallback_steps'] or 0 for r in rr],color=COL[1]);ax[0].set(xlabel='Completed physical time (ms)',ylabel='Measured task wall time (s)',title='Shared CPU/GPU resources, not exclusive benchmark');ax[0].legend();ax[1].set(xlabel='Completed physical time (ms)',ylabel='Original HF fallback steps',title='Near-contact fallback cost retained');manifest.append(save(f,9,'Long-time integration and cost',['SCIENTIFIC_RESULTS.json','LF_LADDER_RESULTS.json']))
 f,ax=plt.subplots(3,2,figsize=(12,11))
 with h5py.File(BASE/'ROTATION_STATIC.h5') as raw:
  for j,a in enumerate(['ROT_PLUS','ROT_MINUS']):
   rr=[r for r in rot['rows'] if r['action']==a and r['z_m']==0];radii=[r['radius_m']*1e6 for r in rr];radial=[]
   for r in rr:
    pos=raw[r['id']+'/positions_m'][:];rad=pos.copy();rad[:,2]=0;rad/=np.linalg.norm(rad,axis=1)[:,None];radial.append(np.sum(raw[r['id']+'/primary_force_N'][:]*rad,axis=1).mean()*1e12)
   ax[0,0].plot(radii,[r['mean_F1_tangent_N']*1e12 for r in rr],'o-',color=COL[j],label=a);ax[0,1].plot(radii,radial,'o-',color=COL[j],label=a)
   ax[1,0].plot(radii,[r['mean_vtheta_m_s']*1e6 for r in rr],'o-',color=COL[j],label=a+' tangent');ax[1,0].plot(radii,[r['mean_vradial_m_s']*1e6 for r in rr],'--',color=COL[j],label=a+' radial')
   z=series('HF','ring50_'+a);ax[1,1].plot(z['omega_time_s'][::100]*1000,z['omega_s'][::100],color=COL[j],label=a);ax[2,0].plot(z['angle_time_s']*1000,z['angle_deg'],color=COL[j],label=a);ax[2,1].plot(z['omega_time_s'][::100]*1000,z['rigid_fraction'][::100],color=COL[j],label=a)
 for a in ax.ravel():a.legend(fontsize=8);a.axhline(0,color='#aaaaaa',lw=.5)
 ax[0,0].set(xlabel='Prescribed ring radius (µm)',ylabel='Mean tangential F1 (pN)');ax[0,1].set(xlabel='Prescribed ring radius (µm)',ylabel='Mean radial F1 (pN)');ax[1,0].set(xlabel='Prescribed ring radius (µm)',ylabel='Constrained velocity (µm/s)');ax[1,1].set(xlabel='Physical time (ms)',ylabel='Group angular velocity (rad/s)');ax[2,0].set(xlabel='Physical time (ms)',ylabel='Integrated group angle (deg)',title='Prescribed ring, full particle dynamics');ax[2,1].axhline(.8,color='black',ls=':');ax[2,1].set(xlabel='Physical time (ms)',ylabel='Rigid rotation fraction',title='Strong criterion: 0.8');manifest.append(save(f,10,'Rotation force and velocity',['ROTATION_STATIC.json','ROTATION_STATIC.h5','series/HF_ring50_*.npz']))
 m=read(BASE/'MISSING_PHYSICS_PREFLIGHT.json');f,ax=plt.subplots(1,2,figsize=(12,5));ax[0].plot(np.array(m['time_s'])*1000,np.array(m['displacement_m'])*1e6,'o-',color=COL[0]);ax[0].set(xlabel='Time (ms)',ylabel='Unbounded buoyancy displacement scale (µm)',title='Theory only; not added to trajectories');ax[1].axis('off');ax[1].text(.02,.95,f'R = 1.78 µm; µ = 1 mPa s\nStokes buoyancy speed ≈ {m["unbounded_immobilized_buoyancy_speed_m_s"]*1e6:.2f} µm/s\nViscous layer at 2.5 MHz ≈ {m["boundary_layer_m"]["2500000"]*1e6:.3f} µm\nBuoyancy / acoustic |Fz| near Z=58 µm ≈ 3.11\n\nStreaming: NOT_MODELED\nBuoyancy: NOT_MODELED\nNear-wall buoyant mobility: NOT_MODELED\nPolydisperse responses: NOT_VALIDATED\n\nNext: field calibration + independent\nstreaming / OFF-motion benchmarks',va='top',fontsize=11);manifest.append(save(f,11,'Missing-physics preflight',['MISSING_PHYSICS_PREFLIGHT.json']))
 f,ax=plt.subplots(1,2,figsize=(12,5));ax[0].bar(['HF complete','HF failed','LF complete','LF failed'],[s['HF_complete'],s['HF_attempts']-s['HF_complete'],s['LF_complete'],s['LF_attempts']-s['LF_complete']],color=[COL[0],COL[1],COL[2],COL[1]]);ax[0].set(ylabel='Actual trajectory attempts',title='Failures retained');ax[0].tick_params(axis='x',rotation=20);ax[1].axis('off');ax[1].text(.02,.95,f'Original F1/F2/contact: preserved\nIndependent accepted-step audit: {"PASS" if s["audit_pass"] else "FAIL"}\nLongest full HF: {s["max_HF_time_s"]*1000:.0f} ms\nLongest LF: {s["max_LF_time_s"]*1000:.0f} ms\nHF/LF agreeing window: {s["max_paired_admitted_s"]*1000:.0f} ms\nStep + cache refined pair: {time_windows(s)['max_new_paired_and_refined_s']*1000:.0f} ms\n\nRandom-cloud stripes: NOT_CONFIRMED\nSustained rigid rotation: NOT_CONFIRMED\nReal PDMS / camera validation: NOT_RUN\nElectrical calibration: NOT_RUN\ntraining_ready = false',va='top',fontsize=11);manifest.append(save(f,12,'Numerical and scientific admission',['SCIENTIFIC_RESULTS.json','FINAL_GATE.json']))
 for r in manifest:r['png_sha256']=sha(OUT/r['png']);r['pdf_sha256']=sha(OUT/r['pdf'])
 write(OUT/'MANIFEST.json',dict(figures=manifest,wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
