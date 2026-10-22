"""Twelve publication figures drawn only from retained scientific records."""
import time,numpy as np,h5py,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .common import *
COLORS=['#7653a6','#247bb2','#b36f35']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelcolor':'#20252a','text.color':'#20252a','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':200})
def loadall():return read(BASE/'ALL_RUNS.json') if (BASE/'ALL_RUNS.json').exists() else [read(p) for p in (BASE/'runs').glob('*/SUMMARY.json')]
def traces(r):
 with h5py.File(BASE/r['trajectory_path']) as h:return h['time_s'][:],h['positions_m'][:]
def metrics(r):return dict(np.load(BASE/'runs'/r['id']/'METRICS.npz'))
def main():
 begin=time.perf_counter();out=BASE/'figures';out.mkdir(exist_ok=True);manifest=[];rows=loadall();byid={r['id']:r for r in rows};valid=[r for r in rows if r['status']=='PASS'];comparison=read(BASE/'COMPARISON_V2.json');validation=read(BASE/'VALIDATION_COMBINED.json');selected=read(BASE/'ELITE_VALIDATION_PREREGISTRATION.json')['selected'];long=read(AD/'LONG.json');forces=read(AD/'FORCE_AUDIT.json')
 def save(fig,num,title,sources):
  fig.suptitle(title,fontsize=13,fontweight='bold');fig.text(.01,.01,'EXPLORATORY | Homogeneous water | Uncalibrated shell / contact / source amplitude',fontsize=8,color='#424b54');fig.tight_layout(rect=[0,.035,1,.95])
  stem=f'Fig_{num:02d}_'+title.lower().replace(' ','_').replace('/','_')
  for ext in ['png','pdf']:fig.savefig(out/(stem+'.'+ext),bbox_inches='tight')
  manifest.append(dict(number=num,title=title,png=stem+'.png',pdf=stem+'.pdf',sources=sources,sha256={e:sha(out/(stem+'.'+e)) for e in ['png','pdf']}));plt.close(fig)
 # 1
 fig,axes=plt.subplots(2,3,figsize=(12,7));ax=axes[0];reg=read(BASE/'REGRESSION.json');ax[0].axis('off');ax[0].text(.05,.85,'2.5 MHz historical regression\n\n'+ '\n'.join(r['id']+': '+('PASS' if r['PASS'] else 'FAIL') for r in reg)+'\n\nPosition, pressure, F1, F2:\nbitwise identical in tested cases',va='top')
 for i,r in enumerate(selected):
  t,x=traces(r);ax[1].scatter(x[0,:,0]*1e6,x[0,:,1]*1e6,s=10,label='Search initial state' if i==0 else None)
  if i==0:break
 ax[1].set(xlabel='X (um)',ylabel='Y (um)',title='Fixed-ID N=64 initial geometry',aspect='equal');p=read(BASE/'PILOT.json')['parallel'];ax[2].bar([str(r['workers']) for r in p],[r['wall_s'] for r in p],color=COLORS);ax[2].set(xlabel='Workers',ylabel='Wall time (s)',title='Four 0.25 ms trajectories');fg=read(AD/'2500000/FIELD_GATE.json')['checks'];keys=['reference_pressure','coarse_pressure','fine_pressure'];axes[1,0].bar(['Reference','3 um grid','1.5 um grid'],[fg[k] for k in keys],color=COLORS);axes[1,0].set(yscale='log',ylabel='Pressure error / source',title='Rebuilt 2.5 MHz field check')
 rr=read(AD/'2500000/RESPONSE_1.0.json')['regression']['points'];axes[1,1].semilogx([r['pressure_Pa'] for r in rr],[r['error']*100 for r in rr],'o-');axes[1,1].axhline(2,color='#a94545',ls='--');axes[1,1].set(xlabel='Actual pressure (Pa)',ylabel='Complex response difference (%)',title='New APECSS vs historical table')
 axes[1,2].axis('off');axes[1,2].text(.01,.9,'Historical trajectory path is unchanged.\n\nRebuilt field is converged;\ncache nodes are not bitwise identical.\n\nHistorical 2.5 MHz response table is\nretained for production trajectories.\nNew radial table is a separate check.',va='top',fontsize=9)
 save(fig,1,'Baseline and reproducibility' ,['REGRESSION.json','PILOT.json','INITIAL_STATES.h5'])
 # 2
 fig,ax=plt.subplots(1,3,figsize=(12,3.8))
 for i,(name,v) in enumerate(comparison['methods'].items()):
  h=v['history'];ax[0].plot([r['evaluations'] for r in h],[r['cells'] for r in h],'-o',ms=3,label=name,color=COLORS[i]);ax[1].plot([r['elapsed_s'] for r in h],[r['qd_score'] for r in h],label=name,color=COLORS[i]);ax[2].bar(i,v['optimizer_wall_s'],color=COLORS[i])
 ax[0].set(xlabel='Charged evaluations',ylabel='Occupied cells');ax[0].legend(fontsize=8);ax[1].set(xlabel='Charged wall time (s)',ylabel='QD score');ax[2].set(xticks=range(3),xticklabels=list(comparison['methods']),ylabel='Optimizer overhead (s)');save(fig,2,'QD algorithm comparison',['COMPARISON_V2.json','search_v2/*/RESULT.json'])
 # 3
 fig,ax=plt.subplots(1,3,figsize=(12,3.8));ranges=np.array(read(CONFIG)['ranges'])
 for a,(name,v) in zip(ax,comparison['methods'].items()):
  e=v['archive']['elites'];a.set(xlim=ranges[0],ylim=ranges[1],xlabel='Delta log Rg (3D)',ylabel='Delta log aspect (3D)',title=name)
  for z in np.linspace(*ranges[0],9):a.axvline(z,c='#dce1e5',lw=.6)
  for z in np.linspace(*ranges[1],9):a.axhline(z,c='#dce1e5',lw=.6)
  if len(e['objective']):m=np.array(e['measures']);im=a.scatter(m[:,0],m[:,1],c=e['objective'],vmin=0,vmax=100,cmap='Purples',s=85,edgecolors='#333333')
 fig.colorbar(im,ax=ax[-1],label='Quality (COM penalty)');save(fig,3,'Shape diversity archives',['COMPARISON_V2.json'])
 # 4
 fig,ax=plt.subplots(1,3,figsize=(12,3.8))
 search_cloud=[]
 for name in comparison['methods']:
  search_cloud += [r['descriptor'] for r in read(BASE/'search_v2'/name/'RESULT.json')['rows'] if r['descriptor'] is not None]
 if search_cloud:
  cloud=np.array(search_cloud);ax[2].scatter(cloud[:,0],cloud[:,1],s=13,c='#aab4bd',alpha=.6,label='All completed search results')
 ax[2].legend(fontsize=7)
 for i,r in enumerate(selected):
  ax[0].plot(range(4),r['task']['amplitudes'],'-o',label=f'Elite {i}');ax[1].plot(range(4),np.degrees(r['task']['phases']),'-o');d=r['descriptor'];ax[2].scatter(d[0],d[1],s=80);ax[2].annotate(f'E{i}',d)
 for a in ax[:2]:a.set(xticks=range(4),xticklabels=['Y-','Y+','X-','X+'])
 ax[0].set(ylabel='Source amplitude (Pa-eq)');ax[0].legend();ax[1].set(ylabel='Phase (deg)');ax[2].set(xlabel='Delta log Rg',ylabel='Delta log aspect');save(fig,4,'Actions and selected modes',['ELITE_VALIDATION_PREREGISTRATION.json'])
 # 5: Actual fixed-ID center geometry plus quantitative time histories.
 fig=plt.figure(figsize=(12,7.6));ax=[]
 clouds=[traces(r)[1] for r in selected];limit=max(float(abs(x[[0,-1]]).max())*1e6 for x in clouds)+3
 for i,(r,x) in enumerate(zip(selected,clouds)):
  a=fig.add_subplot(2,3,i+1,projection='3d');a.scatter(*x[0].T*1e6,s=12,c='#b9c0c5',alpha=.5,label='Initial centers');a.scatter(*x[-1].T*1e6,s=18,c=COLORS[i],depthshade=True,label='Final centers')
  a.set(xlabel='X (um)',ylabel='Y (um)',zlabel='Z (um)',xlim=(-limit,limit),ylim=(-limit,limit),zlim=(-limit,limit),title=f'Elite {i} | {r["task"].get("frequency",2500000.)/1e6:g} MHz');a.set_box_aspect((1,1,1));a.view_init(27,-60);a.tick_params(labelsize=8)
  if i==0:a.legend(fontsize=7,loc='upper left')
 for i in range(3):ax.append(fig.add_subplot(2,3,i+4))
 for i,r in enumerate(selected):
  m=metrics(r);label=f'E{i}, {r["task"].get("frequency",2500000.)/1e6:g} MHz';ax[0].plot(m['time_s']*1e3,m['Rg_3d_m']*1e6,label=label);ax[1].plot(m['time_s']*1e3,m['aspect_3d']);ax[2].step(m['time_s']*1e3,m['proximity_components'],where='post')
 for a,y in zip(ax,['Rg (um)','3D aspect ratio','Geometric components']):a.set(xlabel='Physical time (ms)',ylabel=y)
 ax[0].legend(fontsize=8);save(fig,5,'Aggregation splitting and deformation',[r['trajectory_path'] for r in selected])
 # 6
 fig,ax=plt.subplots(1,3,figsize=(12,3.8))
 for i,g in enumerate(validation['groups']):
  rr=g['replicates'];x=np.arange(len(rr))+.12*(i-1);d=np.array([r['descriptor'] if r['descriptor'] else [np.nan,np.nan] for r in rr]);noise=np.array([r['noise_bound'] for r in rr]);ax[0].errorbar(x,d[:,0],5*noise[:,0],fmt='o-',label=str(g['elite'])+(' (not admitted)' if not g['signed_descriptor_repeated'] else ''));ax[1].errorbar(x,d[:,1],5*noise[:,1],fmt='o-');ax[2].plot(x,[r['components'] for r in rr],'o-')
 for a,y in zip(ax,['Delta log Rg','Delta log aspect','Components at endpoint']):a.set(xticks=range(4),xticklabels=['Search','New 1','New 2','New 3'],ylabel=y)
 ax[0].legend();save(fig,6,'Independent seed generalization',['VALIDATION_COMBINED.json','OUTSIDE_VALIDATION.json'])
 # 7: Same pressure scale; four separate complex sources and consistent gradients.
 fig,ax=plt.subplots(3,3,figsize=(12,10.3));freqs=[2400000,2500000,2600000];source_names=['Y-','Y+','X-','X+']
 field_vmax=max(float(abs(np.load(AD/str(f)/'FIELD_DIAGNOSTICS.npz')['xy_pressure'][0]).max())*700 for f in freqs)
 for col,f in enumerate(freqs):
  d=np.load(AD/str(f)/'FIELD_DIAGNOSTICS.npz');a=ax[0,col];im=a.imshow(abs(d['xy_pressure'][0])*700,origin='lower',extent=[-60,60,-60,60],cmap='magma',aspect='equal',vmin=0,vmax=field_vmax);a.set(title=f'Y- only | {f/1e6:g} MHz',xlabel='Y (um)',ylabel='X (um)');fig.colorbar(im,ax=a,label='Pressure magnitude (Pa)')
  for j,name in enumerate(source_names):
   style='-' if j%2==0 else '--';ax[1,col].plot(d['line_y_m']*1e6,abs(d['line_pressure'][j])*700,style,label=name)
   ax[2,col].plot(d['line_y_m']*1e6,np.linalg.norm(d['line_gradient'][j],axis=-1)*700,style,label=name)
  ax[1,col].set(xlabel='Y (um), X = Z = 0',ylabel='Single-source |P| (Pa)',title='Each disc driven at 700 Pa-eq');ax[1,col].legend(ncol=2,fontsize=8)
  ax[2,col].set(xlabel='Y (um), X = Z = 0',ylabel='Complex gradient norm (Pa/m)')
 save(fig,7,'Multifrequency pressure fields',['../f10ad_multifrequency/*/FIELD_DIAGNOSTICS.npz'])
 # 8
 fig,axes=plt.subplots(2,3,figsize=(12,7));ax=axes[0];freqs=[2400000,2500000,2600000]
 for shell in [.8,1.,1.2]:
  ns=[]
  for f in freqs:
   d=read(AD/str(f)/f'RESPONSE_{shell:.1f}.json');ns.append(next(n for n in d['nodes'] if n['pressure_Pa']==1500))
  ax[0].plot(np.array(freqs)/1e6,[n['metrics']['max_R_over_R0'] for n in ns],'o-',label=f'Elasticity x{shell:g}');ax[1].plot(np.array(freqs)/1e6,[np.degrees(n['metrics']['phase_rad']) for n in ns],'o-');ax[2].plot(np.array(freqs)/1e6,[n['metrics']['mean_dot_volume_squared_m6_s2'] for n in ns],'o-')
 for a,y in zip(ax,['Maximum R/R0','Volume phase (deg)','Mean volume-rate squared (m6/s2)']):a.set(xlabel='Frequency (MHz)',ylabel=y)
 ax[0].legend(fontsize=8);ax[1].set_title('Frequency curves at 1500 Pa')
 for f in freqs:
  d=read(AD/str(f)/'RESPONSE_1.0.json');ns=d['nodes'];pp=[n['pressure_Pa'] for n in ns]
  for a,key in zip(axes[1],['max_R_over_R0','phase_rad','maximum_wall_mach']):a.semilogx(pp,[np.degrees(n['metrics'][key]) if key=='phase_rad' else n['metrics'][key] for n in ns],'-o',ms=3,label=f'{f/1e6:g} MHz')
 for a,y in zip(axes[1],['Maximum R/R0','Volume phase (deg)','Maximum wall Mach']):a.set(xlabel='Actual bubble pressure (Pa)',ylabel=y)
 axes[1,0].legend(fontsize=8);axes[1,1].set_title('Shell x1; verified 0.001-12000 Pa')
 save(fig,8,'APECSS frequency and shell responses',['../f10ad_multifrequency/*/RESPONSE_*.json'])
 # 9
 fig,ax=plt.subplots(1,3,figsize=(12,3.8))
 for name in ['X180','Y90','FOUR']:
  rr=sorted([r for r in forces if r['PASS'] and r['shell']==1 and r['offset_y_m']==0 and r['action']==name],key=lambda r:r['frequency']);f=np.array([r['frequency'] for r in rr])/1e6
  for a,key,factor in zip(ax,['primary_max_N','secondary_max_N','relative_velocity_rms_m_s'],[1e12,1e12,1e6]):a.plot(f,[r[key]*factor for r in rr],'o-',label=name)
 for a,y in zip(ax,['Max primary force (pN)','Max secondary force (pN)','Free relative velocity RMS (um/s)']):a.set(xlabel='Frequency (MHz)',ylabel=y)
 ax[0].legend();save(fig,9,'Acoustic force versus frequency',['../f10ad_multifrequency/FORCE_AUDIT.json'])
 # 10
 fig,ax=plt.subplots(1,3,figsize=(12,3.8))
 for name in ['X180','Y90','FOUR']:
  rr=sorted([r for r in long if r['status']=='PASS' and r['id'].endswith(name)],key=lambda r:r['task']['frequency']);f=[r['task']['frequency']/1e6 for r in rr]
  for a,key in zip(ax,['Rg_3d_m','aspect_3d','similarity_rmsd_3d_m']):a.plot(f,[(r['final'][key]/r['initial'][key]-1)*100 if key!='similarity_rmsd_3d_m' else r['final'][key]*1e9 for r in rr],'o-',label=name)
 for a,y in zip(ax,['Rg change (%)','3D aspect change (%)','Nonrigid RMS (nm)']):a.set(xlabel='Frequency (MHz)',ylabel=y)
 ax[0].legend();save(fig,10,'Frequency dependent shape response',['../f10ad_multifrequency/LONG.json'])
 # 11
 fig,axes=plt.subplots(2,3,figsize=(12,7));ax=axes[0];cc=validation['comparisons']+read(AD/'LONG_CONVERGENCE.json');errs=[r['errors']['fixed_id_m']*1e9 for r in cc if 'errors' in r];ax[0].plot(errs,'.',color=COLORS[0]);ax[0].axhline(50,c='#b74949',ls='--',label='Original 50 nm limit');ax[0].set(xlabel='Comparison index',ylabel='Max fixed-ID difference (nm)');ax[0].legend(fontsize=8)
 robust=read(AD/'ROBUSTNESS.json')
 for f in freqs:
  rr=[r for r in robust if r['status']=='PASS' and r['task']['frequency']==f and 'shell' in r['task']];base=next((r for r in long if r['task']['frequency']==f and r['id'].endswith('X180') and r['status']=='PASS'),None)
  if base:rr.append(base)
  rr=sorted(rr,key=lambda r:r['task'].get('shell',1));ax[1].plot([r['task'].get('shell',1) for r in rr],[r['descriptor'][1] for r in rr],'o-',label=f'{f/1e6:g} MHz')
  yy=[r for r in forces if r['PASS'] and r['frequency']==f and r['shell']==1 and r['action']=='Y90'];ax[2].plot([r['offset_y_m']*1e6 for r in yy],[r['primary_mean_N'][1]*1e12 for r in yy],'o-')
 ax[1].set(xlabel='Shell elasticity multiplier',ylabel='Delta log aspect');ax[1].legend(fontsize=8);ax[2].set(xlabel='Initial Y offset (um)',ylabel='Mean primary force Y (pN)');ar=[r['audit'] for r in read(BASE/'AUDIT.json')['rows'] if r['scientific_status']=='VALID_COMPLETE_TRAJECTORY'];axes[1,0].plot([r['force_gate_fraction_max'] for r in ar],'.',ms=3);axes[1,0].axhline(1,c='#b74949',ls='--');axes[1,0].set(xlabel='Complete trajectory index',ylabel='Maximum force residual / gate')
 axes[1,1].plot([r['overlap_m']*1e9 for r in ar],'.',ms=3);axes[1,1].set(xlabel='Complete trajectory index',ylabel='Maximum overlap (nm)')
 for tag,color in [('DT05','#7653a6'),('FINE','#247bb2')]:
  subset=[r for r in cc if 'errors' in r and tag in r.get('candidate',r.get('comparison',''))];axes[1,2].plot([r['errors']['fixed_id_m']*1e9 for r in subset],'.',label=tag,color=color)
 axes[1,2].axhline(50,c='#b74949',ls='--');axes[1,2].set(xlabel='Comparison within group',ylabel='Fixed-ID difference (nm)');axes[1,2].legend(fontsize=8)
 save(fig,11,'Numerical and model sensitivity' ,['VALIDATION_COMBINED.json','../f10ad_multifrequency/ROBUSTNESS.json','../f10ad_multifrequency/FORCE_AUDIT.json'])
 # 12
 fig,axes=plt.subplots(2,3,figsize=(12,7.5));ax=axes[0];names=list(comparison['methods']);ax[0].bar(names,[comparison['methods'][m]['total_charged_wall_s'] for m in names],color=COLORS);ax[0].set(ylabel='Charged method wall (s)');status={k:sum(r['scientific_status']==k for r in rows) for k in sorted(set(r['scientific_status'] for r in rows))};ax[1].bar(range(len(status)),list(status.values()),color=COLORS[:len(status)]);ax[1].set(xticks=range(len(status)),xticklabels=[s.replace('_','\n') for s in status],ylabel='Actual attempts');ax[1].tick_params(axis='x',labelsize=7);ax[2].axis('off');ax[2].text(.03,.9,f'Completed trajectories: {len(valid)} / {len(rows)}\n\nIndependent seeds per elite: 3\n\nSigned response repeated: {sum(g["signed_descriptor_repeated"] for g in validation["groups"])} / {len(validation["groups"])}\n\nArchive cells are not independent modes.\nNo experimental or training admission.',va='top');profile=read(BASE/'evidence/FIELD_TRANSFER_PROFILE.json');axes[1,0].bar(range(5),np.array(profile['mean_s'])*1e6,color='#247bb2');axes[1,0].set(xticks=range(5),xticklabels=['H2D','F1','D2H','F2','D2H'],ylabel='Profile mean (us)',title='N64 GPU transfer and kernels')
 axes[1,1].bar(names,[comparison['methods'][m]['archive']['cells']/comparison['methods'][m]['total_charged_wall_s']*3600 for m in names],color=COLORS);axes[1,1].set(ylabel='Archive cells / charged hour',title='Coverage throughput; not physical modes')
 axes[1,2].axis('off');axes[1,2].text(.02,.9,'Confirmed response sign categories: 1\n\nNo resolved absolute principal extension.\n\nPer-method independent-mode yield:\nNOT_RUN (selected-candidate validation).\n\nStrong outside-archive response:\nthree new-seed dt checks FAIL.',va='top',fontsize=9)
 save(fig,12,'Compute cost and scientific yield' ,['ALL_RUNS.json','COMPARISON_V2.json','VALIDATION_COMBINED.json','evidence/FIELD_TRANSFER_PROFILE.json'])
 write(out/'MANIFEST.json',dict(figures=manifest,wall_s=time.perf_counter()-begin,identity=identity()))
if __name__=='__main__':main()
