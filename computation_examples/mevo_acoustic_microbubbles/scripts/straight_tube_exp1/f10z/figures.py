"""Paper figures from actual stored particle data. No fabricated success curves."""
import argparse
import numpy as np,h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from .common import *

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelcolor':'#222222','xtick.color':'#222222','ytick.color':'#222222','axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white','pdf.fonttype':42})
COLORS=['#244c71','#4477aa','#66a5ad','#bb8844','#aa3377','#228833']
OUT=BASE/'figures'
def save(fig,name):
 OUT.mkdir(exist_ok=True);fig.savefig(OUT/(name+'.png'),dpi=220);fig.savefig(OUT/(name+'.pdf'));plt.close(fig)
def data(name):
 with h5py.File(BASE/'dataset'/(name+'.h5')) as h:return {k:h[k][:] for k in h if isinstance(h[k],h5py.Dataset)}
def key(initial,action='PAIR',duration=2000,variant='DT1'):return f'{initial}__{action}__{duration}us__{variant}'
def balls(ax,q,color='#3e84bc',alpha=.8):
 for j in np.argsort(q[:,2]):ax.add_patch(Circle(q[j,:2]*1e6,R*1e6,facecolor=color,edgecolor='#163b5a',linewidth=.35,alpha=alpha))
 ax.set_aspect('equal');ax.autoscale_view();ax.set_xlabel('X (µm)');ax.set_ylabel('Y (µm)')
def initial_figure():
 cfg=read(CONFIG);audit=read(BASE/'INITIAL_STATE_AUDIT.json');fig=plt.figure(figsize=(15,11),layout='constrained')
 with h5py.File(BASE/'INITIAL_STATES.h5') as h:
  for i,name in enumerate(cfg['primary_states']):
   x=h[name][:]*1e6;row=i//2;col=2*(i%2);a=fig.add_subplot(3,4,row*4+col+1,projection='3d');a.scatter(*x.T,c='#3479b1',s=20,depthshade=True);span=max(np.ptp(x,axis=0))*1.15
   for dim in ['x','y','z']:getattr(a,'set_'+dim+'lim')(-span/2,span/2)
   a.set_box_aspect((1,1,1));a.set_xlabel('X (µm)',labelpad=0);a.set_ylabel('Y (µm)',labelpad=0);a.set_zlabel('Z (µm)',labelpad=0);a.tick_params(labelsize=8);a.view_init(22,-54);a.set_title(name.replace('_',' ')+' | 3D centres',fontsize=10)
   b=fig.add_subplot(3,4,row*4+col+2);balls(b,h[name][:]);b.set_title('XY projection | min d = %.2f µm\n%d proximity components'%(audit[name]['min_distance_m']*1e6,audit[name]['components']),fontsize=10);b.grid(alpha=.12)
 fig.suptitle('F10Z | Preregistered initial configurations\n32 bubbles, R = 1.78 µm; legal geometry; prescribed near-contact states are not acoustically formed',fontsize=14)
 save(fig,'Fig_F10Z_1_initial_configurations')

def main():
 p=argparse.ArgumentParser();p.add_argument('--initial-only',action='store_true');a=p.parse_args();initial_figure()
 if a.initial_only:return
 runs=read(BASE/'analysis/RUN_ANALYSIS.json');checks=read(BASE/'analysis/NUMERICAL_CHECKS.json');response=read(BASE/'analysis/ACTION_RESPONSE.json');ready=read(BASE/'analysis/WORLD_MODEL_READINESS.json');cfg=read(CONFIG);states=cfg['primary_states'];m=data(key('medium_train'))
 # 2: common force arrow scaling, distinct velocity units.
 with h5py.File(BASE/'dataset'/(key('medium_train')+'.h5')) as h:f={k:h['force_samples'][k][0] for k in ['primary_force_N','secondary_force_N','total_force_N']}
 fig,ax=plt.subplots(1,4,figsize=(15,4.7),layout='constrained');q=m['positions_m'][0];scale=max(np.linalg.norm(v,axis=1).max() for v in f.values())
 for axis,(name,v) in zip(ax[:3],f.items()):
  balls(axis,q,alpha=.35);axis.quiver(q[:,0]*1e6,q[:,1]*1e6,v[:,0]/scale*7,v[:,1]/scale*7,angles='xy',scale_units='xy',scale=1,color='#b65c34',width=.005);axis.set_title(name.replace('_force_N','').replace('_',' ').title());axis.text(.02,.02,'Mean |F| = %.3g pN'%(np.linalg.norm(v,axis=1).mean()*1e12),transform=axis.transAxes,fontsize=9)
 v=m['relative_velocity_m_s'][0];vscale=np.linalg.norm(v,axis=1).max();balls(ax[3],q,alpha=.35);ax[3].quiver(q[:,0]*1e6,q[:,1]*1e6,v[:,0]/vscale*7,v[:,1]/vscale*7,angles='xy',scale_units='xy',scale=1,color='#336655');ax[3].set_title('Measured relative velocity');ax[3].text(.02,.02,'Max |v relative| = %.2f µm/s'%(vscale*1e6),transform=ax[3].transAxes,fontsize=9)
 fig.suptitle('F10Z | X pair, medium dispersed initial state, t = 0\nForce panels share scale: 7 µm arrow = %.3g pN; velocity has a separate scale. Homogeneous-water proxy.'%(scale*1e12),fontsize=12);save(fig,'Fig_F10Z_2_force_and_relative_motion')
 # 3: actual full histories; no endpoint-only inference.
 fig,ax=plt.subplots(2,3,figsize=(13,8),layout='constrained');spec=[('rg_m','ΔRg (nm)',1e9),('nn_m','Δmedian nearest neighbour (nm)',1e9),('components','Proximity components',1),('edge_count','Proximity edges (d ≤ 2.1R)',1),('radial_mean_m_s','Mean relative radial speed (µm/s)',1e6),('pair_quantiles_m','Median pair-distance change (nm)',1e9)]
 for color,state in zip(COLORS,states):
  d=data(key(state));t=d['time_s']*1e3
  for axis,(metric,label,factor) in zip(ax.flat,spec):
   y=d[metric];xx=d['interval_time_s']*1e3 if metric=='radial_mean_m_s' else t
   if metric=='pair_quantiles_m':y=y[:,1]
   if metric in ['rg_m','nn_m','pair_quantiles_m']:y=y-y[0]
   axis.plot(xx,y*factor,label=state.replace('_',' '),color=color,lw=1.4);axis.set_xlabel('Physical time (ms)');axis.set_ylabel(label);axis.grid(alpha=.18)
 for axis,(metric,_,factor) in zip(ax.flat,spec):
  d=data(key('medium_train','OFF'));y=d[metric];tt=d['interval_time_s'] if metric=='radial_mean_m_s' else d['time_s']
  if metric=='pair_quantiles_m':y=y[:,1]
  if metric in ['rg_m','nn_m','pair_quantiles_m']:y=y-y[0]
  axis.plot(tt*1e3,y*factor,'k--',lw=1,label='medium OFF')
 ax[0,0].legend(fontsize=8);fig.suptitle('F10Z | Relative convergence and proximity networks\nActual X-pair trajectories; OFF reference; geometric proximity is not fusion',fontsize=14);save(fig,'Fig_F10Z_3_dispersed_convergence')
 # 4: translational and centered state changes separated.
 fig,ax=plt.subplots(2,6,figsize=(21,8),layout='constrained')
 for row,state in enumerate(['original','medium_train']):
  d=data(key(state));t=d['time_s']*1e3
  ax[row,0].plot(t,(d['centroid_m']-d['centroid_m'][0])*1e6);ax[row,0].set_ylabel(state.replace('_',' ')+'\nCOM displacement (µm)');ax[row,0].legend(['X','Y','Z'],fontsize=8)
  balls(ax[row,1],d['centered_positions_m'][0],color='#c0c6cd',alpha=.5);balls(ax[row,1],d['centered_positions_m'][-1],alpha=.65);ax[row,1].set_title('Centred XY: grey initial, blue final',fontsize=9)
  ax[row,2].plot(t,(d['rg_m']-d['rg_m'][0])*1e9);ax[row,2].set_ylabel('ΔRg (nm)')
  Q=d['tensor_m2'];ax[row,3].plot(t,(Q[:,0,0]/Q[0,0,0]-1)*100,label='Qxx');ax[row,3].plot(t,(Q[:,1,1]/Q[0,1,1]-1)*100,label='Qyy');ax[row,3].plot(t,(Q[:,2,2]/Q[0,2,2]-1)*100,label='Qzz');ax[row,3].legend(fontsize=8);ax[row,3].set_ylabel('Diagonal covariance change (%)')
  eig,vec=np.linalg.eigh(Q);aspect=np.sqrt(eig[:,-1]/eig[:,0]);ax[row,4].plot(t,(aspect/aspect[0]-1)*100);ax[row,4].set_ylabel('Aspect ratio change (%)')
  valid=(eig[:,-1]-eig[:,-2])>1e-3*eig[:,-1]
  if valid[0]:
   angle=np.degrees(np.arccos(np.clip(abs(vec[:,:,-1]@vec[0,:,-1]),0,1)));angle[~valid]=np.nan;ax[row,5].plot(t,angle);ax[row,5].set_ylabel('Principal-axis change (degree)')
  else:ax[row,5].text(.08,.5,'Axis NOT IDENTIFIABLE\nDegenerate initial eigenvalues',transform=ax[row,5].transAxes,fontsize=9)
  for j in [0,2,3,4,5]:ax[row,j].set_xlabel('Physical time (ms)');ax[row,j].set_xlim(t[0],t[-1]);ax[row,j].grid(alpha=.2)
 fig.suptitle('F10Z | Translation is separated from internal deformation\nTrue coordinates and radius; XY projection can overlap distinct Z layers',fontsize=14);save(fig,'Fig_F10Z_4_translation_vs_deformation')
 # 5: discrete action contrasts and finite differences.
 actions=['OFF','XM','XP','PAIR','YM','YP'];fig,ax=plt.subplots(2,3,figsize=(13,7.5),layout='constrained')
 for column,(idx,title) in enumerate([(0,'ΔCOM X (nm)'),(1,'ΔCOM Y (nm)'),(2,'ΔRg XY (nm)')]):
  values=np.array([[runs.get(key(state,act),{}).get('endpoint_state_delta',[np.nan]*6)[idx]*1e9 for act in actions] for state in ['original','medium_train']]);vmax=max(np.nanmax(abs(values)),1e-12);im=ax[0,column].imshow(values,cmap='RdBu_r',vmin=-vmax,vmax=vmax,aspect='auto');ax[0,column].set_xticks(range(6),actions);ax[0,column].set_yticks([0,1],['original','medium']);ax[0,column].set_title(title)
  for i in range(2):
   for j in range(6):
    rr=runs.get(key(['original','medium_train'][i],actions[j]));label='NOT\nRUN'
    if rr:
     noise=rr.get('state_feature_numerical_difference');suffix='\nU' if noise is None else '\n*' if abs(values[i,j])<=5*noise[idx]*1e9 else '';v=values[i,j];number=f'{v:.0f}' if abs(v)>=1 else f'{v:.2f}' if abs(v)>=.01 else f'{v:.1g}'.replace('e-0','e-');label=number+suffix
    ax[0,column].text(j,i,label,ha='center',va='center',fontsize=7,color='white' if np.isfinite(values[i,j]) and abs(values[i,j])>.6*vmax else '#111111')
  fig.colorbar(im,ax=ax[0,column],shrink=.7)
 cols=response['columns'];labels=[c['low']+' / '+c['high'] for c in cols]
 for j,idx in enumerate([0,2]):
  ax[1,j].bar(range(len(cols)),[c['derivative_SI_per_unit'][idx]*1e9 for c in cols],color=COLORS[:len(cols)]);ax[1,j].set_xticks(range(len(cols)),labels,rotation=20,fontsize=8);ax[1,j].set_ylabel(('COM X' if idx==0 else 'Rg XY')+' sensitivity (nm / parameter unit)');ax[1,j].set_title('Amplitude: per Pa; phase: per degree',fontsize=10)
 ax[1,2].bar(range(len(response['singular_values'])),response['singular_values'],color='#4477aa');ax[1,2].set_title('Normalized local sensitivity');ax[1,2].set_ylabel('Singular value');ax[1,2].set_ylim(0,max(response['singular_values'])*1.7);ax[1,2].text(.05,.76,'Two refined X-pair parameters\nLocal rank = '+str(response['local_rank'])+'\nNot a global 4-PZT control proof',transform=ax[1,2].transAxes,fontsize=9)
 fig.suptitle('F10Z | Action effects and bounded local sensitivity\nY-only is proxy specific | * response ≤ 5× measured numerical difference; U = no matched refinement',fontsize=13);save(fig,'Fig_F10Z_5_action_response_matrix')
 # 6: actual time extension, errors and costs.
 fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained');rr=[r for r in runs.values() if r['initial']=='medium_train' and r['action']=='PAIR' and r['name'].endswith('DT1')];rr=sorted(rr,key=lambda r:r['requested_s'])
 ax[0,0].bar([str(r['requested_s']*1e3) for r in rr],[r['completed_s']*1e3 for r in rr],color=['#4477aa' if r['status']=='PASS' else '#bb5533' for r in rr]);ax[0,0].set_xlabel('Requested window (ms)');ax[0,0].set_ylabel('Actually completed time (ms)')
 for r in rr:
  k=r['name'][:-5];g=checks.get(k)
  if not g or not g.get('admitted'):continue
  ax[0,1].plot([0,1,2],[max(c[a] for c in g['comparisons'])*1e9 for a in ['COM_m','Rg_m','fixed_id_m']],marker='o',label=f"{r['requested_s']*1e3:g} ms")
 ax[0,1].set_xticks([0,1,2],['COM','Rg','Fixed ID']);ax[0,1].set_ylabel('Maximum dt/cache difference (nm)');ax[0,1].legend(fontsize=8)
 ax[1,0].bar([str(r['requested_s']*1e3) for r in rr],[r['wall_s'] for r in rr],color='#4477aa');ax[1,0].set_xlabel('Physical window (ms)');ax[1,0].set_ylabel('Actual single-run wall time (s)')
 longest=max(rr,key=lambda r:r['completed_s']);d=data(longest['name']);st=d['step_metrics'];tt=d['step_time_duration_segment'][:,0]*1e3;ax[1,1].plot(tt,st[:,2]*1e9,label='Overlap');ax[1,1].plot(tt,st[:,5]*1e9,label='Wall correction');ax[1,1].axhline(1,color='#bb5533',ls='--',label='1 nm hard geometry limit');ax[1,1].set_xlabel('Physical time (ms)');ax[1,1].set_ylabel('Geometric diagnostics (nm)');ax[1,1].legend(fontsize=8)
 fig.suptitle('F10Z | Time extension and unchanged numerical gates\nMissing longer windows are NOT RUN; no interpolated or stitched trajectories',fontsize=14);save(fig,'Fig_F10Z_6_time_extension_convergence')
 # 7: model truth to observation semantics.
 fig=plt.figure(figsize=(12,8),layout='constrained');a=fig.add_subplot(221,projection='3d');a.scatter(*m['positions_m'][0].T*1e6,c='#b1b7be',s=18);a.scatter(*m['positions_m'][-1].T*1e6,c='#337ab0',s=18);a.set_title('SIMULATION_3D_TRUTH');a.set_xlabel('X (µm)');a.set_ylabel('Y (µm)');a.set_zlabel('Z (µm)');a.set_box_aspect((1,1,1));b=fig.add_subplot(222);balls(b,m['positions_m'][-1]);b.set_title('SIMULATION_XY_PROJECTION')
 c=fig.add_subplot(223);c.plot(m['time_s']*1e3,m['xy_rg_m']/m['xy_rg_m'][0],label='Rg XY / initial');c.plot(m['time_s']*1e3,m['xy_nn_m']/m['xy_nn_m'][0],label='NN XY / initial');c.set_xlabel('Physical time (ms)');c.set_ylabel('Normalized projected distribution');c.legend()
 d=fig.add_subplot(224);d.axis('off');d.text(0,.95,'Shared definitions\nCentroid, covariance, Rg, median NN\nRelative motion and fixed spatial support\n\nSimulation: 32 fixed material IDs\nVision: variable visible candidates, no IDs\n\nNo calibrated PSF, depth or occlusion model\nNo synthetic camera image generated\nNo real Chronos/PZT validation',va='top',linespacing=1.7,fontsize=11)
 fig.suptitle('F10Z | Alignment with F10Y observation semantics\nProjected model positions are not measured camera detections',fontsize=14);save(fig,'Fig_F10Z_7_simulation_F10Y_state')
 # 8: coverage and honest model readiness.
 fig,ax=plt.subplots(2,2,figsize=(12,8),layout='constrained');coverage=np.zeros((len(states),len(actions)))
 for i,s in enumerate(states):
  for j,act in enumerate(actions):
   r=runs.get(key(s,act));coverage[i,j]=0 if r is None else -1 if r['status']!='PASS' else 2 if r.get('refinement_admitted') else 1
 im=ax[0,0].imshow(coverage,cmap='viridis',vmin=-1,vmax=2,aspect='auto');ax[0,0].set_xticks(range(6),actions);ax[0,0].set_yticks(range(len(states)),[s.replace('_',' ') for s in states],fontsize=8);ax[0,0].set_title('Independent 2 ms initial/action coverage');fig.colorbar(im,ax=ax[0,0],ticks=[-1,0,1,2],label='-1 failed; 0 not run; 1 complete; 2 refined')
 rr=[r for r in runs.values() if r['name'].endswith('DT1') and r['action']=='PAIR' and r['requested_s']==.002];ax[0,1].bar(range(len(rr)),[-r['relative_Rg_change']*100 for r in rr],color='#4477aa');ax[0,1].set_xticks(range(len(rr)),[r['initial'].replace('_','\n') for r in rr],fontsize=7);ax[0,1].set_ylabel('Rg reduction (%)');ax[0,1].set_title('Measured contraction is not arbitrary shape control')
 bs=ready['baselines'];methods=['state_hold','constant_velocity_first_0p5ms','training_action_mean'];xx=np.arange(len(bs))
 for i,method in enumerate(methods):ax[1,0].bar(xx+(i-1)*.25,[np.linalg.norm(b['errors_m'][method])*1e6 for b in bs],.25,label=method.replace('_',' '))
 ax[1,0].set_xticks(xx,[b['name'].replace('medium_','').replace('__2000us__DT1','') for b in bs],rotation=20,fontsize=7);ax[1,0].set_ylabel('Endpoint [cx, cy, Rg] error norm (µm)');ax[1,0].legend(fontsize=7)
 ax[1,1].axis('off');ax[1,1].text(0,.95,ready['status']+'\n\nIndependent train groups: '+str(ready['independent_training_trajectories'])+'\nHeld-out initial groups: '+str(ready['heldout_initial_trajectories'])+'\nRefined X-pair local rank: '+str(ready['local_action_rank'])+'\n\nWorld model: NOT TRAINED\nMPPI: NOT RUN\nPhysical calibration: NOT VERIFIED\ntraining_ready = false',va='top',linespacing=1.7,fontsize=10)
 fig.suptitle('F10Z | Dataset coverage, simple prediction and scientific admission',fontsize=14);save(fig,'Fig_F10Z_8_world_model_readiness')
 write(OUT/'FIGURE_INPUTS.json',dict(input_analysis_sha256=sha(BASE/'analysis/RUN_ANALYSIS.json'),source_identity=identity(),figures=[p.name for p in sorted(OUT.glob('*.png'))],data='Saved accepted trajectories only; scientific negative and unresolved results retained'))
if __name__=='__main__':main()
