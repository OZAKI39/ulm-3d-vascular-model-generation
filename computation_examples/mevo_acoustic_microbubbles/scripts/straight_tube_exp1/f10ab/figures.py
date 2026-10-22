"""Publication figures from saved states, numerical evidence and explicit NOT_RUN."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np,h5py
from .common import *
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
COLORS=['#1f77b4','#ae4968','#347b67','#ad7a26','#6b57a3']

def all_rows():
 return [r for p in sorted((BASE/'analysis').glob('*/SUMMARY.json')) for r in read(p)]
def load(row):
 with h5py.File(BASE/row['dataset']) as h:return {k:h[k][:] for k in h if isinstance(h[k],h5py.Dataset)}
def note(ax,text):ax.axis('off');ax.text(.04,.95,text,transform=ax.transAxes,va='top',wrap=True,fontsize=11)
def save(fig,name,inputs,manifest):
 fig.text(.01,.006,'EXPLORATORY | Homogeneous-water proxy | Physics and electrical calibration unavailable',fontsize=8,color='#555555')
 fig.tight_layout(rect=(0,.035,1,.96));out=BASE/'figures';out.mkdir(exist_ok=True)
 for ext in ['png','pdf']:fig.savefig(out/(name+'.'+ext),dpi=220,bbox_inches='tight',pad_inches=.12)
 plt.close(fig);manifest.append(dict(name=name,inputs=inputs,png_sha256=sha(out/(name+'.png')),pdf_sha256=sha(out/(name+'.pdf'))))

def main():
 rows=all_rows();formal=[r for r in rows if r['stage']=='discovery' and r['status']=='PASS'];manifest=[]
 choices=[r for r in formal if r['action']!='A0_OFF'] or [r for r in rows if r['status']=='PASS']
 chosen=next(r for r in choices if r['state_id']=='S2_r0' and r['action']=='A5_X180');m=load(chosen)
 fig=plt.figure(figsize=(20,7));ax=np.empty((2,7),object)
 for k in range(7):
  ax[0,k]=fig.add_subplot(2,7,k+1,projection='3d');ax[1,k]=fig.add_subplot(2,7,k+8)
 with h5py.File(BASE/'INITIAL_STATES.h5') as h:
  for k in range(7):
   x=h[f'S{k}_r0'][:]*1e6;a=read(BASE/'INITIAL_STATE_AUDIT.json')[f'S{k}_r0']
   ax[0,k].scatter(x[:,0],x[:,1],x[:,2],s=12,c=x[:,2],cmap='viridis');ax[0,k].set_box_aspect(np.maximum(np.ptp(x,axis=0),1));ax[0,k].tick_params(labelsize=6);ax[0,k].set_zlabel('Z (um)',fontsize=7);ax[0,k].set_title(f'S{k}: N={len(x)}\ncomponents={a["proximity_components"]}')
   ax[1,k].scatter(x[:,0],x[:,1],s=12,color=COLORS[0]);ax[1,k].set_aspect('equal')
   ax[0,k].set_xlabel('X (um)');ax[0,k].set_ylabel('Y (um)');ax[1,k].set_xlabel('X (um)');ax[1,k].set_ylabel('Y (um)')
 fig.suptitle('Initial shape library | Prescribed legal configurations, not acoustic formation')
 save(fig,'Fig_F10AB_1_initial_shapes',['INITIAL_STATES.h5'],manifest)
 static=read(BASE/'static/SUMMARY.json')['records'];rs=[r for r in static if r['state_id']=='S2_r0' and r['status']=='PASS']
 fig,ax=plt.subplots(2,2,figsize=(12,8));names=[r['action'].split('_',1)[1] for r in rs];xx=np.arange(len(rs))
 for a,key,title in [(ax[0,0],'pressure_max_Pa','Maximum local pressure (Pa)'),(ax[0,1],'force_relative_rms_N','Relative force RMS (N)'),(ax[1,0],'free_deviatoric_rate_s1','Unconstrained deviatoric rate (1/s)')]:
  a.bar(xx,[r[key] for r in rs],color=COLORS[0]);a.set_xticks(xx,names,rotation=70);a.set_title(title)
 for key,col in [('center_primary_force_N',COLORS[0]),('center_secondary_force_N',COLORS[1]),('center_force_N',COLORS[2])]:ax[1,1].plot(xx,[np.linalg.norm(r[key]) for r in rs],'.-',label=key,color=col)
 ax[1,1].set_xticks(xx,names,rotation=70);ax[1,1].set_title('Mean force magnitude (N)');ax[1,1].legend(fontsize=7)
 fig.suptitle('Fixed-state force diagnostic | Not a dynamic shape-control claim')
 save(fig,'Fig_F10AB_2_acoustic_mechanism',['static/SUMMARY.json'],manifest)
 fig,ax=plt.subplots(2,2,figsize=(11,8));t=m['time_s']*1e3
 for key,col in [('unregistered_rms_3d_m',COLORS[0]),('translation_removed_rms_3d_m',COLORS[1]),('kabsch_rmsd_3d_m',COLORS[2]),('similarity_rmsd_3d_m',COLORS[3])]:ax[0,0].plot(t,m[key]*1e9,label=key.replace('_3d_m',''),color=col)
 ax[0,0].legend(fontsize=8);ax[0,0].set_ylabel('Fixed-ID RMS (nm)')
 with h5py.File(BASE/chosen['trajectory_path']) as h:x0=h['positions_m'][0];xf=h['positions_m'][-1]
 for a,x,title in [(ax[0,1],xf,'Original endpoint'),(ax[1,0],xf-xf.mean(0),'Translation removed'),(ax[1,1],m['aligned_positions_m'][-1],'Rigid motion removed')]:
  a.scatter((x0-x0.mean(0))[:,0]*1e6,(x0-x0.mean(0))[:,1]*1e6,s=20,facecolors='none',edgecolors='#777777',label='Initial centered')
  a.scatter(x[:,0]*1e6,x[:,1]*1e6,s=10,c=COLORS[0],label='Final');a.set_aspect('equal');a.set_title(title);a.set_xlabel('X (um)');a.set_ylabel('Y (um)')
 fig.suptitle('Translation, rotation and internal deformation | '+chosen['run_id'])
 save(fig,'Fig_F10AB_3_rigid_shape_separation',[chosen['trajectory_path']],manifest)
 fig,ax=plt.subplots(3,3,figsize=(13,10))
 for k in range(3):ax[0,0].plot(t,m['eigenvalues_3d_m2'][:,k]*1e12,label=f'lambda {k+1}')
 ax[0,0].set_title('3D covariance eigenvalues (um^2)');ax[0,0].legend()
 for k in range(2):ax[0,1].plot(t,m['eigenvalues_xy_m2'][:,k]*1e12)
 ax[0,1].set_title('XY covariance eigenvalues (um^2)')
 for key in ['Rg_3d_m','Rg_xy_m']:ax[0,2].plot(t,m[key]*1e6,label=key)
 ax[0,2].legend();ax[0,2].set_title('Rg (um)')
 for a,key in [(ax[1,0],'aspect_3d'),(ax[1,1],'deviatoric_stretch_3d'),(ax[1,2],'nonaffine_residual_3d_m')]:a.plot(t,m[key]);a.set_title(key)
 for k in range(3):ax[2,0].plot(t,m['principal_rms_lengths_3d_m'][:,k]*1e6,label=f'Axis {k+1}')
 ax[2,0].set_title('Principal RMS lengths (um)');ax[2,0].legend()
 ax[2,1].plot(t,m['affine_valid_3d'].astype(int),label='Affine valid');ax[2,1].plot(t,m['orientation_valid_3d'][:,0].astype(int),'--',label='Major axis identifiable');ax[2,1].set_ylim(-.05,1.2);ax[2,1].legend(fontsize=8)
 ax[2,2].plot(t,m['shear_norm_3d'],label='Symmetric shear');ax[2,2].plot(t,m['directional_xy_stretch_3d'],label='Uxx - Uyy');ax[2,2].legend(fontsize=8)
 for a in ax.ravel():a.set_xlabel('Physical time (ms)')
 fig.suptitle('Covariance, anisotropy and non-affine evolution')
 save(fig,'Fig_F10AB_4_shape_covariance',[chosen['trajectory_path']],manifest)
 fig,ax=plt.subplots(figsize=(12,max(5,len(formal)*.19)))
 for r in formal:
  mm=load(r);r['metrics'].update(com_x_nm=r['metrics']['com_um'][0]*1000,com_y_nm=r['metrics']['com_um'][1]*1000,
   rotation_deg=float(np.degrees(np.arccos(np.clip((np.trace(mm['rigid_rotation_3d'][-1])-1)/2,-1,1)))))
 keys=['com_x_nm','com_y_nm','rg_fraction','deviatoric','signed_xy','shear','rotation_deg','similarity_residual_nm'];arr=np.array([[np.nan if r['metrics'].get(k) is None else r['metrics'][k] for k in keys] for r in formal])
 scale=np.maximum(np.nanmax(abs(arr),axis=0),1e-12);im=ax.imshow(arr/scale,aspect='auto',cmap='coolwarm',vmin=-1,vmax=1)
 ax.set_yticks(np.arange(len(formal)),[r['state_id']+' '+r['action'] for r in formal],fontsize=7);ax.set_xticks(range(len(keys)),keys,rotation=20);fig.colorbar(im,ax=ax,label='Signed response / column maximum')
 fig.suptitle('Actual action-to-shape responses | White cells: invalid global affine fit')
 save(fig,'Fig_F10AB_5_action_response',[r['trajectory_path'] for r in formal],manifest)
 fig,ax=plt.subplots(1,3,figsize=(14,4));repeated=[r for r in rows if r['stage'] in ['discovery','confirmation'] and r['status']=='PASS' and r['variant']=='DT1']
 for j,key in enumerate(['signed_xy','shear','similarity_residual_nm']):
  for c,action in enumerate(sorted(set(r['action'] for r in repeated))):
   rr=[r for r in repeated if r['action']==action];ax[j].scatter([int(r['state_id'][1])+.06*c for r in rr],[r['metrics'][key] for r in rr],s=18,label=action,color=COLORS[c%len(COLORS)])
  ax[j].set_xlabel('Initial shape S0–S6');ax[j].set_title(key)
 ax[0].legend(fontsize=7);fig.suptitle('Initial-shape dependence | Each point is one trajectory, not one frame')
 save(fig,'Fig_F10AB_6_initial_shape_dependence',[r['trajectory_path'] for r in repeated],manifest)
 fig,ax=plt.subplots(2,3,figsize=(15,9));path=BASE/'analysis/SENSITIVITY.json'
 cases=[c for c in read(path).get('cases',[]) if 'normalized_jacobian' in c] if path.exists() else []
 if cases:
  for i,state in enumerate(['S2_r0','S4_r0']):
   case=next(c for c in cases if c['state_id']==state and c['time_s']==.002);j=np.asarray(case['normalized_jacobian']);im=ax[i,0].imshow(j,aspect='auto',cmap='coolwarm');fig.colorbar(im,ax=ax[i,0]);ax[i,0].set_xticks([0,1],['Amplitude','Relative phase']);ax[i,0].set_yticks(range(len(case['outputs'])),case['outputs'],fontsize=8);ax[i,0].set_title(state+' normalized response')
   for c in [c for c in cases if c['state_id']==state]:
    ax[i,1].semilogy([1,2],c['shape_singular_values'],'o-',label=f'{c["time_s"]*1000:g} ms: rank {c["shape_effective_rank"]}');ax[i,1].axhline(c['rank_threshold'],ls='--',alpha=.6)
   ax[i,1].set_xticks([1,2]);ax[i,1].set_title('Shape singular values; dashed: error bound');ax[i,1].legend(fontsize=8)
   note(ax[i,2],f'{state}: {case["status"]}\n\nTwo inputs only: common X-pair amplitude and relative phase.\n\nShape rows exclude COM and uniform Rg scale.\n\nFinite-difference discrepancy: {case["stencil_discrepancy"]:.3g}\nThreshold: {case["rank_threshold"]:.3g}\nNo global four-PZT controllability claim.')
 else:
  for a in ax.ravel():note(a,'NOT_RUN: insufficient sensitivity evidence')
 fig.suptitle('Local action sensitivity and identifiability')
 save(fig,'Fig_F10AB_7_shape_controllability',['analysis/SENSITIVITY.json'] if path.exists() else [],manifest)
 fig,ax=plt.subplots(1,3,figsize=(13,4))
 for c,state in enumerate(['S0_r0','S2_r0','S6_r0']):
  rr=next((r for r in formal if r['state_id']==state and r['action']=='A7_X210'),None)
  if not rr:continue
  mm=load(rr);tt=mm['time_s']*1e3
  for a,key in zip(ax,['Rg_3d_m','similarity_rmsd_3d_m','proximity_components']):a.plot(tt,mm[key],label=state,color=COLORS[c]);a.set_title(key);a.set_xlabel('Physical time (ms)')
 ax[0].legend();fig.suptitle('Aggregation versus deformation | Fixed members; proximity is not fusion')
 save(fig,'Fig_F10AB_8_aggregation_vs_deformation',[r['trajectory_path'] for r in formal if r['action']=='A7_X210'],manifest)
 fig,axes=plt.subplots(2,3,figsize=(14,8));ax=axes[0];path=BASE/'analysis/NUMERICAL_COMPARISONS.json';comp=read(path) if path.exists() else []
 comparison_inputs=[str(path.relative_to(BASE))]
 for name in ['CONFIRMATION_NUMERICAL','SENSITIVITY_NUMERICAL','EXTENSION_NUMERICAL']:
  p=BASE/'analysis'/(name+'.json')
  if p.exists():comp+=read(p);comparison_inputs.append(str(p.relative_to(BASE)))
 if (BASE/'analysis/TARGET.json').exists():
  td=read(BASE/'analysis/TARGET.json');comp+=td['comparisons']+[td['independent_recompute']];comparison_inputs.append('analysis/TARGET.json')
 write(BASE/'figures/CONVERGENCE_COLUMN_ORDER.json',[dict(index=i,reference=c['reference'],comparison=c['comparison']) for i,c in enumerate(comp)])
 if comp:
  for k,(key,label,fac) in enumerate([('fixed_id_m','Maximum fixed-ID difference (nm)',1e9),('Rg_m','Maximum Rg difference (nm)',1e9)]):
   vals=[r['original'].get('errors',{}).get(key,np.nan) for r in comp];ax[k].bar(range(len(vals)),np.array(vals)*fac,color=[COLORS[0] if 'DT05' in r['comparison'] else COLORS[1] for r in comp]);ax[k].set_title(label);ax[k].set_xlabel('Matched comparison (blue: half dt; red: cache/repeat)')
 else:
  note(ax[0],'NOT_RUN: numerical comparisons unavailable');note(ax[1],'No extrapolated error bounds')
 audits=[r['independent_audit'] for r in rows if r['status']=='PASS'];ax[2].plot([r.get('force_gate_fraction_max',np.nan) for r in audits],'.');ax[2].axhline(1,color='r');ax[2].set_title('Independent force residual / original limit');ax[2].set_xlabel('Trajectory')
 for a,key,title,limit in [(axes[1,0],'overlap_m','True geometric overlap (nm)',policy()['hard_geometry_m']),(axes[1,1],'actual_complementarity_m','True-gap complementarity (nm)',policy()['actual_geometry_complementarity_m'])]:
  a.plot([r.get(key,np.nan)*1e9 for r in audits],'.');a.axhline(limit*1e9,color='r',ls='--');a.set_title(title);a.set_xlabel('Completed trajectory')
 axes[1,2].scatter(range(len(rows)),[r['wall_s'] for r in rows],c=[COLORS[0] if r['status']=='PASS' else COLORS[1] for r in rows],s=10);axes[1,2].set_title('Actual trajectory wall time (s)');axes[1,2].set_xlabel('Red: stopped at original rejection gate')
 fig.suptitle('Numerical convergence and original acceptance limits')
 save(fig,'Fig_F10AB_9_numerical_convergence',comparison_inputs if comp else [],manifest)
 fig,ax=plt.subplots(2,2,figsize=(12,9));target=BASE/'analysis/TARGET.json'
 if target.exists():
  d=read(target)
  for k,c in enumerate(d['cases']):
   rr=next(r for r in rows if r['run_id']==c['run_id']);mm=load(rr);tt=mm['time_s']*1e3;label=c['state_id']+(' PASS' if c['admitted'] else ' FAIL')
   ax[0,1].plot(tt,100*(mm['aspect_3d']/mm['aspect_3d'][0]-1),label=label,color=COLORS[k]);ax[1,0].plot(tt,np.linalg.norm(mm['center_3d_m']-mm['center_3d_m'][0],axis=1)*1e9,label=c['state_id'],color=COLORS[k]);ax[1,1].plot(tt,mm['proximity_components'],label=c['state_id'],color=COLORS[k])
   if k==0:
    with h5py.File(BASE/rr['trajectory_path']) as h:x0=h['positions_m'][0];xf=h['positions_m'][-1]
    ax[0,0].scatter(x0[:,0]*1e6,x0[:,1]*1e6,s=40,facecolors='none',edgecolors='#888888',label='Initial');ax[0,0].scatter(xf[:,0]*1e6,xf[:,1]*1e6,s=14,color=COLORS[0],label='Actual endpoint');ax[0,0].set_aspect('equal');ax[0,0].set_xlabel('X (um)');ax[0,0].set_ylabel('Y (um)');ax[0,0].legend();ax[0,0].set_title('S2_r0 fixed-ID XY distribution')
  ax[0,1].axhspan(2.5,3.5,color='gray',alpha=.2,label='Target: +3% ±0.5 percentage point');ax[0,1].set_ylabel('3D aspect ratio change (%)');ax[0,1].legend(fontsize=8)
  ax[1,0].axhline(20,ls='--',color='gray',label='20 nm limit');ax[1,0].set_ylabel('COM drift (nm)');ax[1,0].legend(fontsize=8)
  ax[1,1].set_ylabel('Proximity components');ax[1,1].set_title('Connectivity is part of target admission');ax[1,1].legend(fontsize=8)
  for a in [ax[0,1],ax[1,0],ax[1,1]]:a.set_xlabel('Physical time (ms)')
 else:
  for a in ax.ravel():note(a,'TARGETED_SHAPE_CONTROL_NOT_RUN\nNo target success is inferred from Rg contraction alone.')
 fig.suptitle('Targeted transformation | X pair: 700 Pa-eq each, 180 degrees | Actual limitations')
 save(fig,'Fig_F10AB_10_target_or_limitations',[str(target.relative_to(BASE))] if target.exists() else [],manifest)
 write(BASE/'figures/MANIFEST.json',manifest)
if __name__=='__main__':main()
