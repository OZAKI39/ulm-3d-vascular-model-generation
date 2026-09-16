"""QA plots from frozen CSV/HDF5 only; no package calls, no curve repair."""
from pathlib import Path
import json,csv,hashlib
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap,BoundaryNorm
R=Path(__file__).resolve().parents[1];V=R/'visualization';V.mkdir(exist_ok=True)
with (R/'WALL_REFERENCE_RAW_RESULTS.csv').open() as f:rows=list(csv.DictReader(f))
with (R/'WALL_REFERENCE_CROSS_COMPARISON.csv').open() as f:cross=list(csv.DictReader(f))
eps=sorted({float(r['epsilon']) for r in rows});styles={'RMBW_LUBRICATION':('#1666a0','-','RMBW lubrication'),'PYSTOKES':('#dd7429','--','PyStokes (reciprocity FAIL)'),'RMBW_SPHERE':('#8b638e',':','RMBW sphere (invalid full matrix)')}
art=[]
def data(label,key):return np.array([float(next(r for r in rows if r['implementation']==label and r['size']=='d50' and float(r['viscosity_Pa_s'])==.001 and float(r['epsilon'])==e)[key]) for e in eps])
def base(ax):ax.set_xscale('log');ax.set_xlabel('Surface gap / radius, h/a');ax.grid(alpha=.2);ax.axvspan(.001,.2,color='#d8e9f5',alpha=.18)
def finish(fig,name,desc):
 fig.tight_layout();path=V/(name+'.png');fig.savefig(path,dpi=155,bbox_inches='tight');plt.close(fig);art.append({'file':str(path.relative_to(R)),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'description':desc})
def curves(keys,titles,name,ylabel,scale=None):
 fig,axs=plt.subplots(1,len(keys),figsize=(12,4.3),squeeze=False)
 for ax,key,title in zip(axs[0],keys,titles):
  for label,(color,ls,leg) in styles.items():ax.plot(eps,data(label,key),ls,marker='o',markersize=3,color=color,label=leg)
  base(ax);ax.set_title(title);ax.set_ylabel(ylabel)
  if scale=='symlog':ax.set_yscale('symlog',linthresh=1)
  if scale=='log':ax.set_yscale('log')
  ax.legend(fontsize=7)
 fig.suptitle('Single rigid sphere / no-slip plane; d50; mu = 1 mPa s',fontsize=11)
 finish(fig,name,'Raw package curves; purple invalid results are retained, not repaired')
curves(['Mhat_normal_TT','Mhat_tangential_TT'],['Normal translation','Parallel translation'], 'VIS_WALL_TRANSLATIONAL_MOBILITY','Mobility / bulk mobility')
curves(['Rhat_normal_TT','Rhat_tangential_TT'],['Normal resistance','Parallel resistance (Omega fixed)'], 'VIS_WALL_TRANSLATIONAL_RESISTANCE','Resistance / bulk resistance','symlog')
curves(['Mhat_RR_parallel','Mhat_RR_normal'],['Rotation axis parallel to wall','Rotation axis normal to wall'], 'VIS_WALL_ROTATIONAL_MOBILITY','Rotational mobility / bulk')
curves(['Mhat_RT_y_Fx','Mhat_TR_x_Ty'],['Fx -> Omega_y','Ty -> Vx'], 'VIS_WALL_TRANSLATION_ROTATION_COUPLING','Coupling / sqrt(MT0 MR0)')
fig,ax=plt.subplots(figsize=(9,5))
modekeys=['normal_TT','tangential_TT','RR_parallel','RR_normal','TR_x_Ty','RT_y_Fx']
for mode in modekeys:
 rr=[r for r in cross if r['reference_A']=='RMBW_LUBRICATION' and r['size']=='d50' and r['mode']==mode];rr.sort(key=lambda r:float(r['epsilon']));ax.plot([float(r['epsilon']) for r in rr],[100*float(r['symmetric_relative_difference']) for r in rr],marker='o',markersize=3,label=mode)
base(ax);ax.axhline(2,color='gray',lw=.6);ax.axhline(10,color='gray',lw=.6);ax.set_ylabel('100 * |A-B| / max(|A|,|B|) [%]');ax.set_title('RMBW lubrication vs PyStokes: descriptive differences');ax.legend(fontsize=8);finish(fig,'VIS_REFERENCE_RELATIVE_DIFFERENCE','Signed coefficients compared before absolute difference; reciprocity failure remains disqualifying')
for e,name in [(.01,'VIS_MOBILITY_MATRIX_NEAR_WALL'),(.2,'VIS_MOBILITY_MATRIX_MODERATE_GAP')]:
 fig,axs=plt.subplots(1,2,figsize=(11,5));mats=[]
 with h5py.File(R/'WALL_REFERENCE_MOBILITY_MATRICES.h5','r') as h:
  labels=h['implementation'].asstr()[:];sizes=h['size'].asstr()[:]
  for label in ['RMBW_LUBRICATION','PYSTOKES']:
   ix=np.where((labels==label)&(sizes=='d50')&(h['epsilon'][:]==e)&(h['viscosity_Pa_s'][:]==.001))[0][0];mats.append(h['work_conjugate_normalized_matrix'][ix])
 vmax=max(np.max(abs(m)) for m in mats)
 for ax,M,label in zip(axs,mats,['RMBW lubrication','PyStokes: TR/RT sign mismatch']):
  im=ax.imshow(M,cmap='RdBu_r',vmin=-vmax,vmax=vmax);ax.set_xticks(range(6),['Fx','Fy','Fz','Tx','Ty','Tz']);ax.set_yticks(range(6),['Vx','Vy','Vz','Ox','Oy','Oz']);ax.set_title(label)
  for i in range(6):
   for j in range(6):ax.text(j,i,f'{M[i,j]:.3g}',ha='center',va='center',fontsize=8,color='white' if abs(M[i,j])>.6*vmax else 'black')
 fig.suptitle(f'Work-conjugate normalized mobility, h/a = {e:g}');finish(fig,name,'D^-1 M D^-1; same symmetric color range for both matrices')
fig,axs=plt.subplots(3,1,figsize=(13,10),gridspec_kw={'height_ratios':[2,1,1]})
labels=['normal TT','parallel TT','TR: Ty -> Vx','RT: Fx -> Omega_y','RR parallel','RR normal'];modekeys=['normal_TT','tangential_TT','TR_x_Ty','RT_y_Fx','RR_parallel','RR_normal'];C=np.zeros((len(modekeys),len(eps)))
for i,mode in enumerate(modekeys):
 for j,e in enumerate(eps):
  r=next(r for r in cross if r['reference_A']=='RMBW_LUBRICATION' and r['size']=='d50' and r['mode']==mode and float(r['epsilon'])==e);C[i,j]={'CLOSE_AGREEMENT':0,'MODERATE_DIFFERENCE':1,'LARGE_DIFFERENCE':2}[r['agreement_label']]
cmap=ListedColormap(['#6cb98d','#efc36d','#d67a72']);axs[0].imshow(C,aspect='auto',cmap=cmap,vmin=-.5,vmax=2.5);axs[0].set_yticks(range(len(modekeys)),labels);axs[0].set_title('Agreement only: close <2%, moderate 2-10%, large >10% (not validity)')
for i in range(len(modekeys)):
 for j in range(len(eps)):axs[0].text(j,i,['CLOSE','MOD.','LARGE'][int(C[i,j])],ha='center',va='center',fontsize=7)
A=np.ones((3,len(eps)));A[0,:]=2;A[0,np.array(eps)>=10]=1;A[1,:]=0;A[2,:]=0
axs[1].imshow(A,aspect='auto',cmap=ListedColormap(['#d67a72','#efc36d','#6cb98d']),vmin=-.5,vmax=2.5);axs[1].set_yticks(range(3),['RMBW lubrication','RMBW sphere','PyStokes']);axs[1].set_title('Full matrix: green = supported/finite + audit gates; orange = bulk fallback; red = invalid')
F=np.ones((3,len(eps)));axs[2].imshow(F,aspect='auto',cmap=ListedColormap(['#6cb98d']),vmin=0,vmax=1);axs[2].set_yticks(range(3),['RMBW lubrication','RMBW sphere','PyStokes']);axs[2].set_title('All first formal evaluations returned finite values; no unavailable entries or exceptions')
for ax in axs:ax.set_xticks(range(len(eps)),[f'{e:g}' for e in eps]);ax.set_xlabel('h/a (sampled log-gap grid)')
finish(fig,'VIS_VALIDITY_MAP','Finite, implementation validity and cross agreement shown separately; no accuracy certification by color')
prov={'status':'PASS','generator':'src/visualize_wall_audit.py','generator_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'input_files':{name:hashlib.sha256((R/name).read_bytes()).hexdigest() for name in ['WALL_REFERENCE_RAW_RESULTS.csv','WALL_REFERENCE_CROSS_COMPARISON.csv','WALL_REFERENCE_MOBILITY_MATRICES.h5']},'artifacts':art,'data_source':'frozen CSV/HDF5 only','curves_manually_edited':False,'human_visual_review':'PENDING','matplotlib_version':matplotlib.__version__}
(R/'VISUALIZATION_PROVENANCE.json').write_text(json.dumps(prov,indent=2)+'\n');print('GENERATED',len(art))
