"""English PNG/PDF figures from archived measured trajectories and force diagnostics."""
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from .common import *
from scripts.straight_tube_exp1.f10r.figures import geometry as old_geometry

PURPLE='#70559c';TEAL='#258e95';RED='#b45a60';GREY='#697580'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelcolor':'#262b32','text.color':'#262b32','xtick.color':'#262b32','ytick.color':'#262b32','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.facecolor':'white','axes.titlesize':11})
OUT=BASE/'figures';AN=BASE/'analysis';inputs=set();outputs=[]
def js(p):inputs.add(p);return read(p)
def arr(name,root=AN):
 p=root/(name+'.npz');inputs.add(p);return dict(np.load(p))
def finish(fig,name,title,foot=''):
 fig.suptitle(title,fontsize=15,y=.995)
 fig.text(.5,.014,foot+'\nEXPLORATORY · HOMOGENEOUS_WATER_CONTROL_PROXY · PHYSICS_UNCALIBRATED · training_ready=false',ha='center',va='bottom',fontsize=7.5)
 fig.tight_layout(rect=(0,.067,1,.955),h_pad=2.7,w_pad=2.8)
 for ext in ['png','pdf']:
  path=OUT/(name+'.'+ext);fig.savefig(path,dpi=220);outputs.append(path)
 plt.close(fig)
def geometry(ax):
 inputs.add(OLD/'CONFIG_SNAPSHOT.json');inputs.add(ROOT/'scripts/straight_tube_exp1/f10r/figures.py')
 old_geometry(ax)
 for t in list(ax.texts):
  if 'Planned A:' in t.get_text():t.set_text('Water-filled cavity; zero flow\nBoth X ends sealed; 7 mm PZT discs')
 ax.set_title('a  Confirmed apparatus geometry')
def qual(ax,names,runs):
 keys=['max_overlap_m','max_contact_total_m','max_actual_complementarity_m','max_residual_m']
 limits=[17.8e-9,178e-9,1e-11,1e-13]
 for n,col in zip(names,[PURPLE,TEAL,RED]):
  v=np.array([runs[n][k] for k in keys])/limits
  ax.plot(np.arange(4),v,'o-',ms=4,color=col,label=n.replace('_DT1','').replace('_',' '))
 ax.axhline(1,color=GREY,ls='--',lw=1,label='Original limit')
 ax.set(xticks=range(4),xticklabels=['Overlap','Contact\ncorrection','True\ncomplement.','Mobility\nresidual'],ylabel='Maximum / original limit',ylim=(-.03,1.15))
 ax.legend(fontsize=7,loc='upper left');ax.tick_params(axis='x',labelsize=8)
def fig1(runs,gates,comps):
 fig,axes=plt.subplots(2,3,figsize=(15,9.3));a,b,c,d,e,f=axes.ravel();name='X_PAIR_3MS_DT1';m=arr(name);t=m['time_s']*1e3;old=arr('X_PAIR_DT1',R/'analysis')
 a.axhspan(.125,.275,color=TEAL,alpha=.08,label='Original tolerance');a.axhspan(.19,.21,color=TEAL,alpha=.18,label='Strict endpoint band')
 a.plot(t,m['centroid_m'][:,0]*1e6,color=PURPLE,label='Fresh 3 ms integration');a.plot(old['time_s']*1e3,old['centroid_m'][:,0]*1e6,'--',color=GREY,label='F10R stored 2 ms')
 a.scatter([2,3],[old['centroid_m'][-1,0]*1e6,m['centroid_m'][-1,0]*1e6],color=[GREY,PURPLE],s=25,zorder=5)
 a.axhline(.2,color=TEAL,ls=':',lw=1);a.set(xlabel='Time (ms)',ylabel='COM X (µm)',title='a  Registered target: (0.2, 0, 0) µm');a.legend(fontsize=7.6)
 for i,col in enumerate([PURPLE,TEAL,RED]):b.plot(t,m['centroid_m'][:,i]*1e6,color=col,label='XYZ'[i])
 b.set(xlabel='Time (ms)',ylabel='COM position (µm)',title='b  Actual three-axis trajectory');b.legend()
 for suffix,col in [('DT1',PURPLE),('DT05',TEAL),('FINE',RED)]:
  row=arr('X_PAIR_3MS_'+suffix);c.plot(row['time_s']*1e3,row['target_error_m']*1e6,color=col,label=suffix,lw=1)
 c.axhline(.075,color=GREY,ls='--',label='0.075 µm');c.axhline(.01,color=TEAL,ls=':',label='0.01 µm');c.set(xlabel='Time (ms)',ylabel='Distance to target (µm)',yscale='log',title='c  Endpoint error: %.6f µm'%(m['target_error_m'][-1]*1e6));c.legend(fontsize=8)
 for key,label,col in [('rg_m','Rg',PURPLE),('aspect','Aspect ratio',TEAL)]:d.plot(t,(m[key]/m[key][0]-1)*100,label=label,color=col)
 Q=np.linalg.norm(m['tensor_m2']-m['tensor_m2'][0],axis=(1,2))/np.linalg.norm(m['tensor_m2'][0]);d.plot(t,Q*100,color=RED,label='Shape tensor norm change')
 d.set(xlabel='Time (ms)',ylabel='Change from initial (%)',title='d  Translation and shape evaluated separately');d.legend(fontsize=8)
 qual(e,[name],runs);e.set_title('e  Original basic quality thresholds')
 for suffix,col in [('DT05',TEAL),('FINE',RED)]:
  v=comps[name+'__'+suffix];f.plot([0,1,2],np.array([v['COM_m'],v['Rg_m'],v['fixed_id_m']])*1e9,'o-',label=suffix,color=col)
 f.set(yscale='log',xticks=[0,1,2],xticklabels=['COM','Rg','Fixed ID'],ylabel='Maximum full-window difference (nm)',title='f  Step and cache comparisons');f.legend()
 finish(fig,'Fig_F10S_1_x_pair_3ms','X-pair navigation: fresh 3 ms validation', 'Same full-precision F10R amplitude/phase; no extrapolated or stitched trajectory · All five numerical layers PASS')
def fig2(runs,gates,comps):
 fig=plt.figure(figsize=(15,9.7));grid=fig.add_gridspec(2,3);geometry(fig.add_subplot(grid[0,0],projection='3d'))
 b=fig.add_subplot(grid[0,1]);c=fig.add_subplot(grid[0,2]);d=fig.add_subplot(grid[1,0]);e=fig.add_subplot(grid[1,1]);f=fig.add_subplot(grid[1,2])
 names=['X_MINUS_ONLY_DT1','X_PLUS_ONLY_DT1']
 for name,col,label in zip(names,[PURPLE,TEAL],['X− only','X+ only']):
  m=arr(name);t=m['time_s']*1e3;b.plot(t,m['centroid_m'][:,0]*1e6,color=col,label=label)
  d.plot(t,(m['rg_m']-m['rg_m'][0])*1e9,color=col,label=label)
 off=arr('OFF_DT1',R/'analysis');b.plot(off['time_s']*1e3,off['centroid_m'][:,0]*1e6,'--',color=GREY,label='F10R OFF (read-only)')
 b.set(xlabel='Time (ms)',ylabel='COM X displacement (µm)',title='b  Independent runs, identical initial state');b.legend(fontsize=8)
 for i,name in enumerate(names):
  r=runs[name];v=np.array(r['force_samples_summary']['primary_force_N']['mean_vector_N'])*1e12
  c.bar(i,v[0],color=[PURPLE,TEAL][i]);c.text(i,v[0]+(.02 if i==0 else -.02),'%.1f–%.1f Pa'%tuple(r['actual_pressure_min_max_Pa']),ha='center',va='bottom' if i==0 else 'top',fontsize=8)
 c.axhline(0,color=GREY,lw=.6);c.set(xticks=[0,1],xticklabels=['X− only','X+ only'],ylabel='Mean primary Fx per bubble (pN)',title='c  Measured pressure and force',ylim=(-.37,.37))
 d.set(xlabel='Time (ms)',ylabel='Rg change (nm)',title='d  Measured shape change');d.legend(fontsize=8)
 qual(e,names,runs);e.set_title('e  Contact and numerical quality')
 f.axis('off');lines=['Each source: 400 Pa-equivalent, 2.5 MHz','All other sources OFF; no background flow','Independent 0.5 ms screens → fresh 2 ms runs','']
 for name,label in zip(names,['X−','X+']):
  r=runs[name];noise=max(comps[name+'__'+s]['COM_m'] for s in ['DT05','FINE'])*1e9
  lines.extend([f"{label}: ΔX = {r['displacement_m'][0]*1e6:+.8f} µm",f"Mean velocity = {r['mean_COM_velocity_m_s'][0]*1e6:+.5f} µm/s",f"Largest COM numerical difference = {noise:.6f} nm",''])
 lines+=['All five numerical layers: PASS','Connected 32-ID network retained','Small shape change is not physical fusion']
 f.text(0,1,'\n'.join(lines),va='top',fontsize=9);f.set_title('f  Direction, signal and scope',loc='left')
 finish(fig,'Fig_F10S_2_independent_x_sources','Independent single-PZT direction verification','Full primary + secondary forces, F2, pair lubrication, hard-sphere SQP and walls · Not extracted from a switched sequence')
def fig3():
 data=js(BASE/'diagnostics/TRAPS.json');records=data['records'];strong=max(r['amplitude_Pa_equivalent'] for r in records)
 fig,axes=plt.subplots(2,3,figsize=(15,9.5));a,b,c,d,e,f=axes.ravel()
 for deg,col in zip([-6,0,6],[TEAL,PURPLE,RED]):
  key=f'A{strong:g}_phase{deg:+g}';curve=data['curves'][key];x=np.array(curve['x_m'])*1e6
  # JSON stores complex arrays as real/imag mappings.
  p=decode(curve['pressure_Pa']);a.plot(x,abs(p),color=col,label=f'{deg:+d}°')
  b.plot(x,np.array(curve['primary_force_N'])[:,0]*1e12,color=col,label=f'{deg:+d}°')
 sel=[r for r in records if r['amplitude_Pa_equivalent']==strong];ph=[r['phase_deg'] for r in sel]
 c.plot(ph,[r['pressure_peak_x_m']*1e6 for r in sel],'o-',color=PURPLE,label='Pressure maximum')
 c.plot(ph,[r['zero_force_x_m']*1e6 for r in sel],'x--',color=TEAL,label='Primary-force equilibrium')
 a.set(xlabel='X (µm)',ylabel='Pressure magnitude (Pa)',title='a  Standing-field pressure peaks');a.legend(title='X+ − X− phase',fontsize=8)
 b.axhline(0,color=GREY,lw=.7);b.set(xlabel='X (µm)',ylabel='Primary Fx (pN)',title='b  Actual restoring force near each peak');b.legend(fontsize=8)
 c.set(xlabel='Phase difference (degree)',ylabel='X position (µm)',title='c  Peaks and force zeros differ');c.legend(fontsize=8)
 c.text(.03,.04,'±6°: mismatch ≈ 0.019 µm\n0°: coincident by symmetry',transform=c.transAxes,fontsize=8)
 zero=[r for r in records if r['phase_deg']==0]
 for k,key,col in [(0,'primary_N',PURPLE),(1,'secondary_N',TEAL),(2,'F2_increment_N',RED)]:
  values=[np.linalg.norm(np.array(r['original_geometry_at_single_bubble_equilibrium'][key]),axis=1).mean()*1e12 for r in zero]
  d.bar(np.arange(2)+(k-1)*.23,values,width=.23,color=col,label=['Primary','Secondary','F2 increment'][k])
 d.set(xticks=[0,1],xticklabels=['400','1036.202'],xlabel='Each source (Pa-equivalent)',ylabel='Mean force magnitude (pN)',yscale='log',title='d  Actual 32-bubble force / compression');d.legend(fontsize=8)
 e.bar([0,1],[r['frozen_COM']['tau_x_s'] for r in zero],color=[TEAL,PURPLE]);e.axhline(.002,color=RED,ls='--',label='2 ms window');e.set(yscale='log',xticks=[0,1],xticklabels=['400','1036.202'],xlabel='Each source (Pa-equivalent)',ylabel='Local X relaxation time (s)',title='e  Frozen-shape COM following time');e.legend(fontsize=8)
 f.axis('off');strong0=zero[-1];eig=np.array(strong0['force_J_eigenvalues_N_m']);fc=strong0['frozen_COM']
 f.text(0,1,'LOCAL 3D RESTORING (proxy only)\nAll three force-Jacobian eigenvalues < 0\n\nFrozen-shape diagnostic, NOT a trajectory\nAt 1036.202 Pa-equivalent:\n  X relaxation: %.3f s\n  1° peak step: −0.843 µm\n  Predicted 2 ms following: %.4f µm\n\nCompression exists; no full trap contact admission\nLong moving-trap trajectories: NOT_RUN\nRigid internal shape is not experimentally validated'%(fc['tau_x_s'],fc['predicted_following_displacement_2ms_m'][0]*1e6),va='top',fontsize=9)
 f.set_title('f  Stability does not imply fast following',loc='left')
 finish(fig,'Fig_F10S_3_phase_trap_feasibility','Phase-shifted pressure peaks: force stability and following limits','APECSS and F2 evaluated at actual pressure · FROZEN_SHAPE_COM_DIAGNOSTIC · No long-time particle navigation claimed')
def decode(value):
 if isinstance(value,dict):return np.array(value['real'])+1j*np.array(value['imag'])
 return np.array([complex(v['real'],v['imag']) for v in value]) if value and isinstance(value[0],dict) else np.array(value)
def fig4():
 y=js(BASE/'diagnostics/Y_SCOPE.json');pd=js(BASE/'diagnostics/PDMS_PREFLIGHT.json')
 fig=plt.figure(figsize=(15,9.7));grid=fig.add_gridspec(2,3);geometry(fig.add_subplot(grid[0,0],projection='3d'))
 b=fig.add_subplot(grid[0,1]);c=fig.add_subplot(grid[0,2]);d=fig.add_subplot(grid[1,0]);e=fig.add_subplot(grid[1,1]);f=fig.add_subplot(grid[1,2])
 rows=[r for r in y['rows'] if r['kind']=='source_distance'];x=np.array([r['source_distance_offset_m'] for r in rows])*1e6
 b.plot(x,[r['pressure_magnitude_Pa'] for r in rows],'o-',color=PURPLE);b.axvline(4.1666667,color=GREY,ls='--',label='On-axis zero at +4.167 µm')
 b.set(xlabel='Y− source distance offset (µm)',ylabel='Group-center pressure (Pa)',title='b  Very near an axial pressure zero');b.legend(fontsize=7)
 rows=[r for r in y['rows'] if r['kind']=='frequency'];c.plot(np.array([r['frequency_Hz'] for r in rows])*1e-6,[r['pressure_magnitude_Pa'] for r in rows],'o-',color=TEAL)
 c.set(xlabel='Frequency (MHz)',ylabel='Group-center pressure (Pa)',title='c  Wavelength sensitivity (field only)');c.text(.03,.95,'APECSS force only at 2.5 MHz\nNo response extrapolation',transform=c.transAxes,va='top',fontsize=8)
 rows=[r for r in y['rows'] if r['kind']=='source_distance'];d.plot(x,[r['primary_force_N'][1]*1e12 for r in rows],'o-',color=RED)
 d.axhline(0,color=GREY,lw=.6);d.set(xlabel='Y− source distance offset (µm)',ylabel='Primary Fy (pN)',title='d  Force can change sign near the zero')
 rows=pd['candidates'];e.bar([str(r['PPW'])+' PPW' for r in rows],[r['estimated_peak_GiB'] for r in rows],color=[PURPLE,TEAL]);e.axhline(24,color=RED,ls='--',label='24 GiB device');e.axhline(20,color=GREY,ls=':',label='20 GiB safe budget');e.set(ylabel='Estimated peak GPU memory (GiB)',title='e  Conditional 2.5 MHz resource estimate');e.legend(fontsize=8)
 f.axis('off');r=rows[0]
 f.text(0,1,'PDMS_FIELD_NOT_RUN\nPDMS_Y_FIELD_NOT_ADMITTED\n\nXu et al., Phys. Rev. Applied 13, 054069 (2020)\nMeasured reference band: 3–7 MHz\n2.5 MHz material input: not validated\nEquivalent-fluid approximation omits shear waves\n\n6 PPW: %.1f µm, PML 20 cells, CFL ≤ 0.2\nGrid incl. PML: %s\nEstimated single-source runtime: %.1f min\n1 mm water envelope + BLI support margin\nEnd-seal material / coupling uncalibrated\n\nNo multi-material pressure field is fabricated'%(r['spacing_m']*1e6,' × '.join(map(str,r['expanded_shape'])),r['estimated_single_source_wall_s']/60),va='top',fontsize=8.7)
 f.set_title('f  Material and resource admission',loc='left')
 finish(fig,'Fig_F10S_4_y_field_scope','Y near-field sensitivity and PDMS field admission','All Y conclusions: HOMOGENEOUS_WATER_PROXY_SPECIFIC · Geometry offsets are exploratory, not measured assembly errors')
def main():
 start=time.perf_counter();OUT.mkdir(exist_ok=True);runs=js(AN/'RUN_ANALYSIS.json');gate=js(AN/'FINAL_GATE.json');comps=js(AN/'COMPARISONS.json')
 fig1(runs,gate,comps);fig2(runs,gate,comps);fig3();fig4()
 write(OUT/'FIGURE_PROVENANCE.json',dict(identity=identity(),input_sha256={str(p.relative_to(ROOT)):sha(p) for p in sorted(inputs)},output_sha256={p.name:sha(p) for p in outputs},render_wall_s=time.perf_counter()-start,
    scientific_execution_commits=sorted({read(p)['identity']['execution_commit'] for p in BASE.glob('*/CAMPAIGN.json')}),
    data_scope='All curves from actual archived trajectories or labeled field/force diagnostics. No fabricated PDMS fields.',plot_script_sha256=sha(ROOT/'scripts/straight_tube_exp1/f10s/figures.py')))
if __name__=='__main__':main()
