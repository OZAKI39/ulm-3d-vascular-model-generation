from pathlib import Path
import csv,json,hashlib
import numpy as np,h5py
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from analyze_matrices import R,reference,C
plt.rcParams.update({'font.size':9,'figure.dpi':130})
def read(name):
 with (R/name).open() as f:r=list(csv.DictReader(f))
 for row in r:
  for k,v in row.items():
   try:row[k]=float(v)
   except (ValueError,TypeError):pass
 return r
def save(fig,name):fig.tight_layout();fig.savefig(R/(name+'.png'));plt.close(fig)
def placeholder(name,title,detail):
 fig,ax=plt.subplots(figsize=(9,4));ax.axis('off');ax.text(.5,.7,title,ha='center',va='center',fontsize=16);ax.text(.5,.43,'NOT RUN: PLANAR GATE FAILED',ha='center',color='firebrick',fontsize=14);ax.text(.5,.2,detail,ha='center',wrap=True);save(fig,name)
rows=read('PLANAR_CONVERGENCE.csv');phase=read('LATERAL_PHASE_SENSITIVITY.csv');anis=read('ANISOTROPY.csv');rough=read('FIRST_CONTACT.csv');leaks=read('LEAKAGE_SCREENING.csv');perf=read('PERFORMANCE_SCALING.csv')
base=[r for r in rows if r['category']=='PLANAR'];betas=sorted(set(r['beta'] for r in base),reverse=True);colors={b:plt.cm.tab10(i) for i,b in enumerate(betas)}
with h5py.File(R/'FIXED_MULTIBLOB_WALL_AUDIT.h5') as f:
 fig,axes=plt.subplots(3,4,figsize=(14,10))
 for col,la in enumerate(['HEX','SQUARE']):
  for ly in [1,2,3]:
   g=f[f'RAW_planar_beta0.5/PLANAR_{la}_L{ly}_R4_S1.01'];x=g['coords'][:];a=g['radii'][:]
   top,side=axes[ly-1,2*col:2*col+2]
   for p,aa in zip(x,a):
    top.add_patch(Circle(p[:2],aa,fill=False,alpha=.3));side.add_patch(Circle(p[[0,2]],aa,fill=False,alpha=.2))
   top.set(xlim=(-2,2),ylim=(-2,2),aspect='equal',title=f'{la}, {ly} layer(s): top',xlabel='x/a',ylabel='y/a');side.set(xlim=(-2,2),ylim=(-3.2,2.2),aspect='equal',title='side; a_w/a=0.5, s/(2a_w)=1.01',xlabel='x/a',ylabel='z/a');side.axhline(0,color='k',ls='--');side.add_patch(Circle((0,1.01),1,fill=False,color='red'));side.text(.15,1.01,'a',color='red')
 save(fig,'VIS_WALL_DISCRETIZATION')
chosen='BEAD'
fig,axes=plt.subplots(1,2,figsize=(12,4.5))
for ax,la in zip(axes,['HEX','SQUARE']):
 e=np.array(C['epsilon']);ax.loglog(e,[reference(v)[2,2] for v in e],'k-',lw=2,label='RMBW')
 for b in betas:
  for ly in [1,2,3]:
   rr=sorted([r for r in base if r['layout']==la and r['beta']==b and r['layers']==ly and r['phase']==chosen],key=lambda x:x['epsilon'])
   if not rr:continue
   ax.loglog([r['epsilon'] for r in rr],[r['normal'] for r in rr],ls=['--','-.',':'][ly-1],color=colors[b],label=f'beta={b}, L={ly}')
 ax.set(title=la+'; bead phase (all phases in CSV)',xlabel='h/a',ylabel='Rnormal / (6 pi mu a)');ax.legend(fontsize=7)
save(fig,'VIS_PLANAR_NORMAL_CONVERGENCE')
fig,axes=plt.subplots(2,2,figsize=(12,8))
for ax,key,ind in zip(axes.flat,['tangent','RR_parallel','RR_normal','TR'],[(0,0),(3,3),(5,5),(1,3)]):
 e=np.array(C['epsilon']);ax.semilogx(e,[reference(v)[ind] for v in e],'k-',lw=2,label='RMBW')
 for b in betas:
  for ly in [1,2,3]:
   rr=sorted([r for r in base if r['layout']=='HEX' and r['beta']==b and r['layers']==ly and r['phase']==chosen],key=lambda r:r['epsilon'])
   if not rr:continue
   ax.semilogx([r['epsilon'] for r in rr],[r[key] for r in rr],ls=['--','-.',':'][ly-1],color=colors[b],label=f'beta={b}, L={ly}')
 ax.set(title=key+'; HEX bead phase',xlabel='h/a',ylabel='scaled resistance');ax.legend(fontsize=7)
save(fig,'VIS_PLANAR_TANGENTIAL_ROTATION')
fig,ax=plt.subplots(figsize=(9,5))
for b in betas:
 for la in ['HEX','SQUARE']:
  for ly in [1,2,3]:
   rr=[r for r in base if r['beta']==b and r['layout']==la and r['layers']==ly]
   if not rr:continue
   e=sorted(set(r['epsilon'] for r in rr));val=[max(r['matrix_action_error'] for r in rr if r['epsilon']==v) for v in e];ax.semilogx(e,val,label=f'{la}, beta={b}, L={ly}')
for gate in [.05,.1]:ax.axhline(gate,ls='--',color='gray');ax.text(3,gate,f'{100*gate:g}%')
ax.set(xlabel='h/a',ylabel='max action relative error (256 q; 19 phases)',title='All required phases retained');ax.legend(fontsize=7,ncol=3);save(fig,'VIS_MATRIX_ACTION_ERROR')
fig,axes=plt.subplots(len(betas),3,figsize=(13,3.5*len(betas)))
for axrow,b in zip(axes,betas):
 rr=[r for r in base if r['beta']==b and r['layout']=='HEX' and r['layers']==1 and r['epsilon']==.01]
 for ax,key in zip(axrow,['normal_error','tangent_error','TR_error']):
  im=ax.scatter([r['x'] for r in rr],[r['y'] for r in rr],c=[r[key] for r in rr],s=100,cmap='magma');
  from mpl_toolkits.axes_grid1 import make_axes_locatable
  cax=make_axes_locatable(ax).append_axes('right',size='5%',pad=.07);fig.colorbar(im,cax=cax,label='relative error')
  step=2*b*C['spacing_ratio_primary'];v1=np.array([step,0.]);v2=np.array([step/2,np.sqrt(3)*step/2]);poly=np.array([np.zeros(2),v1,v1+v2,v2,np.zeros(2)]);ax.plot(poly[:,0],poly[:,1],color='gray',lw=.7)
  ax.set(title=f'{key}; beta={b}\nHEX L1, h/a=0.01',xlabel='x/a',ylabel='y/a',aspect='equal')
save(fig,'VIS_LATERAL_PHASE_SENSITIVITY')
fig,axes=plt.subplots(1,3,figsize=(13,4.5))
for ax,key in zip(axes,['TT','RR','TR']):
 for b in betas:
  for ph in ['BEAD','BRIDGE','PORE']:
   rr=[r for r in anis if r['beta']==b and r['layout']=='HEX' and r['layers']==1 and r['epsilon']==.01 and r['phase']==ph][0];ax.plot([0,45,90],[rr[f'{key}_{a}'] for a in [0,45,90]],'o-',label=f'beta={b}, {ph}')
 ax.set(title=f'{key}; HEX L1 h/a=.01',xlabel='direction (degree)',ylabel='scaled response');ax.legend(fontsize=7)
save(fig,'VIS_ARTIFICIAL_ANISOTROPY')
fig,axes=plt.subplots(1,2,figsize=(12,5))
for b in betas:
 for la in ['HEX','SQUARE']:
  rr=[r for r in rough if r['beta']==b and r['layout']==la and r['layers']==1 and r['phase'] in ['BEAD','BRIDGE','PORE']];axes[0].plot(range(3),[r['first_contact_nominal_h_over_a'] for r in rr],'o-',label=f'{la}, beta={b}')
axes[0].set(xticks=range(3),xticklabels=['BEAD','BRIDGE','PORE'],ylabel='first contact nominal h/a',title='Hard-sphere first contact; no fitted offset');axes[0].axhline(0,color='k',ls='--');axes[0].legend()
b=.5;s=2*b*1.01
for x in [-s,0,s]:axes[1].add_patch(Circle((x,-b),b,fill=False,color='steelblue'))
for i,x in enumerate([0,s/2]):
 z=max(-b+np.sqrt((1+b)**2-(x-xx)**2) for xx in [-s,0,s] if abs(x-xx)<=1+b);axes[1].add_patch(Circle((x,z),1,fill=False,color=['red','green'][i],ls=['-',':'][i]))
axes[1].axhline(0,color='k',ls='--');axes[1].set(xlim=(-2,2),ylim=(-1.2,2.3),aspect='equal',xlabel='x/a',ylabel='z/a',title='Side section: bead vs bridge contact');save(fig,'VIS_WALL_ROUGHNESS_AND_CONTACT')
fig,axes=plt.subplots(1,2,figsize=(12,4.5))
for ax,la in zip(axes,['HEX','SQUARE']):
 for b in betas:
  for sp in [1.01,1.05,1.1]:
   rr=sorted([r for r in leaks if r['layout']==la and r['beta']==b and r['spacing_ratio']==sp],key=lambda r:r['layers']);ax.plot([r['layers'] for r in rr],[r['transmission'] for r in rr],'o-',label=f'beta={b}, spacing={sp}')
 ax.axhline(.05,color='k',ls='--',label='frozen amplitude gate');ax.axhline(-.05,color='k',ls='--');ax.set_yscale('symlog',linthresh=1e-4);ax.set(title=la,xlabel='layers',ylabel='signed lower Vz / no-wall lower Vz',xticks=[1,2,3]);ax.legend(fontsize=7)
save(fig,'VIS_LEAKAGE_SCREENING')
fig,axes=plt.subplots(1,3,figsize=(13,4))
for ax,key in zip(axes,['reciprocity','min_eigenvalue','condition']):
 v=np.array([r[key] for r in rows]);ax.plot(np.arange(len(v)),np.maximum(v,1e-20),'.',ms=1);ax.set(yscale='log',xlabel='completed static query',ylabel=key)
save(fig,'VIS_MATRIX_PROPERTIES')
placeholder('VIS_CYLINDER_CURVATURE','Cylinder curvature test','No cylinder discretization or hydrodynamic result is represented as tested.')
placeholder('VIS_REAL_STL_BLOB_PATCHES','Real STL patch audit','Original STL is retained byte-for-byte as input; no patch or bubble site was tested.')
fig,axes=plt.subplots(1,3,figsize=(13,4))
valid=[r for r in perf if isinstance(r.get('median_query_seconds'),float)]
for ax,key,label in zip(axes,['wall_assembly_seconds','median_query_seconds','max_rss_bytes'],['wall assembly seconds','reused query seconds','process peak RSS GiB']):
 v=valid if key=='median_query_seconds' else [r for r in perf if isinstance(r.get('wall_assembly_seconds'),float) and r.get('max_rss_bytes',0)>0]
 v=sorted(v,key=lambda r:r['N_wall'])
 ax.plot([r['N_wall'] for r in v],[r[key]/(2**30 if key=='max_rss_bytes' else 1) for r in v],'o-');ax.set(xscale='log',yscale='log',xlabel='N_wall',ylabel=label)
fig.suptitle('CPU, 4 BLAS threads; N=2000 query failed; N=5000/10000 blocked by memory');save(fig,'VIS_PERFORMANCE_SCALING')
# Point-only VTP: glyph radius is explicit point data; no fabricated curved-wall results.
import vtk
from vtk.util.numpy_support import numpy_to_vtk

def vtp(name,coords,radii,layers,status):
 poly=vtk.vtkPolyData();pts=vtk.vtkPoints();pts.SetData(numpy_to_vtk(np.asarray(coords,dtype=float).reshape(-1,3),deep=True));poly.SetPoints(pts);verts=vtk.vtkCellArray()
 for i in range(len(coords)):verts.InsertNextCell(1);verts.InsertCellPoint(i)
 poly.SetVerts(verts);n=len(coords)
 for key,val in dict(object_type=np.zeros(n),blob_radius=np.asarray(radii),layer=np.asarray(layers),patch_id=np.zeros(n),bubble_id=np.full(n,-1),gap=np.full(n,np.nan),matrix_error=np.full(n,np.nan),validity_flags=np.zeros(n)).items():a=numpy_to_vtk(val.astype(float),deep=True);a.SetName(key);poly.GetPointData().AddArray(a)
 txt=vtk.vtkStringArray();txt.SetName('STATUS');txt.InsertNextValue(status);poly.GetFieldData().AddArray(txt);unit=vtk.vtkStringArray();unit.SetName('COORDINATE_UNITS');unit.InsertNextValue('meters; frozen d50 radius used for display');poly.GetFieldData().AddArray(unit)
 writer=vtk.vtkXMLPolyDataWriter();writer.SetFileName(str(R/name));writer.SetInputData(poly);assert writer.Write()==1
with h5py.File(R/'FIXED_MULTIBLOB_WALL_AUDIT.h5') as f:
 b=min(betas);g=f[f'RAW_planar_beta{b:g}/PLANAR_HEX_L1_R4_S1.01'];ar=C['radii_m']['d50'];vtp('PLANAR_WALL_BLOBS.vtp',g['coords'][:]*ar,g['radii'][:]*ar,g['layer'][:],'REPRESENTATIVE_FAILED_PLANAR_TEST; not selected production geometry')
for n in ['CYLINDER_WALL_BLOBS.vtp','REAL_STL_WALL_BLOBS.vtp','REAL_STL_SAMPLE_BUBBLES.vtp']:vtp(n,np.empty((0,3)),np.empty(0),np.empty(0),'NOT_RUN_PLANAR_GATE_FAILED')
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
p={'source_script_sha256':sha(Path(__file__)),'data_sha256':{p.name:sha(p) for p in [R/'FIXED_MULTIBLOB_WALL_AUDIT.h5']+list(R.glob('*.csv'))},'PNG_sha256':{p.name:sha(p) for p in R.glob('VIS_*.png')},'VTP_sha256':{p.name:sha(p) for p in R.glob('*.vtp')},'HUMAN_VISUAL_REVIEW':'PENDING','gate_blocked_visuals':'two explicitly marked PNG placeholders and three empty labelled VTP; not simulated evidence'};(R/'VISUALIZATION_PROVENANCE.json').write_text(json.dumps(p,indent=2)+'\n');print('PNG',len(p['PNG_sha256']),'VTP',len(p['VTP_sha256']))
