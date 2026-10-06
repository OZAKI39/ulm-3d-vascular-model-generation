"""Display existing local/author comparisons only; never execute either solver."""
import csv
import hashlib
import json
from itertools import product
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.text import Text
from PIL import Image

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];OUT=HERE/'output'
B='verification/pfpd_baqp_finite_domain_boundary_consistency_v1/'
A=B+'author_L075/'
SOURCES={'state':B+'L075/AFFINE_STATE.npz','case':B+'L075/SUMMARY.json','comparison':A+'COMPARISON.json',
         'X':A+'input_X.csv','x':A+'input_x.csv','volume':A+'input_volume.csv','parameters':A+'input_parameters.csv',
         'ids':A+'julia_ids.csv','force':A+'julia_force.csv','run':A+'RUN.log',
         'python_adapter':'scripts/finite_domain_author_oracle.py','julia_adapter':'scripts/finite_domain_author_oracle.jl',
         'local_kernel':'pd_clot/pfpd_candidate/kernel.py','local_geometry':'pd_clot/pfpd_candidate/boundary_force_diagnostics.py',
         'status':B+'SUMMARY.json'}
BLUE='#0072B2';ORANGE='#D55E00';GREEN='#009E73';GRAY='#a9b6c0';INK='#233746';MAG=50
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':11,
 'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,
 'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':600,
 'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(k):return json.loads((ROOT/SOURCES[k]).read_text())
def arr(k,**kwargs):return np.loadtxt(ROOT/SOURCES[k],delimiter=',',**kwargs)
def table(name,rows):
 with (OUT/'tables'/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def load():
 q=read('case');comparison=read('comparison')
 with np.load(ROOT/SOURCES['state'],allow_pickle=False) as z:
  X=z['X'];x=z['x'];volume=z['volume'];local=z['force'];core=z['full_support']
 author_X=arr('X');author_x=arr('x');author_V=arr('volume');params=arr('parameters')
 input_checks={name:bool(np.array_equal(a,b)) for name,a,b in [('X',X,author_X),('x',x,author_x),('volume',volume,author_V)]}
 input_checks['horizon']=bool(params[0]==q['delta_m'])
 input_checks['mu']=bool(params[1]==q['material']['mu_Pa'])
 input_checks['bulk']=bool(params[2]==q['material']['bulk_Pa'])
 assert all(input_checks.values())
 ids=arr('ids',dtype=int)-1;raw=arr('force');N=len(X)
 assert X.shape==x.shape==local.shape==raw.shape==(comparison['particles'],3)
 assert np.array_equal(np.sort(ids),np.arange(N))
 author=np.empty_like(raw);author[ids]=raw
 delta=local-author;dn=np.linalg.norm(delta,axis=1)
 ln=np.linalg.norm(local,axis=1);an=np.linalg.norm(author,axis=1)
 relative=float(np.linalg.norm(delta)/np.linalg.norm(local))
 core_relative=float(np.linalg.norm(delta[core])/np.linalg.norm(local[core]))
 assert np.isclose(relative,comparison['force_relative_error'],rtol=1e-14,atol=0)
 assert np.isclose(core_relative,comparison['full_support_core_force_relative_error'],rtol=1e-14,atol=0)
 assert comparison['passed'] and comparison['phase']=='DISABLED' and comparison['integrations_run']==0
 assert np.all(dn>0), 'Log view requires positive archived differences; do not insert an artificial floor.'
 assert comparison['strain']==q['epsilon'] and N==q['particles']
 assert read('status')['status']['ROOT_CAUSE_PRIMARY']=='INCONCLUSIVE'
 # The wrapper exports nodal forces by multiplying the author's force density by nodal volumes.
 wrapper=(ROOT/SOURCES['julia_adapter']).read_text()
 assert 'storage.b_int.*transpose(system.volume)' in wrapper
 assert 'quadratic_kernel(delta,length)=(1-length/delta)^2' in wrapper
 assert 'monomial=:C1,epsilon=0' in wrapper
 inverse=np.empty_like(ids);inverse[ids]=np.arange(N)
 display=X+MAG*(x-X)
 rows=[]
 for i in range(N):
  r=dict(particle_ID=i,author_row_0based=int(inverse[i]),X_m=X[i,0],Y_m=X[i,1],Z_m=X[i,2],
         x_m=x[i,0],y_m=x[i,1],z_m=x[i,2],nodal_volume_m3=volume[i],full_family=bool(core[i]))
  for j,axis in enumerate('xyz'):
   r[f'local_f{axis}_N']=local[i,j];r[f'author_f{axis}_N']=author[i,j];r[f'difference_f{axis}_N']=delta[i,j]
  r.update(local_force_norm_N=ln[i],author_force_norm_N=an[i],absolute_vector_difference_N=dn[i])
  rows.append(r)
 table('all_particle_force_comparisons.csv',rows)
 table('common_input_identity.csv',[dict(quantity=k,exact_match=v) for k,v in input_checks.items()])
 metrics=dict(particles=N,components_compared=N*3,relative_vector_force_error=relative,core_relative_vector_force_error=core_relative,
              maximum_absolute_particle_difference_N=float(dn.max()),maximum_local_particle_force_N=float(ln.max()),
              maximum_difference_particle_ID=int(dn.argmax()),full_family_particles=int(core.sum()),author_commit=comparison['commit'])
 table('comparison_metrics.csv',[metrics])
 return dict(q=q,comparison=comparison,X=X,x=x,volume=volume,display=display,local=local,author=author,
             core=core,difference=delta,difference_norm=dn,local_norm=ln,author_norm=an,metrics=metrics,input_checks=input_checks)

def edges(points):
 lo=points.min(0);hi=points.max(0)
 for axis in range(3):
  other=[j for j in range(3) if j!=axis]
  for choices in product((0,1),repeat=2):
   e=np.zeros((2,3));e[:,axis]=[lo[axis],hi[axis]]
   for j,choice in zip(other,choices):e[:,j]=hi[j] if choice else lo[j]
   yield e

def cloud(ax,d,color=BLUE):
 X=d['X']*1e6;x=d['display']*1e6;half=d['q']['L_m']*1e6/2
 ax.scatter(x[:,0],x[:,1],x[:,2],s=5,color=color,alpha=.58,depthshade=False,edgecolors='none')
 for e in edges(X):ax.plot(*e.T,color=GRAY,lw=.9,ls='--',alpha=.9)
 for e in edges(x):ax.plot(*e.T,color=color,lw=.7,alpha=.65)
 top=x[:,2].max();ax.set(xlim=(-half,half),ylim=(-half,half),zlim=(-half,top+30),
        xlabel=r'$x$ ($\mu$m)',ylabel=r'$y$ ($\mu$m)',zlabel=r'$z$ ($\mu$m)')
 ax.set_xticks([-half,0]);ax.set_yticks([0,half]);ax.set_zticks([-half,0,half])
 ax.set_box_aspect((2*half,2*half,top+half+30));ax.view_init(elev=21,azim=-52);ax.set_proj_type('ortho')
 for axis in (ax.xaxis,ax.yaxis,ax.zaxis):axis.pane.fill=False;axis.pane.set_edgecolor('#eeeeee');axis.labelpad=1
 ax.grid(False);ax.tick_params(pad=0,labelsize=7)
 ax.text2D(1,.99,f'$u$ × {MAG}\n(display only)',transform=ax.transAxes,ha='right',va='top',fontsize=8,color=ORANGE)

def box(fig,x,y,w,h,heading,detail,color):
 fig.add_artist(FancyBboxPatch((x,y),w,h,transform=fig.transFigure,boxstyle='round,pad=0.005,rounding_size=.008',
                             facecolor='#f7f9fa',edgecolor=color,linewidth=1))
 fig.text(x+w/2,y+h*.66,heading,ha='center',va='center',fontsize=10,color=color,weight='bold')
 fig.text(x+w/2,y+h*.27,detail,ha='center',va='center',fontsize=8,color=INK)

def arrow(fig,a,b):
 fig.add_artist(FancyArrowPatch(a,b,transform=fig.transFigure,arrowstyle='->',color=GRAY,mutation_scale=12,linewidth=1))

def parity(ax,d,compact=False):
 # Norm parity is paired with the full vector difference panel to retain directional information.
 x=d['author_norm']*1e9;y=d['local_norm']*1e9
 hi=float(max(x.max(),y.max()))*1.055
 ax.plot([0,hi],[0,hi],ls='--',color=GRAY,lw=1,label='Expected equality')
 ax.scatter(x,y,s=12 if compact else 18,facecolors='none',edgecolors=BLUE,linewidths=.6,alpha=.65,label='All particle IDs')
 ax.set(xlim=(-.4,hi),ylim=(-.4,hi),xlabel='Author force magnitude (nN)',ylabel='Local magnitude (nN)')
 ax.set_xticks([0,5,10,15]);ax.set_yticks([0,5,10,15]);ax.grid(alpha=.17)
 ax.legend(frameon=False,loc='upper left',fontsize=6.5 if compact else 8)

def differences(ax,d,compact=False):
 N=d['metrics']['particles'];v=d['difference_norm']
 ax.scatter(np.arange(N),v,s=3 if compact else 5,color=ORANGE,alpha=.7,edgecolors='none')
 ax.set(xlim=(-35,N+35),yscale='log',xlabel='Particle ID',ylabel=r'$\|\mathbf{f}_{local}-\mathbf{f}_{author}\|$ (N)')
 ax.set_xticks([0,500,1000,1500,2000]);ax.grid(axis='y',alpha=.17)
 ax.tick_params(labelsize=7 if compact else 8)


def save(fig,stem):
 fig.canvas.draw();texts=[x.get_text() for x in fig.findobj(Text) if x.get_text()]
 assert not any(any('\u3400'<=c<='\u9fff' for c in s) for s in texts)
 for ext in ('png','pdf','svg'):fig.savefig(OUT/'figures'/f'{stem}.{ext}',dpi=600,facecolor='white')
 im=Image.open(OUT/'figures'/f'{stem}.png');assert min(im.info['dpi'])>599
 dimensions=im.size;im.thumbnail((1700,1200));im.save(OUT/'figures'/f'{stem}_preview.png');plt.close(fig)
 return dict(stem=stem,dimensions_px=dimensions,png_dpi=600,figure_text=texts)

def main_figure(d):
 m=d['metrics'];q=d['q'];f=plt.figure(figsize=(12.8,7.2))
 f.suptitle("3. Compare with the Original Authors' Code",x=.045,ha='left',y=.966,fontsize=20,weight='bold')
 f.text(.045,.897,f"Existing L075 comparison  |  Peridynamics.jl / RKCMaterial  |  Commit {m['author_commit'][:12]}",fontsize=11,color='#586b78')
 f.text(.26,.826,'(a) Same specimen, same inputs',ha='center',fontsize=12)
 f.text(.755,.826,"(b) Compare every particle's force",ha='center',fontsize=12)
 f.text(.26,.775,f"{m['particles']:,} particles  |  Strain {q['epsilon']*100:g}%  |  Horizon {q['delta_m']*1e6:g} µm",ha='center',fontsize=9)
 ax=f.add_axes([.105,.365,.30,.37],projection='3d');cloud(ax,d)
 arrow(f,(.26,.33),(.14,.294));arrow(f,(.26,.33),(.375,.294))
 box(f,.045,.205,.19,.083,'Local code','Python / BAQP',BLUE)
 box(f,.275,.205,.21,.083,'Original author code','Julia / RKCMaterial',ORANGE)
 f.text(.26,.156,'Identical coordinates, volumes, horizon and material.',ha='center',fontsize=9,color='#586b78')
 a=f.add_axes([.598,.575,.355,.205]);parity(a,d,compact=True)
 a=f.add_axes([.598,.282,.355,.17]);differences(a,d,compact=True)
 f.text(.765,.187,f"Vector force relative difference: {m['relative_vector_force_error']:.2e}",ha='center',fontsize=10.5,weight='bold')
 f.text(.765,.145,f"All {m['particles']:,} particles; maximum absolute difference {m['maximum_absolute_particle_difference_N']:.2e} N",ha='center',fontsize=8.5,color='#586b78')
 f.text(.045,.076,'Particle IDs are aligned before comparing all three force components. Coincident points overlap; none are subsampled.',fontsize=9)
 f.text(.045,.034,'Agreement checks the implementation; it does not resolve the physical boundary-loading question. Phase and time integration were disabled.',fontsize=9,color='#586b78')
 return f

def input_figure(d):
 m=d['metrics'];q=d['q'];f=plt.figure(figsize=(10.4,6.5))
 f.suptitle('Same 3D specimen in both codes',fontsize=18,weight='bold',y=.97)
 f.text(.5,.909,f"{m['particles']:,} particles  |  Same affine state, strain {q['epsilon']*100:g}%  |  $u$ × {MAG} for display",ha='center',fontsize=11,color='#586b78')
 for left,color,heading,sub in [(.055,BLUE,'Local implementation','Python / BAQP'),(.535,ORANGE,'Author implementation','Peridynamics.jl / RKCMaterial')]:
  ax=f.add_axes([left,.292,.39,.535],projection='3d');cloud(ax,d,color)
  f.text(left+.195,.835,heading,ha='center',fontsize=12,weight='bold',color=color)
  f.text(left+.195,.249,sub,ha='center',fontsize=11)
 f.text(.5,.172,'Input arrays match exactly: reference X, deformed x and nodal volumes.',ha='center',fontsize=11,weight='bold')
 f.text(.5,.121,f"Same horizon {q['delta_m']*1e6:g} µm; quadratic kernel; SVK μ = {q['material']['mu_Pa']:g} Pa, λ = {q['material']['lambda_Pa']:g} Pa.",ha='center',fontsize=10)
 f.text(.5,.070,f"Author package pinned at {m['author_commit']}",ha='center',fontsize=8,color='#586b78')
 f.text(.5,.029,'Existing intact affine evaluation; no fracture evolution or time integration. Numerical coupon, not anatomical clot geometry.',ha='center',fontsize=8,color='#586b78')
 return f

def force_figure(d):
 m=d['metrics'];f=plt.figure(figsize=(10.4,6.5))
 f.suptitle('Force on every particle compared',fontsize=18,weight='bold',y=.97)
 f.text(.5,.908,f"{m['particles']:,} aligned particle IDs  |  {m['components_compared']:,} force components  |  Existing local and author outputs",ha='center',fontsize=11,color='#586b78')
 ax=f.add_axes([.085,.34,.37,.47]);parity(ax,d);ax.set_title('(a) Force-magnitude agreement',pad=14)
 ax=f.add_axes([.60,.34,.36,.47]);differences(ax,d);ax.set_title('(b) Full-vector difference by particle',pad=14)
 f.text(.5,.210,r'$\|\mathbf{f}_{local}-\mathbf{f}_{author}\|_2/\|\mathbf{f}_{local}\|_2$'+f" = {m['relative_vector_force_error']:.3e}",ha='center',fontsize=14,weight='bold')
 f.text(.5,.149,f"Maximum per-particle absolute difference: {m['maximum_absolute_particle_difference_N']:.3e} N",ha='center',fontsize=11)
 f.text(.5,.102,f"Full-family core ({m['full_family_particles']} particles): vector relative difference {m['core_relative_vector_force_error']:.3e}",ha='center',fontsize=10,color='#586b78')
 f.text(.5,.044,'All particles retained. Differences use complete 3D vectors; no per-particle division by a near-zero force.',ha='center',fontsize=9,color='#586b78')
 return f

def interactive(d):
 import plotly.graph_objects as go
 x=d['display']*1e6;m=d['metrics']
 custom=np.column_stack([np.arange(len(x)),d['local'],d['author'],d['difference_norm'],d['x'][:,2]*1e6])
 mag=max(abs(d['local'][:,2]).max(),abs(d['author'][:,2]).max())*1e9
 colors=d['local'][:,2]*1e9
 fig=go.Figure(go.Scatter3d(x=x[:,0],y=x[:,1],z=x[:,2],mode='markers',customdata=custom,
  marker=dict(size=3,color=colors,colorscale='RdBu',cmin=-mag,cmax=mag,colorbar=dict(title='Local f_z (nN)')),
  hovertemplate='Particle %{customdata[0]:.0f}<br>Local f: (%{customdata[1]:.4e}, %{customdata[2]:.4e}, %{customdata[3]:.4e}) N<br>Author f: (%{customdata[4]:.4e}, %{customdata[5]:.4e}, %{customdata[6]:.4e}) N<br>Vector difference: %{customdata[7]:.4e} N<br>Actual z: %{customdata[8]:.5f} µm<extra></extra>'))
 buttons=[]
 for label,val,scale,lo,hi,bar in [('Local axial force',d['local'][:,2]*1e9,'RdBu',-mag,mag,'Local f_z (nN)'),
                                ('Author axial force',d['author'][:,2]*1e9,'RdBu',-mag,mag,'Author f_z (nN)'),
                                ('Vector difference',d['difference_norm'],'Viridis',0,float(d['difference_norm'].max()),'Difference (N)')]:
  buttons.append(dict(label=label,method='restyle',args=[{'marker.color':[val.tolist()],'marker.colorscale':scale,'marker.cmin':lo,'marker.cmax':hi,'marker.colorbar.title.text':bar},[0]]))
 fig.update_layout(title=dict(text='Compare the force on every particle | Existing L075 data',x=.035),
  font=dict(family='DejaVu Sans, Arial',color=INK,size=13),margin=dict(l=10,r=70,t=110,b=70),paper_bgcolor='white',
  scene=dict(xaxis_title='x (µm)',yaxis_title='y (µm)',zaxis_title='Displayed z (µm); u ×50',aspectmode='data',camera=dict(eye=dict(x=1.5,y=-1.8,z=1.2))),
  updatemenus=[dict(type='buttons',direction='right',buttons=buttons,x=.03,y=1.11)],
  annotations=[dict(text=f"All {m['particles']:,} particle IDs; hover for both force vectors. Display displacement ×{MAG}; actual strain 0.3%. No new simulation.",x=.5,y=-.065,xref='paper',yref='paper',showarrow=False)],uirevision='preserve-camera')
 p=OUT/'figures/particle_force_comparison_interactive.html'
 fig.write_html(p,include_plotlyjs=True,config=dict(responsive=True,displaylogo=False,scrollZoom=True))
 assert len(fig.data[0].x)==m['particles'] and len(buttons)==3
 return str(p.relative_to(OUT))

def main():
 for name in ('figures','tables','manifests'):(OUT/name).mkdir(parents=True,exist_ok=True)
 before={k:sha(ROOT/v) for k,v in SOURCES.items()};d=load()
 figures=[save(main_figure(d),'author_code_comparison_two_points'),save(input_figure(d),'01_same_specimen_two_codes'),save(force_figure(d),'02_every_particle_force_comparison')]
 html=interactive(d)
 assert before=={k:sha(ROOT/v) for k,v in SOURCES.items()}
 ledger=dict(source_files=[dict(role=k,path=v,sha256=before[k],bytes=(ROOT/v).stat().st_size) for k,v in SOURCES.items()],
  metrics=d['metrics'],exact_input_checks=d['input_checks'],source_hashes_unchanged=True,scientific_simulations_run=0,
  author_package_rerun=False,phase_disabled_in_source=True,particle_ID_mapping='julia_ids.csv is one-based; subtract one and scatter Julia rows to the original particle IDs.',
  force_quantity='Nodal force in N, not force density; archived Julia wrapper already multiplied b_int by system.volume.',
  displacement_magnification=MAG,plotting_scope='All particle IDs, no subsampling, smoothing, interpolation or artificial log floors.',
  vector_relative_norm='Euclidean norm over the flattened full three-component force field, denominator local-force norm.',
  interpretation='Agreement with the pinned author implementation; physical boundary-loading interpretation remains INCONCLUSIVE.',
  interactive_html=html,figures=figures)
 (OUT/'manifests/SOURCE_AND_VALIDATION.json').write_text(json.dumps(ledger,indent=2)+'\n')
 print(json.dumps({'metrics':d['metrics'],'input_checks':d['input_checks'],'sources_unchanged':True,'new_simulations':0},indent=2))

if __name__=='__main__':main()
