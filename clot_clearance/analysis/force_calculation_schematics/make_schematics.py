"""Archived-data figures for three force-calculation checks; no solver calls."""
import ast
import csv
import hashlib
import json
from pathlib import Path
from itertools import product
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.text import Text
from matplotlib.lines import Line2D
from PIL import Image

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=HERE/'output'
B='verification/pfpd_baqp_finite_domain_boundary_consistency_v1/'
S='verification/pfpd_at1_spectral_split_v1/'
PATHS={
 'state':B+'L075/AFFINE_STATE.npz',
 'summary':B+'L075/SUMMARY.json',
 'variation':B+'FIRST_VARIATION.json',
 'audit_source':'scripts/finite_domain_boundary_audit.py',
 'split':S+'constitutive/VERIFICATION.json',
 'fixtures':S+'constitutive/BAQP_SMALL_FIXTURES.json',
 'tensor_source':'scripts/spectral_split_audit.py',
 'split_source':'pd_clot/pfpd_candidate/svk_spectral_split.py',
 'physics':'verification/pfpd_at1_evolution_v1/provenance/PREDECLARED_CONFIG.json',
 'boundary_status':B+'SUMMARY.json'}
BLUE='#0072B2';ORANGE='#D55E00';GREEN='#009E73';GRAY='#a9b6c0';INK='#233746'
DISPLAY_DISPLACEMENT_SCALE=50
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,
 'axes.titlesize':11,'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
 'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':600,'text.color':INK,
 'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(key):return json.loads((ROOT/PATHS[key]).read_text())
def table(name,rows):
 with (OUT/'tables'/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def existing_tensor_inputs():
 # Read only the small existing input-definition function. Do not import the audit.
 tree=ast.parse((ROOT/PATHS['tensor_source']).read_text())
 fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='tensor_cases')
 assert {n.id for n in ast.walk(fn) if isinstance(n,ast.Name) and isinstance(n.ctx,ast.Load)} <= {'dict','np'}
 ns={'np':np};exec(compile(ast.Module(body=[fn],type_ignores=[]),PATHS['tensor_source'],'exec'),ns)
 return ns['tensor_cases']()

def existing_fd_tolerance():
 tree=ast.parse((ROOT/PATHS['audit_source']).read_text())
 fn=next(n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name=='fd_probe')
 checks=[n for n in ast.walk(fn) if isinstance(n,ast.Compare) and isinstance(n.left,ast.Call)
         and isinstance(n.left.func,ast.Name) and n.left.func.id=='min']
 assert len(checks)==1 and isinstance(checks[0].comparators[0],ast.Constant)
 return float(checks[0].comparators[0].value)

def load():
 q=read('summary');F=np.asarray(q['continuum']['F']);physics=read('physics')
 assert q['material']['mu_Pa']==physics['material']['mu'] and q['material']['lambda_Pa']==0
 with np.load(ROOT/PATHS['state']) as z:
  X=z['X'];x=z['x'];Fbar=z['Fbar']
 assert X.shape==x.shape==(q['particles'],3)
 expected=X@F.T;expected[:,2]+=q['epsilon']*q['L_m']/2
 assert np.allclose(x,expected,rtol=0,atol=1e-18)
 ferr=float(abs(Fbar-F).max());assert ferr==q['affine_Fbar_max_error']
 ids=np.arange(len(X))
 display=X+DISPLAY_DISPLACEMENT_SCALE*(x-X)
 extension=float(np.ptp(x[:,2])-np.ptp(X[:,2]));assert np.isclose(extension,q['epsilon']*q['L_m'],rtol=1e-12)
 table('affine_particles_3d.csv',[dict(ID=int(i),X_m=X[i,0],Y_m=X[i,1],Z_m=X[i,2],x_m=x[i,0],y_m=x[i,1],z_m=x[i,2],
       actual_u_z_m=x[i,2]-X[i,2],display_x_um=display[i,0]*1e6,display_y_um=display[i,1]*1e6,
       display_z_um=display[i,2]*1e6,display_displacement_multiplier=DISPLAY_DISPLACEMENT_SCALE) for i in ids])
 variation=read('variation');assert variation['passed'];vr=[];tol=existing_fd_tolerance()
 for probe in variation['probes']:
  assert probe['passed']
  for r in probe['probes']:
   # Stored normalization has a small energy floor; do not replace it with a new convention.
   assert np.isfinite([r['finite_difference_J'],r['force_work_J'],r['relative_error']]).all()
   vr.append(dict(seed=probe['seed'],kind=probe['kind'],**r,acceptance_tolerance=tol))
 assert len(vr)==24 and len(variation['probes'])==8
 maxerr=max(r['relative_error'] for r in vr);assert maxerr<tol
 table('energy_force_all_probes.csv',vr)
 tensors=existing_tensor_inputs();split=read('split');assert split['gate1_passed'] and split['gate2_passed']
 tr=[]
 for case,label in [('hydrostatic_tension','Pure tension'),('hydrostatic_compression','Pure compression')]:
  r=next(r for r in split['tensors'] if r['case']==case and r['lambda_Pa']==q['material']['lambda_Pa'])
  E=tensors[case];assert np.allclose(E,np.eye(3)*E[0,0])
  stretch=float(np.sqrt(1+2*E[0,0]))
  # Algebraic cross-check of existing scalar records, not a fresh constitutive solve.
  ep=np.maximum(np.diag(E),0);em=np.minimum(np.diag(E),0);mu=q['material']['mu_Pa']
  assert np.isclose(r['Psi_plus'],mu*np.dot(ep,ep),atol=1e-13)
  assert np.isclose(r['Psi_minus'],mu*np.dot(em,em),atol=1e-13)
  tr.append(dict(label=label,case=case,lambda_Pa=r['lambda_Pa'],mu_Pa=mu,Green_E11=float(E[0,0]),Green_E22=float(E[1,1]),
                 Green_E33=float(E[2,2]),principal_stretch=stretch,Psi_plus_Pa=r['Psi_plus'],Psi_minus_Pa=r['Psi_minus'],phase_driver_Pa=r['Psi_plus']))
 assert tr[0]['Psi_minus_Pa']==0 and tr[1]['Psi_plus_Pa']==0
 fixture=next(r for r in read('fixtures')['fixtures'] if r['case']=='compression')
 assert fixture['driver_max']==0 and fixture['maximum_phase_increase']==0
 table('spectral_material_point_cases.csv',tr)
 status=read('boundary_status');assert status['status']['ROOT_CAUSE_PRIMARY']=='INCONCLUSIVE'
 return dict(q=q,F=F,X=X,x=x,ids=ids,display=display,Ferror=ferr,extension_m=extension,
             variation=vr,variation_max_error=maxerr,tolerance=tol,tensors=tr,compression_fixture=fixture)

def box_edges(points):
 lo=points.min(axis=0);hi=points.max(axis=0)
 for axis in range(3):
  other=[j for j in range(3) if j!=axis]
  for choices in product((0,1),repeat=2):
   edge=np.zeros((2,3));edge[:,axis]=[lo[axis],hi[axis]]
   for j,choice in zip(other,choices):edge[:,j]=hi[j] if choice else lo[j]
   yield edge

def affine(ax,d):
 X=d['X']*1e6;x=d['display']*1e6
 # Every archived particle is included; the reference is an outline, not a slice.
 ax.scatter(x[:,0],x[:,1],x[:,2],s=5,c=BLUE,alpha=.58,depthshade=False,edgecolors='none')
 for edge in box_edges(X):ax.plot(*edge.T,color=GRAY,lw=1.0,ls='--',alpha=.95)
 for edge in box_edges(x):ax.plot(*edge.T,color=BLUE,lw=.75,alpha=.65)
 top=x[:,2].max();half=d['q']['L_m']*1e6/2
 for xx,yy in [(-half,-half),(half,half),(0,0)]:
  ax.quiver(xx,yy,top+12,0,0,70,color=ORANGE,linewidth=1.2,arrow_length_ratio=.25)
 ax.set(xlim=(-half,half),ylim=(-half,half),zlim=(-half,top+95),
        xlabel=r'$x$ ($\mu$m)',ylabel=r'$y$ ($\mu$m)',zlabel=r'$z$ ($\mu$m)')
 ax.set_xticks([-half,0]);ax.set_yticks([0,half]);ax.set_zticks([-half,0,half])
 ax.set_box_aspect((2*half,2*half,top+95+half));ax.view_init(elev=21,azim=-52);ax.set_proj_type('ortho')
 for axis in (ax.xaxis,ax.yaxis,ax.zaxis):axis.pane.fill=False;axis.pane.set_edgecolor('#eeeeee');axis.labelpad=1
 ax.grid(False);ax.tick_params(pad=0,labelsize=7)
 ax.legend(handles=[Line2D([],[],color=GRAY,ls='--',lw=1,label='Reference outline'),
                    Line2D([],[],color=BLUE,marker='o',ls='none',ms=3,label='Stretched particles')],
           loc='upper left',bbox_to_anchor=(-.01,1.07),frameon=False,fontsize=7)
 ax.text2D(1,.99,f'$u$ × {DISPLAY_DISPLACEMENT_SCALE}\n(display only)',transform=ax.transAxes,ha='right',va='top',fontsize=8,color=ORANGE)
 ax.text2D(.5,-.065,f"All {len(X):,} particles in 3D",transform=ax.transAxes,ha='center',fontsize=8,color='#586b78')

def interactive_affine(d):
 import plotly.graph_objects as go
 X=d['X']*1e6;x=d['x']*1e6;display=d['display']*1e6
 def outline(points):
  vertices=[]
  for edge in box_edges(points):vertices.extend([edge[0].tolist(),edge[1].tolist(),[None]*3])
  return list(zip(*vertices))
 ref=outline(X);actual=outline(x);magnified=outline(display)
 custom=np.column_stack([np.arange(len(X)),x[:,2],(x[:,2]-X[:,2])])
 fig=go.Figure([
  go.Scatter3d(x=display[:,0],y=display[:,1],z=display[:,2],mode='markers',name='All stretched particles',
   marker=dict(size=2.3,color=BLUE,opacity=.65),customdata=custom,
   hovertemplate='Particle %{customdata[0]:.0f}<br>x: %{x:.3f} µm<br>y: %{y:.3f} µm<br>Displayed z: %{z:.3f} µm<br>Actual z: %{customdata[1]:.5f} µm<br>Actual u_z: %{customdata[2]:.5f} µm<extra></extra>'),
  go.Scatter3d(x=ref[0],y=ref[1],z=ref[2],mode='lines',name='Reference outline',line=dict(color=GRAY,width=4,dash='dash'),hoverinfo='skip'),
  go.Scatter3d(x=magnified[0],y=magnified[1],z=magnified[2],mode='lines',name='Stretched outline',line=dict(color=BLUE,width=3),hoverinfo='skip')])
 buttons=[]
 for label,coords,edges,mult in [('Actual displacement ×1',x,actual,1),('Display displacement ×50',display,magnified,DISPLAY_DISPLACEMENT_SCALE)]:
  buttons.append(dict(label=label,method='update',args=[
   {'x':[coords[:,0].tolist(),list(edges[0])],'y':[coords[:,1].tolist(),list(edges[1])],'z':[coords[:,2].tolist(),list(edges[2])]},
   {'title.text':f'Uniform stretch in 3D | Actual strain 0.3% | Display displacement ×{mult}'},[0,2]]))
 fig.update_layout(title=dict(text='Uniform stretch in 3D | Actual strain 0.3% | Display displacement ×50',x=.03),
  font=dict(family='DejaVu Sans, Arial',size=13,color=INK),paper_bgcolor='white',margin=dict(l=10,r=10,t=120,b=65),
  scene=dict(xaxis_title='x (µm)',yaxis_title='y (µm)',zaxis_title='Displayed z (µm)',
   xaxis=dict(range=[-400,400]),yaxis=dict(range=[-400,400]),zaxis=dict(range=[-400,600]),
   aspectmode='manual',aspectratio=dict(x=1,y=1,z=1.25),camera=dict(eye=dict(x=1.5,y=-1.8,z=1.2))),
  legend=dict(x=.01,y=.92,bgcolor='rgba(255,255,255,0.85)'),
  updatemenus=[dict(type='buttons',direction='right',active=1,buttons=buttons,x=.02,y=1.10)],
  annotations=[dict(text='Drag to rotate; scroll to zoom. All 2,197 archived particles. Actual end extension: 2.25 µm. No new simulation.',
                    x=.5,y=-.065,xref='paper',yref='paper',showarrow=False)],uirevision='preserve-camera')
 destination=OUT/'figures/01_uniform_stretch_interactive.html'
 fig.write_html(destination,include_plotlyjs=True,full_html=True,config={'responsive':True,'displaylogo':False,'scrollZoom':True})
 payload=fig.to_plotly_json();assert len(payload['data'])==3 and len(fig.data[0].x)==len(X)
 return dict(path=str(destination.relative_to(OUT)),self_contained=True,particle_count=len(X),display_options=[1,DISPLAY_DISPLACEMENT_SCALE])

def energy_force(ax,d):
 ax.axis('off')
 # Symbolic perturbation positions are explanatory glyphs, not new trajectory data.
 xx=[.19,.50,.81];yy=.90
 ax.scatter(xx,[yy]*3,s=[40,65,40],facecolors=['white',ORANGE,'white'],edgecolors=[GRAY,ORANGE,GRAY],transform=ax.transAxes)
 ax.annotate('',xy=(.26,yy),xytext=(.43,yy),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color=BLUE,lw=1.2))
 ax.annotate('',xy=(.74,yy),xytext=(.57,yy),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color=BLUE,lw=1.2))
 for xpos,label in zip(xx,[r'$\mathbf{x}-\eta\mathbf{v}$',r'$\mathbf{x}$',r'$\mathbf{x}+\eta\mathbf{v}$']):
  ax.text(xpos,.98,label,ha='center',va='top',fontsize=11,transform=ax.transAxes)
 ax.text(.5,.795,r'$\mathrm{d}U/\mathrm{d}\eta=-\sum_i\mathbf{f}_i\cdot\mathbf{v}_i$',ha='center',va='center',transform=ax.transAxes,fontsize=11)
 p=ax.inset_axes([.14,.08,.80,.56])
 rows=d['variation'];allx=np.array([r['force_work_J'] for r in rows])*1e12
 lo=float(allx.min());hi=float(allx.max());pad=(hi-lo)*.10
 p.plot([lo-pad,hi+pad],[lo-pad,hi+pad],ls='--',color=GRAY,lw=.8,label='Expected identity')
 for kind,label,color,marker in [('core_node','Core node',ORANGE,'s'),('core_random','Core directions',BLUE,'o'),('collar_random','Collar directions',GREEN,'^')]:
  rr=[r for r in rows if r['kind']==kind]
  p.scatter([r['force_work_J']*1e12 for r in rr],[r['finite_difference_J']*1e12 for r in rr],s=27,
            facecolors='none',edgecolors=color,marker=marker,lw=1,label=label,zorder=3)
 p.set(xlim=(lo-pad,hi+pad),ylim=(lo-pad,hi+pad),xlabel='From nodal forces (pJ)',ylabel='From energy difference (pJ)')
 p.tick_params(labelsize=7);p.set_xticks([-.4,-.2,0]);p.set_yticks([-.4,-.2,0]);p.grid(alpha=.15)
 p.legend(frameon=False,fontsize=6.5,loc='upper left',labelspacing=.3,handletextpad=.4)
 ax.text(.54,-.115,'All 24 saved comparisons; coincident points overlap.',transform=ax.transAxes,ha='center',fontsize=7,color='#586b78')

def material_split(ax,d):
 ax.axis('off')
 for j,(r,col) in enumerate(zip(d['tensors'],(ORANGE,BLUE))):
  gx=.04+j*.50;g=ax.inset_axes([gx,.60,.42,.36]);g.set(xlim=(-1,1),ylim=(-.86,.95));g.set_aspect('equal');g.axis('off')
  stretch=r['principal_stretch'];a=.47;f=a*stretch
  g.add_patch(Rectangle((-a,-a),2*a,2*a,fill=False,edgecolor=GRAY,ls=':',lw=.8))
  g.add_patch(Rectangle((-f,-f),2*f,2*f,fill=False,edgecolor=col,lw=1.6))
  # Arrow direction communicates the prescribed tension/compression tensor.
  for axis in (0,1):
   for sign in (-1,1):
    near=np.zeros(2);far=np.zeros(2);near[axis]=sign*.59;far[axis]=sign*.87
    tail,head=(near,far) if j==0 else (far,near)
    g.annotate('',xy=head,xytext=tail,arrowprops=dict(arrowstyle='->',lw=1.3,color=col))
  g.text(0,.94,r['label'],ha='center',va='bottom',fontsize=9,color=col,weight='bold')
  ax.text(gx+.21,.55,f"$E={r['Green_E11']:+.2f}I$",transform=ax.transAxes,ha='center',fontsize=9)
 p=ax.inset_axes([.18,.10,.76,.34]);v=d['tensors'];x=np.arange(2)
 plus=[r['Psi_plus_Pa'] for r in v];minus=[r['Psi_minus_Pa'] for r in v]
 p.bar(x,plus,.50,color=ORANGE,label=r'Tensile $\Psi^+$')
 p.bar(x,minus,.50,bottom=plus,color=BLUE,label=r'Compressive $\Psi^-$')
 for i,r in enumerate(v):p.text(i,plus[i]+minus[i]+.065,f"{plus[i]+minus[i]:.2f}",ha='center',fontsize=8)
 p.set(xticks=x,xticklabels=['Tension','Compression'],ylim=(0,2.05),ylabel='Energy density (Pa)')
 p.set_yticks([0,.6,1.2,1.8]);p.tick_params(labelsize=8);p.legend(frameon=False,fontsize=7,loc='upper left',ncols=1)
 ax.text(.56,-.055,r'Phase driver: $Y=\Psi^+$',transform=ax.transAxes,ha='center',fontsize=10)
 ax.text(.56,-.12,f"Tension: {v[0]['phase_driver_Pa']:.2f} Pa    Compression: {v[1]['phase_driver_Pa']:.0f} Pa",transform=ax.transAxes,ha='center',fontsize=8,color='#586b78')

DRAW=[affine,energy_force,material_split]
HEADS=['(a) Uniform stretch in 3D','(b) Energy–force check','(c) Stretching versus squeezing']
STEMS=['01_uniform_stretch','02_energy_force_consistency','03_tension_compression_split']

def save(fig,stem):
 fig.canvas.draw();texts=[t.get_text() for t in fig.findobj(Text) if t.get_text()]
 assert not any(any('\u3400'<=c<='\u9fff' for c in s) for s in texts)
 for ext in ('png','pdf','svg'):fig.savefig(OUT/'figures'/f'{stem}.{ext}',dpi=600,facecolor='white')
 im=Image.open(OUT/'figures'/f'{stem}.png');dpi=im.info['dpi'];assert min(dpi)>599
 size=im.size;im.thumbnail((1600,1200));im.save(OUT/'figures'/f'{stem}_preview.png');plt.close(fig)
 return dict(stem=stem,dimensions_px=size,png_dpi=dpi,text=texts)

def main():
 for sub in ('figures','tables','manifests'):(OUT/sub).mkdir(parents=True,exist_ok=True)
 hashes={k:sha(ROOT/p) for k,p in PATHS.items()};d=load();reports=[]
 lines=[('Affine deformation reproduced',f"Actual strain: {100*d['q']['epsilon']:g}%; extension: {d['extension_m']*1e6:.2f} µm\nMax gradient error: {d['Ferror']:.2e}"),
        ('Energy and forces agree',f"8 directions × 3 steps; max error {d['variation_max_error']:.2e}"),
        ('Only tensile energy drives phase','Compression still stores elastic energy')]
 fig=plt.figure(figsize=(12.8,7.2));fig.suptitle('2. Test the Force Calculation',x=.045,ha='left',y=.966,fontsize=20,weight='bold')
 fig.text(.045,.899,'Archived L075 affine / force checks and spectral material-point tests',fontsize=11,color='#586b78')
 boxes=[[.027,.31,.30,.49],[.376,.33,.26,.465],[.705,.33,.265,.465]]
 for i,(box,draw) in enumerate(zip(boxes,DRAW)):
  ax=fig.add_axes(box,projection='3d' if i==0 else None);draw(ax,d);center=box[0]+box[2]/2
  fig.text(center,.828,HEADS[i],ha='center',fontsize=12)
  fig.text(center,.197,lines[i][0],ha='center',va='top',fontsize=11,weight='bold')
  fig.text(center,.151,lines[i][1],ha='center',va='top',fontsize=8.5,color='#586b78')
 fig.text(.045,.078,'Affine accuracy does not establish force equilibrium. Perturbation checks cover the saved directions, not every particle independently.',fontsize=9)
 fig.text(.045,.035,'Tensile phase driving is a constitutive rule; these verification results do not establish physical clot tearing.',fontsize=9,color='#586b78')
 reports.append(save(fig,'force_calculation_three_points'))
 for i,draw in enumerate(DRAW):
  fig=plt.figure(figsize=(6.2,6.5));fig.suptitle(HEADS[i][4:],fontsize=16,weight='bold',y=.966)
  fig.text(.5,.91,'Existing L075 affine state' if i==0 else 'Archived first-variation probes' if i==1 else 'Archived spectral material-point tests',ha='center',fontsize=10,color='#586b78')
  box=[.08,.28,.82,.59] if i==0 else [.18,.30,.66,.55]
  ax=fig.add_axes(box,projection='3d' if i==0 else None);draw(ax,d)
  fig.text(.5,.172,lines[i][0],ha='center',fontsize=12,weight='bold')
  fig.text(.5,.138,lines[i][1],ha='center',va='top',fontsize=9.5,color='#586b78')
  notes=[f"Actual strain = {d['q']['epsilon']*100:g}%; only the displayed displacement is multiplied by {DISPLAY_DISPLACEMENT_SCALE}.\nAccurate deformation reconstruction does not imply affine force equilibrium.",
         f"All {len(d['variation'])} saved comparisons are plotted; markers can overlap.\nParticle-position glyphs are schematic, not new simulation states.",
         'Pure tensile and compressive Green-strain inputs; the same material is used.\nMixed/shear states can also contain tensile energy. This is not a rupture test.']
  fig.text(.5,.022,notes[i],ha='center',va='bottom',fontsize=8,linespacing=1.5,color='#586b78')
  reports.append(save(fig,STEMS[i]))
 interactive=interactive_affine(d)
 assert hashes=={k:sha(ROOT/p) for k,p in PATHS.items()}
 ledger=dict(source_files=[dict(role=k,path=p,sha256=hashes[k],bytes=(ROOT/p).stat().st_size) for k,p in PATHS.items()],
             affine_case='L075',affine_particle_count=d['q']['particles'],affine_displayed_particle_count=len(d['ids']),
             affine_view='Full three-dimensional particle cloud; reference outline; no slicing or particle subsampling',
             interactive_affine=interactive,
             actual_engineering_strain=d['q']['epsilon'],actual_end_extension_m=d['extension_m'],
             affine_Fbar_error=d['Ferror'],archived_Ftilde_error=d['q']['affine_Ftilde_max_error'],
             display_displacement_multiplier=DISPLAY_DISPLACEMENT_SCALE,
             displayed_position_rule='X + 50*(x - X), using the archived translated affine x with fixed bottom plane',
             variation_directions=8,variation_steps_per_direction=3,variation_max_relative_error=d['variation_max_error'],
             original_variation_gate='For each direction, min(error over three steps) < tolerance; all 24 saved comparisons also lie below this tolerance.',
             original_variation_tolerance=d['tolerance'],spectral_cases=d['tensors'],
             independent_compression_fixture=d['compression_fixture'],
             scientific_simulations_run=0,source_hashes_unchanged=True,all_figure_text_English=True,
             limited_claims=['Affine kinematics is reproduced, not affine force equilibrium.',
              'Saved directional derivative checks, not an exhaustive per-particle coordinate test.',
              'Tensile spectral energy drives phase; no physical tearing claim.',
              'Geometric cell glyphs encode archived strain inputs and are not specimen simulations.'],figures=reports)
 (OUT/'manifests/SOURCE_AND_VALIDATION.json').write_text(json.dumps(ledger,indent=2)+'\n')
 print(json.dumps({k:v for k,v in ledger.items() if k not in ('source_files','figures','spectral_cases','independent_compression_fixture')},indent=2))

if __name__=='__main__':main()
