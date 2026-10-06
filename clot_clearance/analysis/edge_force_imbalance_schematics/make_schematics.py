"""Archived affine-force / coefficient evidence. No mechanics or evolution solve."""
import csv
import argparse
import hashlib
import json
from itertools import product
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyBboxPatch
from matplotlib.text import Text
from matplotlib.lines import Line2D
from PIL import Image

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];OUT=HERE/'output'
B='verification/pfpd_baqp_finite_domain_boundary_consistency_v1/'
T='results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/'
SOURCES={'state':B+'L075/AFFINE_STATE.npz','case':B+'L075/SUMMARY.json','status':B+'SUMMARY.json','config':B+'CONFIG.json',
 'particles':T+'L075_particles.csv','layers':T+'L075_layers.csv','channels':T+'L075_coefficient_channels.csv',
 'normalization':B+'NORMALIZATION_CHANNELS.json','first_variation':B+'FIRST_VARIATION.json',
 'author_comparison':B+'author_L075/COMPARISON.json','force_source':'pd_clot/pfpd_candidate/boundary_force_diagnostics.py',
 'channel_source':'scripts/finite_domain_normalization.py'}
BLUE='#0072B2';ORANGE='#D55E00';GREEN='#009E73';GRAY='#a9b6c0';INK='#233746'
MAG=50;FORCE_ARROW_UM_PER_NN=8
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,'axes.titlesize':11,
 'xtick.labelsize':8,'ytick.labelsize':8,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,
 'axes.linewidth':.7,'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':600,
 'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(k):return json.loads((ROOT/SOURCES[k]).read_text())
def csvread(k):return list(csv.DictReader((ROOT/SOURCES[k]).open()))
def table(name,rows):
 with (OUT/'tables'/name).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def load():
 q=read('case');status=read('status');config=read('config');points=csvread('particles');layers=csvread('layers');channels=csvread('channels')
 with np.load(ROOT/SOURCES['state'],allow_pickle=False) as z:
  X=z['X'];x=z['x'];f=z['force'];full=z['full_support']
 N=len(X);h=q['h_m'];L=q['L_m'];delta=q['delta_m'];R=q['continuum']['reaction_N']
 assert X.shape==x.shape==f.shape==(q['particles'],3)
 assert np.array_equal(X,np.array([[float(r[k]) for k in ('X_m','Y_m','Z_m')] for r in points]))
 assert np.array_equal(f,np.array([[float(r[k]) for k in ('force_x_N','force_y_N','force_z_N')] for r in points]))
 dmin=L/2-abs(X).max(1);assert np.array_equal(full,dmin>=delta)
 fixed=np.isclose(abs(X[:,2]),L/2,rtol=0,atol=h*1e-9);free=~fixed
 assert free.sum()==q['free_force_stats']['count'] and full.sum()==q['full_support_count']
 mag=np.linalg.norm(f,axis=1);assert np.isclose(mag[free].max(),q['free_force_stats']['max_N'],rtol=1e-13,atol=0)
 assert np.isclose(mag[full].max()/R,q['bulk']['max_over_R'],rtol=1e-13,atol=0)
 assert all(r['all_neighbors_full']=='False' and r['all_Z_normalization_inputs_full']=='False' for r in points)
 assert q['two_hop_closed_count']==q['three_hop_closed_count']==0
 for r in layers:
  mask=np.isclose(X[:,2],float(r['z_m']),rtol=0,atol=h*1e-9)
  assert abs(f[mask,2].sum()/R-float(r['layer_fz_over_R']))<1e-13
 lower=[r for r in layers if float(r['z_m'])<=h*1e-9]
 core_rows=[r for r in channels if r['full_support']=='True']
 row=max(core_rows,key=lambda r:abs(float(r['force_z_N'])));i=int(row['ID'])
 assert full[i] and row['direct_self_signed_fraction']!=''
 keys=['direct_self','direct_neighbor','Z_self','Z_neighbor']
 values=np.array([float(row[k+'_z_N']) for k in keys]);fractions=np.array([float(row[k+'_signed_fraction']) for k in keys])
 assert np.isclose(values.sum(),f[i,2],rtol=1e-12,atol=0) and abs(fractions.sum()-1)<1e-12
 # An actual short coefficient-dependency path. These links are not per-link force attributions.
 target_j=X[i].copy();target_j[2]=-L/2+h
 target_k=X[i].copy();target_k[2]=-L/2
 j=int(np.linalg.norm(X-target_j,axis=1).argmin());k=int(np.linalg.norm(X-target_k,axis=1).argmin())
 assert np.linalg.norm(X[j]-target_j)<1e-14 and np.linalg.norm(X[k]-target_k)<1e-14
 path=[k,j,i];pathrows=[]
 for role,pid in zip(['surface k','neighbor j','full-family i'],path):
  neighbor_count=int(np.sum((np.linalg.norm(X-X[pid],axis=1)>0)&(np.linalg.norm(X-X[pid],axis=1)<=delta)))
  assert neighbor_count==int(points[pid]['neighbor_count'])
  pathrows.append(dict(role=role,ID=pid,X_m=X[pid,0],Y_m=X[pid,1],Z_m=X[pid,2],depth_from_bottom_over_h=(X[pid,2]+L/2)/h,
                       neighbors=neighbor_count,full_family=bool(full[pid]),force_z_N=f[pid,2]))
 lengths=np.array([np.linalg.norm(X[j]-X[k]),np.linalg.norm(X[i]-X[j])])
 assert np.all(lengths<=delta) and np.linalg.norm(X[i]-X[k])>delta
 assert q['bulk_patch_test']=='FAIL' and status['status']['ROOT_CAUSE_PRIMARY']=='INCONCLUSIVE'
 assert status['status']['BOUNDARY_CONTROL_RUN']=='NO' and not status['force_stencil_isolated_core_available']
 assert read('first_variation')['passed'] and read('author_comparison')['passed']
 # Cross-check the selected coefficient row against its JSON copy.
 saved=next(c for c in read('normalization')['cases'] if c['case']=='L075')
 saved=next(r for r in saved['rows'] if r['ID']==i)
 assert float(row['Z_neighbor_signed_fraction'])==saved['Z_neighbor_signed_fraction']
 table('particles_and_internal_force.csv',[dict(ID=n,X_m=X[n,0],Y_m=X[n,1],Z_m=X[n,2],x_m=x[n,0],y_m=x[n,1],z_m=x[n,2],
   f_x_N=f[n,0],f_y_N=f[n,1],f_z_N=f[n,2],force_norm_N=mag[n],fixed_loading_plane=bool(fixed[n]),full_family=bool(full[n]),
   nearest_surface_distance_over_delta=dmin[n]/delta) for n in range(N)])
 table('lower_half_signed_layers.csv',[dict(depth_from_bottom_over_h=(float(r['z_m'])+L/2)/h,**r) for r in lower])
 table('actual_dependency_path.csv',pathrows)
 table('selected_core_force_channels.csv',[dict(ID=i,channel=key,force_z_N=v,signed_fraction=frac) for key,v,frac in zip(keys,values,fractions)])
 metrics=dict(particles=N,free_particles=int(free.sum()),fixed_particles=int(fixed.sum()),maximum_free_force_N=float(mag[free].max()),
              full_family_particles=int(full.sum()),two_hop_closed=q['two_hop_closed_count'],three_hop_closed=q['three_hop_closed_count'],
              maximum_core_force_N=float(mag[full].max()),core_max_over_R=q['bulk']['max_over_R'],core_collar_ratio=q['bulk_collar_max_ratio'],
              core_patch_status=q['bulk_patch_test'],noise_gate=config['bulk_noise_normalized_tolerance'],ratio_gate=config['bulk_collar_ratio_max'],
              selected_core_ID=i,selected_core_signed_force_N=float(f[i,2]),neighbor_Z_signed_fraction=float(fractions[3]),
              dependency_path_IDs=path,link_lengths_m=lengths.tolist(),root_cause='INCONCLUSIVE')
 table('evidence_summary.csv',[metrics])
 return dict(q=q,X=X,x=x,display=(X+MAG*(x-X))*1e6,f=f,full=full,free=free,fixed=fixed,layers=lower,channels=values,
             fractions=fractions,path=path,pathrows=pathrows,path_lengths=lengths,metrics=metrics,config=config)

def edges(points):
 lo=points.min(0);hi=points.max(0)
 for axis in range(3):
  other=[j for j in range(3) if j!=axis]
  for cc in product((0,1),repeat=2):
   e=np.zeros((2,3));e[:,axis]=[lo[axis],hi[axis]]
   for j,c in zip(other,cc):e[:,j]=hi[j] if c else lo[j]
   yield e

def format3d(ax,d,labels=True):
 p=d['display'];half=d['q']['L_m']*1e6/2;top=p[:,2].max()
 for e in edges(p):ax.plot(*e.T,color=GRAY,lw=.6,alpha=.7)
 ax.set(xlim=(-half,half),ylim=(-half,half),zlim=(-half,top+25));ax.set_box_aspect((2*half,2*half,top+half+25))
 ax.set_proj_type('ortho');ax.view_init(elev=22,azim=-52);ax.grid(False)
 if labels:
  ax.set(xlabel=r'$x$ ($\mu$m)',ylabel=r'$y$ ($\mu$m)',zlabel=r'$z$ ($\mu$m)')
  ax.set_xticks([-half,0]);ax.set_yticks([0,half]);ax.set_zticks([-half,0,half]);ax.tick_params(pad=0,labelsize=7)
  for a in (ax.xaxis,ax.yaxis,ax.zaxis):a.pane.fill=False;a.pane.set_edgecolor('#eeeeee');a.labelpad=1
 else:ax.set_axis_off()

def force_cloud(ax,d,arrows=False):
 p=d['display'];free=d['free'];fixed=d['fixed'];v=d['f'][:,2]*1e9;lim=float(abs(v[free]).max())
 sc=ax.scatter(*p[free].T,c=v[free],s=7,cmap='RdBu_r',norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim),alpha=.8,depthshade=False,edgecolors='none')
 ax.scatter(*p[fixed].T,c=GRAY,s=7,alpha=.28,depthshade=False,marker='s',edgecolors='none')
 if arrows:
  X=d['X'];mask=np.isclose(X[:,0],0,atol=1e-14)&np.isclose(X[:,1],0,atol=1e-14)&free
  vectors=d['f'][mask]*1e9*FORCE_ARROW_UM_PER_NN
  ax.quiver(*p[mask].T,*vectors.T,color=INK,linewidth=1.1,arrow_length_ratio=.25,normalize=False)
 format3d(ax,d);ax.text2D(.98,.99,f'$u$ × {MAG}\n(display only)',transform=ax.transAxes,ha='right',va='top',color=ORANGE,fontsize=8)
 return sc

def force_colorbar(fig,sc,rect):
 c=fig.colorbar(sc,cax=fig.add_axes(rect),orientation='horizontal');c.set_label('Free-node internal $f_z$ (nN)',fontsize=8);c.ax.tick_params(labelsize=7);c.set_ticks([-15,0,15])

def dependency(ax,d):
 ax.set_axis_off();path=d['path'];positions=[.10,.49,.89];colors=[GRAY,ORANGE,GREEN]
 for pos,col,pid,role in zip(positions,colors,path,['Surface k','Neighbor j','Full-family i']):
  ax.scatter([pos],[.53],s=135,color=col,edgecolors='white',linewidths=.8,transform=ax.transAxes,zorder=3)
  ax.text(pos,.89,role,ha='center',va='top',transform=ax.transAxes,fontsize=8.5,weight='bold',color=INK)
  ax.text(pos,.19,f'ID {pid}',ha='center',va='top',transform=ax.transAxes,fontsize=8)
 for a,b,L in zip(positions[:-1],positions[1:],d['path_lengths']):
  ax.annotate('',xy=(b-.07,.53),xytext=(a+.07,.53),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color=INK,lw=1.2))
  ax.text((a+b)/2,.66,f"{L/d['q']['h_m']:.0f}h",ha='center',transform=ax.transAxes,fontsize=8)
 ax.text(.5,-.13,'One real dependency path; layout not to scale',ha='center',transform=ax.transAxes,fontsize=7,color='#586b78')

def layer_plot(ax,d):
 rows=d['layers'];q=d['q'];x=np.array([(float(r['z_m'])+q['L_m']/2)/q['h_m'] for r in rows]);y=np.array([float(r['layer_fz_over_R']) for r in rows])
 ax.plot(x,y,color=BLUE,lw=.9,zorder=2);ax.scatter(x[1:],y[1:],c=BLUE,s=23,label='Free layers',zorder=3)
 ax.scatter(x[0],y[0],c=GRAY,marker='s',s=27,label='Loaded face',zorder=3)
 ax.axhline(0,color=GRAY,lw=.75);ax.axvline(q['delta_m']/q['h_m'],color=ORANGE,ls='--',lw=.85,label=r'One horizon $\delta$')
 ax.set(xlabel='Layers inward from bottom face (d/h)',ylabel=r'Layer $\sum f_z/R_{cont}$',xlim=(-.3,6.3),ylim=(-.18,.75),xticks=range(7))
 ax.grid(axis='y',alpha=.17);ax.legend(frameon=False,fontsize=6.7,loc='upper right')

def core_cloud(ax,d,path=False,labels=False):
 p=d['display'];full=d['full'];ax.scatter(*p[~full].T,s=5,c=GRAY,alpha=.12,depthshade=False,edgecolors='none')
 if not path:ax.scatter(*p[full].T,s=26,color=GREEN,depthshade=False,edgecolors='white',linewidths=.35)
 else:
  ax.scatter(*p[full].T,s=9,color=GREEN,alpha=.25,depthshade=False,edgecolors='none')
  for pid,col in zip(d['path'],[INK,ORANGE,GREEN]):ax.scatter(*p[pid],s=70,color=col,depthshade=False,edgecolors='white',linewidths=.7)
  pathp=p[d['path']];ax.plot(*pathp.T,color=INK,lw=2)
  for pid,label in zip(d['path'],['k','j','i']):
   a=p[pid];ax.text(a[0]+30,a[1]-15,a[2]+12,f'{label}: {pid}',fontsize=8,color=INK)
 format3d(ax,d,labels=labels)
 ax.text2D(.98,.99,f'$u$ × {MAG}\n(display only)',transform=ax.transAxes,ha='right',va='top',color=ORANGE,fontsize=8 if labels else 7)

def open_question(fig,rect):
 x,y,w,h=rect;fig.add_artist(FancyBboxPatch((x,y),w,h,transform=fig.transFigure,
  boxstyle='round,pad=0.005,rounding_size=.008',facecolor='#f6f8fa',edgecolor=GRAY,linewidth=.8))
 fig.text(x+w/2,y+h*.72,'Boundary influence or bulk defect?',ha='center',va='center',fontsize=9.5,weight='bold')
 fig.text(x+w/2,y+h*.29,'Physical cause: UNRESOLVED',ha='center',va='center',fontsize=9,color=ORANGE)


def combined(d):
 m=d['metrics'];q=d['q'];f=plt.figure(figsize=(12.8,7.2))
 f.suptitle('3. Trace the Force Imbalance Near Edges',x=.045,ha='left',y=.966,fontsize=20,weight='bold')
 f.text(.045,.897,f"Existing L075 affine state  |  Strain {q['epsilon']*100:g}%  |  {m['particles']:,} particles  |  Phase disabled",fontsize=11,color='#586b78')
 heads=['(a) Uniform stretch, nonzero free-node force','(b) Neighbor dependence reaches inward','(c) Boundary effect or bulk defect?']
 for pos,head in zip([.177,.512,.835],heads):f.text(pos,.825,head,ha='center',fontsize=10.5)
 ax=f.add_axes([.025,.375,.30,.40],projection='3d');sc=force_cloud(ax,d)
 force_colorbar(f,sc,[.075,.311,.205,.015])
 f.text(.177,.253,'Gray end faces: constrained nodes',ha='center',fontsize=8,color='#586b78')
 f.text(.177,.199,f"Max free-node force: {m['maximum_free_force_N']*1e9:.2f} nN",ha='center',fontsize=10.5,weight='bold')
 f.text(.177,.152,'Accurate deformation is not force equilibrium.',ha='center',fontsize=8,color='#586b78')
 ax=f.add_axes([.367,.638,.28,.13]);dependency(ax,d)
 ax=f.add_axes([.417,.329,.215,.22]);layer_plot(ax,d);ax.tick_params(labelsize=7);ax.xaxis.label.set_size(8);ax.yaxis.label.set_size(8)
 f.text(.516,.208,f"Selected core particle {m['selected_core_ID']}:",ha='center',fontsize=9)
 f.text(.516,.160,f"Neighbor Z term: {100*m['neighbor_Z_signed_fraction']:.1f}% of signed $f_z$",ha='center',fontsize=9.5,weight='bold')
 ax=f.add_axes([.715,.535,.24,.26],projection='3d');core_cloud(ax,d)
 f.text(.835,.488,f"{m['full_family_particles']} geometrically full families",ha='center',fontsize=10,color=GREEN,weight='bold')
 f.text(.835,.440,'0 closed two-hop / three-hop cores',ha='center',fontsize=9)
 f.text(.835,.393,f"Core max |f| / R = {m['core_max_over_R']:.2e}",ha='center',fontsize=9)
 f.text(.835,.348,f"Noise criterion: {m['noise_gate']:.0e}  →  FAIL",ha='center',fontsize=9,color=ORANGE)
 open_question(f,[.712,.166,.25,.118])
 f.text(.045,.076,'Static internal-force evidence. Dependency arrows do not show a time-evolving propagation front or per-link force attribution.',fontsize=9)
 f.text(.045,.034,'Constraint reactions must be distinguished from unbalanced free-node forces. Geometric full families do not isolate the complete force stencil.',fontsize=9,color='#586b78')
 return f

def first_figure(d):
 m=d['metrics'];q=d['q'];f=plt.figure(figsize=(9.6,6.5))
 f.suptitle('Uniform stretch still leaves free-node forces',fontsize=18,weight='bold',y=.967)
 f.text(.5,.909,'Archived L075 intact affine state; no additional external force on free nodes',ha='center',fontsize=10,color='#586b78')
 ax=f.add_axes([.025,.237,.52,.61],projection='3d');sc=force_cloud(ax,d,arrows=True)
 force_colorbar(f,sc,[.10,.190,.34,.020])
 ax=f.add_axes([.653,.485,.28,.293]);layer_plot(ax,d)
 f.text(.797,.336,f"{m['free_particles']:,} free particles\nMax |f| = {m['maximum_free_force_N']*1e9:.3f} nN",ha='center',fontsize=11,weight='bold',linespacing=1.6)
 f.text(.797,.240,'Surface internal force can be balanced\nby a constraint reaction.\nFree-node residuals are the imbalance.',ha='center',fontsize=9,linespacing=1.55,color='#586b78')
 f.text(.5,.104,f"Actual strain {q['epsilon']*100:g}%; geometry displacement ×{MAG}. Gray squares: fixed top/bottom particles.",ha='center',fontsize=9)
 f.text(.5,.053,f'Black arrows: saved centerline internal forces, drawn at {FORCE_ARROW_UM_PER_NN} µm per nN. Colors use actual forces.',ha='center',fontsize=8.5,color='#586b78')
 return f

def second_figure(d):
 from mpl_toolkits.mplot3d import proj3d
 m=d['metrics'];q=d['q'];f=plt.figure(figsize=(14,9.5))
 f.suptitle('Neighbor dependence and inward force distribution',fontsize=22,weight='bold',y=.967)
 f.text(.5,.916,'Actual particle links, saved layer forces, and exact force decomposition',ha='center',fontsize=12,color='#586b78')
 for x,y,title in [(.045,.856,'(a) Actual particles in three dimensions'),(.565,.856,'(b) Signed lower-half layer forces'),
                   (.045,.414,'(c) One actual dependency path'),(.565,.414,f"(d) Force decomposition at particle {m['selected_core_ID']}")]:
  f.text(x,y,title,fontsize=13,weight='bold',va='center')
 for xx,yy in [([.045,.955],[.452,.452]),([.53,.53],[.137,.875]),([.045,.955],[.106,.106])]:
  f.add_artist(Line2D(xx,yy,transform=f.transFigure,lw=.7,color='#e1e7eb'))

 # Reserve a separate label column: no particle names lie on the 3-D cloud.
 ax=f.add_axes([.025,.500,.315,.308],projection='3d')
 core_cloud(ax,d,path=True,labels=True)
 for label in list(ax.texts):label.remove()
 ax.set_zlabel('')
 f.text(.045,.817,f'Geometry displacement ×{MAG}; coordinates in µm',fontsize=10,color='#586b78')
 f.canvas.draw()
 zx,zy,_=proj3d.proj_transform(q['L_m']*5e5,q['L_m']*5e5,ax.get_zlim()[1],ax.get_proj())
 ax.annotate(r'$z$ ($\mu$m)',xy=(zx,zy),xycoords='data',xytext=(4,8),textcoords='offset points',ha='left',va='bottom',fontsize=9)
 for pid,role,color,y in zip(d['path'],['Surface k','Neighbor j','Full family i'],[INK,ORANGE,GREEN],[.516,.600,.692]):
  point=d['display'][pid];px,py,_=proj3d.proj_transform(*point,ax.get_proj())
  ax.annotate(f'{role}\nID {pid}',xy=(px,py),xycoords='data',xytext=(.295,y),textcoords=f.transFigure,
              ha='left',va='center',fontsize=10.5,weight='bold',color=color,linespacing=1.6,
              bbox=dict(boxstyle='round,pad=.35',fc='white',ec='none'),
              arrowprops=dict(arrowstyle='-',color=color,lw=.9,shrinkA=6,shrinkB=5,connectionstyle='arc3,rad=0'))

 ax=f.add_axes([.602,.564,.343,.221]);layer_plot(ax,d)
 ax.legend(frameon=False,loc='lower left',bbox_to_anchor=(-.01,1.035),ncol=3,fontsize=9,
           borderaxespad=0,handlelength=1.5,columnspacing=1.2)
 ax.tick_params(labelsize=10,pad=5);ax.xaxis.label.set_size(10.5);ax.yaxis.label.set_size(10.5)
 ax.xaxis.labelpad=9;ax.yaxis.labelpad=10
 layer5=d['layers'][5]
 f.text(.565,.479,f"Layer 5: signed sum / R = {float(layer5['layer_fz_over_R']):.5f}; depth = 5h > δ",fontsize=10,color='#586b78')

 # The reference-distance diagram occupies its own panel beneath the cloud.
 ax=f.add_axes([.045,.196,.45,.181]);ax.set_axis_off()
 positions=[.09,.50,.91];colors=[INK,ORANGE,GREEN]
 for x,pid,color,role,row in zip(positions,d['path'],colors,['Surface k','Neighbor j','Full family i'],d['pathrows']):
  ax.text(x,.93,role,transform=ax.transAxes,ha='center',va='top',fontsize=11,weight='bold')
  ax.scatter([x],[.60],s=200,color=color,edgecolors='white',linewidths=.8,transform=ax.transAxes,zorder=3)
  ax.text(x,.34,f'ID {pid}',transform=ax.transAxes,ha='center',va='top',fontsize=10.5)
  ax.text(x,.12,f"{row['neighbors']} neighbors",transform=ax.transAxes,ha='center',va='top',fontsize=9.5,color='#586b78')
 for a,b,length in zip(positions[:-1],positions[1:],d['path_lengths']):
  ax.annotate('',xy=(b-.08,.60),xytext=(a+.08,.60),xycoords='axes fraction',arrowprops=dict(arrowstyle='->',color=INK,lw=1.3))
  ax.text((a+b)/2,.74,f"{length/q['h_m']:.0f}h = {length*1e6:g} µm",transform=ax.transAxes,ha='center',fontsize=10)
 f.text(.045,.145,'Reference link lengths; diagram spacing is schematic.',fontsize=10,color='#586b78')

 f.text(.565,.379,f"Saved axial force: {m['selected_core_signed_force_N']*1e9:.6f} nN",fontsize=10,color='#586b78')
 ax=f.add_axes([.686,.216,.258,.137]);fractions=d['fractions']*100
 ax.barh(np.arange(4),fractions,color=[GRAY,BLUE,GRAY,ORANGE],height=.52)
 ax.set(yticks=range(4),yticklabels=['Self direct','Neighbor direct','Self Z','Neighbor Z'],
        xlabel='Signed share of particle force (%)',xlim=(-3,114))
 ax.invert_yaxis();ax.axvline(0,color=GRAY,lw=.6);ax.tick_params(labelsize=10,pad=6);ax.xaxis.label.set_size(10)
 ax.xaxis.labelpad=8
 for index,value in enumerate(fractions):
  label=f'{value:.2f}%' if abs(value)>1e-7 else f'{value:.1e}%'
  ax.annotate(label,(max(value,0),index),xytext=(4,0),textcoords='offset points',va='center',ha='left',fontsize=9.5)
 f.text(.565,.139,'Neighbor Z sums over all neighbors, not just this path.',fontsize=10,color='#586b78')
 f.text(.045,.068,f"L075  |  {m['particles']:,} particles  |  Strain {q['epsilon']*100:g}%  |  h = {q['h_m']*1e6:g} µm  |  δ = {q['delta_m']*1e6:g} µm",fontsize=10.5)
 f.text(.045,.035,'Arrows indicate coefficient dependence; these static data do not establish a time-evolving propagation front.',fontsize=10.5,color='#586b78')
 return f

def third_figure(d):
 m=d['metrics'];f=plt.figure(figsize=(9.6,6.5))
 f.suptitle('Boundary influence or a bulk operator defect?',fontsize=18,weight='bold',y=.967)
 f.text(.5,.908,'The existing test has not isolated a boundary-independent force stencil',ha='center',fontsize=11,color='#586b78')
 ax=f.add_axes([.022,.30,.44,.54],projection='3d');core_cloud(ax,d,labels=True)
 f.text(.245,.226,f"Green: {m['full_family_particles']} full-family particles\nGray: surrounding finite-domain particles",ha='center',fontsize=10,linespacing=1.6)
 ax=f.add_axes([.555,.48,.395,.32]);ax.axis('off')
 rows=[['Geometrically full families',str(m['full_family_particles'])],['Closed two-hop cores',str(m['two_hop_closed'])],['Closed three-hop cores',str(m['three_hop_closed'])],
       ['Core max |f| / R',f"{m['core_max_over_R']:.3e}"],['Declared noise criterion',f"{m['noise_gate']:.0e}"],
       ['Core / collar maximum force',f"{m['core_collar_ratio']:.5f}"],['Declared ratio criterion',f"{m['ratio_gate']:.0e}"]]
 tab=ax.table(cellText=rows,colWidths=[.72,.28],cellLoc='left',bbox=[0,0,1,1]);tab.auto_set_font_size(False);tab.set_fontsize(9)
 for (r,c),cell in tab.get_celld().items():cell.set_edgecolor('#dde2e6');cell.set_linewidth(.6);cell.set_facecolor('#f7f9fa' if r%2==0 else 'white')
 f.text(.753,.414,'Geometric core test: FAIL',ha='center',fontsize=12,weight='bold',color=ORANGE)
 open_question(f,[.56,.247,.39,.105])
 f.text(.5,.137,'Energy derivatives and the pinned reference implementation reproduce these forces.',ha='center',fontsize=10)
 f.text(.5,.078,'Next discriminating test: close the full coefficient-dependency stencil. This test has not been run.',ha='center',fontsize=9,color='#586b78')
 f.text(.5,.033,'A possible 3δ dependence bound is not a validated boundary-layer thickness or a boundary-condition prescription.',ha='center',fontsize=8.5,color='#586b78')
 return f

def save(fig,stem):
 fig.canvas.draw();texts=[a.get_text() for a in fig.findobj(Text) if a.get_text()]
 assert not any(any('\u3400'<=c<='\u9fff' for c in a) for a in texts)
 for ext in ('png','pdf','svg'):fig.savefig(OUT/'figures'/f'{stem}.{ext}',dpi=600,facecolor='white')
 im=Image.open(OUT/'figures'/f'{stem}.png');assert min(im.info['dpi'])>599;size=im.size
 im.thumbnail((1700,1200));im.save(OUT/'figures'/f'{stem}_preview.png');plt.close(fig)
 return dict(stem=stem,dimensions_px=size,png_dpi=600,figure_text=texts)

def interactive(d):
 import plotly.graph_objects as go
 p=d['display'];N=len(p);free=d['free'];fixed=d['fixed'];full=d['full'];fz=d['f'][:,2]*1e9;lim=float(abs(fz[free]).max())
 custom=np.column_stack([np.arange(N),d['f'],full,d['fixed'],d['X'][:,2]*1e6])
 hover='ID %{customdata[0]:.0f}<br>Internal f: (%{customdata[1]:.4e}, %{customdata[2]:.4e}, %{customdata[3]:.4e}) N<br>Full family: %{customdata[4]:.0f}<br>Fixed end plane: %{customdata[5]:.0f}<br>Reference z: %{customdata[6]:.3f} µm<extra></extra>'
 fig=go.Figure()
 fig.add_trace(go.Scatter3d(x=p[free,0],y=p[free,1],z=p[free,2],customdata=custom[free],hovertemplate=hover,mode='markers',name='Free-node internal force',marker=dict(size=3,color=fz[free],colorscale='RdBu',reversescale=True,cmin=-lim,cmax=lim,colorbar=dict(title='Internal f_z (nN)'))))
 fig.add_trace(go.Scatter3d(x=p[fixed,0],y=p[fixed,1],z=p[fixed,2],customdata=custom[fixed],hovertemplate=hover,mode='markers',name='Fixed end planes',marker=dict(size=2,color=GRAY,opacity=.3)))
 fig.add_trace(go.Scatter3d(x=p[full,0],y=p[full,1],z=p[full,2],customdata=custom[full],hovertemplate=hover,mode='markers',name='Full-family core',visible=False,marker=dict(size=5,color=GREEN)))
 fig.add_trace(go.Scatter3d(x=p[:,0],y=p[:,1],z=p[:,2],mode='markers',name='Finite-domain context',visible=False,marker=dict(size=2,color=GRAY,opacity=.12),hoverinfo='skip'))
 ids=d['path'];pp=p[ids]
 fig.add_trace(go.Scatter3d(x=pp[:,0],y=pp[:,1],z=pp[:,2],mode='lines+markers+text',customdata=custom[ids],hovertemplate=hover,
   text=[f'{name}: ID {i}' for name,i in zip(['Surface k','Neighbor j','Full-family i'],ids)],textposition='middle right',name='One coefficient path',visible=False,
   marker=dict(size=7,color=[INK,ORANGE,GREEN]),line=dict(color=INK,width=6)))
 buttons=[dict(label=lab,method='restyle',args=[{'visible':v}]) for lab,v in [
  ('Free-node forces',[True,True,False,False,False]),('Full families',[False,False,True,True,False]),('Actual neighbor path',[False,False,False,True,True])]]
 fig.update_layout(title=dict(text='Trace finite-domain force imbalance | Existing L075 data',x=.03),
  font=dict(family='DejaVu Sans, Arial',color=INK,size=13),margin=dict(l=10,r=70,t=115,b=70),paper_bgcolor='white',
  scene=dict(xaxis_title='x (µm)',yaxis_title='y (µm)',zaxis_title='Displayed z (µm); u ×50',aspectmode='data',camera=dict(eye=dict(x=1.6,y=-1.8,z=1.2))),
  updatemenus=[dict(type='buttons',direction='right',buttons=buttons,x=.02,y=1.11)],legend=dict(x=.01,y=.88),
  annotations=[dict(text='All saved particles retained. Geometry displacement ×50; forces unscaled. Dependency links are not a dynamic propagation front.',
                    x=.5,y=-.065,xref='paper',yref='paper',showarrow=False)],uirevision='keep-camera')
 destination=OUT/'figures/edge_force_imbalance_interactive.html';fig.write_html(destination,include_plotlyjs=True,config=dict(responsive=True,displaylogo=False,scrollZoom=True))
 assert len(fig.data[0].x)+len(fig.data[1].x)==N and len(fig.data[2].x)==27
 return str(destination.relative_to(OUT))

def main():
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('--neighbor-only',action='store_true',help='Redraw only figure 02 and update its existing export metadata.')
 args=parser.parse_args()
 for sub in ('figures','tables','manifests'):(OUT/sub).mkdir(parents=True,exist_ok=True)
 hashes={k:sha(ROOT/p) for k,p in SOURCES.items()};d=load()
 if args.neighbor_only:
  path=OUT/'manifests/SOURCE_AND_VALIDATION.json';manifest=json.loads(path.read_text())
  assert all(hashes[item['role']]==item['sha256'] for item in manifest['source_files'])
  assert manifest['metrics']==d['metrics']
  revised=save(second_figure(d),'02_neighbor_dependence_and_layers')
  revised['layout_revision']='V3: panel (a) particle labels moved closer to the cloud with shorter leaders; panels (b-d) unchanged'
  manifest['figures']=[revised if item['stem']==revised['stem'] else item for item in manifest['figures']]
  assert hashes=={k:sha(ROOT/p) for k,p in SOURCES.items()}
  path.write_text(json.dumps(manifest,indent=2)+'\n')
  print(json.dumps(dict(figure=revised['stem'],layout_revision=revised['layout_revision'],sources_unchanged=True,metrics_unchanged=True),indent=2))
  return
 figures=[save(combined(d),'edge_force_imbalance_three_points'),save(first_figure(d),'01_affine_free_node_imbalance'),
          save(second_figure(d),'02_neighbor_dependence_and_layers'),save(third_figure(d),'03_boundary_or_bulk_unresolved')]
 html=interactive(d)
 assert hashes=={k:sha(ROOT/p) for k,p in SOURCES.items()}
 manifest=dict(source_files=[dict(role=k,path=p,sha256=hashes[k],bytes=(ROOT/p).stat().st_size) for k,p in SOURCES.items()],
  metrics=d['metrics'],sources_unchanged=True,new_scientific_simulations=0,geometry_displacement_multiplier=MAG,
  arrow_force_scale_um_per_nN=FORCE_ARROW_UM_PER_NN,force_quantity='Archived internal nodal force in N; surface constraint reaction is not an unbalanced free-node force.',
  dependency_interpretation='Actual i-j-k neighbor links illustrate coefficient support only. Neighbor Z fraction is the sum over all neighbors, not attributed to the pictured links.',
  force_profile='All seven saved bottom-half/midplane layer sums, signed and normalized by R_cont. Segments guide the eye; no interpolation data or propagation dynamics.',
  closed_stencil_test_run=False,root_cause='INCONCLUSIVE',interactive_html=html,figures=figures)
 (OUT/'manifests/SOURCE_AND_VALIDATION.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(json.dumps(d['metrics'],indent=2))

if __name__=='__main__':main()
