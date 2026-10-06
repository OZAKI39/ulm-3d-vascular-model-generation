"""Three model-construction schematics from archived L075 coordinates and counts.
Only geometric neighbor extraction is performed. No mechanics or phase solve.
"""
import csv
import hashlib
import json
from pathlib import Path
from itertools import product
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.collections import LineCollection
from matplotlib.text import Text
from PIL import Image

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
OUT=HERE/'output'
BLUE='#0072B2';ORANGE='#D55E00';GRAY='#b8c3ca';INK='#233746'
SOURCES=[
 'results/pfpd_baqp_finite_domain_boundary_consistency_v1/tables/L075_particles.csv',
 'verification/pfpd_baqp_finite_domain_boundary_consistency_v1/L075/SUMMARY.json',
 'pd_clot/pfpd_candidate/boundary_force_diagnostics.py']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':9,'axes.labelsize':9,
 'axes.titlesize':11,'xtick.labelsize':8,'ytick.labelsize':8,
 'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.7,
 'pdf.fonttype':42,'svg.fonttype':'none','savefig.dpi':600,
 'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK})

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write_csv(p,rows):
 with p.open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def load():
 rows=list(csv.DictReader((ROOT/SOURCES[0]).open()))
 q=json.loads((ROOT/SOURCES[1]).read_text())
 X=np.array([[float(r[k]) for k in ('X_m','Y_m','Z_m')] for r in rows])
 ids=np.array([int(r['ID']) for r in rows]);assert np.array_equal(ids,np.arange(len(X)))
 tree=cKDTree(X);pairs=tree.query_pairs(q['delta_m'],output_type='ndarray')
 counts=np.bincount(pairs.ravel(),minlength=len(X))
 assert np.array_equal(counts,np.array([int(r['neighbor_count']) for r in rows]))
 assert len(pairs)==q['bonds'] and len(X)==q['particles']
 selected={};edges=[]
 for key,target in [('interior',[0,0,0]),('surface',[-q['L_m']/2,0,0])]:
  i=int(np.linalg.norm(X-target,axis=1).argmin());assert np.linalg.norm(X[i]-target)<1e-14
  nn=np.array(sorted(j for j in tree.query_ball_point(X[i],q['delta_m']) if j!=i))
  plane=np.isclose(X[nn,1],X[i,1],rtol=0,atol=1e-14);shown=nn[plane]
  assert len(nn)==counts[i]
  selected[key]=dict(ID=i,coordinates_m=X[i].tolist(),neighbors_3D=len(nn),neighbors_in_slice=len(shown),full_family=rows[i]['full_support']=='True',neighbors=nn,shown=shown)
  for j in nn:
   edges.append(dict(selected_class=key,selected_ID=i,neighbor_ID=int(j),length_m=float(np.linalg.norm(X[j]-X[i])),
                     neighbor_X_m=X[j,0],neighbor_Y_m=X[j,1],neighbor_Z_m=X[j,2],
                     shown_in_midplane=bool(abs(X[j,1]-X[i,1])<1e-14)))
 assert selected['interior']['full_family'] and not selected['surface']['full_family']
 assert selected['surface']['neighbors_3D']<selected['interior']['neighbors_3D']
 write_csv(OUT/'tables/particles.csv',[dict(ID=int(ids[i]),X_m=X[i,0],Y_m=X[i,1],Z_m=X[i,2],neighbor_count=int(counts[i]),full_family=rows[i]['full_support']) for i in range(len(X))])
 write_csv(OUT/'tables/selected_neighbor_links.csv',edges)
 write_csv(OUT/'tables/selected_particles.csv',[dict(selection=k,**{n:v for n,v in r.items() if n not in ('neighbors','shown','coordinates_m')},X_m=r['coordinates_m'][0],Y_m=r['coordinates_m'][1],Z_m=r['coordinates_m'][2]) for k,r in selected.items()])
 return X,q,selected

def cloud(ax,X,q):
 x=X*1e6;half=q['L_m']*1e6/2
 # All 2,197 particles are drawn; no artificial jitter or point sampling.
 ax.scatter(x[:,0],x[:,1],x[:,2],s=6,c=BLUE,alpha=.55,depthshade=False,edgecolors='none')
 for axis in range(3):
  others=[i for i in range(3) if i!=axis]
  for signs in product((-half,half),repeat=2):
   ends=np.zeros((2,3));ends[:,axis]=[-half,half];ends[:,others]=signs
   ax.plot(*ends.T,color='#7c8c98',lw=.8,alpha=.65)
 ax.set(xlabel=r'$x$ ($\mu$m)',ylabel=r'$y$ ($\mu$m)',zlabel=r'$z$ ($\mu$m)')
 ax.set(xlim=(-half,half),ylim=(-half,half),zlim=(-half,half))
 for axis in (ax.xaxis,ax.yaxis,ax.zaxis):
  axis.set_ticks([-half,0,half]);axis.pane.fill=False;axis.pane.set_edgecolor('#eeeeee')
 # Omit the two coincident front-corner tick labels in the chosen projection.
 ax.set_xticks([-half,0]);ax.set_yticks([0,half])
 ax.tick_params(pad=0,labelsize=7);ax.xaxis.labelpad=1;ax.yaxis.labelpad=1;ax.zaxis.labelpad=1
 ax.set_box_aspect((1,1,1));ax.view_init(elev=22,azim=-48);ax.grid(False)
 ax.set_proj_type('ortho')

def neighborhood(ax,X,q,sel,kind):
 entry=sel[kind];pid=entry['ID'];origin=X[pid];xy=(X[:,[0,2]]-origin[[0,2]])*1e6
 delta=q['delta_m']*1e6;half=q['L_m']*1e6/2;lim=delta*1.23
 plane=np.isclose(X[:,1],origin[1],rtol=0,atol=1e-14)
 visible=plane&(np.max(abs(xy),axis=1)<=lim)
 shown=entry['shown'];outside=visible.copy();outside[shown]=False;outside[pid]=False
 if kind=='surface':
  ax.axvspan(-lim,0,facecolor='#f1f2f3',zorder=0)
  ax.axvline(0,color='#465563',lw=1.3,zorder=1)
  ax.text(-lim*.92,-lim*.82,'Outside\nthe specimen',ha='left',va='center',fontsize=8,color='#65747e')
  ax.text(.045,.98,rf'Surface at $x={origin[0]*1e6:g}$ $\mu$m',transform=ax.transAxes,va='top',fontsize=8)
 ax.add_patch(Circle((0,0),delta,fill=False,edgecolor=ORANGE,lw=1.2,ls=(0,(4,3)),zorder=2))
 segments=np.stack([np.zeros((len(shown),2)),xy[shown]],axis=1)
 ax.add_collection(LineCollection(segments,colors=BLUE,linewidths=.7,alpha=.53,zorder=2))
 ax.scatter(xy[outside,0],xy[outside,1],s=16,c=GRAY,alpha=.7,edgecolors='none',zorder=1)
 ax.scatter(xy[shown,0],xy[shown,1],s=23,c=BLUE,edgecolors='white',linewidths=.3,zorder=3)
 ax.scatter([0],[0],s=75,c=ORANGE,edgecolors='white',linewidths=.8,zorder=4)
 ax.annotate(r'$\delta$',xy=(delta*.707,delta*.707),xytext=(delta*.88,delta*.88),color=ORANGE,fontsize=12,
             arrowprops=dict(arrowstyle='-',color=ORANGE,lw=.8))
 ax.set(xlim=(-lim,lim),ylim=(-lim,lim),xticks=[-250,0,250],yticks=[-250,0,250],
        xlabel=r'Local $x$ ($\mu$m)',ylabel=r'Local $z$ ($\mu$m)')
 ax.set_aspect('equal');ax.spines[['top','right']].set_visible(False)
 return entry

def decorate(fig,ax,head,lines,kind):
 # Captions use figure coordinates for identical alignment below 3D and 2D axes.
 box=ax.get_position();xx=(box.x0+box.x1)/2
 fig.text(xx,.82,head,ha='center',va='bottom',fontsize=11)
 fig.text(xx,.205,lines[0],ha='center',va='top',fontsize=11,weight='bold',color=INK)
 fig.text(xx,.162,lines[1],ha='center',va='top',fontsize=9,color='#586b78')

def save(fig,stem):
 fig.canvas.draw()
 texts=[t.get_text() for t in fig.findobj(Text) if t.get_text()]
 assert not any(any('\u3400'<=c<='\u9fff' for c in s) for s in texts)
 for ext in ('png','pdf','svg'):fig.savefig(OUT/'figures'/f'{stem}.{ext}',dpi=600,facecolor='white')
 image=Image.open(OUT/'figures'/f'{stem}.png');assert min(image.info['dpi'])>=599
 # Small derivative for chat preview; publication originals stay at 600 dpi.
 image.thumbnail((1500,1100));image.save(OUT/'figures'/f'{stem}_preview.png')
 plt.close(fig)
 return dict(stem=stem,text=texts,png_dpi=600)

def main():
 before={p:sha(ROOT/p) for p in SOURCES};X,q,s=load();reports=[]
 interior=s['interior'];surface=s['surface'];lost=1-surface['neighbors_3D']/interior['neighbors_3D']
 labels=[('(a) Model split into particles',[f"{q['particles']:,} particles",f"Actual {q['n']} × {q['n']} × {q['n']} discretization"]),
         ('(b) Each particle linked to neighbors',[f"Interior: {interior['neighbors_3D']} neighbors in 3D",f"{interior['neighbors_in_slice']} in-plane connections shown"]),
         ('(c) Surface particles have fewer neighbors',[f"Surface: {surface['neighbors_3D']} neighbors in 3D",f"{surface['neighbors_in_slice']} in-plane connections shown"])]
 fig=plt.figure(figsize=(12.8,7.2))
 fig.suptitle('1. Check How the Model Is Built',x=.045,ha='left',y=.966,fontsize=20,weight='bold')
 fig.text(.045,.896,f"Existing L075 numerical specimen  |  Spacing h = {q['h_m']*1e6:g} µm  |  Interaction radius δ = {q['delta_m']*1e6:g} µm",fontsize=11,color='#586b78')
 axes=[fig.add_axes([.020,.32,.27,.46],projection='3d'),fig.add_axes([.37,.32,.265,.46]),fig.add_axes([.705,.32,.265,.46])]
 cloud(axes[0],X,q);neighborhood(axes[1],X,q,s,'interior');neighborhood(axes[2],X,q,s,'surface')
 for ax,(head,lines) in zip(axes,labels):decorate(fig,ax,head,lines,None)
 fig.text(.045,.082,'Blue: particles and connections    Orange point: selected particle    Dashed circle: interaction radius',fontsize=9)
 fig.text(.045,.038,'Connection panels show the actual y = 0 slice; counts refer to all 3D neighbors. This specimen is not an anatomical clot reconstruction.',fontsize=9,color='#586b78')
 reports.append(save(fig,'model_build_three_points'))
 for index,stem in enumerate(('01_particle_discretization','02_interior_neighborhood','03_surface_neighborhood')):
  fig=plt.figure(figsize=(5.8,6.0));head,lines=labels[index]
  ax=fig.add_axes([.16,.28,.69,.57],projection='3d' if index==0 else None)
  if index==0:cloud(ax,X,q)
  else:neighborhood(ax,X,q,s,'interior' if index==1 else 'surface')
  fig.suptitle(head[4:],y=.969,fontsize=13,weight='bold')
  fig.text(.5,.914,f"L075  |  h = {q['h_m']*1e6:g} µm  |  δ = {q['delta_m']*1e6:g} µm",ha='center',fontsize=10,color='#586b78')
  fig.text(.5,.18,lines[0],ha='center',fontsize=13,weight='bold')
  fig.text(.5,.13,lines[1],ha='center',fontsize=10)
  text='All archived particles are shown.\nNumerical specimen, not anatomical clot geometry.' if index==0 else 'Actual y = 0 slice; local coordinates relative to the selected particle.\nBlue lines: real in-plane links. Dashed circle: interaction radius.'
  fig.text(.5,.045,text,ha='center',fontsize=8,color='#586b78',linespacing=1.5)
  reports.append(save(fig,stem))
 assert before=={p:sha(ROOT/p) for p in SOURCES}
 result=dict(source_files=[dict(path=p,sha256=h,bytes=(ROOT/p).stat().st_size) for p,h in before.items()],
             case='L075',particles=len(X),total_bonds=q['bonds'],spacing_m=q['h_m'],horizon_m=q['delta_m'],
             selected={k:{n:v for n,v in r.items() if n not in ('neighbors','shown')} for k,r in s.items()},
             surface_neighbor_reduction_percent=lost*100,all_archived_neighbor_counts_match=True,
             original_sources_unchanged=True,new_scientific_simulations=0,
             extraction='cKDTree.query_pairs(delta) on archived X coordinates, same geometric rule as existing make_cube; validate every neighbor count and total bonds.',
             display='Panel a draws all particles. Panels b/c draw only actual y=0 neighbors, without subsampling within that slice; no invented outside particles.',
             figures=reports)
 (OUT/'manifests/SOURCE_AND_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ('source_files','figures')},indent=2))

if __name__=='__main__':main()
