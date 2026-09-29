"""Scientific plots of derived audit values; raw fields are never edited."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pyvista as pv
from audit_field import js
ap=argparse.ArgumentParser();ap.add_argument('--population-root',type=Path,required=True);a=ap.parse_args();root=a.population_root;out=root/'cfd_wall_audit';figdir=out/'figures';figdir.mkdir(exist_ok=True)
f=pd.read_csv(out/'all_cells.csv');z=np.load(out/'recomputed_arrays.npz');cl=np.load(out/'cluster_labels.npz');vis=pd.read_csv(out/'particle_overlap.csv');ident=json.loads((out/'AUDIT_IDENTITY.json').read_text())
m=pv.read(ident['input_paths']['mesh']);m.points=m.points*1000;wall=pv.read(out/'wall_audit.vtp');wall.points=wall.points*1000
plt.rcParams.update({'font.size':11,'axes.titlesize':13,'figure.dpi':120,'savefig.dpi':300})
figures=[]
def save(fig,name,caption):
 fig.text(.015,.015,caption,fontsize=9,va='bottom');fig.savefig(figdir/(name+'.png'),dpi=300,bbox_inches='tight');fig.savefig(figdir/(name+'.pdf'),dpi=300,bbox_inches='tight');plt.close(fig);figures.append(dict(file=name,caption=caption));print(name,flush=True)
def view(mesh,scalar=None,clim=(0,1),context=False,zoom=1.,edges=False,cmap='viridis',vector=None):
 p=pv.Plotter(off_screen=True,window_size=(1400,1050));p.set_background('white')
 if context:p.add_mesh(wall,color='#c7d0d8',opacity=.08)
 if scalar:p.add_mesh(mesh,scalars=scalar,preference='cell',clim=clim,cmap=cmap,show_edges=edges,edge_color='#333333',smooth_shading=False,show_scalar_bar=False)
 else:p.add_mesh(mesh,color='#df8843',show_edges=edges,edge_color='#111111',smooth_shading=False)
 if vector is not None:
  point,vec=vector;p.add_arrows(point.reshape(1,3),vec.reshape(1,3),mag=1,color='#c02b32')
 p.view_isometric();p.enable_parallel_projection()
 pts=wall.points if context else mesh.points
 cam=p.camera;direction=np.asarray(cam.position)-np.asarray(cam.focal_point);direction/=np.linalg.norm(direction)
 right=np.cross(np.asarray(cam.up),direction);right/=np.linalg.norm(right);up=np.cross(direction,right)
 xp=pts@right;yp=pts@up;fp=np.asarray(cam.focal_point);fp=fp+((xp.max()+xp.min())/2-fp@right)*right+((yp.max()+yp.min())/2-fp@up)*up
 cam.focal_point=fp;cam.position=fp+direction*max(np.ptp(pts,axis=0))*4
 cam.parallel_scale=max(np.ptp(yp)/2,np.ptp(xp)/(2*1400/1050))*1.10/zoom
 if context:p.add_point_labels(np.mean(mesh.points,axis=0).reshape(1,3),['Cluster'],font_size=24,point_size=12,always_visible=True)
 p.add_axes();im=p.screenshot(return_img=True);p.close();return im
def cmapbar(fig,axes,lo,hi,label,cmap='viridis'):
 fig.colorbar(plt.cm.ScalarMappable(norm=plt.Normalize(lo,hi),cmap=cmap),ax=axes,shrink=.72,pad=.015,label=label)
def clean(ax,title):ax.set_axis_off();ax.set_title(title)
fig,axs=plt.subplots(1,3,figsize=(15,6))
for k,ax in enumerate(axs):
 sel=m.extract_cells(np.flatnonzero(f.layer.values==k));sel.cell_data['eta']=f.loc[f.layer==k,'eta_div'].values;ax.imshow(view(sel,'eta',(0,1.05)));clean(ax,f'W{k}: {sel.n_cells:,} tetrahedra')
cmapbar(fig,list(axs),0,1.05,r'$\eta=|\nabla\cdot u|/||G||_F$ (dimensionless)')
save(fig,'01_all_wall_layers_eta','Disjoint face-neighbor layers. Visible faces of selected tetrahedra; original cell values, no smoothing. Geometry: mm.')
for idx,scalar,title,label,clim in [(2,'eta_div','W0 values on the official WALL',r'$\eta$ (dimensionless)',(0,1.05)),(3,'normal_fraction','W0 interior-vertex normal fraction',r'$|u_n|/|u|$ (dimensionless)',(0,1))]:
 fig,ax=plt.subplots(figsize=(10,8));ax.imshow(view(wall,scalar,clim));clean(ax,title);cmapbar(fig,ax,*clim,label)
 save(fig,f'{idx:02d}_'+('W0_eta' if idx==2 else 'wall_normal_fraction'),'Each WALL triangle inherits its owning tetrahedron. Normal fraction is evaluated INSIDE, not on zero-velocity WALL nodes. Geometry: mm.')
for idx,key,label in [(4,'abs_div_s_inv',r'$|\nabla\cdot u|$ (s$^{-1}$)'),(5,'eta_div',r'$\eta$ (dimensionless)')]:
 fig,axs=plt.subplots(1,2,figsize=(12,5))
 for k in range(3):
  v=f.loc[f.layer==k,key].values;pos=v[v>0];bins=np.geomspace(max(pos.min(),1e-9),pos.max(),90) if idx==4 else np.linspace(0,1.75,90)
  axs[0].hist(v,bins=bins,histtype='step',density=True,label=f'W{k}');sv=np.sort(v);axs[1].plot(sv,np.arange(1,len(v)+1)/len(v),label=f'W{k}')
 for ax in axs:ax.set_xlabel(label);ax.legend();ax.grid(alpha=.2)
 axs[0].set_ylabel('Probability density');axs[1].set_ylabel('Cell-count ECDF')
 if idx==4:
  for ax in axs:ax.set_xscale('symlog',linthresh=1e-3)
 caption='All cells in ECDF; positive |div| in logarithmic histogram. Exact-zero counts W0/W1/W2: '+str([int(np.sum(f.loc[f.layer==k,key]==0)) for k in range(3)]) if idx==4 else 'All cells of each layer included; count weighted. Screening thresholds are not physical CFD acceptance criteria.'
 save(fig,f'{idx:02d}_'+('divergence_distributions' if idx==4 else 'eta_distributions'),caption)
for idx,key,label in [(6,'h_max_m','Maximum edge (mm)'),(7,'mean_ratio_quality','Mean-ratio shape quality (regular tetra = 1)'),(8,'normal_fraction','Interior centroid normal fraction')]:
 fig,axs=plt.subplots(1,3,figsize=(15,4.7))
 for k,ax in enumerate(axs):
  d=f[f.layer==k];xx=d[key].values*(1000 if idx==6 else 1);hb=ax.hexbin(xx,d.eta_div,gridsize=60,mincnt=1,bins='log',cmap='magma');ax.set_xlabel(label);ax.set_ylabel(r'$\eta$');ax.set_title(f'W{k}');fig.colorbar(hb,ax=ax,label='Cell count (log color)')
 save(fig,f'{idx:02d}_eta_vs_'+{6:'size',7:'quality',8:'normal_fraction'}[idx],'Every cell contributes; hexagons summarize density, no random point selection. W0 normals use largest owned WALL face.')
ids=np.flatnonzero(cl['eta_ge_0p5']==0);mesh=m.extract_cells(ids);mesh.cell_data['eta']=f.eta_div.values[ids]
fig,axs=plt.subplots(1,2,figsize=(13,6));axs[0].imshow(view(mesh,'eta',(.5,1.05),True));clean(axs[0],'Largest high-eta component: vessel context');axs[1].imshow(view(mesh,'eta',(.5,1.05),False,edges=True));clean(axs[1],f'{len(ids):,} face-connected tetrahedra');cmapbar(fig,list(axs),.5,1.05,r'$\eta$')
save(fig,'09_largest_eta_cluster','Largest component by cell count within W0-W2, eta >= 0.5. This is a screening cluster, not a diagnosed physical defect. Geometry: mm.')
for idx,tid in [(10,388769),(11,138584)]:
 reg=json.loads((out/f'known_regression_{tid}.json').read_text());adj=z['adjacency'];nb=set([tid]);front=[tid]
 for _ in range(2):
  front=[int(j) for i in front for j in adj[i] if j>=0 and j not in nb];nb.update(front)
 ids=np.array(sorted(nb));mesh=m.extract_cells(ids);mesh.cell_data['eta']=f.eta_div.values[ids]
 p=pv.Plotter(off_screen=True,window_size=(1600,1200));p.set_background('white');p.add_mesh(mesh,scalars='eta',clim=(0,1.05),show_edges=True,opacity=.6,show_scalar_bar=False)
 target=m.extract_cells([tid]);p.add_mesh(target,style='wireframe',color='red',line_width=5)
 if tid==138584:
  zeroid=reg['zero_component_id'];zero=np.flatnonzero(cl['all_zero_velocity']==zeroid);p.add_mesh(m.extract_cells(zero),color='limegreen',opacity=.65,show_edges=True)
 p.add_point_labels(np.mean(z['xyz'][z['connectivity'][tid]],axis=0).reshape(1,3)*1000,[str(tid)],font_size=22,point_size=6,always_visible=True)
 p.view_isometric();p.add_axes();im=p.screenshot(return_img=True);p.close()
 fig,axs=plt.subplots(1,2,figsize=(13,6),gridspec_kw={'width_ratios':[2,1]});axs[0].imshow(im);clean(axs[0],f'Tetra {tid}: two face-neighbor rings')
 text=f"Canonical zero-based ID: {tid}\n\n"+f"div(u) = {reg['metrics']['div_s_inv']:.6g} s^-1\n||G|| = {reg['metrics']['G_norm_s_inv']:.6g} s^-1\neta = {reg['metrics']['eta_div']:.6g}\n\n"+f"Node speeds (m/s):\n{np.linalg.norm(reg['velocity_m_s'],axis=1)}\n\nPressure (Pa):\n{np.round(reg['pressure_Pa'],3)}\n\nWALL distance = {reg['metrics']['d_wall_m']*1e6:.3f} um\nZero component cells: {reg['zero_component_cell_count']}"
 axs[1].text(.02,.95,text,va='top',fontsize=11,transform=axs[1].transAxes);clean(axs[1],'Raw-node independent reconstruction')
 cmapbar(fig,axs[0],0,1.05,r'$\eta$');save(fig,f'{idx:02d}_tetra_{tid}','Red edges: target tetra. Green: connected all-four-node-zero patch (figure 11 only). No geometry or field alteration. Geometry: mm.')
owner=z['owner'];wall.cell_data['log_visit']=np.log10(1+vis.state_count.values[owner]);fig,axs=plt.subplots(1,2,figsize=(13,6))
axs[0].imshow(view(wall,'log_visit',(0,float(wall['log_visit'].max()))));clean(axs[0],'Actual particle states in W0 owners');cmapbar(fig,axs[0],0,float(wall['log_visit'].max()),'log10(1 + saved-state count)')
mask=(f.layer<=2)&(vis.state_count>0);hb=axs[1].hexbin(f.eta_div[mask],vis.state_count[mask],gridsize=60,yscale='log',bins='log',mincnt=1,cmap='magma');axs[1].set_xlabel('eta (dimensionless)');axs[1].set_ylabel('Saved particle states per visited cell');axs[1].set_title('Visited W0-W2 cells');fig.colorbar(hb,ax=axs[1],label='Cell count')
save(fig,'12_particle_visitation_overlap','Historical realization_0000_contactfix only. Accepted-time state counts are not time weights; adaptive steps can inflate counts.')
js(out/'figures/FIGURE_INDEX.json',figures)
