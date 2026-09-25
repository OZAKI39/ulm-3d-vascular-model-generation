#!/usr/bin/env python3
from pathlib import Path
import csv,json,collections
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,FancyArrowPatch,Circle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
S=Path(__file__).resolve().parents[1];V=S/'visualization';V.mkdir(exist_ok=True);OLD=S.parent/'wallslide_20260916_233923';plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':130})
D='NOT EXPERIMENTAL CONCENTRATION'
def save(fig,name):
 fig.text(.5,.012,D+' | Geometry / kinematics only; wall hydrodynamics NOT IMPLEMENTED',ha='center',fontsize=8,color='#555555');fig.savefig(V/name,dpi=180,bbox_inches='tight');plt.close(fig)
def read(p):
 with p.open() as f:return list(csv.DictReader(f))
def array(rows,names):return np.array([[float(r[k]) for k in names] for r in rows])
root=json.loads((S/'validation/OLD_AMBIGUITY_ROOT_CAUSE.json').read_text());q=root['reference'];x=np.array(root['old']['query_point']);origin=np.array(q['clusters'][0]['point']);tri=np.array([c['triangle_coordinates'] for c in q['candidates']]);pseudonormal=np.array(q['clusters'][0]['normal']);off=x-origin;off/=np.linalg.norm(off)
for filename,zoom in [('old_long_ambiguity_geometry.png',False),('triangle_candidates_and_surface_clusters.png',True)]:
 fig=plt.figure(figsize=(9,6));ax=fig.add_subplot(111,projection='3d');colors=['#337eae','#53abc2'];local=(tri-origin)*1e6
 ax.add_collection3d(Poly3DCollection(local,facecolors=colors,edgecolors='#16465f',alpha=.7))
 for i,t in enumerate(local):ax.text(*t.mean(axis=0),f"triangle {q['candidates'][i]['id']}\ncluster 11572",fontsize=9)
 if not zoom:
  center=(x-origin)*1e6;a=root['old']['radius']*1e6;u=np.linspace(0,2*np.pi,36);v=np.linspace(0,np.pi,24);xx=center[0]+a*np.outer(np.cos(u),np.sin(v));yy=center[1]+a*np.outer(np.sin(u),np.sin(v));zz=center[2]+a*np.outer(np.ones_like(u),np.cos(v));ax.plot_wireframe(xx,yy,zz,color='#d8974b',alpha=.12,rstride=3,cstride=3);ax.scatter(*center,color='#b87920',s=25);ax.plot([0,center[0]],[0,center[1]],[0,center[2]],color='#b87920');ax.text(*center,'bubble center',fontsize=10);bound=1.7
 else:bound=.23
 for vec,col,label in [(pseudonormal,'#70379a','mesh pseudonormal'),(off,'#b87920','offset contact normal')]:ax.quiver(0,0,0,*(vec*.17),color=col,linewidth=2,label=label)
 ax.set(xlim=(-bound,bound),ylim=(-bound,bound),zlim=(-bound,bound),xlabel='local x (µm)',ylabel='local y (µm)',zlabel='local z (µm)');ax.set_box_aspect((1,1,1));ax.view_init(22,-45);ax.legend(loc='upper left',fontsize=8)
 ax.set_title('Old LONG failure: two triangles, ONE smooth patch\nshared manifold edge · dihedral 2.227° · same biological region');save(fig,filename)
# Smooth-edge normal explanation, including the finite-radius distinction.
fig,ax=plt.subplots(figsize=(9,5));ax.plot([-2,0,2],[.25,0,.25],color='#276786',lw=5);ax.fill_between([-2,0,2],[-.6,-.6,-.6],[.25,0,.25],color='#cbdde5');center=np.array([.4,1.45]);ax.add_patch(Circle(center,.52,fill=False,color='#b87920',lw=2));ax.scatter(0,0,c='k',s=25)
for p,v,col,label in [(np.array([-.9,.11]),np.array([.08,.7]),'#4a91ae','face n₁'),(np.array([.9,.11]),np.array([-.08,.7]),'#4a91ae','face n₂'),(np.zeros(2),np.array([0,.9]),'#70379a','weighted mesh pseudonormal'),(np.zeros(2),center*.65,'#b87920','sphere offset distance gradient')]:ax.annotate('',xy=p+v,xytext=p,arrowprops=dict(arrowstyle='->',lw=2,color=col));ax.text(*(p+v),label,fontsize=9,color=col)
ax.text(-1.9,1.8,'One smooth patch\nweighted face/one-ring normals orient the mesh\ncenter-to-feature gradient constrains a finite sphere',fontsize=11);ax.set(xlim=(-2.1,2.3),ylim=(-.65,2.25));ax.set_aspect('equal');ax.axis('off');save(fig,'smooth_edge_pseudonormal_scheme.png')
# True multi-surface example.
fig,ax=plt.subplots(figsize=(8,5));ax.plot([0,0,2],[2,0,0],lw=5,color='#315a74');ax.add_patch(Circle((.65,.65),.5,fill=False,color='#b87920',lw=2));ax.annotate('n₁',xy=(.8,.65),xytext=(.08,.65),arrowprops=dict(arrowstyle='->',color='#326b9d',lw=2));ax.annotate('n₂',xy=(.65,.8),xytext=(.65,.08),arrowprops=dict(arrowstyle='->',color='#bb6230',lw=2));ax.text(1,1.55,'Distinct patches stay distinct\nmin ½‖V − Vraw‖²\nsubject to N V ≥ 0\nactive-set rank ≤ 3',fontsize=12);ax.text(.6,-.25,'No averaging across a sharp feature');ax.set(xlim=(-.3,3.1),ylim=(-.4,2.3));ax.set_aspect('equal');ax.axis('off');save(fig,'true_multi_surface_scheme.png')
fig,axes=plt.subplots(1,2,figsize=(12,5))
for ax in axes:ax.set_xlim(0,10);ax.set_ylim(0,8);ax.axis('off')
ax=axes[0];verts=np.array([[.5,4],[3,4.5],[5.5,4],[.5,5.5],[3,6],[5.5,5.5]])
for i,ids in enumerate([[0,1,3],[1,3,4],[1,2,4],[2,4,5]]):ax.add_patch(Polygon(verts[ids],facecolor=['#95c8de','#b8dae7'][i%2],edgecolor='#21607d'));ax.text(*verts[ids].mean(axis=0),f'triangle {chr(65+i)}',ha='center',fontsize=9)
ax.text(3,3.2,'↓',ha='center',fontsize=25);ax.text(3,2.5,'same continuous local patch',ha='center',fontsize=12);ax.text(3,1.9,'↓',ha='center',fontsize=25);ax.text(3,1,'one stable pseudonormal\n→ one sphere contact constraint',ha='center',fontsize=12);ax.set_title('Mesh discretization tie',fontsize=15)
ax=axes[1];ax.add_patch(Polygon([[1,4],[1,6.5],[3.5,5.5],[3.5,3]],facecolor='#95c8de',edgecolor='#21607d'));ax.add_patch(Polygon([[3.5,3],[3.5,5.5],[6,6.5],[6,4]],facecolor='#efc898',edgecolor='#a0611f'));ax.text(1.4,5,'patch 1');ax.text(4.1,5,'patch 2');ax.text(3.5,2.5,'↓',ha='center',fontsize=25);ax.text(3.5,1.65,'n₁ + n₂: independent constraints',ha='center',fontsize=12);ax.text(3.5,.8,'minimum-change multi-normal projection',ha='center',fontsize=11);ax.set_title('True distinct surfaces',fontsize=15)
fig.suptitle('TRIANGLES ≠ PHYSICAL WALLS',fontsize=20,fontweight='bold');save(fig,'continuous_surface_contact_3d_scheme.png')
# Real topology distribution.
a=np.load(S/'raw/TOPOLOGY_ARRAYS.npz')['manifold_edges'][:,-1];fig,(ax,bx)=plt.subplots(1,2,figsize=(11,4));ax.hist(a,bins=np.linspace(0,30,61),color='#4c9db6');ax.set(xlabel='adjacent-face dihedral (degrees)',ylabel='manifold edges',title='Smooth continuation range');
for d,col in [(10,'#d69b38'),(15,'#924aa4'),(25,'#b43d3d')]:ax.axvline(d,color=col,label=str(d)+'°')
ax.legend();bx.hist(a,bins=np.linspace(0,180,73),color='#4c9db6');bx.set_yscale('log');bx.set(xlabel='dihedral (degrees)',ylabel='edge count (log)',title='Full geometry: sharp features retained');save(fig,'dihedral_angle_distribution.png')
new=read(S/'runs/LONG_TRANSPORT_mpi1/TRAJECTORIES.csv');old=read(OLD/'runs/LONG_TRANSPORT_release_mpi1/TRAJECTORIES.csv');contacts=read(S/'runs/LONG_TRANSPORT_mpi1/WALL_SURFACE_CONTACT_EVENTS.csv');accepted=[r for r in contacts if r['accepted']=='1'];hist=read(S/'runs/LONG_TRANSPORT_mpi1/SOLVER_HISTORY.csv');state=json.loads((S/'runs/LONG_TRANSPORT_mpi1/RUN_STATE.json').read_text())
fig,ax=plt.subplots(figsize=(8,4));counts=collections.Counter(int(r['surface_cluster_count']) for r in accepted);ax.bar([str(k) for k in counts],[counts[k] for k in counts],color='#4c9db6');ax.set(xlabel='local surface cluster count',ylabel='accepted RK2 particle stages',title='Real LONG: triangle ties remain a single physical patch');ax.text(.5,.85,f"{len(accepted):,} independently checked stages\nTrue multi-surface solver qualified separately on wedges/corner",transform=ax.transAxes,ha='center');save(fig,'surface_cluster_count_distribution.png')
fig=plt.figure(figsize=(10,6));ax=fig.add_subplot(111,projection='3d')
for dataset,style,prefix in [(old,'--','previous'),(new,'-','new')]:
 for i,col in [(1,'#276eab'),(2,'#c07820')]:
  rr=[r for r in dataset if int(r['particle_id'])==i];xyz=array(rr,['x_m','y_m','z_m'])*1e6;ax.plot(*xyz.T,style,color=col,label=f'{prefix}, particle {i}',lw=1.5);ax.scatter(*xyz[-1],color=col,s=20)
ax.set(xlabel='x (µm)',ylabel='y (µm)',zlabel='z (µm)',title=f"LONG: old ambiguity 0.494123 s → new pair blocker {state['time_s']:.6f} s");ax.legend(fontsize=8);save(fig,'LONG_old_vs_new_trajectory.png')
fig,ax=plt.subplots(figsize=(10,4))
for i in [1,2]:
 rr=[r for r in new if int(r['particle_id'])==i];ax.semilogy([float(r['time_s']) for r in rr],[float(r['surface_gap_m'])*1e9 for r in rr],label=f'particle {i}')
ax.axhline(.1,color='k',ls=':',label='inherited 0.1 nm margin');ax.axvline(.494122765648,color='#963f83',ls='--',label='old ambiguity');ax.set(xlabel='physical time (s)',ylabel='wall surface gap (nm)',title='Accepted LONG states maintain wall clearance');ax.legend();save(fig,'LONG_gap_vs_time.png')
fig,ax=plt.subplots(figsize=(10,4));modes=['NONE','SINGLE_FACE','SMOOTH_EDGE','SMOOTH_VERTEX','SINGLE_SURFACE_PATCH','TRUE_MULTI_SURFACE','BOUNDED_POSITION_CORRECTION']
for i,col in [(1,'#276eab'),(2,'#c07820')]:
 rr=[r for r in accepted if int(r['particle_id'])==i and r['stage']=='1'];ax.scatter([float(r['time']) for r in rr],[modes.index(r['constraint_mode'])+.08*(i-1) for r in rr],s=3,color=col,label=f'particle {i}')
ax.set_yticks(range(len(modes)),modes);ax.axvline(.494122765648,color='#963f83',ls='--');ax.set(xlabel='physical time (s)',title='Accepted final RK2 stage constraint mode');ax.legend();save(fig,'LONG_constraint_mode_vs_time.png')
fig,ax=plt.subplots(figsize=(10,4));ax.semilogy([float(r['time_s']) for r in hist],[float(r['dt_s']) for r in hist],color='#246b93',lw=.8);ax.axvline(.494122765648,color='#963f83',ls='--',label='old ambiguity');ax.axvline(state['time_s'],color='#bb493d',ls=':',label='pair safety blocker');ax.set(xlabel='physical time (s)',ylabel='accepted dt (s)',title='Physical-time refinement; rejected trials do not advance time');ax.legend();save(fig,'LONG_dt_history.png')
synthetic=json.loads((S/'raw/SYNTHETIC_TRAJECTORIES.json').read_text());fig,axes=plt.subplots(1,3,figsize=(12,4))
for ax,k in zip(axes,range(3)):
 data=np.load(S/'raw/synthetic'/f'plane{k}.npz');v=data['vertices']*1e6;ax.triplot(v[:,0],v[:,1],data['faces'],color='#aac3cf',lw=.5);a=next(r for r in synthetic if r['name']==f'plane{k}' and abs(r['time']-.01)<1e-12);tr=np.array(a['trajectory']);ax.plot(tr[:,1]*1e6,tr[:,2]*1e6,color='#c17426',lw=3);ax.set(xlim=(-2,2),ylim=(-2,2),title=f'Triangulation {k+1}',xlabel='x (µm)',ylabel='y (µm)');ax.set_aspect('equal')
fig.suptitle('Same plane, same constrained trajectory; triangle reordering also invariant');save(fig,'triangle_tessellation_invariance.png')
sens=json.loads((S/'validation/REAL_LONG_TOLERANCE_SENSITIVITY.json').read_text());fig,axes=plt.subplots(1,3,figsize=(12,4));keys=['max_position_difference_m','max_used_velocity_difference_m_s','mode_changes'];titles=['max trajectory difference (m)','max used-velocity difference (m/s)','constraint-mode differences']
for ax,key,title in zip(axes,keys,titles):
 z=np.array([[next(r[key] for r in sens['cases'] if r['tie_factor']==f and r['smooth_deg']==a) for a in [10,15,25]] for f in [.5,1,2]]);im=ax.imshow(z,cmap='Blues',vmin=0,vmax=max(z.max(),1e-30));ax.set_xticks(range(3),['10°','15°','25°']);ax.set_yticks(range(3),['0.5×','1×','2×']);ax.set(xlabel='smooth dihedral threshold',ylabel='tie tolerance factor',title=title)
 for i in range(3):
  for j in range(3):ax.text(j,i,f'{z[i,j]:.1g}',ha='center',va='center')
fig.suptitle('Nine actual LONG integrations from t=0 to 0.51 s');save(fig,'tolerance_sensitivity.png')
print('VISUALIZATIONS',len(list(V.glob('*.png'))))
