"""Figures from exact accepted coordinates and event topology, with no surface synthesis."""
from pathlib import Path
import sys,json,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from design import js,save
from metrics import Topology
RED='#a93040';BLUE='#2375a3';TEAL='#208e84';GOLD='#c58a2c';GRAY='#777777';PURPLE='#90549d'
COLORS=[TEAL,BLUE,GOLD,RED]
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.titlesize':13,'axes.labelsize':12,'legend.fontsize':9,'axes.grid':True,'grid.color':'#989898','grid.alpha':.5,'svg.fonttype':'none'})

def view(ax,title,limits=None):
 ax.view_init(elev=14,azim=-56);ax.set(xlabel='x (mm)',ylabel='y (mm)',zlabel='z (mm)',title=title,xticks=[-.3,0,.3],yticks=[-.3,0,.3]);ax.tick_params(labelsize=10,pad=1)
 if limits is None:limits=[[-.34,.34],[-.34,.34],[-.54,.54]]
 ax.set(xlim=limits[0],ylim=limits[1],zlim=limits[2]);ax.set_box_aspect([q[1]-q[0] for q in limits])
 for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.set_facecolor((.96,.96,.96,1));axis._axinfo['grid'].update(color=(.5,.5,.5,.55),linewidth=.65)
def notch(ax):
 ax.plot([-.3,-.1,-.1,-.3,-.3],[-.3,-.3,.3,.3,-.3],[0]*5,c='forestgreen',lw=1.6);ax.plot([-.1,-.1],[-.3,.3],[0,0],c='limegreen',lw=2.3)
def figsave(fig,out,name,svg=False):
 fig.savefig(out/(name+'.png'),dpi=200,bbox_inches='tight',facecolor='white')
 if svg:fig.savefig(out/(name+'.svg'),bbox_inches='tight',facecolor='white')
 plt.close(fig)
def state(base,name,key):return dict(np.load(base/'analysis'/name/'states'/(key+'.npz')))
def deformed(ax,g,z,title,lim=None,components=False):
 X=g['X'];b=g['bonds'];x=z['x_m']*1e3;alive=z['alive'].astype(bool);mask=g['prescribed_mask'][:,2]
 ax.scatter(*x.T,c=z['component_labels'] if components else RED,cmap='tab10' if components else None,s=7,alpha=.75,depthshade=False)
 if not components:ax.scatter(*x[mask].T,c=BLUE,s=13,alpha=.9,depthshade=False)
 cut=alive&(abs(X[b[:,0],1])<55e-6)&(abs(X[b[:,1],1])<55e-6);ax.add_collection3d(Line3DCollection(x[b[cut]],colors='#68615c',linewidths=.32,alpha=.38));view(ax,title,lim)

def main(base,figs):
 figs.mkdir(exist_ok=True,parents=True);Q=js(base/'CRACK_PROPAGATION_AUDIT.json');B=Q['M5B'];C=Q['M5C'];area=js(base/'CRACK_AREA_ATTEMPT.json');kin=js(base/'GRIP_KINEMATIC_AUDIT.json')['kinematics'];g={n:dict(np.load(base/'analysis'/n/'geometry_events.npz')) for n in ['M5B','M5C']};h={n:dict(np.load(base/'analysis'/n/'review_curves.npz')) for n in g};bd={n:dict(np.load(base/'analysis'/n/'boundary_review.npz')) for n in g};k=dict(np.load(base/'KINEMATIC_BOND_SOURCE.npz'));X=k['X']*1e3;b=k['bonds'];tend=C['end_time_s']*1e3
 # Fig1: real group membership and every dual-axial-prescribed bond; no assumed fixed xy.
 fig=plt.figure(figsize=(16,11));fig.subplots_adjust(top=.9,bottom=.14,hspace=.3,wspace=.28)
 for j,n in enumerate(g):
  ax=fig.add_subplot(2,3,j+1,projection='3d');m=g[n]['prescribed_mask'][:,2];ax.scatter(*X[~m].T,c=RED,s=6,alpha=.45,depthshade=False);ax.scatter(*X[m].T,c=BLUE,s=13,alpha=.9,depthshade=False);notch(ax)
  for zz in [-.35,.35]:ax.plot([-.3,.3,.3,-.3,-.3],[-.3,-.3,.3,.3,-.3],[zz]*5,c=GRAY,ls='--',lw=1)
  view(ax,f'{n}: {m.sum()} axial-controlled nodes\n0 fully prescribed nodes; all x/y dynamic')
  ax=fig.add_subplot(2,3,j+4,projection='3d');pp=k[n+'_pp'];v=k[n+'_prescribed_axial_lower_bound_max_s'];ax.add_collection3d(Line3DCollection(X[b[pp]],colors=plt.cm.viridis(Normalize(-1,.8)(v[pp])),linewidths=.32,alpha=.6));notch(ax);view(ax,f'{n}: {pp.sum():,} dual-controlled bonds\n{kin[n]["guaranteed_failure_pairs"]:,} axially forced to fail')
 ax=fig.add_subplot(2,3,3)
 for n,c in [('M5B',GOLD),('M5C',BLUE)]:v=k[n+'_prescribed_axial_lower_bound_max_s'][k[n+'_pp']];ax.hist(v,bins=np.linspace(-1,1,41),alpha=.6,label=n,color=c)
 ax.axvline(.15,color='black',ls='--',label='s00=0.15');ax.set(xlabel='Maximum axial lower bound on stretch',ylabel='Dual-controlled bond count',title='Full 24 ms kinematic audit');ax.legend()
 ax=fig.add_subplot(2,3,6);t=np.linspace(0,24,300);d=.00068*.5*(1-np.cos(np.pi*t/24));ax.plot(t,d/.00085,c=GOLD,label='M5B vertical-bond lower bound');ax.axhline(.15,c='black',ls='--',label='Failure threshold');ax.axvline(kin['M5B']['threshold_s']*1e3,c=GRAY,ls=':');ax.plot(t,np.zeros_like(t),c=BLUE,label='M5C fixed-reference-x/y proxy');ax.set(xlabel='Planned physical time (ms)',ylabel='Stretch',title='M5C actual x/y response is NOT prescribed');ax.legend(fontsize=8)
 cb=fig.colorbar(plt.cm.ScalarMappable(norm=Normalize(-1,.8),cmap='viridis'),cax=fig.add_axes([.15,.045,.4,.015]),orientation='horizontal');cb.set_label('Bond colors: maximum axial lower bound on actual stretch')
 fig.suptitle('1 | Real point lattice and prescribed axial kinematics\nBlue: z-controlled; red: all DOFs free; green: original pre-cut; dashed: inherited grip planes',fontsize=16);figsave(fig,figs,'Fig1_grips_and_kinematic_bounds')
 # Fig2: actual pre-first utilization, actual first markers and first coherent cluster.
 fig=plt.figure(figsize=(15,11));fig.subplots_adjust(top=.86,bottom=.08,hspace=.35,wspace=.24)
 for j,n in enumerate(g):
  q=Q[n];z=state(base,n,'pre_first_break');gg=g[n];bb=gg['bonds'];u=z['utilization'];node=np.zeros(len(X));eligible=z['alive'].astype(bool);np.maximum.at(node,bb[eligible].ravel(),np.repeat(u[eligible],2));ax=fig.add_subplot(2,2,j+1,projection='3d',computed_zorder=False);sc=ax.scatter(*X.T,c=node,cmap='magma',vmin=0,vmax=1.01,s=8,depthshade=False,zorder=1);notch(ax);ids=np.array([e['bond_index_zero_based'] for e in q['first_events']]);ax.add_collection3d(Line3DCollection(X[bb[ids]],colors='cyan',linewidths=3,zorder=20));mid=X[bb[ids]].mean(1);ax.scatter(*mid.T,c='cyan',marker='*',s=170,edgecolors='black',depthshade=False,zorder=21);view(ax,f'{n} accepted pre-first: {float(z["time_s"])*1e3:.4f} ms\nCyan: next actual first bond(s) at {q["first_time_s"]*1e3:.4f} ms')
 fig.colorbar(sc,ax=fig.axes[:2],shrink=.65,pad=.01,label='Max adjacent live-bond stretch / s00')
 ax=fig.add_subplot(2,2,3);rb=dict(np.load(base/'analysis/M5B/risk_review.npz'))
 for vals,timekey,n,c,ls in [(rb,'time_s','M5B',GOLD,'--'),(bd['M5C'],'time_s','M5C',BLUE,'-')]:
  for region,alpha in [('notch',1),('grip',.55)]:key=region+('_including_newbreak_max' if n=='M5B' else '_pre_event_util');keep=vals[timekey]*1e3<=tend;ax.plot(vals[timekey][keep]*1e3,vals[key][keep],ls,c=c,alpha=alpha,label=n+' '+region)
 ax.axhline(1,c='black',ls=':');ax.set(xlabel='Physical time (ms)',ylabel='Alive + newly broken utilization',title='Competition under unchanged region definitions');ax.legend(ncol=2)
 ax=fig.add_subplot(2,2,4,projection='3d');z=state(base,'M5C','first_coherent_cluster');dead=~z['alive'].astype(bool);mid=X[b].mean(1);ax.scatter(*X.T,c='#dddddd',s=3,alpha=.2);ax.scatter(*mid[dead].T,c=np.array(COLORS)[g['M5C']['region'][dead]],s=23,depthshade=False);notch(ax);view(ax,f'M5C first coherent event: {float(z["time_s"])*1e3:.4f} ms\nReference broken midpoints; cluster is near grip')
 fig.suptitle('2 | First-event evidence and competing damage clusters\nAccepted stretch shown on reference coordinates; cyan markers overlaid for visibility',fontsize=16);figsave(fig,figs,'Fig2_first_break_and_competition')
 # Fig3: exact full fields, common axes, no displacement magnification.
 keys=['initial','pre_first_break','first_break','first_coherent_cluster','peak_tensile','end'];keys.sort(key=lambda s:C['states'][s]['step']);zz={key:state(base,'M5C',key) for key in keys};pts=np.concatenate([z['x_m'] for z in zz.values()])*1e3;lo=pts.min(0)-.05;hi=pts.max(0)+.05;lim=list(zip(lo,hi));fig=plt.figure(figsize=(17,11));fig.subplots_adjust(top=.89,bottom=.09,wspace=.14,hspace=.25)
 for j,key in enumerate(keys):
  z=zz[key];ax=fig.add_subplot(2,3,j+1,projection='3d');deformed(ax,g['M5C'],z,f'{key.replace("_"," ")} | {float(z["time_s"])*1e3:.4f} ms\n{np.sum(~z["alive"].astype(bool)):,} broken bonds',lim)
 fig.legend(handles=[Line2D([],[],color=RED,marker='o',ls='',label='Actual material points'),Line2D([],[],color=BLUE,marker='o',ls='',label='Axial-controlled end nodes'),Line2D([],[],color=GRAY,label='Actual surviving bonds in central y slab')],loc='lower center',ncol=3);fig.suptitle('3 | M5C accepted three-dimensional states; deformation scale = 1\nNo surface reconstruction; all states retained through the frozen stop at 7.05 ms',fontsize=16);figsave(fig,figs,'Fig3_exact_accepted_states')
 # Fig4: defined boundary observables, energy and dynamics.
 fig,axs=plt.subplots(2,2,figsize=(14,10),layout='constrained')
 for n,col in [('M5B',GOLD),('M5C',BLUE)]:
  a=h[n];keep=a['time_s']*1e3<=tend;axs[0,0].plot(a['time_s'][keep]*1e3,a['displacement_m'][keep]*1e6,color=col,label=n);axs[0,1].plot(a['displacement_m'][keep]*1e6,a['stress_Pa'][keep],color=col,label=n);q=bd[n];keep=q['time_s']*1e3<=tend;Qkey='generalized_external_Q_N' if n=='M5B' else 'generalized_Q_N';axs[1,0].plot(q['time_s'][keep]*1e3,q['top_Rz_N'][keep]*1e6,c=col,label=n+' top R');axs[1,0].plot(q['time_s'][keep]*1e3,q[Qkey][keep]*1e6,c=col,ls='--',label=n+' Q')
 axs[0,0].set(xlabel='Physical time (ms)',ylabel='Mean controlled-boundary opening (µm)',title='Same mean opening; H changes from 850 to 950 µm');axs[0,1].set(xlabel='Actual mean boundary opening (µm)',ylabel='Top reaction / A0 (Pa)',title='Specimen response, not a constitutive material curve');axs[1,0].set(xlabel='Physical time (ms)',ylabel='Force (µN)',title='Q includes all prescribed DOFs')
 a=h['M5C']
 for key,label,col in [('W_ext_J','Boundary work',BLUE),('U_bulk_J','Bulk U',TEAL),('K_free_J','Free K',PURPLE),('D_num_J','Numeric damping',GRAY),('D_break_independent_J','Signed topology term',RED),('R_energy_J','Balance residual',GOLD)]:axs[1,1].plot(a['time_s']*1e3,a[key]*1e9,label=label,c=col)
 axs[1,1].set(xlabel='Physical time (ms)',ylabel='Energy (nJ)',title='M5C independent ledger: no clamping of negative releases')
 for ax in axs.flat:ax.legend(fontsize=9)
 fig.suptitle('4a | Matched physical time and mean opening; transmitted work is integrated over all controlled DOFs',fontsize=15);figsave(fig,figs,'Fig4a_boundary_force_work',True)
 fig,axs=plt.subplots(2,2,figsize=(14,10),layout='constrained')
 for j,n in enumerate(['M5B','M5C']):
  a=h[n]
  for key,label,col in [('W_ext_J','W',BLUE),('U_bulk_J','U bulk',TEAL),('K_free_J','K free',PURPLE),('D_num_J','D numerical',GRAY),('D_break_independent_J','Signed topology',RED),('R_energy_J','Residual',GOLD)]:axs[0,j].plot(a['time_s']*1e3,a[key]*1e9,c=col,label=label)
  axs[0,j].set(xlabel='Physical time (ms)',ylabel='Energy (nJ)',title=n+' complete retained history');axs[0,j].axvline(tend,c='black',ls=':');axs[0,j].legend(ncol=2)
  axs[1,0].plot(a['time_s']*1e3,a['KE_over_U']*100,c=GOLD if n=='M5B' else BLUE,label=n)
 a=h['M5C'];inc=np.diff(a['D_break_independent_J'],prepend=0)*1e9;axs[1,1].plot(a['time_s']*1e3,inc,c=RED,lw=.65,label='Signed event increment');axs[1,1].plot(a['time_s']*1e3,a['R_energy_J']*1e9,c=TEAL,label='Cumulative balance residual');axs[1,1].set(xlabel='Physical time (ms)',ylabel='Energy (nJ)',title='Negative topology increments are preserved');axs[1,1].legend()
 axs[1,0].set(xlabel='Physical time (ms)',ylabel='Free kinetic / bulk energy (%)',title='Startup and post-break inertia both disclosed',yscale='symlog');axs[1,0].axhline(5,c='black',ls='--');axs[1,0].legend();fig.suptitle('4b | Energy closure does not identify irreversible fracture energy',fontsize=16);figsave(fig,figs,'Fig4b_complete_energy_and_dynamics',True)
 # Fig5: reference damage region, surviving crossing map, strips and actual components.
 z=state(base,'M5C','end');gg=g['M5C'];T=Topology(gg,gg['prescribed_mask']);dead=~z['alive'].astype(bool);mid=T.mid*1e3
 fig=plt.figure(figsize=(16,11));fig.subplots_adjust(top=.88,bottom=.08,wspace=.32,hspace=.35)
 ax=fig.add_subplot(2,3,1,projection='3d');ax.scatter(*X.T,c='#ddd',s=2,alpha=.15);ax.scatter(*mid[dead].T,c=np.array(COLORS)[T.reg[dead]],s=9,alpha=.7,depthshade=False);notch(ax);view(ax,'All new breaks by inherited region\nReference midpoints, no surface');ax.legend(handles=[Line2D([],[],c=c,marker='o',ls='',label=n) for n,c in zip(['Notch','Grip','Surface','Interior'],COLORS)],fontsize=8,loc='upper left')
 for j,(vals,title) in enumerate([(z['initial_crossings'],'Initial surviving cross-z0 bonds'),(z['remaining_crossings'],'Surviving cross-z0 bonds at stop')]):
  ax=fig.add_subplot(2,3,j+2);im=ax.imshow(vals.T,origin='lower',extent=[-100,300,-300,300],aspect='auto',cmap='viridis',vmin=0,vmax=z['initial_crossings'].max());ax.set(xlabel='Reference crossing x (µm)',ylabel='Reference crossing y (µm)',title=title);ax.axhspan(-100,100,fill=False,edgecolor='white',lw=1.5);fig.colorbar(im,ax=ax,label='Actual surviving connections')
 ax=fig.add_subplot(2,3,4);y=np.arange(12)*50-275;ax.plot(z['advance_m']*1e6,y,'o-',c=BLUE,label='M5C end');ax.plot(np.array(B['matched_common_time']['front_advance_m'])*1e6,y,'s--',c=GOLD,label='M5B same time');ax.axhspan(-100,100,color=TEAL,alpha=.1);ax.set(xlabel='Contiguous loss prefix from notch (µm)',ylabel='Thickness strip center y (µm)',title='Connectivity advance; never physical crack area',xlim=(-10,410));ax.legend()
 ax=fig.add_subplot(2,3,5,projection='3d');survive=T.cross&~dead;ax.add_collection3d(Line3DCollection(X[b[survive]],colors=TEAL,linewidths=.42,alpha=.35));notch(ax);view(ax,f'{survive.sum():,} residual cross-z0 bridges\nAll genuine survivors; none removed in analysis')
 ax=fig.add_subplot(2,3,6,projection='3d');deformed(ax,gg,z,'Actual connected components at stop\n'+str(C['states']['end']['component_sizes']),lim,components=True)
 fig.suptitle('5 | Gc unavailable / crack area unresolved\nFinite-horizon damage clouds and crossing depletion do not certify a unique new material surface',fontsize=16);figsave(fig,figs,'Fig5_topology_and_area_limits')
 # Compact plot provenance is unrounded exact arrays in analysis/, not an image-derived dataset.
 save(base/'FIGURE_PROVENANCE.json',dict(figures=[p.name for p in figs.glob('*.png')],state_source='analysis/<case>/states/*.npz full-precision accepted coordinates; geometry_events.npz exact native events; review_curves.npz all original scalar samples',deformation_scale=1,synthetic_surfaces=False,trajectory_interpolation=False,bond_display='All nodes and all relevant damage midpoints; for deformed scenes surviving bonds displayed only in the explicitly labelled central reference-y slab for legibility.',area_claim='No reconstructed surface, no certified physical new area, no Gc.'))
 print('saved',len(list(figs.glob('*.png'))),'PNGs')

if __name__=='__main__':main(Path(sys.argv[1]).resolve(),Path(sys.argv[2]).resolve())
