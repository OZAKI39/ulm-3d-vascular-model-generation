"""Figures use actual frozen inputs, independent audits and native accepted states."""
from common import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,LinearSegmentedColormap
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.lines import Line2D
import pandas as pd,time,resource
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.titlesize':14,'axes.labelsize':12,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'figure.facecolor':'white','savefig.facecolor':'white','axes.grid':True,'grid.color':'#788590','grid.alpha':.6,'pdf.fonttype':42})
DAMAGE=LinearSegmentedColormap.from_list('clot_damage',['#37779b','#d6ba78','#c54a54','#54234f'])
NAMES={'uniform':'U · Uniform reference','heterogeneous':'H · Correlated softening'}
COLORS={'uniform':'#3e7aa2','heterogeneous':'#a04e6a'}
graph=module('fig_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')

def save(fig,name):
 p=F/'figures';p.mkdir(exist_ok=True,parents=True);fig.savefig(p/(name+'.png'),dpi=250,bbox_inches='tight');fig.savefig(p/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
def ax3(fig,pos,title,limits):
 ax=fig.add_subplot(pos,projection='3d',computed_zorder=False);ax.view_init(elev=16,azim=-64);ax.set_title(title,pad=18);lo,hi=limits;pad=np.maximum((hi-lo)*.06,20.)
 for j,lim in enumerate([ax.set_xlim,ax.set_ylim,ax.set_zlim]):lim(lo[j]-pad[j],hi[j]+pad[j])
 ax.set_box_aspect(np.maximum(hi-lo,50));ax.set_xlabel('X (µm)',labelpad=8);ax.set_ylabel('Y (µm)',labelpad=8);ax.set_zlabel('Z (µm)',labelpad=8)
 for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=False;axis.set_major_locator(MaxNLocator(4));axis._axinfo['grid'].update(color='#7b8792',linewidth=.7)
 return ax

def damage(ref,d):
 b=ref['bonds'];sel=~ref['initial_cut'];count=np.bincount(b[sel].ravel(),minlength=len(ref['X']));return np.bincount(b[sel].ravel(),weights=np.repeat(d[sel],2),minlength=len(ref['X']))/np.maximum(count,1)
def component_colors(ref,s):
 lab,sz,att=graph.topology(ref,s['d']);n=int(ref['clot_count']);colors=np.tile(matplotlib.colors.to_rgba('#aebbc5'),(n,1))
 for label in np.flatnonzero(~att):
  ids=np.flatnonzero(lab==label);hue=((int(ids.min())+1)*.61803398875)%1;colors[ids]=np.r_[matplotlib.colors.hsv_to_rgb([hue,.65,.8]),1]
 return colors,lab,sz,att

def scatter(ax,ref,s,kind='damage',size=7):
 x=s['x']*1e6;n=int(ref['clot_count']);ax.scatter(*x[n:].T,c='#78818b',s=5,alpha=.25,depthshade=False)
 if kind=='components':ax.scatter(*x[:n].T,c=component_colors(ref,s)[0],s=size,depthshade=False)
 else:ax.scatter(*x[:n].T,c=damage(ref,s['d'])[:n],s=size,cmap=DAMAGE,vmin=0,vmax=1,depthshade=False,edgecolors='none')
def inputs_figures(ref,field):
 n=int(ref['clot_count']);X=ref['X']*1e6;lim=(X.min(0),X.max(0));bulk=ref['bond_types']==1;sc=field['bond_sc'][bulk]
 fig=plt.figure(figsize=(18,11));fig.subplots_adjust(wspace=.18,hspace=.25,bottom=.1)
 for j,(title,val,cmap,vmin,vmax) in enumerate([('A · Original material and anchors',None,None,None,None),('B · Uniform bulk sc = 0.18',np.full(n,.18),'viridis',sc.min(),sc.max()),('C · Frozen reference field q',field['q'][:n],'coolwarm',-1,1),('D · Node softening proxy',field['node_bulk_sc_proxy'][:n],'viridis',sc.min(),sc.max())]):
  ax=ax3(fig,231+j,title,lim);ax.scatter(*X[n:].T,c='#727d87',s=5,alpha=.3,depthshade=False)
  if val is None:ax.scatter(*X[:n].T,c='#b33e51',s=6,depthshade=False)
  else:
   im=ax.scatter(*X[:n].T,c=val,cmap=cmap,vmin=vmin,vmax=vmax,s=6,depthshade=False);pos=ax.get_position();cax=fig.add_axes([pos.x0+.2*pos.width,pos.y0-.005,pos.width*.6,.008]);fig.colorbar(im,cax=cax,orientation='horizontal')
 ax=ax3(fig,235,'E · Original bulk / interface bonds',lim);b=ref['bonds'];near=ref['r0']<=1.01*DX
 for mask,color in [(near&(ref['bond_types']==1),'#4f83a4'),(near&(ref['bond_types']==2),'#cf9a30')]:ax.add_collection3d(Line3DCollection(X[b[mask]],colors=color,linewidths=.4,alpha=.55))
 ax.text2D(.02,.96,'Display: nearest neighbors only; topology: all original bonds',transform=ax.transAxes,fontsize=9)
 ax=fig.add_subplot(236);ax.hist(sc,bins=32,color='#368b85');ax.axvline(.18,color='#293b4a',ls='--',label='Weighted mean = 0.18');ax.set_xlabel('Actual original-bond sc');ax.set_ylabel('Bulk bond count');ax.set_title(f'F · Actual bulk-bond sc\nRange {sc.min():.6f}–{sc.max():.6f}');ax.legend(loc='upper left');ax.set_ylim(0,6500)
 
 for axis in fig.axes:
  if hasattr(axis,'get_zlim'):
   pos=axis.get_position();axis.set_position([pos.x0,pos.y0+.05,pos.width,pos.height-.07])
 fig.suptitle('Reference specimen and spatially correlated softening · exploratory, uncalibrated',fontsize=18)
 fig.text(.06,.02,'D is a node display proxy; mechanics uses symmetric qij and actual pair-specific sc. Original IDs, bonds, weighted volume and precut retained.',fontsize=11);save(fig,'Fig1_Reference_and_heterogeneity')
 fig,axs=plt.subplots(2,3,figsize=(18,10),layout='constrained');q=field['q'][:n]
 axs[0,0].hist(q,bins=30,color='#8c749d');axs[0,0].set(xlabel='q',ylabel='Material node count',title='Bounded q distribution')
 axs[0,1].hist(sc-.15,bins=30,color='#368b85');axs[0,1].set(xlabel='Actual bond Δs',ylabel='Bulk bond count',title='Softening interval distribution')
 corr=pd.read_csv(O/'FIELD_CORRELATION.csv')
 for c,g in corr.groupby('axis'):axs[0,2].plot(g.distance_m*1e6,g.empirical_q_covariance_normalized,label='Finite q: '+c)
 d=np.linspace(0,750,200);axs[0,2].plot(d,np.exp(-(d/300)**2),'k--',label='Nominal Gaussian z');axs[0,2].axvline(300,color='#986b33',ls=':');axs[0,2].set(xlabel='Reference separation (µm)',ylabel='Normalized covariance',title='Finite field ≠ nominal covariance');axs[0,2].legend(fontsize=9)
 for j,c in enumerate(range(3)):
  unique=np.unique(ref['X'][:n,c]);val=unique[np.argmin(abs(unique))];sel=ref['X'][:n,c]==val;other=[k for k in range(3) if k!=c];im=axs[1,j].scatter(X[:n][sel,other[0]],X[:n][sel,other[1]],c=q[sel],cmap='coolwarm',vmin=-1,vmax=1,s=65,marker='s');axs[1,j].set_aspect('equal');axs[1,j].set(xlabel='XYZ'[other[0]]+' (µm)',ylabel='XYZ'[other[1]]+' (µm)',title=f"Reference {'XYZ'[c]} = {val*1e6:.0f} µm slice")
 fig.colorbar(im,ax=axs[1,:],shrink=.6,label='q');fig.suptitle('Scale audit · nominal correlation 300 µm · dx 50 µm · horizon 150.5 µm',fontsize=17);save(fig,'Fig2_Field_statistics_and_scale')
 path=pd.read_csv(O/'PRESCRIBED_PATH_RESPONSE.csv');audit=json.loads((O/'VARIATIONAL_ENERGY_AUDIT.json').read_text());fig,axs=plt.subplots(2,2,figsize=(14,10),layout='constrained')
 for name,g in path.groupby('case',sort=False):
  axs[0,0].plot(g.opening_m*1e6,g.reaction_N*1e6,label=name);axs[0,1].plot(g.opening_m*1e6,g.reaction_N*1e6,label=name);axs[1,0].plot(g.opening_m*1e6,g.D_J*1e9,label=name)
 axs[0,0].set(xlim=(0,60),xlabel='Prescribed opening (µm)',ylabel='Reaction (µN)',title='Full discrete opening response');axs[0,0].legend(fontsize=9)
 axs[0,1].set(xlim=(0,12),xlabel='Prescribed opening (µm)',ylabel='Reaction (µN)',title='Common elastic response and onset')
 axs[1,0].set(xlabel='Prescribed opening (µm)',ylabel='Path damage dissipation (nJ)',title='Offline prescribed-path energy · NOT Gc')
 rows=audit['checks']['variational'];axs[1,1].bar(['Intact','Near onset','Softened'],[r['relative_gradient_error'] for r in rows],color='#6a8c9f');axs[1,1].set_yscale('log');axs[1,1].set(ylabel='Relative force-gradient error',title='Native force vs independent ∂U/∂x | d');axs[1,1].axhline(1e-5,color='#ac4b57',ls='--',label='Check tolerance');axs[1,1].legend()
 fig.suptitle('Material-path verification · same complete endpoint-star energy · immutable pair properties',fontsize=17);save(fig,'Fig3_Material_energy_validation')
def results_figures(ref,field):
 n=int(ref['clot_count']);b=ref['bonds'];cut=ref['initial_cut'];near=ref['r0']<=1.01*DX
 data=json.loads((O/'COMPARISON.json').read_text());cases=data['cases'];states={m:{p.stem:load_npz(p) for p in (O/'cases'/m/'keys').glob('*.npz') if p.stem!='failed_unaccepted'} for m in CASES};hist={m:pd.read_csv(O/'cases'/m/'history.csv') for m in CASES}
 allx=np.concatenate([s['x'] for d in states.values() for s in d.values()])*1e6;global_limits=(allx.min(0),allx.max(0));X=ref['X']*1e6;ref_limits=(X.min(0),X.max(0))
 fig=plt.figure(figsize=(18,13));fig.subplots_adjust(hspace=.48,wspace=.15,bottom=.12)
 for row,m in enumerate(CASES):
  for col,key in enumerate(['first_softening','first_failure','final']):
   s=states[m].get(key);ax=ax3(fig,231+3*row+col,NAMES[m]+'\n'+key.replace('_',' '),ref_limits)
   ax.scatter(*X[:n].T,c=field['q'][:n],cmap='coolwarm',vmin=-1,vmax=1,s=3,alpha=.2,depthshade=False)
   if s is None:ax.text2D(.3,.5,'NOT_OBSERVED',transform=ax.transAxes);continue
   # Reference positions reveal event location without hiding final deformation.
   event=cases[m]['native']['events'].get(key,{})
   fail=(s['d']==1)&~cut;shown=near&~cut
   if np.any(shown&(s['d']<1)):ax.add_collection3d(Line3DCollection(X[b[shown&(s['d']<1)]],colors='#55947e',linewidths=.3,alpha=.12))
   if np.any(fail&near):ax.add_collection3d(Line3DCollection(X[b[fail&near]],colors='#a52e45',linewidths=.65,alpha=.65,linestyles='dashed'))
   ids=np.array(event.get('original_bond_ids',[]),int)-1
   if len(ids):
    ax.add_collection3d(Line3DCollection(X[b[ids]],colors='#e68a13',linewidths=2.8));nodes=np.unique(b[ids]);ax.scatter(*X[nodes].T,color='#e68a13',s=25,depthshade=False)
   ax.text2D(.02,.94,f't = {float(s["time"])*1e3:.4f} ms',transform=ax.transAxes)
 fig.suptitle('True original-bond events mapped to reference q · reference coordinates',fontsize=17)
 fig.text(.05,.035,'Green: displayed live nearest bonds; dashed red: naturally failed nearest bonds; orange: actual event bonds (all lengths). Connectivity always uses all bonds.',fontsize=10);save(fig,'Fig4_Original_bond_fracture_locations')
 for m in CASES:
  late=max((v for k,v in states[m].items() if k.startswith('frame_')),key=lambda s:float(s['time'])<=.002 and float(s['time']) or -1)
  keys=['initial','first_softening','first_failure',None,'max_multifragment','final'];labels=['Initial','First softening','First complete failure','Crack propagation','Maximum 20–200 count','Final accepted']
  fig=plt.figure(figsize=(18,11));fig.subplots_adjust(wspace=.17,hspace=.28,bottom=.12)
  for j,(key,label) in enumerate(zip(keys,labels)):
   s=late if key is None else states[m].get(key)
   if key=='max_multifragment' and cases[m]['native']['max_medium_components']==0:s=None
   ax=ax3(fig,231+j,label+(' · NOT_OBSERVED' if s is None else f'\n{float(s["time"])*1e3:.4f} ms'),global_limits)
   if s is not None:scatter(ax,ref,s,size=5)
  fig.suptitle(NAMES[m]+' · native accepted states · displacement ×1 · common world bounds',fontsize=17);cax=fig.add_axes([.35,.045,.3,.018]);fig.colorbar(ScalarMappable(norm=Normalize(0,1),cmap=DAMAGE),cax=cax,orientation='horizontal',label='Mean incident damage (original precut excluded)');save(fig,'Fig5_Damage_keyframes_'+m)
 fig=plt.figure(figsize=(16,12));fig.subplots_adjust(hspace=.35,wspace=.17)
 dist=pd.read_csv(O/'FRAGMENT_SIZE_DISTRIBUTION.csv')
 for j,m in enumerate(CASES):
  s=states[m]['final'];ax=ax3(fig,221+j,NAMES[m]+' · true final components',global_limits);scatter(ax,ref,s,'components',size=9)
  ax=fig.add_subplot(223+j);r=dist[dist['case']==m];bars=ax.bar(r['bin'],r['count'],color='#4d8d9c');ax.set(title=NAMES[m]+' · free sizes',xlabel='Original nodes per component',ylabel='Component count');ax.yaxis.set_major_locator(MaxNLocator(integer=True));ax.bar_label(bars,fmt='%d');ax.set_ylim(0,max(1,r['count'].max()*1.2));ax2=ax.twinx();ax2.plot(r['bin'],100*r.fraction_initial,'o--',color='#a74b62');ax2.set_ylabel('Mass / initial clot (%)',color='#a74b62');ax2.grid(False)
 fig.suptitle('True surviving-graph fragments · gray attached body / anchors · independent free colors',fontsize=17);save(fig,'Fig6_True_fragments_and_sizes')
 # Separate actual-coordinate local view: no fictitious shell or displacement scale.
 fig=plt.figure(figsize=(14,7))
 for j,m in enumerate(CASES):
  s=states[m]['final'];_,lab,sz,att=component_colors(ref,s);labels=np.flatnonzero(~att);ids=np.flatnonzero(lab==labels[np.argmax(sz[labels])]) if len(labels) else np.array([],int)
  center=s['x'][ids].mean(0)*1e6 if len(ids) else np.zeros(3);ax=ax3(fig,121+j,NAMES[m]+f'\nLargest free component: {len(ids)} nodes',(center-150,center+150))
  ax.set_xlim(center[0]-150,center[0]+150);ax.set_ylim(center[1]-150,center[1]+150);ax.set_zlim(center[2]-150,center[2]+150)
  if len(ids):
   ax.scatter(*s['x'][ids].T*1e6,s=90,c='#b34c60',edgecolors='#273b47',depthshade=False);sel=np.isin(b,ids).all(1)&(s['d']<1);ax.add_collection3d(Line3DCollection(s['x'][b[sel]]*1e6,colors='#325d7b',linewidths=2))
   ax.text2D(.02,.94,f'{sel.sum()} live internal bonds',transform=ax.transAxes)
   for i in ids:ax.text(*(s['x'][i]*1e6+[3,3,4]),f'ID {i+1}',fontsize=10)
 fig.suptitle('Independent local detail · identical 300 µm view width · actual world coordinates ×1',fontsize=16);save(fig,'Fig6b_Largest_free_component_detail')
 for m in CASES:
  h=hist[m];t=h.time_s*1e3;fig,axs=plt.subplots(2,3,figsize=(18,10),layout='constrained')
  axs[0,0].plot(t,h.shear_peak_Pa/1000,label='Actual peak facet shear');axs[0,0].plot(t,h.nominal_amplitude_Pa/1000,ls='--',label='Frozen amplitude before envelope');axs[0,0].set(ylabel='Shear / amplitude (kPa)',title='Actual response and fixed drive');axs[0,0].legend(fontsize=9)
  axs[0,1].plot(t,h.bulk_failed,label='Bulk failed');axs[0,1].plot(t,h.interface_failed,label='Interface failed');axs[0,1].set(ylabel='Natural failure count',title='Original bond failures');axs[0,1].legend()
  for k,l in [('detached_nodes','Detached nodes'),('largest_nodes','Largest free'),('N20_200','20–200 components')]:axs[0,2].plot(t,h[k],label=l)
  axs[0,2].set(ylabel='Count',title='True material topology');axs[0,2].legend()
  for k,l in [('W_J','Net W'),('W_drive_J','Drive'),('W_drag_J','Passive drag')]:axs[1,0].plot(t,h[k]*1e9,label=l)
  axs[1,0].set(ylabel='Actual work (nJ)',title='Signed work channels');axs[1,0].legend()
  for k,l in [('U_J','Stored U'),('K_J','Kinetic K'),('D_J','Damage D'),('Uc_J','Contact Uc'),('Jc_J','Contact creation Jc'),('R_J','Numerical R')]:axs[1,1].plot(t,h[k]*1e9,label=l)
  axs[1,1].set(ylabel='Energy (nJ)',title='Full energy ledger');axs[1,1].legend(fontsize=9,ncol=2)
  axs[1,2].plot(t,h.relative_R*100);axs[1,2].axhline(1,color='#a74b62',ls='--');axs[1,2].set(ylabel='|R| / historical energy scale (%)',title='Numerical admission ≤1%')
  for ax in axs.ravel():ax.set_xlabel('Time (ms)');ax.axvline(OFF*1e3,color='#808080',ls=':',lw=1)
  fig.suptitle(NAMES[m]+' · accepted loading, damage and complete energy ledger',fontsize=17);save(fig,'Fig7_Loading_damage_energy_'+m)
 fig,axs=plt.subplots(3,3,figsize=(16,12),layout='constrained')
 metrics=[('Target count 20–200',lambda q:q['objective']['n_fragments_20_200']),('Target mass / initial (%)',lambda q:100*q['objective']['mass_fraction_target_initial']),('Largest free (nodes)',lambda q:max(q['final_free_sizes'],default=0)),('Fines / detached (%)',lambda q:100*(q['objective']['fine_debris_mass_fraction_detached'] or 0)),('Damage dissipation (nJ)',lambda q:q['native']['last_row']['D_J']*1e9),('Actual net work (nJ)',lambda q:q['work']['actual_net_J']*1e9),('Naturally failed bulk bonds',lambda q:q['native']['last_row']['bulk_failed']),('Damage midpoint spread (µm)',None),('Native wall time (s)',lambda q:q['native']['wall_s'])]
 for i,(ax,(title,fn)) in enumerate(zip(axs.ravel(),metrics)):
  vals=[fn(cases[m]) if fn else data['spatial'][m]['failed_bond_midpoint_rms_spread_m']*1e6 for m in CASES];bars=ax.bar(['Uniform','Heterogeneous'],vals,color=[COLORS[m] for m in CASES]);ax.bar_label(bars,fmt='%.4g');ax.set_title(title);ax.set_ylim(0,max(1,max(vals)*1.2));ax.grid(axis='x',visible=False)
  if i in [0,2,6]:ax.yaxis.set_major_locator(MaxNLocator(integer=True))
 fig.suptitle('Same frozen shear drive · only bulk softening interval differs · one spatial realization',fontsize=17);save(fig,'Fig8_Core_scientific_comparison')
 # Supporting association figure, with all eligible bulk bonds in denominators.
 fig,axs=plt.subplots(1,2,figsize=(14,5),layout='constrained')
 for m in CASES:
  r=data['spatial'][m]['q_quartile_damage'];axs[0].plot(range(1,5),[100*x['failed_fraction'] for x in r],'o-',label=NAMES[m]);axs[1].plot(range(1,5),[x['mean_final_damage'] for x in r],'o-',label=NAMES[m])
 axs[0].set(ylabel='Failed bulk bonds / eligible (%)',title='Local q quartile vs complete failure');axs[1].set(ylabel='Mean bulk damage d',title='Local q quartile vs damage')
 for ax in axs:ax.set_xlabel('Reference q quartile (low → high)');ax.set_xticks(range(1,5));ax.legend()
 fig.suptitle('Descriptive spatial associations · single field · no independent-sample or causal claim',fontsize=15);save(fig,'Fig8b_Local_property_damage_association')
 dump(O/'VISUALIZATION_GEOMETRY.json',dict(common_global_bounds_um=np.array(global_limits),camera=dict(elevation_deg=16,azimuth_deg=-64),displacement_scale=1,connectivity='all original live bonds',display_bond_cutoff_m=1.01*DX,node_damage='mean incident d, original precut excluded',local_detail_width_um=300))
def main():
 start=time.perf_counter();cpu=time.process_time();ref=load_npz(O/'provenance/reference.npz');field=load_npz(O/'HETEROGENEITY_FIELD.npz');inputs_figures(ref,field);results_figures(ref,field)
 dump(O/'FIGURE_RESOURCES.json',dict(wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024));print('FIGURES_DONE',flush=True)
if __name__=='__main__':main()
