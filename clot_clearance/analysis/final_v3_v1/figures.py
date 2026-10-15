from common import *
import numpy as np,json,time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from pd_clot.lammps_fragmentation.early_flow_coupling import DTYPE
from pd_clot.lammps_fragmentation.detached_fragment_metrics import graph_components
start=time.perf_counter();F.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':13,'axes.titlesize':15,'axes.labelsize':13,'legend.fontsize':12,'svg.fonttype':'none','axes.grid':True,'grid.color':'#8393a2','grid.alpha':.7,'figure.facecolor':'white','axes.spines.top':False,'axes.spines.right':False})
BLUE='#3579a6';ORANGE='#cf633e';GREEN='#399777';PURPLE='#8c58a3';GRAY='#a4aeb8'
ref=np.load(P/'reference.npz');X=ref['X'];b=ref['bonds'];typ=ref['bond_types'];m=ref['reference_m'];n=1927;N=len(X);V=50e-6**3
s=np.genfromtxt(P/'accepted_step_scalars.csv',delimiter=',',names=True);ht=np.genfromtxt(O/'HISTORICAL_TOPOLOGY_EVERY_STEP.csv',delimiter=',',names=True);comp=np.genfromtxt(O/'V3_OLD_NEW_COMPARISON.csv',delimiter=',',names=True);q=json.loads((O/'SCIENTIFIC_SUMMARY.json').read_text());mi=json.loads((O/'MPI_SORT_AND_MIGRATION_EVENTS.json').read_text())
def save(fig,name,svg=True):
 fig.savefig(F/(name+'.png'),dpi=210,bbox_inches='tight')
 if svg:fig.savefig(F/(name+'.svg'),bbox_inches='tight')
 plt.close(fig)
def setup3(ax,limits=None):
 ax.view_init(elev=17,azim=-55);ax.set(xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)');ax.tick_params(labelsize=10,pad=1)
 if limits is not None:
  ax.set_xlim(limits[0]);ax.set_ylim(limits[1]);ax.set_zlim(limits[2]);ax.set_box_aspect(np.maximum(np.diff(limits,axis=1).ravel(),.25))
 for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.set_major_locator(MaxNLocator(4));axis._axinfo['grid'].update(color='#8c9ba8',linewidth=.65);axis.pane.set_facecolor('#f3f5f6');axis.pane.set_alpha(.55)
def scatter(ax,x,c,s=13,**kwargs):return ax.scatter(*(x*1e3).T,c=c,s=s,depthshade=False,edgecolors='none',rasterized=True,**kwargs)
def limits_for(xs):
 xx=np.concatenate(xs)*1e3;lo=xx.min(0);hi=xx.max(0);pad=np.maximum((hi-lo)*.055,.02);return np.c_[lo-pad,hi+pad]
def clot(ax,x,alive,title,limits,anchors=True):
 lab,sz,at=graph_components(n,b,typ,alive);det=~at[lab];tiny=det&(sz[lab]<=2);resolved=det&~tiny
 eb=b[alive&(typ==1)];ax.add_collection3d(Line3DCollection(x[eb]*1e3,colors=BLUE,linewidths=.17,alpha=.035,rasterized=True))
 if anchors:scatter(ax,x[n:],GRAY,s=2.0,alpha=.16)
 scatter(ax,x[:n][~det],BLUE,s=11);scatter(ax,x[:n][tiny],ORANGE,s=17);scatter(ax,x[:n][resolved],PURPLE,s=23)
 if resolved.any():
  eb=b[alive&(typ==1)&det[b[:,0].clip(max=n-1)]&det[b[:,1].clip(max=n-1)]];ax.add_collection3d(Line3DCollection(x[eb]*1e3,colors=PURPLE,linewidths=1.8,alpha=.9))
 setup3(ax,limits);ax.set_title(title,pad=8);return sz,at
legend=[Line2D([0],[0],marker='o',ls='',color=c,label=l,markersize=7) for c,l in [(BLUE,'Attached material'),(GRAY,'Fixed wall anchors'),(ORANGE,'Detached 1–2 node debris'),(PURPLE,'Detached components ≥3 nodes')]]
# 1: exact reference / owner / sorter data.
z=np.load(P/'states/step_00455.npz');fig=plt.figure(figsize=(19,11));lim=limits_for([X[:n]])
ax=fig.add_subplot(231,projection='3d');a=scatter(ax,X[:n],m[:n]/1e-15,s=14,cmap='viridis');setup3(ax,lim);ax.set_title('Original V3: fixed reference m');fig.colorbar(a,ax=ax,shrink=.72,pad=.17,location='bottom',label='m (10$^{-15}$ m$^4$)')
ax=fig.add_subplot(232,projection='3d');scatter(ax,z['position_m'][:n],np.where(z['owner'][:n]==0,BLUE,ORANGE),s=14);setup3(ax,lim);ax.set_title('Actual owners at step 455');ax.legend(handles=[Line2D([0],[0],marker='o',ls='',color=BLUE,label='Rank 0'),Line2D([0],[0],marker='o',ls='',color=ORANGE,label='Rank 1')],loc='upper left')
ax=fig.add_subplot(233,projection='3d');allmig=np.array(mi['ownership_unique_atom_ids'])-1;allmig=allmig[allmig<n];evs=[json.loads(v) for r in [0,1] for v in (P/f'identity_events_rank{r}.jsonl').read_text().splitlines()];sortids=np.unique([i for e in evs if e['sort_called'] for i in e['sorted_ids']])-1;sortids=sortids[sortids<n];scatter(ax,X[:n],'#c4d1d9',s=10,alpha=.4);scatter(ax,X[sortids],GREEN,s=13,alpha=.6);scatter(ax,X[allmig],PURPLE,s=26);setup3(ax,lim);ax.set_title('IDs exercised over all 4660 steps');ax.legend(handles=[Line2D([0],[0],marker='o',ls='',color=GREEN,label='Actually permuted'),Line2D([0],[0],marker='o',ls='',color=PURPLE,label='Changed owner')],loc='upper left')
ax=fig.add_subplot(234);ax.scatter(np.arange(1,N+1),abs(z['m']-m)/m,s=4,color=BLUE);ax.set(xlabel='Immutable atom ID',ylabel='|m native − m reference| / m reference',title='Per-ID error at first migration');ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
ax=fig.add_subplot(235);ax.plot(s['time_s']*1e3,s['m_max_relative'],color=GREEN,lw=2,label='Maximum owned + ghost error');ax.axhline(2e-12,c='black',ls='--',label='Frozen gate: 2×10$^{-12}$');ax.set_yscale('log');ax.set_ylim(3e-16,5e-12)
for eventstep in sorted({e['step'] for e in mi['actual_sort_events']}):ax.axvline(eventstep*.0005,c=BLUE,alpha=.4,lw=1)
ax.set(xlabel='Actual physical time (ms)',ylabel='Maximum per-ID relative m error',title='Full run: 4661 audited states');ax.legend(loc='upper right')
ax=fig.add_subplot(236);ax.step(s['time_s']*1e3,s['ownership_changes'],where='post',color=PURPLE,lw=1.6,label='IDs changing rank');ay=ax.twinx();ay.step(s['time_s']*1e3,s['sort_nonidentity_nodes'],where='post',color=GREEN,alpha=.65,label='Nonidentity sort: nodes');ay.set_ylabel('Actually sorted nodes',color=GREEN);ay.grid(False);ax.set(xlabel='Actual physical time (ms)',ylabel='IDs changing owner per step',title='Migration and sorting are distinct');ax.set_ylim(0,3);ax.set_yticks([0,1,2,3]);ax.legend(loc='upper left');ay.legend(loc='upper right')
fig.suptitle('1 | Final two-fix library: native identity on the real V3 clot',fontsize=20);fig.tight_layout(rect=[0,0,1,.95]);save(fig,'Fig1_native_identity')
# 2: actual exact-time old/full-precision states. Historical sparse force file reused read-only.
oldraw=np.memmap(OLD/'saved_frame_force_audit.bin',dtype=DTYPE,mode='r').reshape(-1,n);oldlookup={int(st):i for i,st in enumerate(oldraw['step'][:,0])};olde=np.load(OLD/'broken_bonds.npz');oes=olde['events'][:,1].astype(int);oei=olde['initial_bond_indices'];steps=[0,241,455,593,2000,4660];states=[]
for st in steps:
 z=np.load(P/'states'/f'step_{st:05d}.npz');ox=X.copy();ox[:n]=oldraw[oldlookup[st]]['x_m'];oa=np.ones(len(b),bool);oa[oei[oes<=st]]=False;states.append((st,ox,oa,z['position_m'],z['alive']))
lim=limits_for([x[:n] for _,ox,_,nx,_ in states for x in [ox,nx]])
fig=plt.figure(figsize=(21,22))
for k,(st,ox,oa,nx,na) in enumerate(states):
 group=k//3;col=k%3
 for row,(name,x,al) in enumerate([('Historical',ox,oa),('Final two-fix',nx,na)]):
  ax=fig.add_subplot(4,3,(2*group+row)*3+col+1,projection='3d');sz,at=clot(ax,x,al,f'{name} | step {st} · {st*.0005:.4f} ms',lim);ax.text2D(.01,.89,f'{int(sz[~at].sum())} detached nodes · {int((~at).sum())} components',transform=ax.transAxes,fontsize=11)
fig.suptitle('2 | Same accepted step, same IDs, same camera and physical scale',fontsize=23,y=.995);fig.legend(handles=legend,loc='lower center',ncol=4,bbox_to_anchor=(.5,.005),fontsize=14);fig.subplots_adjust(left=.01,right=.98,bottom=.055,top=.955,hspace=.17,wspace=.07);save(fig,'Fig2_old_new_evolution',False)
# 3: audit curves, separate force coverage.
fig,axs=plt.subplots(2,3,figsize=(19,10.5));t=s['time_s']*1e3;htime=ht['time_s']*1e3
panels=[('broken_bulk','Bulk broken bonds'),('broken_interface','Interface broken bonds'),('detached_nodes','Detached mass (% of initial clot)'),('detached_components','Detached connected components'),('largest_detached_nodes','Largest detached equivalent diameter (µm)')]
for ax,(key,label) in zip(axs.ravel(),panels):
 yy=s[key];oy=ht[key]
 if key=='detached_nodes':yy=yy/n*100;oy=oy/n*100
 if key=='largest_detached_nodes':yy=(6*yy*(50e-6)**3/np.pi)**(1/3)*1e6;oy=(6*oy*(50e-6)**3/np.pi)**(1/3)*1e6
 ax.plot(htime,oy,c=BLUE,lw=2.5,label='Historical');ax.plot(t,yy,c=ORANGE,ls='--',lw=1.8,label='Final two-fix');ax.axvline(.2275,c='gray',ls=':',lw=1.4);ax.set(xlabel='Time (ms)',ylabel=label)
 if key=='largest_detached_nodes':ax.axhline(200,color='black',ls=':',lw=1);ax.set_ylim(0,215)
axs[0,0].legend();ax=axs[1,2];ax.plot(comp['time_s']*1e3,comp['old_vmax_saved_frame_m_s'],c=BLUE,label='Old: every 20-step saved frames');ax.plot(t,s['vmax_m_s'],c=ORANGE,lw=1.6,label='Final: every accepted step');ax.axvline(.2275,c='gray',ls=':',lw=1.4);ax.set(xlabel='Time (ms)',ylabel='Maximum material-point speed (m/s)');ax.legend(fontsize=10)
fig.suptitle('3 | Actual fracture, free topology and dynamics · dotted line = first saved force difference',fontsize=18);fig.tight_layout(rect=[0,0,1,.95]);save(fig,'Fig3_fracture_dynamics')
# 4: final real graph and auditable tiny/resolved counts.
z=np.load(P/'states/step_04660.npz');x=z['position_m'];al=z['alive'];lab=z['component_label'];sz=z['component_sizes'];at=z['component_attached'];nf=q['new_final'];of=q['old_final']
fig=plt.figure(figsize=(20,11));ax=fig.add_subplot(231,projection='3d');clot(ax,x,al,'Full actual V3 at 2.33 ms',limits_for([x]),anchors=True)
ax=fig.add_subplot(232,projection='3d');det=~at[lab];tiny=det&(sz[lab]<=2);res=det&~tiny;scatter(ax,x[:n][tiny],ORANGE,s=22);scatter(ax,x[:n][res],PURPLE,s=32);eb=b[al&(typ==1)&det[b[:,0].clip(max=n-1)]&det[b[:,1].clip(max=n-1)]];ax.add_collection3d(Line3DCollection(x[eb]*1e3,colors=PURPLE,linewidths=2.6));setup3(ax,limits_for([x]));ax.set_title('Detached points: same coordinates and scale')
ax=fig.add_subplot(233);countvals=sorted(set(map(int,nf['detached_component_size_histogram']))|set(map(int,of['detached_component_size_histogram'])));pos=np.arange(len(countvals));ax.bar(pos-.18,[of['detached_component_size_histogram'].get(str(v),0) for v in countvals],.36,color=BLUE,label='Historical');ax.bar(pos+.18,[nf['detached_component_size_histogram'].get(str(v),0) for v in countvals],.36,color=ORANGE,label='Final two-fix');ax.set_xticks(pos,countvals);ax.set(xlabel='Nodes in a detached component',ylabel='Component count',title='True detached size distribution');ax.legend()
ax=fig.add_subplot(234);ax.plot(t,s['tiny_detached_nodes']/n*100,c=ORANGE,label='1–2 node detached debris');ax.plot(t,s['detached_nodes']/n*100,c=PURPLE,ls='--',label='All detached');ax.set(xlabel='Time (ms)',ylabel='% initial material mass',title='Unresolved debris remains dominant');ax.legend()
ax=fig.add_subplot(235);ax.plot(t,(s['attached_nodes']+s['detached_nodes']),c=GREEN,lw=2,label='Attached + detached');ax.axhline(n,c='black',ls='--',label='1927 initial material points');ax.set_ylim(1926,1928);ax.set_yticks([1926,1927,1928]);ax.set(xlabel='Time (ms)',ylabel='Material node count',title='Conserved independently of wall anchors');ax.legend()
ax=fig.add_subplot(236);ax.axis('off');ax.text(0,.98,f"FINAL PATCHED DIAGNOSTIC LIBRARY\n\n4660 accepted steps · both ranks exit 0\n{nf['detached_nodes']} detached nodes / {nf['detached_component_count']} components\n{nf['tiny_detached_nodes']} nodes in 1–2 node debris\nDebris = {100*nf['tiny_fraction_of_detached_mass']:.2f}% of detached mass\nLargest final free component: {nf['largest_detached_nodes']} nodes\nEquivalent diameter: {nf['largest_detached_diameter_um']:.2f} µm\n≥200 µm: {nf['N_detached_ge200']} · ≥300 µm: {nf['N_detached_ge300']}\n\nNo experimental fracture validation.\nExplicit theta time layer retained.\nProduction library unchanged.",va='top',fontsize=14,linespacing=1.55)
fig.suptitle('4 | Final attachment-aware fragment review · no shells or merged particles',fontsize=20);fig.legend(handles=legend,loc='lower center',ncol=4);fig.tight_layout(rect=[0,.055,1,.95]);save(fig,'Fig4_final_fragments')
dump(O/'FIGURE_PROVENANCE.json',{'figures':[p.name for p in F.glob('*.png')],'physical_data':'C1_final reference, true accepted states, alive initial bond IDs and components; old saved force records reused by exact step','camera':'elev17 azim−55, fixed within each comparison','illustrated_steps':steps,'no_interpolation':True,'no_synthetic_shell':True,'runtime_s':time.perf_counter()-start})
print('FIGURES_DONE',time.perf_counter()-start)
