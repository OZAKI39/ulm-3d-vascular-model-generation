from common import *
from graph_ops import topology
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.lines import Line2D
import html,time
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':15,'axes.labelsize':14,'xtick.labelsize':12,'ytick.labelsize':12,'legend.fontsize':12,'figure.facecolor':'white','axes.grid':True,'grid.color':'#8e959c','grid.alpha':.65,'grid.linewidth':.7})
BLUE='#1676a1';RED='#ce3550';GOLD='#df971a';GRAY='#9aa4ad';PURPLE='#7258a5'
CAPTIONS={}
def save(fig,name,caption):
 fig.savefig(F/(name+'.png'),dpi=240,bbox_inches='tight');plt.close(fig);CAPTIONS[name]=caption;print('FIG',name,flush=True)
def limits(arrays,origin):
 points=np.concatenate(arrays)-origin;lo=points.min(0)*1e3;hi=points.max(0)*1e3;center=(lo+hi)/2;span=(hi-lo).max()*1.08;return center,span
def axes3(ax,bounds):
 c,span=bounds
 for j,a in enumerate('xyz'):getattr(ax,'set_'+a+'lim')(c[j]-span/2,c[j]+span/2);getattr(ax,'set_'+a+'label')(a.upper()+' [mm]',labelpad=9)
 ax.set_box_aspect((1,1,1));ax.view_init(elev=18,azim=-58);ax.tick_params(pad=1,labelsize=10)
 for a in [ax.xaxis,ax.yaxis,ax.zaxis]:a.pane.set_facecolor((.97,.98,.99,1));a._axinfo['grid']['color']=(.47,.51,.55,.7)
def draw(ax,ref,s,origin,bounds,candidate=None,local=None):
 n=int(ref['clot_count']);x=(s['x']-origin)*1e3;labels,sizes,att=topology(ref,s['d']);free=~att[labels];mask=np.ones(len(x),bool) if local is None else local
 edges=ref['bonds'][(ref['bond_types']==1)&(s['d']<1)];edges=edges[mask[edges].all(1)]
 ax.add_collection3d(Line3DCollection(x[edges],colors=BLUE,linewidths=.2,alpha=.12,rasterized=True))
 for ids,col,size in [(np.flatnonzero(mask[n:])+n,GRAY,3),(np.flatnonzero(mask[:n]&~free),BLUE,7),(np.flatnonzero(mask[:n]&free),RED,16)]:
  if len(ids):ax.scatter(*x[ids].T,s=size,c=col,depthshade=False,edgecolors='none')
 if candidate is not None:ax.scatter(*x[candidate].T,s=65,c=GOLD,edgecolors='black',linewidths=.6,depthshade=False)
 axes3(ax,bounds)
def main():
 start=time.perf_counter();F.mkdir(parents=True,exist_ok=True);ref=reference();n=int(ref['clot_count']);states={int(s['step']):s for _,s in keyframes()};origin=ref['X'][:n].mean(0);bound=limits([s['x'] for s in states.values()],origin);cand=json.loads((O/'V3_CHUNK_CANDIDATES.json').read_text());ids=np.asarray(cand['candidates'][0]['atom_ids'])-1
 fig=plt.figure(figsize=(18,6.4));handles=[Line2D([],[],marker='o',linestyle='',color=c,label=l) for c,l in [(BLUE,'Attached material'),(RED,'ACTUAL DETACHED'),(GOLD,'TETHERED CANDIDATE'),(GRAY,'Fixed anchor')]]
 for i,step in enumerate([1588,5182,10502]):
  ax=fig.add_subplot(1,3,i+1,projection='3d');s=states[step];draw(ax,ref,s,origin,bound,candidate=ids if step==5182 else None);ax.set_title(f"{s['time']*1e3:.3f} ms | actual points + live bulk bonds")
 fig.legend(handles=handles,loc='lower center',ncol=4,bbox_to_anchor=(.5,-.02));fig.suptitle('V3 accepted states — displacement ×1; common physical scale / camera',fontsize=18)
 save(fig,'Fig01_V3_actual_states','ACTUAL DETACHED and TETHERED CANDIDATE are graph labels. All surviving bulk bonds shown; scaffold/interface lines omitted for readability. Axes use mm relative to the same reference centroid. No interpolated fracture surface.')
 center=states[5182]['x'][ids].mean(0);rad=.22e-3;local=np.linalg.norm(states[5182]['x']-center,axis=1)<rad;local[ids]=True;bridge=90571-1;local[ref['bonds'][bridge]]=True;bb=((center-origin)*1e3,rad*2e3)
 fig=plt.figure(figsize=(13,6))
 for i,step in enumerate([5182,10502]):
  s=states[step];ax=fig.add_subplot(1,2,i+1,projection='3d');draw(ax,ref,s,origin,bb,candidate=ids,local=local)
  if s['d'][bridge]<1:
   xx=(s['x'][ref['bonds'][bridge]]-origin)*1e3;ax.plot(*xx.T,c='black',lw=2.5);ax.text(*xx.mean(0),'bond 90571',fontsize=12)
  internal=np.isin(ref['bonds'][:,0],ids)&np.isin(ref['bonds'][:,1],ids)&(s['d']<1)
  ax.add_collection3d(Line3DCollection((s['x'][ref['bonds'][internal]]-origin)*1e3,colors=GOLD,linewidths=2,alpha=1))
  title='TETHERED CANDIDATE (3 nodes)' if step==5182 else 'ACTUAL DETACHED — same IDs, now split'
  ax.set_title(f"{title}\n{s['time']*1e3:.3f} ms; IDs 1625 / 1632 / 1649")
 save(fig,'Fig02_candidate_bridge','Only the three-node candidate passes the frozen graph screen at 1 ms. Orange keeps the same identities at final time; it is not an assertion that the candidate criteria still pass. Black bridge is drawn only while actually alive. Large candidates ≥34/114 nodes were not found. Candidate local shear is NOT MEASURED because the star is rank deficient.')
 h=np.genfromtxt(B/'V3_NATIVE_PILOT_HISTORY.csv',delimiter=',',names=True);t=h['time_s']*1e3;c=json.loads((O/'COHESIVE_ADHESIVE_COMPETITION.json').read_text());tr=c['historical_maximum_fragment']['trajectory']
 fig,ax=plt.subplots(1,3,figsize=(17,4.6));ax[0].plot(t,h['bulk_failed']/82270*100,label='Bulk',color=BLUE);ax[0].plot(t,h['interface_failed']/22040*100,label='Interface',color=PURPLE);ax[0].set(xlabel='Time [ms]',ylabel='Failed fraction [%]',title='Actual typed failure histories');ax[0].legend();ax[1].plot(t,h['detached_nodes'],label='Detached nodes',color=RED);ax[1].plot(t,h['largest_detached_nodes'],label='Largest component',color=BLUE);ax[1].set(xlabel='Time [ms]',ylabel='Nodes',title='ACTUAL DETACHED');ax[1].legend();sel=[v for v in tr if v['step']>=3281];ax[2].step([v['time_s']*1e3 for v in sel],[max(v['pieces']) for v in sel],where='post',label='Largest daughter',color=RED);ax[2].step([v['time_s']*1e3 for v in sel],[v['internal_live'] for v in sel],where='post',label='Surviving internal bonds',color=BLUE);ax[2].set(xlabel='Time [ms]',ylabel='Count',title='Exact 13-node identity set');ax[2].legend();fig.tight_layout();save(fig,'Fig03_competition_and_breakup','Every accepted state is checked against original-bond connectivity. The 13-node component forms at 0.554125 ms with only 13/42 original internal bonds remaining, then splits at 0.5745 ms. Step plots represent event changes, not interpolated node states.')
 ct=json.loads((O/'CONTACT_0376MS_AUDIT.json').read_text());peak=load_npz(O/'derived/recovered_accepted/before_002383.npz');near=np.asarray(ct['reference_neighborhood_atom_ids'])-1;center=peak['x'][[573,1428]].mean(0);local=np.zeros(len(ref['X']),bool);local[near]=True
 fig=plt.figure(figsize=(14,5.5));ax=fig.add_subplot(121,projection='3d');draw(ax,ref,peak,origin,((center-origin)*1e3,.6),local=local);pair=np.array([573,1428]);xx=(peak['x'][pair]-origin)*1e3;ax.scatter(*xx.T,c=[GOLD,RED],s=80,edgecolor='black',depthshade=False);ff=np.array(ct['states'][1]['pair_force_on_first_N']);ax.quiver(*xx[0],*(ff*1e4),color='black',linewidth=2,arrow_length_ratio=.2);ax.quiver(*xx[1],*(-ff*1e4),color='black',linewidth=2,arrow_length_ratio=.2);ax.set_title('0.376 ms accepted state\nPair repulsion: 12.64 µN (arrow scale 10 µm/µN)')
 ax=fig.add_subplot(122);win=ct['local_bulk_failure_windows'];ax.bar(np.arange(len(win)),[v['local_bulk_new_failures'] for v in win],color=[BLUE,BLUE,GOLD,RED,RED]);ax.set_xticks(np.arange(len(win)),['0–.2305','.2305–.37575','.37575–.37625','.37625–.45','.45–2.33'],rotation=25);ax.set(xlabel='Actual failure time interval [ms]',ylabel='Local bulk failure events',title='DIAGNOSTIC ONLY: time ordering, not causality');fig.tight_layout();save(fig,'Fig04_contact_peak_audit','Positions and forces use SHA-exact recovered accepted state 1943. Arrows show computed original nonbonded pair repulsion only, with stated force scale. Neighborhood is fixed in reference coordinates within delta of either contact atom; widths of histogram time intervals differ. Damage predates peak contact.')
 su=c['surface_horizon_audit'];bins=su['distance_bins'];fig,ax=plt.subplots(1,2,figsize=(12,4.8));xx=np.arange(4);ax[0].bar(xx-.17,[v['early_failure_node_fraction'] or 0 for v in bins],width=.34,color=RED,label='Early failed-bond endpoints');ax[0].bar(xx+.17,[v['load_area_node_fraction'] or 0 for v in bins],width=.34,color=BLUE,label='Applied-load nodes');ax[0].set_xticks(xx,['0–0.5','0.5–1','1–2','>2']);ax[0].set(xlabel='Distance to original exposed surface / δ',ylabel='Fraction within distance bin',title='DIAGNOSTIC ONLY: confounded surface association');ax[0].legend(fontsize=10);groups=['all_material','early_bulk_failed_nodes','ever_bulk_failed_nodes'];ax[1].bar(np.arange(3),[su[k]['m_ratio_mean'] for k in groups],color=[GRAY,RED,BLUE]);ax[1].set_xticks(np.arange(3),['All material','Early failure','Any failure']);ax[1].set(ylabel='Mean original m / complete-lattice m',title='Immutable original neighborhoods');fig.tight_layout();save(fig,'Fig05_surface_horizon','Real exposed_surface.vtp closest-triangle distance, not distance to the 453 loaded nodes. Early means full bulk failure before the first interface softening. Correlation is confounded by applied surface load and geometry; not a demonstrated surface-correction error.')
 fig,ax=plt.subplots(1,2,figsize=(13,5.2))
 for a,step in zip(ax,[1588,10502]):
  q=load_npz(O/'derived'/f'mixed_step_{step:06d}.npz');s=states[step];mask=(ref['bond_types']<3)&(s['d']<1)&q['valid'];ii=np.flatnonzero(mask);colour=np.log10(np.maximum(q['virtual_release_J'][ii],1e-18));p=a.scatter(q['stretch'][ii],q['shear'][ii],c=colour,s=6,cmap='viridis',vmin=-18,vmax=-10,rasterized=True);a.axvline(.15,color=RED,lw=1);a.axhline(.15,color=GRAY,lw=1);a.set(xlabel='Actual relative stretch',ylabel='Objective transverse diagnostic',title=f"DIAGNOSTIC ONLY | {s['time']*1e3:.3f} ms");a.set_xlim(-.35,.2);a.set_ylim(0,1.5)
 fig.colorbar(p,ax=ax,label='log10(ΔU_virtual [J]); not Gc',shrink=.8);save(fig,'Fig06_mixed_mode_quantitative','All valid live bulk/interface connections, no sampling. The 0.15 shear line is a diagnostic scale, not a failure criterion. Rank-deficient local stars omitted and counted in JSON. Energy comes from full endpoint-star virtual removal. Display limits crop extreme kinematic values; original FP64 values retained in NPZ.')
 fig=plt.figure(figsize=(13,5.5));s=states[1588];q=load_npz(O/'derived/mixed_step_001588.npz');valid=q['valid']&(ref['bond_types']<3)&(s['d']<1);a=valid&(s['kappa']<S0[ref['bond_types']])&(q['shear']>=.15);e=valid&(s['kappa']>=S0[ref['bond_types']]);mid=.5*(s['x'][ref['bonds'][:,0]]+s['x'][ref['bonds'][:,1]])
 for k,(mask,title,col) in enumerate([(a,'Below onset + shear diagnostic ≥0.15',GOLD),(e,'Stretch-history softening',RED)]):
  ax=fig.add_subplot(1,2,k+1,projection='3d');ax.scatter(*((s['x'][:n]-origin)*1e3).T,c=GRAY,s=3,alpha=.15);ax.add_collection3d(Line3DCollection((s['x'][ref['bonds'][mask]]-origin)*1e3,colors=col,linewidths=.6,alpha=.7));axes3(ax,bound);ax.set_title('DIAGNOSTIC ONLY\n'+title)
 save(fig,'Fig07_mixed_mode_locations','Original live bond segments at the first detachment keyframe, same physical scale/camera as the global V3 figure. Orange connections are kinematic observations and have NOT been virtually or physically deleted.')
 fig,ax=plt.subplots(1,2,figsize=(12,4.5));lab,sz,att=topology(ref,states[10502]['d']);free=sz[~att];u,cc=np.unique(free,return_counts=True);ax[0].bar(u,cc,color=RED);ax[0].set_xticks([1,2,3]);ax[0].set(xlabel='Nodes per component',ylabel='Number of components',title='ACTUAL DETACHED at 2.33 ms');ax[1].bar(u,cc*u*RHO*V*1e9,color=BLUE);ax[1].set_xticks([1,2,3]);ax[1].set(xlabel='Nodes per component',ylabel='Detached mass [µg]',title='Mass by actual component size');fig.tight_layout();save(fig,'Fig08_actual_fragment_distribution','Mass counts material only, no anchors. One- and two-node debris are unresolved particle-scale products; not labeled clinical microemboli. Discrete equivalent diameters are not converged physical fragment diameters.')
 if (O/'SMALL_SPECIMEN_VALIDATION.json').exists():benchmark_figures()
 dump(F/'FIGURE_PROVENANCE.json',dict(captions=CAPTIONS,origin_m=origin,common_V3_bounds=bound,displacement_scale=1,camera=dict(elevation_deg=18,azimuth_deg=-58),source='accepted FP64 states/original IDs; no artificial surfaces or fragments',wall_s=time.perf_counter()-start))
 (F/'index.html').write_text('<!doctype html><html lang="zh"><meta charset="utf-8"><title>V3 chunk mechanism audit</title><style>body{max-width:1500px;margin:auto;padding:30px;font:18px/1.6 sans-serif;background:#f7f9fb;color:#172b3a}img{width:100%;background:white}article{margin:35px 0}h1{font-size:28px}</style><h1>V3_CHUNK_FRAGMENTATION_MECHANISM_V1</h1><p>真实点与存活原始键；位移 ×1。图例明确区分实际脱离、拴住候选、离线诊断及独立合成试样。</p>'+''.join('<article><h2>'+html.escape(k)+'</h2><a href="'+k+'.png"><img src="'+k+'.png"></a><p>'+html.escape(v)+'</p></article>' for k,v in CAPTIONS.items())+'</html>')

def benchmark_figures():
 cases=[]
 for name in ['opening_z','shear_x']:
  meta=json.loads((O/'SMALL_SPECIMEN_VALIDATION.json').read_text());root=O/meta.get('benchmark_root','benchmark')/name;ref=load_npz(root/'provenance/reference.npz');summary=json.loads((root/'SUMMARY.json').read_text());keys=[load_npz(p) for p in sorted((root/'keys').glob('*.npz'))];hist=np.genfromtxt(root/'history.csv',delimiter=',',names=True);cases.append((name,root,ref,summary,keys,hist))
 origin=np.zeros(3);bound=limits([s['x'] for _,_,_,_,kk,_ in cases for s in kk],origin);fig=plt.figure(figsize=(18,11))
 for row,(name,root,ref,summary,keys,hist) in enumerate(cases):
  peakstep=int(hist['step'][np.argmax(abs(hist['reaction_N']))]);peak=min(keys,key=lambda s:abs(int(s['step'])-peakstep))
  for col,s in enumerate([keys[0],peak,keys[-1]]):
   ax=fig.add_subplot(2,3,row*3+col+1,projection='3d');draw(ax,ref,s,origin,bound);ax.set_title(f"SYNTHETIC BENCHMARK | {name}\n{s['time']*1e3:.3f} ms; saved accepted state")
 fig.suptitle('Original material points + surviving bulk bonds; displacement ×1',fontsize=18);save(fig,'Fig09_benchmark_actual_states','Common scale/camera for both cases. Intermediate state is the saved keyframe nearest the reaction peak, not an interpolated peak state. Red denotes no surviving path to fixed anchors; before 2 ms a red component may still contain prescribed grip nodes and is not freely moving. Free-fragment counts exclude gripped components. The same grip-release protocol is used in both cases.')
 fig,axes=plt.subplots(2,3,figsize=(17,9))
 for row,(name,root,ref,summary,keys,h) in enumerate(cases):
  t=h['time_s']*1e3;ax=axes[row,0]
  for field,label,col in [('W_J','External work','black'),('U_J','Stored U',BLUE),('K_J','Kinetic K',GOLD),('D_J','Damage D',RED),('Uc_J','Contact Uc',PURPLE)]:ax.plot(t,h[field]*1e9,label=label,color=col)
  ax.set(title=name+' | energy',xlabel='Time [ms]',ylabel='Energy [nJ]');ax.legend(fontsize=9);ax=axes[row,1];ax.plot(t,h['reaction_N']*1e6,color=BLUE,label='Top reaction');ax.plot(t,h['anchor_reaction_N']*1e6,color=GRAY,label='Anchor reaction');ax.set(xlabel='Time [ms]',ylabel='Reaction [µN]',title='Actual dynamics (not quasistatic)');ax.legend(fontsize=10);ax=axes[row,2];ax.plot(t,h['bulk_failed'],color=BLUE,label='Bulk failed');ax.plot(t,h['interface_failed'],color=PURPLE,label='Interface failed');ax.plot(t,h['largest_detached_nodes'],color=RED,label='Largest free component');ax.set(xlabel='Time [ms]',ylabel='Count',title='Natural failures + free topology');ax.legend(fontsize=10)
  for ax in axes[row]:ax.axvline(2,color='black',ls='--',lw=.8)
 fig.tight_layout();save(fig,'Fig10_benchmark_energy_and_topology','Only loading direction differs. Dashed line: top displacement constraint released at zero prescribed terminal velocity; no material bond changes at release. Before release free topology excludes all gripped components. Full residual, original IDs, reactions and inertia are retained; U, D and residual are not calibrated Gc.')
if __name__=='__main__':main()
