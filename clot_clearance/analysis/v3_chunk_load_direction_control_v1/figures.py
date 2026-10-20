from common import *
from analyze import Archive
from graph_ops import topology
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import html,time
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':13,'axes.titlesize':15,'axes.labelsize':14,'xtick.labelsize':11,'ytick.labelsize':11,'legend.fontsize':11,'axes.grid':True,'grid.color':'#808a94','grid.alpha':.65,'figure.facecolor':'white'})
BLUE='#1676a1';RED='#ce3550';GOLD='#df971a';GRAY='#9aa4ad';PURPLE='#7258a5';CAP={}
def save(fig,name,caption):
 for ext in ['png','svg']:fig.savefig(F/(name+'.'+ext),dpi=240,bbox_inches='tight')
 plt.close(fig);CAP[name]=caption;print('FIG',name,flush=True)
def bounds(arrays,origin):
 q=(np.concatenate(arrays)-origin)*1e3;lo=q.min(0);hi=q.max(0);return (lo+hi)/2,float((hi-lo).max()*1.10)
def axes3(ax,bound):
 center,span=bound
 for k,a in enumerate('xyz'):getattr(ax,'set_'+a+'lim')(center[k]-span/2,center[k]+span/2);getattr(ax,'set_'+a+'label')(a.upper()+' [mm]',labelpad=8)
 ax.set_box_aspect((1,1,1));ax.view_init(elev=18,azim=-58);ax.tick_params(pad=1,labelsize=10)
 for a in [ax.xaxis,ax.yaxis,ax.zaxis]:a.set_major_locator(MaxNLocator(nbins=4));a.pane.set_facecolor((.97,.98,.99,1));a._axinfo['grid']['color']=(.47,.51,.55,.7)
def draw(ax,ref,s,origin,bound,highlight=None,mask=None):
 x=(s['x']-origin)*1e3;n=int(ref['clot_count']);lab,sz,att=topology(ref,s['d']);free=~att[lab];mask=np.ones(len(x),bool) if mask is None else mask;b=ref['bonds'];live=s['d']<1;edge=b[(ref['bond_types']==1)&live];edge=edge[mask[edge].all(1)]
 ax.add_collection3d(Line3DCollection(x[edge],colors=BLUE,linewidths=.22,alpha=.13,rasterized=True))
 for ids,c,z in [(np.flatnonzero(mask[n:])+n,GRAY,2),(np.flatnonzero(mask[:n]&~free),BLUE,7),(np.flatnonzero(mask[:n]&free),RED,15)]:
  if len(ids):ax.scatter(*x[ids].T,c=c,s=z,depthshade=False,edgecolors='none',rasterized=True)
 if highlight is not None:
  ax.scatter(*x[highlight].T,c=GOLD,s=25,depthshade=False,edgecolor='black',linewidth=.3,rasterized=True);inside=np.isin(b[:,0],highlight)&np.isin(b[:,1],highlight);outside=np.isin(b[:,0],highlight)^np.isin(b[:,1],highlight)
  for ids,c,lw,alpha in [(np.flatnonzero(inside&live),GOLD,1.,.6),(np.flatnonzero(outside&live),'black',1.1,.8),(np.flatnonzero(outside&~live),RED,.3,.12)]:
   if len(ids):ax.add_collection3d(Line3DCollection(x[b[ids]],colors=c,linewidths=lw,alpha=alpha,rasterized=True))
 axes3(ax,bound)
def mapping_plot(ref):
 F.mkdir(parents=True,exist_ok=True);X=ref['X'];n=int(ref['clot_count']);origin=X[:n].mean(0);nf=load_npz(O/'REFERENCE_NORMAL_FIELD.npz');ids=nf['atom_ids']-1;req=json.loads((O/'provenance/request.json').read_text());card=json.loads((O/'provenance/LOAD_PROVIDER.json').read_text());g=np.exp(-np.sum((X[ids]-req['load_centers_m'][0])**2,axis=1)/(2*req['gaussian_sigma_m']**2));magnitude=ref['areas'][ids]*card['scale_Pa']*g*np.linalg.norm(ref['tangent'][ids],axis=1)
 # Fixed farthest-point subset only for arrow legibility; all 453 IDs are plotted.
 chosen=[int(np.argmax(g))];distance=np.full(len(ids),np.inf)
 while len(chosen)<48:distance=np.minimum(distance,np.linalg.norm(X[ids]-X[ids[chosen[-1]]],axis=1));chosen.append(int(np.argmax(distance)))
 chosen=np.array(chosen);arrow_scale=1.2e5;directions=[ref['tangent'][ids]/np.linalg.norm(ref['tangent'][ids],axis=1)[:,None],nf['normal_ref']];bb=bounds([X],origin);bb=(bb[0],bb[1]*1.25);fig=plt.figure(figsize=(15,7))
 for i,(v,title) in enumerate(zip(directions,['HISTORICAL C1 TANGENT','REFERENCE OUTWARD NORMAL'])):
  ax=fig.add_subplot(1,2,i+1,projection='3d');x=(X-origin)*1e3;ax.scatter(*x[n:].T,c=GRAY,s=3,depthshade=False,rasterized=True);ax.scatter(*x[:n].T,c=BLUE,s=4,alpha=.3,depthshade=False,rasterized=True);p=ax.scatter(*x[ids].T,c=g,cmap='viridis',vmin=0,vmax=1,s=13,depthshade=False,rasterized=True);vectors=v[chosen]*magnitude[chosen,None]*arrow_scale;ax.quiver(*x[ids[chosen]].T,*vectors.T,color='black',linewidth=1.2,arrow_length_ratio=.25);axes3(ax,bb);ax.set_title(title+'\n453 loaded IDs | same areas')
 fig.colorbar(p,ax=fig.axes[:2],shrink=.58,pad=.03,label='Gaussian at reference coordinates');fig.legend(handles=[Line2D([],[],marker='o',ls='',c=BLUE,label='Material (1,927)'),Line2D([],[],marker='o',ls='',c=GRAY,label='Fixed anchors (1,561)')],loc='lower center',ncol=2);fig.suptitle('Synthetic direction control | displacement ×1 | arrow scale: 120 µm / µN',fontsize=17)
 save(fig,'Fig01_reference_load_directions','1,927 material points and 1,561 fixed anchors; all 453 loaded IDs colored. Forty-eight deterministic farthest-point arrows show forces at the reference geometry and unit pulse amplitude, 120 micrometres per micronewton. Direction is fixed in reference space; the actual Gaussian uses moving coordinates. Geometry is not scaled. Normal orientation and same-input force magnitudes are independently checked.')
 return origin,ids[chosen]+1,arrow_scale
def main():
 started=time.perf_counter();old=Archive(B);new=Archive(O);ref=new.ref;X=ref['X'];n=int(ref['clot_count']);origin,arrow_ids,arrow_scale=mapping_plot(ref)
 matches=json.loads((O/'VISUAL_STATE_MATCHING.json').read_text())['matches'];states=[[load_npz(O/v[label]['state_path']) for v in matches] for label in ['historical_tangent','reference_normal']];bound=bounds([s['x'] for row in states for s in row],origin)
 for part,start in enumerate(range(0,len(matches),3)):
  subset=matches[start:start+3];fig=plt.figure(figsize=(18,10))
  for row,label in enumerate(['historical_tangent','reference_normal']):
   for col,match in enumerate(subset):
    s=states[row][start+col];ax=fig.add_subplot(2,len(subset),row*len(subset)+col+1,projection='3d');draw(ax,ref,s,origin,bound);ax.set_title(('Historical tangent' if row==0 else 'Reference normal')+f"\nactual {s['time']*1e3:.6f} ms")
  fig.subplots_adjust(top=.88,bottom=.10,hspace=.30);fig.legend(handles=[Line2D([],[],marker='o',ls='',c=c,label=l) for c,l in [(BLUE,'Attached'),(RED,'ACTUAL DETACHED'),(GRAY,'Fixed anchor')]],loc='lower center',ncol=3);fig.suptitle('Actual accepted states | same physical scale / camera | displacement ×1',fontsize=18)
  save(fig,f'Fig02{chr(97+part)}_matched_V3_states','Requested matching times [ms]: '+', '.join(f"{v['requested_time_s']*1e3:.6f}" for v in subset)+'. Each panel labels its actual accepted time. Nearest saved state matching is recorded with exact offsets in VISUAL_STATE_MATCHING.json; no interpolated nodes or fracture events. All live bulk bonds are drawn; interface/scaffold lines omitted for readability.')
 fig=plt.figure(figsize=(15,6.5));first_state_records=[]
 for i,(ar,label) in enumerate([(old,'Historical tangent'),(new,'Reference normal')]):
  q=json.loads((ar.root/'FIRST_EVENTS.json').read_text())['bulk_failure'];s=ar.state(q['time_s']);ax=fig.add_subplot(1,2,i+1,projection='3d');draw(ax,ref,s,origin,bound);ax.set_title(label+f"\nfirst bulk failure {q['time_s']*1e3:.6f} ms\nactual displayed state {s['time']*1e3:.6f} ms");first_state_records.append(dict(case=label,event_time_s=q['time_s'],displayed_time_s=float(s['time']),displayed_step=int(s['step'])))
 fig.subplots_adjust(top=.76,bottom=.12);fig.suptitle("Each case's first natural bulk failure | different event times | displacement ×1",fontsize=16)
 save(fig,'Fig02c_first_failure_events','Complement to the matched-time panels: each case at its own first natural bulk failure. Both event time and actual saved display time are explicit; these two panels are not a simultaneous comparison. Common camera and physical bounds are retained.')
 topo=json.loads((O/'V3_DIRECTION_CONTROL_TOPOLOGY.json').read_text());cand=json.loads((O/'V3_CHUNK_CANDIDATES.json').read_text());group=topo['reference_normal']['largest_historical_group'];candidates=[c for c in cand['candidates'] if c['nodes']>=14]
 if candidates:highlight=np.array(max(candidates,key=lambda c:c['nodes'])['atom_ids'])-1;t=max(candidates,key=lambda c:c['nodes'])['time_s'];label='TETHERED CANDIDATE'
 elif group:highlight=np.array(group['atom_ids'])-1;t=group['first_separation_time_s'];label='ACTUAL DETACHED'
 else:highlight=None;t=.0008;label='NO ACTUAL DETACHED COMPONENT'
 initial=new.state(t);final=new.keys[-1];mask=None
 if highlight is not None:
  near=np.linalg.norm(ref['X'][:,None,:]-ref['X'][highlight][None,:,:],axis=2).min(1)<HORIZON;mask=near;bb=bounds([initial['x'][near],final['x'][near]],origin)
 else:bb=bound
 fig=plt.figure(figsize=(15,6.5))
 for i,s in enumerate([initial,final]):
  ax=fig.add_subplot(1,2,i+1,projection='3d');draw(ax,ref,s,origin,bb,highlight,mask);ax.set_title(f"{label if i==0 else 'SAME ORIGINAL IDENTITIES'}\nactual {s['time']*1e3:.6f} ms")
 selected='none' if highlight is None else ', '.join(map(str,highlight+1));diameter=0 if highlight is None else (6*len(highlight)*V/np.pi)**(1/3)*1e6
 fig.subplots_adjust(top=.75,bottom=.16);fig.suptitle(f"Normal case: ≥14-node candidates {len(candidates)}; largest free {topo['reference_normal']['history_maxima']['largest_detached_nodes']} node(s)\nSelected atom ID(s): {selected}; reference-volume equivalent diameter {diameter:.2f} µm",fontsize=15)
 fig.legend(handles=[Line2D([],[],marker='o',ls='',c=GOLD,label='Selected original atom(s)'),Line2D([],[],c='black',label='Live external bond'),Line2D([],[],c=RED,label='Failed external bond')],loc='lower center',ncol=3)
 save(fig,'Fig03_local_cohesion','Orange highlights the selected original identity set, not an artificially merged fragment. Gold lines: surviving internal bonds; black: surviving external bonds; faint red: originally external bonds that actually failed. Exact formation event topology is in JSON; a sampled display time is labeled explicitly if no full nodal state exists at that event. Negative candidate results use the unchanged prior graph criteria. IDs and internal survival are retained in CHUNK_COHESION_AUDIT.json.')
 fig=plt.figure(figsize=(15,6.5));bb=bounds([X],origin)
 for i,(label,title) in enumerate([('historical_tangent','Historical tangent'),('reference_normal','Reference normal')]):
  fields=load_npz(O/'derived'/f'{label}_event_fields.npz');spatial=load_npz(O/'derived'/f'{label}_spatial_fields.npz');mask=np.isfinite(fields['fail_time_s'])&(ref['bond_types']<3);ax=fig.add_subplot(1,2,i+1,projection='3d');ax.scatter(*((X[:n]-origin)*1e3).T,s=2,c=GRAY,alpha=.15,rasterized=True);pts=(spatial['midpoints_m'][mask]-origin)*1e3;p=ax.scatter(*pts.T,c=fields['fail_time_s'][mask]*1e3,s=7,cmap='plasma',vmin=0,vmax=END*1e3,depthshade=False,rasterized=True);axes3(ax,bb);ax.set_title(title+'\nFailed bond midpoints (diagnostic)')
 fig.colorbar(p,ax=fig.axes[:2],shrink=.65,label='Actual accepted failure time [ms]');save(fig,'Fig03b_damage_organization','Original reference bond midpoints colored by actual full-failure time. Common reference coordinates distinguish changes in where damage initiates and spreads. A connected damage cloud is not a reconstructed continuous crack surface; bulk and interface regions are quantified separately in BULK_INTERFACE_FAILURE_AUDIT.json.')
 fig,axes=plt.subplots(2,3,figsize=(17,9));labels=['Historical tangent','Reference normal']
 for series,(ar,label,col) in enumerate(zip([old,new],labels,[BLUE,RED])):
  h=ar.hist;t=h['time_s']*1e3;ll,sz,att=topology(ref,ar.keys[-1]['d']);free=sz[~att];u,cc=np.unique(free,return_counts=True)
  axes[0,0].plot(t,h['detached_mass_kg']*1e9,label=label,c=col);axes[0,1].plot(t,h['detached_components'],label=label,c=col);axes[0,2].plot(t,h['largest_detached_nodes'],label=label,c=col);axes[1,0].plot(t,h['tiny_nodes'],label=label,c=col)
  bins=np.arange(1,max(topo['historical_tangent']['final']['largest_detached_nodes'],topo['reference_normal']['final']['largest_detached_nodes'])+1);axes[1,1].bar(bins+(series-.5)*.34,[np.sum(free==k) for k in bins],width=.34,label=label,color=col);axes[1,1].set_xticks(bins)
 for ax,y,title in zip(axes.ravel()[:5],['Mass [µg]','Components','Nodes','Nodes','Components'],['Actual detached mass','Actual free component count','Historical maximum visible in time series','Unresolved 1–2-node debris','Final size distribution']):ax.set(ylabel=y,title=title,xlabel='Time [ms]');ax.legend()
 axes[1,1].set_xlabel('Nodes per component');ax=axes[1,2];xx=np.arange(3)
 for k,(key,label,col) in enumerate(zip(['historical_tangent','reference_normal'],labels,[BLUE,RED])):
  q=topo[key];ax.bar(xx+(k-.5)*.34,[q['history_maxima']['N_ge'+str(n)] for n in [14,34,114]],width=.34,color=col,label=label)
 ax.set_ylim(0,1);ax.set_yticks([0,1]);ax.text(.5,.45,'0 in both cases at all accepted times',transform=ax.transAxes,ha='center',fontsize=11);ax.set_xticks(xx,['≥14','≥34 / ≥200µm','≥114 / ≥300µm']);ax.set(ylabel='Maximum concurrent count',title='Actual components: discrete volume thresholds');ax.legend();fig.tight_layout();save(fig,'Fig04_actual_fragment_sizes','Historical maxima and final sizes are distinguished. Component membership uses live original bulk/interface paths to fixed anchors, never contact proximity. The 200/300 micrometre labels are equivalent reference-volume thresholds, not converged physical size or clinical safety targets.')
 fig,axes=plt.subplots(3,3,figsize=(18,13))
 panels=[('load_sum_norm_N',1e6,'Sum of point-force magnitudes [µN]'),('W_load_J',1e9,'Actual artificial-load work [nJ]'),('U_J',1e9,'Recoverable energy U [nJ]'),('K_J',1e9,'Kinetic energy K [nJ]'),('D_J',1e9,'Irreversible damage D [nJ]'),('W_drag_J',1e9,'Background COM drag work [nJ]'),('W_numeric_J',1e9,'Attached numerical-damping work [nJ]'),('max_contact_penetration_dx',1.,'Particle overlap / dx'),('relative_R',100.,'Independent normalized residual [%]')]
 for ax,(field,scale,title) in zip(axes.ravel(),panels):
  for ar,label,col in zip([old,new],labels,[BLUE,RED]):ax.plot(ar.hist['time_s']*1e3,ar.hist[field]*scale,label=label,c=col)
  ax.set(xlabel='Time [ms]',title=title);ax.legend()
 fig.tight_layout();save(fig,'Fig05a_forces_work_energy','Same scalar load rule is not equal actual work. Exact force-displacement work and independent residual are preserved; the residual is not damage dissipation. Force magnitude is computed at each case’s own evolving coordinates. Full signed force vector and all channels are in the FP64 history.')
 fig,axes=plt.subplots(2,3,figsize=(17,8.5))
 for ax,(field,scale,title) in zip(axes.ravel(),[('Uc_J',1e9,'Particle contact storage [nJ]'),('Uw_J',1e9,'Wall contact storage [nJ]'),('Jc_J',1e9,'Contact topology term [nJ]'),('W_wall_damping_J',1e9,'Wall damping work [nJ]'),('bulk_failed',1,'Natural bulk failures'),('interface_failed',1,'Natural interface failures')]):
  for ar,label,col in zip([old,new],labels,[BLUE,RED]):ax.plot(ar.hist['time_s']*1e3,ar.hist[field]*scale,label=label,c=col)
  ax.set(xlabel='Time [ms]',title=title);ax.legend()
 fig.tight_layout();save(fig,'Fig05b_contact_and_typed_failures','Complementary channels in W_load+W_drag+W_numeric+W_wall_damping=ΔU+ΔK+D+ΔUc+ΔUw−Jc+R. Wall elastic work is recorded separately and not counted again as external input when wall storage Uw is included.')
 fig,axes=plt.subplots(1,3,figsize=(17,4.8))
 for ar,label,col in zip([old,new],labels,[BLUE,RED]):
  h=ar.hist;t=h['time_s']*1e3;axes[0].plot(t,h['load_amp'],c=col,label=label)
  for component,style in zip('xyz',['-','--',':']):axes[1].plot(t,h['load_'+component+'_N']*1e6,c=col,ls=style,label=label+' '+component)
  axes[2].plot(t,h['W_wall_elastic_J']*1e9,c=col,label=label)
 for ax,title,y in zip(axes,['Identical prescribed pulse','Actual resultant load vector','Wall elastic work (diagnostic)'],['Amplitude','Force [µN]','Work [nJ]']):ax.set(xlabel='Time [ms]',ylabel=y,title=title);ax.legend(fontsize=9)
 fig.tight_layout();save(fig,'Fig05c_pulse_and_resultant','The original pulse is inherited without override. Resultants differ from the sum of individual force magnitudes because the outward normal varies across the curved surface. Wall elastic work is a diagnostic channel and is not double counted in the primary ledger.')
 fig,ax=plt.subplots(figsize=(15,8));ax.axis('off');audit=json.loads((O/'NUMERICAL_AUDIT.json').read_text());cost=json.loads((O/'RUN_RESOURCES.json').read_text());normal=topo['reference_normal'];oldtop=topo['historical_tangent'];energy=json.loads((O/'DYNAMIC_ENERGY_LEDGER.json').read_text());first=json.loads((O/'FIRST_EVENTS.json').read_text())
 records=[['Reference normals','453 / 453 outward into lumen; original area and scalar rule retained'],['Completed time',f"{new.hist['time_s'][-1]*1e3:.6f} ms; one new native physical case"],['First bulk soft / failure',' / '.join('not observed' if first[k] is None else f"{first[k]['time_s']*1e3:.6f} ms" for k in ['bulk_softening','bulk_failure'])],['Largest free: history / final',f"normal {normal['history_maxima']['largest_detached_nodes']} / {normal['final']['largest_detached_nodes']} nodes; tangent {oldtop['history_maxima']['largest_detached_nodes']} / {oldtop['final']['largest_detached_nodes']}"],['Tethered candidates ≥14 / ≥34',f"{len(candidates)} / {cand['candidate_ge200']} (unchanged diagnostic criteria)"],['Load-work ratio: normal / tangent',f"{energy['normal_over_tangent_load_work']:.6g}; equal-energy efficiency NOT TESTED"],['Numerical admission',audit['status']],['Native cost',f"wall {cost['main_wall_s']:.2f} s; CPU {cost['main_CPU_s']:.2f} s; peak RSS {cost['peak_RSS_bytes']/1024**3:.3f} GiB"],['Scientific limits','Synthetic force; interface tail / Gc / size convergence uncalibrated']]
 records.insert(3,['First failure reference midpoint','Normal: (−6.000, 0.125, 0.700) mm; tangent: (−6.175, −0.050, 0.750) mm']);records.insert(6,['Secondary breakup','Normal: not applicable (all singletons); historical 13 nodes split 9 / 3 / 1'])
 table=ax.table(cellText=records,colLabels=['Measure','Actual result / limitation'],colWidths=[.3,.7],cellLoc='left',loc='center');table.auto_set_font_size(False);table.set_fontsize(11.5);table.scale(1,2.3);ax.set_title('Direction-control evidence and measured cost',fontsize=19,pad=20);save(fig,'Fig06_evidence_and_cost','Summary of actual accepted numerical results; neither normal traction nor the retained artificial C1 loading is a verified microbubble force. A single direction comparison cannot identify equal-input-energy efficiency or calibrate material Gc.')
 dump(F/'FIGURE_PROVENANCE.json',dict(captions=CAP,displacement_scale=1,camera=dict(elevation_deg=18,azimuth_deg=-58),origin_m=origin,V3_common_bounds=bound,first_failure_event_displays=first_state_records,arrow_scale_mm_per_N=arrow_scale,arrow_atom_ids=arrow_ids,wall_s=time.perf_counter()-started,source='Original accepted FP64 positions and live original bond IDs; no invented surfaces'))
 page='<!doctype html><html lang="zh"><meta charset="utf-8"><title>V3 载荷方向对照</title><style>body{max-width:1500px;margin:auto;padding:30px;background:#f7f9fb;color:#172b3a;font:18px/1.6 sans-serif}img{width:100%;background:white}article{margin:40px 0}</style><h1>V3 载荷方向对照：历史切向 / 参考外法向</h1><p>真实接受状态，位移 ×1。材料与Gaussian标量规则相同；实际输入功可能不同。法向为合成力学对照。</p>'+''.join('<article><h2>'+html.escape(k)+'</h2><a href="'+k+'.png"><img src="'+k+'.png"></a><p><a href="'+k+'.svg">SVG</a> · '+html.escape(v)+'</p></article>' for k,v in CAP.items())+'</html>'
 (F/'OPEN_RESULTS.html').write_text(page);(O/'OPEN_RESULTS.html').write_text('<!doctype html><meta charset="utf-8"><a href="../../../results/'+TASK+'/'+RUN+'/OPEN_RESULTS.html">打开本轮全部图表</a>')
if __name__=='__main__':
 import sys
 if '--mapping-only' in sys.argv:mapping_plot(load_npz(O/'provenance/reference.npz'))
 else:main()
