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
from multifragment_shear_objective import shape
import pandas as pd,time,resource
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.titlesize':14,'axes.labelsize':12,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'figure.facecolor':'white','savefig.facecolor':'white','axes.grid':True,'grid.color':'#788590','grid.alpha':.6,'pdf.fonttype':42})
DAMAGE=LinearSegmentedColormap.from_list('clot_damage',['#37779b','#d6ba78','#c54a54','#54234f'])
NAMES={'uniform':'U · Uniform reference','biphasic':'H · Biphasic'}
COLORS={'uniform':'#3e7aa2','biphasic':'#a04e6a'}
graph=module('fig_graph',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')

def save(fig,name):
 p=F/'figures';p.mkdir(exist_ok=True,parents=True);fig.savefig(p/(name+'.png'),dpi=250,bbox_inches='tight');fig.savefig(p/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
def ax3(fig,pos,title,limits):
 ax=fig.add_subplot(*(pos if isinstance(pos,tuple) else (pos,)),projection='3d',computed_zorder=False);ax.view_init(elev=16,azim=-64);ax.set_title(title,pad=18);lo,hi=limits;pad=np.maximum((hi-lo)*.06,20.)
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

def parameter_figures(ref,field):
 audit=json.loads((O/'PREFLIGHT_STRUCTURE_AUDIT.json').read_text());bulk=ref['bond_types']==1
 fig,axs=plt.subplots(2,2,figsize=(13,9),layout='constrained');labels=['WEAK','MIXED','STRONG','UNIFORM'];on=np.array([.10,.175,.25,.15]);off=np.array([.12,.21,.30,.18]);xx=np.arange(4)
 axs[0,0].bar(xx-.15,on,.3,label='s0');axs[0,0].bar(xx+.15,off,.3,label='sc');axs[0,0].set_xticks(xx,labels);axs[0,0].set(title='Immutable pair thresholds',ylabel='Stretch');axs[0,0].legend()
 for j,(key,label) in enumerate([('bond_s0','s0'),('bond_sc','sc')]):
  values,counts=np.unique(field[key][bulk],return_counts=True);axs[0,1].plot(values,counts,'o-',label=label)
 axs[0,1].set(title='Actual original bulk pair properties',xlabel='Threshold',ylabel='Bond count');axs[0,1].legend()
 domain=pd.read_csv(O/'STRONG_DOMAINS.csv');axs[1,0].bar(domain.domain_id,domain.nodes,color=plt.get_cmap('tab10').colors);axs[1,0].set(title='Reference strong regions',xlabel='Domain ID',ylabel='Nodes');axs[1,0].set_xticks(domain.domain_id)
 axs[1,1].axis('off');axs[1,1].text(.03,.95,f'Strong fraction: {audit["strong_volume_fraction"]:.5%}\nStrong domains (6-face grid): {audit["strong_domain_count"]}\nCross-phase bulk bonds: {audit["cross_phase_material_bonds"]:,}\nCross-phase / bulk: {audit["cross_phase_fraction_of_bulk"]:.2%}\nNode mean s0 / sc: {audit["mean_node_s0"]:.6f} / {audit["mean_node_sc"]:.6f}\nWeighted pair mean: {audit["weighted_bulk_s0"]:.6f} / {audit["weighted_bulk_sc"]:.6f}\ndx = 50 µm; horizon = 150.5 µm\nG = 940 Pa; K = 100000/3 Pa; rho = 1000 kg/m³\nNo property-dependent cuts; no energy renormalization',va='top',linespacing=1.7,fontsize=12)
 fig.suptitle('Frozen exploratory biphasic material · threshold ratios 2.5, not assumed energy ratios',fontsize=16);save(fig,'FIG_04_BIPHASIC_PARAMETERS')
 path=pd.read_csv(O/'PRESCRIBED_PATH_RESPONSE.csv');summary=json.loads((O/'PRESCRIBED_PATH_ENERGY.json').read_text());checks=json.loads((O/'VARIATIONAL_ENERGY_AUDIT.json').read_text())['checks']['variational'];fig,axs=plt.subplots(2,2,figsize=(14,9),layout='constrained')
 for name,g in path.groupby('case',sort=False):
  axs[0,0].plot(g.opening_m*1e6,g.reaction_N*1e6,label=name);axs[0,1].plot(g.opening_m*1e6,g.D_J*1e9,label=name)
 axs[0,0].set(title='Same full-model prescribed opening path',xlabel='Opening (µm)',ylabel='Reaction (µN)');axs[0,0].legend()
 axs[0,1].set(title='Damage energy drop summed at fixed coordinates',xlabel='Opening (µm)',ylabel='Model dissipation (nJ)');axs[0,1].legend()
 for j,(name,(v,w)) in enumerate(zip(labels,zip(on,off))):axs[1,0].plot([v,w],[j,j],lw=10,solid_capstyle='butt');axs[1,0].plot([v,w],[j,j],'ko',ms=4)
 axs[1,0].set_yticks(range(4),labels);axs[1,0].set(title='Onset and complete-damage interval',xlabel='Maximum tensile bond stretch κ')
 axs[1,1].bar([str(q['strain']) for q in checks],[q['relative_gradient_error'] for q in checks],color='#598c9d');axs[1,1].set_yscale('log');axs[1,1].axhline(1e-5,ls='--',color='#ae445e');axs[1,1].set(title='Native force vs fixed-damage energy gradient',xlabel='Prescribed verification strain',ylabel='Relative gradient error')
 fig.suptitle(f'Full endpoint-star energy · prescribed-path STRONG/WEAK dissipation = {summary["strong_weak_path_dissipation_ratio"]:.4f}\nDiscrete model-path diagnostic, not calibrated toughness Gc',fontsize=16);save(fig,'FIG_05_MATERIAL_MECHANICS')

def result_figures(ref,field):
 data=json.loads((O/'COMPARISON.json').read_text());cases=data['cases'];n=int(ref['clot_count']);b=ref['bonds'];cut=ref['initial_cut'];near=ref['r0']<=1.01*DX
 states={};hist={}
 for m in CASES:
  states[m]={p.stem:load_npz(p) for p in (O/'cases'/m/'keys').glob('*.npz') if p.stem!='failed_unaccepted'}
  if 'final' not in states[m]:states[m]['final']=states[m]['last_accepted']
  hist[m]=pd.read_csv(O/'cases'/m/'history.csv')
 allx=np.concatenate([s['x'] for d in states.values() for s in d.values()])*1e6;limits=(allx.min(0),allx.max(0));ref_limits=(ref['X'].min(0)*1e6,ref['X'].max(0)*1e6)
 # Five actual accepted states per case. No synthetic 'multisite' event is added.
 fig=plt.figure(figsize=(23,11));fig.subplots_adjust(wspace=.14,hspace=.30,bottom=.09)
 for row,m in enumerate(CASES):
  frames=sorted(states[m].values(),key=lambda s:float(s['time']));event=cases[m]['native']['events'].get('first_softening',{});origin=float(event.get('time_s',0));chosen=[states[m].get('first_softening'),states[m].get('first_failure'),min(frames,key=lambda s:abs(float(s['time'])-(origin+.0001))),min(frames,key=lambda s:abs(float(s['time'])-.0024)) if float(frames[-1]['time'])>=HOLD-1e-12 else None,states[m]['final']]
  for col,(s,title) in enumerate(zip(chosen,['First softening','First complete failure','Early spatial development','Near end of hold','Last accepted'])):
   ax=ax3(fig,(2,5,1+row*5+col),NAMES[m]+'\n'+title,limits)
   if s is None:ax.text2D(.15,.5,'NOT_OBSERVED',transform=ax.transAxes);continue
   scatter(ax,ref,s,size=3);sel=near&~cut&(s['d']==1)
   if sel.any():ax.add_collection3d(Line3DCollection(s['x'][b[sel]]*1e6,colors='#8c233f',alpha=.25,linewidths=.35))
   ax.text2D(.02,.94,f't = {float(s["time"])*1e3:.4f} ms',transform=ax.transAxes,fontsize=10)
   if col<2:
    ev=cases[m]['native']['events'].get('first_softening' if col==0 else 'first_failure',{});ids=np.array(ev.get('original_bond_ids',[]),int)-1
    if len(ids):ax.add_collection3d(Line3DCollection(s['x'][b[ids]]*1e6,colors='#f29a22',linewidths=2.4));ax.text2D(.02,.02,'Bond IDs: '+','.join(map(str,(ids+1)[:6]))+('…' if len(ids)>6 else ''),transform=ax.transAxes,fontsize=9)
 fig.suptitle('Actual accepted-state fracture development · displacement ×1 · identical physical scale\nDamage color excludes the initial notch; dark traces are failed nearest bonds, not surviving material',fontsize=16);save(fig,'FIG_06_TRUE_CRACK_EVOLUTION')
 # Reference-location markers make distinct initiation sites visible even when
 # the final full-panorama bounds are much larger than the initial specimen.
 fig=plt.figure(figsize=(16,9));fig.subplots_adjust(wspace=.18,bottom=.22,top=.84)
 for j,m in enumerate(CASES):
  sites=data['spatial'][m]['early_spatially_separated_softening_sites']
  ax=ax3(fig,(1,2,j+1),NAMES[m]+f' · {len(sites)} separated bulk-softening locations',ref_limits)
  ax.scatter(*ref['X'][:n].T*1e6,s=2,c='#a6b9c4',alpha=.16,depthshade=False)
  ax.scatter(*ref['X'][n:].T*1e6,s=2,c='#6d7884',alpha=.12,depthshade=False)
  text=[]
  for k,event in enumerate(sites):
   xyz=np.array(event['reference_midpoint_m'])*1e6
   ax.scatter(*xyz,s=110,c='#efad35',edgecolors='#503c24',marker='*',depthshade=False,zorder=10)
   ax.text(*(xyz+np.array([0,0,30.])),str(k+1),fontsize=12,color='#5b3545',zorder=11)
   text.append(f"{k+1}: bond {event['original_bond_id']}, t={event['time_s']*1e3:.5f} ms")
  if not sites:ax.text2D(.08,.5,'NOT_OBSERVED',transform=ax.transAxes)
  fig.text(.08+j*.49,.14,'; '.join(text[:3])+'\n'+'; '.join(text[3:6])+('\nAll additional IDs/times: SPATIAL_DAMAGE_ANALYSIS.json' if len(text)>6 else ''),fontsize=10,va='top')
 fig.suptitle('Original reference locations · common camera and physical scale\nFirst bulk onset + 0.1 ms per case; midpoint separation > 301 µm · not an independent continuum-crack count',fontsize=16)
 save(fig,'FIG_06B_EARLY_SOFTENING_LOCATIONS')
 palette=['#aebbc5','#bd4d67','#e5a24b','#428e83','#7654a5'];legend=['Attached body','Free: 1–2','Free: 3–19','Free: 20–200','Free: >200']
 fig=plt.figure(figsize=(16,12));fig.subplots_adjust(wspace=.17,hspace=.30,bottom=.1)
 for j,m in enumerate(CASES):
  s=states[m]['final'];lab,sz,att=graph.topology(ref,s['d']);cat=np.zeros(n,int)
  for label in np.flatnonzero(~att):cat[lab==label]=1 if sz[label]<=2 else 2 if sz[label]<20 else 3 if sz[label]<=200 else 4
  ax=ax3(fig,221+j,NAMES[m]+f' · last accepted\nt = {float(s["time"])*1e3:.5f} ms · {cases[m]["native"]["status"]}',limits);ax.scatter(*s['x'][n:].T*1e6,s=3,c='#747d86',alpha=.25,depthshade=False);ax.scatter(*s['x'][:n].T*1e6,s=5,c=np.array(palette)[cat],depthshade=False)
  nodefree=np.r_[~att[lab],np.zeros(len(ref['X'])-n,bool)];inside=nodefree[b].all(1)&(s['d']<1)&(ref['bond_types']==1)
  if inside.any():ax.add_collection3d(Line3DCollection(s['x'][b[inside]]*1e6,colors='#244e63',alpha=.75,linewidths=.65))
  labels=np.flatnonzero(~att);ids=np.flatnonzero(lab==labels[np.argmax(sz[labels])]) if len(labels) else np.array([],int);center=(s['x'][ids].min(0)+s['x'][ids].max(0))*.5e6 if len(ids) else np.zeros(3);width=max(300.,np.ptp(s['x'][ids]*1e6,axis=0).max()*1.3) if len(ids) else 300.
  ax=ax3(fig,223+j,NAMES[m]+f' · largest free component: {len(ids)} nodes',(center-width/2,center+width/2))
  if len(ids):
   ax.scatter(*s['x'][ids].T*1e6,s=35,c=palette[cat[ids[0]]],depthshade=False);inside=np.isin(b,ids).all(1)&(s['d']<1);ax.add_collection3d(Line3DCollection(s['x'][b[inside]]*1e6,colors='#244e63',linewidths=1.1,alpha=.8));ax.text2D(.02,.94,f'{inside.sum()} live internal bonds; displayed width {width+2*max(.06*width,20.):.0f} µm',transform=ax.transAxes)
 fig.legend(handles=[Line2D([],[],marker='o',ls='',color=c,label=l) for c,l in zip(palette,legend)],loc='lower center',ncol=5);fig.suptitle('Surviving original-bond graph + anchor reachability · no proximity regrouping\nCommon full-panorama scale above; explicitly labeled physical local detail below',fontsize=16);save(fig,'FIG_07_TRUE_COHESIVE_FRAGMENTS')
 if any(cases[m]['native']['max_medium_components']>0 for m in CASES):
  fig=plt.figure(figsize=(16,7))
  for j,m in enumerate(CASES):
   s=states[m]['max_multifragment'];lab,sz,att=graph.topology(ref,s['d']);free=sz[~att];target=int(((free>=20)&(free<=200)).sum());ax=ax3(fig,121+j,NAMES[m]+f' · saved target maximum {target}\nt = {float(s["time"])*1e3:.5f} ms',limits);scatter(ax,ref,s,'components',size=7)
   inside=np.r_[~att[lab],np.zeros(len(ref['X'])-n,bool)];sel=inside[b].all(1)&(s['d']<1)&(ref['bond_types']==1)
   if sel.any():ax.add_collection3d(Line3DCollection(s['x'][b[sel]]*1e6,colors='#254d63',linewidths=.6,alpha=.7))
  fig.suptitle('Actual saved maximum-count states · common camera/scale · each panel shows its own physical time',fontsize=15);save(fig,'FIG_07B_TRANSIENT_TARGET_EVIDENCE')
  # The same material IDs can form a curved strip with full 3D point-cloud rank.
  # Show their reference shape alongside the actual shape, without recentering data.
  target_morphology_figure(ref,states,limits)
 # Separate, explicitly labeled detail panels prevent distant fines from hiding
 # the actual internal connections of target-size components in the panorama.
 target_details=[]
 for m in CASES:
  s=states[m]['final'];lab,sz,att=graph.topology(ref,s['d']);registry=json.loads((O/'cases'/m/'COMPONENT_REGISTRY.json').read_text());lookup={tuple(q['atom_ids']):int(cid) for cid,q in registry.items()}
  for label in np.flatnonzero(~att&(sz>=20)&(sz<=200)):
   ids=np.flatnonzero(lab==label);sel=np.isin(b,ids).all(1)&(s['d']<1);cid=lookup[tuple((ids+1).tolist())]
   target_details.append(dict(case=m,state=s,ids=ids,bonds=b[sel],component_id=cid))
 if target_details:
  width=max(150.,1.3*max(np.ptp(q['state']['x'][q['ids']]*1e6,axis=0).max() for q in target_details));detail_evidence=[]
  for start in range(0,len(target_details),4):
   page=target_details[start:start+4];cols=min(2,len(page));rows=(len(page)+cols-1)//cols;fig=plt.figure(figsize=(5.4*cols,5.3*rows+1));fig.subplots_adjust(wspace=.22,hspace=.32,top=.82 if rows==1 else .88,bottom=.06)
   for j,q in enumerate(page):
    m=q['case'];s=q['state'];ids=q['ids'];center=(s['x'][ids].min(0)+s['x'][ids].max(0))*.5e6;ref_shape=shape(ref['X'][ids]);actual_shape=shape(s['x'][ids])
    ax=ax3(fig,(rows,cols,j+1),NAMES[m]+f' · component {q["component_id"]}\n{len(ids)} nodes; t = {float(s["time"])*1e3:.5f} ms',(center-width/2,center+width/2))
    ax.scatter(*s['x'][ids].T*1e6,s=30,c=COLORS[m],depthshade=False);ax.add_collection3d(Line3DCollection(s['x'][q['bonds']]*1e6,colors='#244e63',linewidths=.7,alpha=.55))
    ax.text2D(.02,.98,f'{len(q["bonds"])} surviving internal bonds\nReference: '+ref_shape['shape_class'].replace('THREE_DIMENSIONAL_POINT_CLOUD','3D point cloud').replace('_',' ').lower(),transform=ax.transAxes,fontsize=10,va='top')
    detail_evidence.append(dict(case=m,component_id=q['component_id'],time_s=float(s['time']),protocol_status=cases[m]['native']['status'],original_atom_ids=ids+1,internal_alive_bonds=len(q['bonds']),reference_shape=ref_shape,current_shape=actual_shape))
   fig.suptitle(f'Every target-band component at the last accepted states · page {start//4+1}\nEqual displayed coordinate width {width+2*max(.06*width,20.):.1f} µm · actual coordinates ×1 · full panoramas retained in FIG_07',fontsize=15);save(fig,f'FIG_07D_TARGET_DETAILS_{start//4+1:02d}')
  dump(O/'LAST_ACCEPTED_TARGET_DETAIL_EVIDENCE.json',dict(common_local_width_um=width+2*max(.06*width,20.),records=detail_evidence,interpretation='Separate local camera panels, never an exploded or artificially separated material view. Coordinates and all surviving internal material bonds are unchanged.'))
 dist=pd.read_csv(O/'FRAGMENT_SIZE_DISTRIBUTION.csv');fig,axs=plt.subplots(2,3,figsize=(17,10),layout='constrained')
 for j,m in enumerate(CASES):
  r=dist[dist['case']==m];xx=np.arange(len(r));axs[0,0].bar(xx+(j-.5)*.36,r['count'],.36,label=NAMES[m],color=COLORS[m]);axs[0,1].bar(xx+(j-.5)*.36,100*r.fraction_initial,.36,label=NAMES[m],color=COLORS[m]);h=hist[m];t=h.time_s*1e3
  for ax,key in [(axs[0,2],'largest_nodes'),(axs[1,0],'N20_200'),(axs[1,1],'mass_fraction20_200'),(axs[1,2],'debris12_detached_fraction')]:ax.plot(t,h[key]*(100 if 'fraction' in key else 1),label=NAMES[m],color=COLORS[m])
 for ax,title in [(axs[0,0],'Last accepted free component counts'),(axs[0,1],'Last accepted mass / initial clot (%)')]:ax.set_xticks(xx,r['bin'],rotation=25);ax.set(title=title,xlabel='Original nodes per component');ax.legend()
 for ax,title,y in [(axs[0,2],'Largest free component','Nodes'),(axs[1,0],'Inclusive 20–200 target count','Components'),(axs[1,1],'Target mass / initial clot','%'),(axs[1,2],'1–2 node fines / detached mass','%')]:ax.set(title=title,xlabel='Time (ms)',ylabel=y);ax.axvline(OFF*1e3,ls=':',color='gray');ax.legend()
 fig.suptitle('Every accepted-state true topology · exact 14–19 / 20–33 split · 20 and 200 included\n'+' | '.join(f'{m}: {cases[m]["native"]["status"]}, {cases[m]["native"]["time_s"]*1e3:.5f} ms' for m in CASES),fontsize=16);save(fig,'FIG_08_FRAGMENT_DISTRIBUTIONS')
 for m in CASES:
  h=hist[m];t=h.time_s*1e3;fig,axs=plt.subplots(2,3,figsize=(17,10),layout='constrained')
  axs[0,0].plot(t,h.shear_peak_Pa,label='Actual slip shear');axs[0,0].plot(t,h.nominal_shear_Pa,'--',label='Specified nominal drive');axs[0,0].set(title='Current-surface traction',ylabel='Pa');axs[0,0].legend()
  for k,label in [('W_drive_J','Drive'),('W_drag_J','Passive drag'),('W_J','Net work')]:axs[0,1].plot(t,h[k]*1e9,label=label)
  axs[0,1].set(title='Signed actual external work',ylabel='nJ');axs[0,1].legend()
  for k,label in [('U_J','Stored U'),('K_J','Kinetic K'),('D_J','Damage D'),('Uc_J','Contact Uc'),('Jc_J','Contact creation Jc')]:axs[0,2].plot(t,h[k]*1e9,label=label)
  axs[0,2].set(title='Complete physical ledger',ylabel='nJ');axs[0,2].legend(fontsize=9,ncol=2)
  axs[1,0].plot(t,h.R_J*1e9,label='Numerical R (nJ)');axs[1,0].set(title='Residual remains separate from damage',ylabel='nJ');other=axs[1,0].twinx();other.plot(t,h.relative_R*100,color='#ad4063');other.set_ylabel('Relative residual (%)');other.grid(False)
  axs[1,1].plot(t,h.bulk_failed,label='Bulk');axs[1,1].plot(t,h.interface_failed,label='Interface');axs[1,1].set(title='Natural complete bond failures',ylabel='Bonds');axs[1,1].legend()
  axs[1,2].plot(t,h.components,label='Free components');axs[1,2].plot(t,h.N20_200,label='20–200 target');axs[1,2].set(title='True graph components',ylabel='Count');axs[1,2].legend()
  for ax in axs.flat:ax.set_xlabel('Time (ms)');ax.axvline(OFF*1e3,ls=':',color='gray')
  fig.suptitle(NAMES[m]+' · '+cases[m]['native']['status']+' · Wdrive + Wdrag = ΔU + ΔK + D + ΔUc − Jc + R',fontsize=16);save(fig,'FIG_09_ENERGY_'+m.upper())
 fig,axs=plt.subplots(3,3,figsize=(18,13),layout='constrained');items=[('Last accepted target count',lambda q:q['objective']['n_fragments_20_200']),('Maximum target count',lambda q:q['native']['max_medium_components']),('Target mass / initial (%)',lambda q:100*q['objective']['mass_fraction_target_initial']),('Total detached / initial (%)',lambda q:100*q['objective']['total_detached_mass_fraction']),('Largest free (nodes)',lambda q:max(q['final_free_sizes'],default=0)),('Fines / detached (%)',lambda q:100*(q['objective']['fine_debris_mass_fraction_detached'] or 0)),('Unchanged target member sets\nsince drive-off',lambda q:q['objective']['drive_off_target_survivors_to_end']),('Native CPU time (s)',lambda q:q['native']['cpu_s'])]
 for ax,(title,fn) in zip(axs.flat,items):
  vals=[fn(cases[m]) for m in CASES];bars=ax.bar(['Uniform','Biphasic'],vals,color=[COLORS[m] for m in CASES]);ax.bar_label(bars,labels=['N/A' if (title=='Fines / detached (%)' and cases[m]['objective']['fine_debris_mass_fraction_detached'] is None) or (title=='Unchanged target member sets\nsince drive-off' and not cases[m]['objective']['drive_off_time_reached']) else f'{v:.4g}' for m,v in zip(CASES,vals)]);ax.set(title=title,ylim=(0,max(1,max(vals)*1.2)));ax.grid(axis='x',visible=False)
 axs[2,2].axis('off');success=cases['biphasic']['objective']['full_target'];axs[2,2].text(0,.95,('INCOMPLETE_PROTOCOL' if not all(cases[m]['native']['status']=='COMPLETED' for m in CASES) else 'TARGET_OBSERVED' if success else 'TARGET_NOT_ACHIEVED')+'\n\nSingle frozen spatial realization\nUncalibrated material thresholds\nPrescribed mechanical shear proxy\nNo validated bubble power\nGraph separation ≠ certified\ncontinuum crack surface',va='top',fontsize=13,linespacing=1.6)
 fig.suptitle('Scientific outcome from actual native trajectories · no grouping or prescribed fragment shapes\n'+' | '.join(f'{m}: {cases[m]["native"]["status"]}, {cases[m]["native"]["time_s"]*1e3:.5f} ms' for m in CASES),fontsize=16);save(fig,'FIG_10_SCIENTIFIC_CONCLUSION')
 fig,axs=plt.subplots(1,3,figsize=(18,6),layout='constrained')
 for m in CASES:
  q=data['spatial'][m];g=q['phase_damage'];axs[0].plot([x['kind'] for x in g],[100*x['failed_fraction'] for x in g],'o-',label=NAMES[m]);r=q['domain_retention'];axs[1].plot([x['strong_domain_id'] for x in r],[x['largest_internal_connected_nodes'] for x in r],'o-',label=NAMES[m]);axs[2].plot([x['strong_domain_id'] for x in r],[x['remaining_cross_phase_bonds'] for x in r],'o-',label=NAMES[m])
 axs[0].set(title='Natural failure / eligible bonds',ylabel='%');axs[0].tick_params(axis='x',rotation=25);axs[1].set(title='Strong-region internal retention',xlabel='Domain ID',ylabel='Largest internally connected node count');axs[2].set(title='Remaining live mixed-phase connections',xlabel='Domain ID',ylabel='Bonds')
 for ax in axs:ax.legend()
 fig.suptitle('Phase and nonlocal bridges · U uses H reference labels only; region labels are not fragments',fontsize=16);save(fig,'FIG_11_PHASE_AND_BRIDGES')
 surface_targets=[]
 for m in CASES:
  key='final' if cases[m]['native']['status']=='COMPLETED' else 'last_accepted'
  surface_targets.extend((m,r) for r in cases[m]['surface']['target_component_surface_support'][key])
 if surface_targets:
  fig,axs=plt.subplots(2,2,figsize=(17,10),layout='constrained');xx=np.arange(len(surface_targets));labels=[('U' if m=='uniform' else 'H')+f" {r['component_id']}" for m,r in surface_targets];colors=[COLORS[m] for m,r in surface_targets]
  axs[0,0].bar(xx,[r['nodes'] for m,r in surface_targets],color=colors);axs[0,0].set(title='Target component size',ylabel='Original material nodes')
  bond_lookup={(r['case'],r['component_id']):len(r['bonds']) for r in target_details}
  axs[0,1].bar(xx,[bond_lookup[(m,r['component_id'])] for m,r in surface_targets],color=colors);axs[0,1].set(title='Actual internal surviving material connections',ylabel='All original live bonds')
  axs[1,0].bar(xx-.18,[r['supported_current_tetra'] for m,r in surface_targets],.36,label='Supported current tetra',color='#7793a5');axs[1,0].bar(xx+.18,[r['active_supported_faces'] for m,r in surface_targets],.36,label='Active load facets',color='#b17785');axs[1,0].legend();axs[1,0].set(title='Inherited surface representation',ylabel='Count')
  axs[1,1].bar(xx,[100*r['node_fraction_on_active_faces'] for m,r in surface_targets],color=colors);axs[1,1].set(title='Nodes on active facets / component nodes',ylabel='% (not continuum area completeness)',ylim=(0,105))
  for ax in axs.flat:ax.set_xticks(xx,labels,rotation=45,ha='right');ax.set_xlabel('Case and exact member-set component ID');ax.grid(axis='x',visible=False)
  fig.suptitle('Last accepted target blocks: cohesive graph vs represented surface loading\nZero active facets require support/accessibility diagnosis; zero modeled force is not evidence of physical zero drag',fontsize=16);save(fig,'FIG_12_TARGET_SURFACE_REPRESENTATION')
 dump(O/'VISUALIZATION_GEOMETRY.json',dict(common_global_bounds_um=np.array(limits),camera=dict(elevation_deg=16,azimuth_deg=-64),displacement_scale=1,connectivity='all original live bonds',display_bond_cutoff_m=1.01*DX,free_internal_bonds='all original surviving internal bonds',node_damage='mean incident d, initial precut excluded'))

def target_morphology_figure(ref,states,limits):
 b=ref['bonds']
 fig=plt.figure(figsize=(16,13));fig.subplots_adjust(wspace=.2,hspace=.32,bottom=.08);morphology=[]
 morphology_widths=[]
 for mode in CASES:
  state=states[mode]['max_multifragment'];ll,ss,aa=graph.topology(ref,state['d']);ee=np.flatnonzero(~aa&(ss>=20)&(ss<=200))
  if len(ee):
   selected=np.flatnonzero(ll==ee[np.argmax(ss[ee])]);morphology_widths.extend(float(np.ptp(pos[selected]*1e6,axis=0).max()) for pos in [ref['X'],state['x']])
 morphology_width=max(300.,1.3*max(morphology_widths,default=0.))
 for row,m in enumerate(CASES):
  s=states[m]['max_multifragment'];lab,sz,att=graph.topology(ref,s['d']);eligible=np.flatnonzero(~att&(sz>=20)&(sz<=200))
  if not len(eligible):
   for col in range(2):
    ax=ax3(fig,(2,2,1+2*row+col),NAMES[m]+' · target morphology',limits);ax.text2D(.1,.5,'TARGET_FRAGMENT_NOT_OBSERVED',transform=ax.transAxes,fontsize=12)
   continue
  label=eligible[np.argmax(sz[eligible])];ids=np.flatnonzero(lab==label);inside=np.isin(b,ids).all(1)&(s['d']<1)
  reference_shape=shape(ref['X'][ids]);actual_shape=shape(s['x'][ids])
  morphology.append(dict(case=m,time_s=float(s['time']),original_atom_ids=ids+1,internal_alive_bonds=int(inside.sum()),reference_shape=reference_shape,actual_shape=actual_shape))
  for col,(positions,description,sh) in enumerate([(ref['X'],'Reference coordinates of the same IDs',reference_shape),(s['x'],'Actual accepted coordinates',actual_shape)]):
   center=(positions[ids].min(0)+positions[ids].max(0))*.5e6
   ax=ax3(fig,(2,2,1+2*row+col),NAMES[m]+f' · {len(ids)} nodes\n'+description,(center-morphology_width/2,center+morphology_width/2))
   ax.scatter(*positions[ids].T*1e6,s=30,c=COLORS[m],depthshade=False);ax.add_collection3d(Line3DCollection(positions[b[inside]]*1e6,colors='#244e63',linewidths=.85,alpha=.65))
   span=' × '.join(f'{v:.1f}' for v in sh['actual_span_m']*1e6)
   display_shape=sh['shape_class'].replace('THREE_DIMENSIONAL_POINT_CLOUD','3D point cloud').replace('_',' ').lower()
   ax.text2D(.02,.98,f't = {float(s["time"])*1e3:.5f} ms; {inside.sum()} live internal bonds\nXYZ span: {span} µm\nPoint-cloud class: {display_shape}',transform=ax.transAxes,fontsize=10,va='top')
 fig.suptitle('Target-count membership does not certify compact 3D block morphology\nSame original IDs and surviving edges · actual coordinates · equal local scale; views centered on each target',fontsize=16);save(fig,'FIG_07C_REFERENCE_AND_CURRENT_TARGET_SHAPE');dump(O/'TARGET_MORPHOLOGY_FIGURE_EVIDENCE.json',morphology)

def main():
 start=time.perf_counter();ref=load_npz(O/'provenance/reference.npz');field=load_npz(O/'BIPHASIC_FIELD.npz');parameter_figures(ref,field);result_figures(ref,field)
 if all((O/'cases'/m/'keys/rise_end.npz').exists() for m in CASES):
  from preview_rise_end import main as rise_end_figure
  rise_end_figure()
 dump(O/'FIGURE_RESOURCES.json',dict(wall_s=time.perf_counter()-start,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024));print('FIGURES_DONE',flush=True)
if __name__=='__main__':main()
