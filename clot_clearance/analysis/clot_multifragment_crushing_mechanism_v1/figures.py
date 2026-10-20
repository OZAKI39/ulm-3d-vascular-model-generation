from common import *
from indenter import *
from analyze import read,histories,graph
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
plt.rcParams.update({'font.size':13,'axes.titlesize':15,'axes.labelsize':14,'legend.fontsize':11,'figure.facecolor':'white','savefig.facecolor':'white','svg.fonttype':'none','axes.grid':True,'grid.color':'#8b949e','grid.alpha':.55})
ref=load_npz(O/'provenance/reference.npz');field=load_npz(O/'HETEROGENEITY_FIELD.npz');X=ref['X'];b=ref['bonds'];n=int(ref['clot_count']);tool=Indenter(ref);names=['uniform','heterogeneous'];titles=['UNIFORM','CORRELATED HETEROGENEITY'];colors=['#326a9c','#bb5a38'];cases={k:read(O/'cases'/k/'ANALYSIS.json') for k in names};summaries={k:read(O/'cases'/k/'SUMMARY.json') for k in names};rows={k:histories(O/'cases'/k) for k in names};figures=[]
ROI=np.array([X[:n].min(0)-[.0002,.0002,.00012],X[:n].max(0)+[.0002,.0002,.00045]])
def save(fig,name,caption):
 fig.text(.99,.006,'KINEMATIC_IMPACT_SURROGATE · EXPLORATORY_UNCALIBRATED',ha='right',fontsize=9,color='#68737d');fig.savefig(F/(name+'.png'),dpi=220,bbox_inches='tight');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig);figures.append((name,caption))
def axis3(ax):
 ax.view_init(elev=16,azim=-63);ax.set(xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)',xlim=ROI[:,0]*1e6,ylim=ROI[:,1]*1e6,zlim=ROI[:,2]*1e6);ax.set_box_aspect(np.ptp(ROI,axis=0));ax.tick_params(labelsize=10,pad=1)
 for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.set_facecolor((.95,.96,.98,.3));axis._axinfo['grid']['color']=(.35,.4,.46,.45)
def sphere(ax,c):
 u=np.linspace(0,2*np.pi,17);v=np.linspace(0,np.pi,11);xx=c[0]+tool.radius*np.outer(np.cos(u),np.sin(v));yy=c[1]+tool.radius*np.outer(np.sin(u),np.sin(v));zz=c[2]+tool.radius*np.outer(np.ones_like(u),np.cos(v));ax.plot_wireframe(xx*1e6,yy*1e6,zz*1e6,color='#515b68',linewidth=.5,alpha=.35)
def state3(ax,s,title,show_indent=True):
 x=s['x'];labels,sizes,attached=graph.topology(ref,s['d']);free=~attached[labels];soft=np.zeros(n,bool);pairs=b[(s['d']>0)&(s['d']<1)&~ref['initial_cut']];ids=pairs[pairs<n];soft[ids]=True;inside=np.all((x[:n]>=ROI[0])&(x[:n]<=ROI[1]),axis=1);col=np.tile(np.array(matplotlib.colors.to_rgba('#b55142')),(n,1));col[soft]=matplotlib.colors.to_rgba('#d4a732');col[free]=matplotlib.colors.to_rgba('#478ba8')
 ax.scatter(*(x[:n][inside]*1e6).T,s=6,c=col[inside],depthshade=False);ax.scatter(*(x[n:]*1e6).T,s=5,color='#8a8c91',alpha=.6,depthshade=False)
 live=np.flatnonzero((s['d']<1)&(ref['bond_types']==1));live=live[::max(1,len(live)//600)];e=x[b[live]];sel=np.all((e>=ROI[0])&(e<=ROI[1]),axis=(1,2));ax.add_collection3d(Line3DCollection(e[sel]*1e6,colors='#404751',linewidths=.25,alpha=.25))
 if show_indent:sphere(ax,tool.center(float(s['time'])))
 axis3(ax);ax.set_title(f'{title}\nt={float(s["time"])*1e3:.4f} ms',pad=1);ax.text2D(.02,.93,f'Outside ROI: {np.sum(~inside)} material nodes',transform=ax.transAxes,fontsize=10,bbox=dict(facecolor='white',edgecolor='none',alpha=.8))
def series(ax,key,scale=1,label=None):
 for k,c,t in zip(names,colors,titles):ax.plot([r['time_s']*1e3 for r in rows[k]],[r[key]*scale for r in rows[k]],color=c,label=t,lw=1.7)
 ax.axvline(CLEAR*1e3,color='#777777',linestyle='--',lw=.8);ax.set_xlabel('Time (ms)');ax.set_ylabel(label or key)
def main():
 F.mkdir(parents=True,exist_ok=True)
 fig=plt.figure(figsize=(12,9));ax=fig.add_subplot(projection='3d');initial=load_npz(O/'cases/uniform/keys/initial.npz');state3(ax,initial,'Original specimen and frozen spherical indenter');cut=np.flatnonzero(ref['initial_cut'])[::24];ax.add_collection3d(Line3DCollection(X[b[cut]]*1e6,colors='black',linewidths=.75,linestyles='--'));c=tool.c0*1e6;ax.quiver(*c,0,0,-120,color='#2a6592',arrow_length_ratio=.2);angle=np.linspace(0,2*np.pi,121);ax.plot((tool.c0[0]+tool.radius*np.cos(angle))*1e6,(tool.c0[1]+tool.radius*np.sin(angle))*1e6,np.full_like(angle,X[:n,2].max())*1e6,color='#2a6592',ls='--',lw=1.5);ax.legend(handles=[Line2D([0],[0],marker='o',color='w',markerfacecolor='#b55142',label='1792 material points'),Line2D([0],[0],marker='o',color='w',markerfacecolor='#8a8c91',label='336 fixed anchors'),Line2D([0],[0],color='black',linestyle='--',label='Original pre-cut bonds (sample)'),Line2D([0],[0],color='#515b68',label='Rigid sphere, R = 182.8 µm'),Line2D([0],[0],color='#2a6592',linestyle='--',label='Reference projected footprint')],loc='upper left',bbox_to_anchor=(0,.95));save(fig,'Fig01_indenter_and_specimen','Original 1792 material points, 336 anchors and original pre-cut bonds. Sphere represents the actual native contact geometry; no new cuts or continuous clot shell. Blue circle is the reference projected footprint, not measured instantaneous contact area. Arrow indicates the inward trajectory, not a force scale.')
 fig,axs=plt.subplots(2,3,figsize=(18,9));spec=[('indent_displacement_m',1e6,'Displacement (µm)'),('indent_speed_m_s',1,'Velocity (m/s)'),('reaction_z_N',1e6,'Tool reaction Z (µN)'),('contact_nodes',1,'Current contact nodes'),('max_overlap_m',1e6,'Maximum overlap (µm)'),('contact_duration_s',1e3,'Accumulated contact duration (ms)')]
 for ax,(key,scale,label) in zip(axs.flat,spec):series(ax,key,scale,label)
 axs[0,1].legend();fig.tight_layout();save(fig,'Fig02_dynamic_contact','Actual imposed motion and native contact response. Dashed line marks scheduled withdrawal completion; any later contact is physical recontact, not persistent reference-node traction. Contact area is represented by node count; it is not a measured continuous area.')
 fig=plt.figure(figsize=(17,8));norm=Normalize(.135,.165);subset=np.flatnonzero(ref['bond_types']==1)[::32]
 for j,k in enumerate(names):
  ax=fig.add_subplot(1,2,j+1,projection='3d');values=field['s0_by_original_bond_id'] if k=='heterogeneous' else S0[ref['bond_types']];lc=Line3DCollection(X[b[subset]]*1e6,cmap='coolwarm',norm=norm,linewidths=.6,alpha=.7);lc.set_array(values[subset]);ax.add_collection3d(lc);ax.scatter(*(X[n:]*1e6).T,s=3,color='#888888');axis3(ax);ax.set_title(titles[j]+' · original bulk s₀');e=summaries[k]['events'].get('first_softening')
  if e:pts=np.array(e['reference_midpoints_m'])*1e6;ax.scatter(*pts.T,marker='*',s=90,c='black',label='First softening');ax.legend(loc='upper left')
  births=read(O/'cases'/k/'CRACK_NETWORK_AUDIT.json')['births']
  if births:
   points=np.array([q['reference_centroid_m'] for q in births])*1e6;ax.scatter(*points.T,marker='x',s=28,c='#343434',label='Failed-network births');ax.legend(loc='upper left',fontsize=10)
 fig.colorbar(lc,ax=fig.axes,shrink=.55,pad=.01,label='Local bulk s₀');fig.suptitle('One frozen field · α = 0.10 bound · nominal lc = 250 µm · sc = 0.18');save(fig,'Fig03_original_bond_onset_field','Same original bulk bond sample, coordinate frame and camera in both panels. Field acts on original bond ID; interface unchanged. Black stars show actual first-softening reference locations. lc is an uncalibrated correlation hypothesis, not a measured clot microstructure.')
 keys=[('initial.npz','Initial'),('first_softening.npz','First softening'),('first_failure.npz','First full failure'),('first_multinode.npz','First multinode detachment'),('max_multifragment.npz','Maximum multifragment score'),('final.npz','End after withdrawal')]
 for k,title in zip(names,titles):
  fig=plt.figure(figsize=(19,14))
  for j,(key,label) in enumerate(keys):
   ax=fig.add_subplot(2,3,j+1,projection='3d');p=O/'cases'/k/'keys'/key
   if p.exists():
    state=load_npz(p);state3(ax,state,label)
    event=summaries[k]['events'].get(key.removesuffix('.npz'))
    if key in ['first_softening.npz','first_failure.npz'] and event:
     ids=np.array(event['original_bond_ids'])-1;points=state['x'][b[ids]].mean(1)*1e6;ax.scatter(*points.T,marker='*',s=80,c='#151515',depthshade=False)
     ax.add_collection3d(Line3DCollection(state['x'][b[ids]]*1e6,colors='#9d206d',linewidths=2))
   else:ax.set_axis_off();ax.text2D(.2,.5,label+'\nNOT OBSERVED',transform=ax.transAxes,fontsize=17)
  fig.suptitle(title+' · accepted coordinates ×1 · common physical ROI and camera');fig.subplots_adjust(wspace=.05,hspace=.3);save(fig,'Fig04_'+k+'_evolution','Actual accepted states at displacement ×1, with fixed physical ROI/camera. Missing events are explicitly absent. Maximum score prioritizes the number and mass of 20–200-node components; it does not invent fragments. Black stars and magenta original bond segments identify the first softening/failure events. Red attached, gold incident to partially softened bonds, blue detached, gray anchors. Out-of-ROI counts are printed; full final extent follows.')
 fig,axs=plt.subplots(1,2,figsize=(15,7))
 for ax,k,t in zip(axs,names,titles):
  s=load_npz(O/'cases'/k/'keys/final.npz');lab,_,at=graph.topology(ref,s['d']);free=~at[lab];ax.scatter(s['x'][:n,0]*1e3,s['x'][:n,2]*1e3,s=5,c=np.where(free,'#478ba8','#b55142'));ax.scatter(s['x'][n:,0]*1e3,s['x'][n:,2]*1e3,s=3,c='#888888');ax.set(xlabel='X (mm)',ylabel='Z (mm)',title=t+' · all final material points');ax.set_aspect('equal',adjustable='box')
 fig.tight_layout();save(fig,'Fig04b_full_final_extent','All final material and anchor coordinates in X–Z projection; no particles are clipped by the main 3D camera ROI.')
 fig,axs=plt.subplots(2,2,figsize=(15,10));bins=[(1,2),(3,13),(14,33),(34,113),(114,100000)];labels=['1–2','3–13','14–33','34–113','≥114']
 for i,(k,t,c) in enumerate(zip(names,titles,colors)):
  sizes=np.array([q['nodes'] for q in cases[k]['final_components']]);counts=[int(np.sum((sizes>=lo)&(sizes<=hi))) for lo,hi in bins];mass=[sizes[(sizes>=lo)&(sizes<=hi)].sum()*RHO*V*1e9 for lo,hi in bins];axs[0,0].bar(np.arange(5)+(i-.5)*.35,counts,width=.35,color=c,label=t);axs[0,1].bar(np.arange(5)+(i-.5)*.35,mass,width=.35,color=c,label=t)
 for ax in axs[0]:ax.set_xticks(range(5),labels);ax.set_xlabel('Nodes per actual free component')
 axs[0,0].set_ylabel('Final component count');axs[0,0].legend();axs[0,1].set_ylabel('Final detached mass (µg)');series(axs[1,0],'N20_200',1,'20–200-node component count');series(axs[1,1],'mass_fraction20_200',100,'Initial material mass in 20–200 nodes (%)');axs[1,0].axhline(5,color='black',linestyle=':',label='Exploratory count target');fig.tight_layout();save(fig,'Fig05_fragment_size_and_mass','Counts use actual live-bond connected components and anchor reachability. Contact pairs never connect fragments. A 20-node component is not a resolved 200 µm continuum fragment. Zero bars mean no such final components.')
 fig,axs=plt.subplots(2,2,figsize=(16,10))
 for col,k in enumerate(names):
  reg={int(cid):v for cid,v in read(O/'cases'/k/'COMPONENT_REGISTRY.json').items()};trace=read(O/'cases'/k/'FRAGMENT_TRACE.json');ax=axs[0,col]
  for cid,q in reg.items():
   if q['nodes']>=2:ax.scatter(q['birth_time_s']*1e3,q['nodes'],s=22,color=colors[col])
  for e in summaries[k]['secondary_events']:
   parent=reg[e['parent_component_id']]
   if parent['nodes']<2:continue
   for child in e['children_component_ids']:ax.plot([parent['birth_time_s']*1e3,e['time_s']*1e3],[parent['nodes'],reg[child]['nodes']],color='#777777',lw=.6,alpha=.6)
  ax.set(title=titles[col],xlabel='Birth / split time (ms)',ylabel='Original nodes in component',yscale='symlog',xlim=(0,END*1e3))
  observed_max=max([q['nodes'] for q in reg.values()]+[0])
  if observed_max<=10:ax.set_yscale('linear');ax.set_ylim(0,max(3,observed_max+1));ax.set_yticks(range(max(3,observed_max+1)+1))
  ax.axvline(CLEAR*1e3,color='#555555',ls='--');selected=sorted(trace,key=lambda cid:reg[int(cid)]['nodes'],reverse=True)[:5]
  if not selected:ax.text(.1,.5,'No ≥2-node free component observed',transform=ax.transAxes)
  for cid in selected:
   ts=trace[cid];axs[1,col].plot([q['time_s']*1e3 for q in ts],[q['internal_survival'] for q in ts],label=f'ID {cid}, n={reg[int(cid)]["nodes"]}')
  if selected:
   axs[1,col].legend()
   axs[1,col].text(.03,.08,'Live external material connections = 0\n2-node pairs are unresolved debris',transform=axs[1,col].transAxes,fontsize=11)
  else:axs[1,col].text(.1,.5,'No ≥2-node free component observed',transform=axs[1,col].transAxes)
  axs[1,col].set(xlabel='Time (ms)',ylabel='Internal original bond survival',ylim=(-.03,1.03),xlim=(0,END*1e3));axs[1,col].axvline(CLEAR*1e3,color='#555555',ls='--')
 fig.tight_layout();save(fig,'Fig06_fragment_lineage_and_cohesion','Births and split links derive from exact original atom membership. Bottom panels show up to five largest observed ≥2-node components; two-node pairs remain unresolved debris, not cohesive continuum blocks. Free components have zero surviving external material connections, independently checked against original graph boundaries. Dashed line marks tool withdrawal completion; stability is limited to the actual observation window.')
 fig,axs=plt.subplots(2,2,figsize=(16,10));terms=[('W_J','External W'),('U_J','Stored U'),('K_J','Kinetic K'),('D_J','Damage D'),('Uc_J','Particle contact Uc'),('U_indenter_J','Indenter Ui'),('R_J','Residual R')]
 for col,k in enumerate(names):
  for key,label in terms:axs[0,col].plot([q['time_s']*1e3 for q in rows[k]],[q[key]*1e9 for q in rows[k]],label=label)
  axs[0,col].set(title=titles[col],xlabel='Time (ms)',ylabel='Energy (nJ)');axs[0,col].legend(ncol=2)
  for key,label in [('R_J','Signed R'),('Jc_J','Pair activation offset Jc')]:axs[1,col].plot([q['time_s']*1e3 for q in rows[k]],[q[key]*1e9 for q in rows[k]],label=label)
  axs[1,col].set(xlabel='Time (ms)',ylabel='Accounting terms (nJ)');axs[1,col].legend()
 fig.tight_layout();save(fig,'Fig07_native_energy_ledger','Native material and contact energies, independently integrated moving-tool work, exact damage energy drop and signed numerical residual. No numerical damping. Pair-contact activation offset is separate; residual is not material dissipation. Equal trajectories do not imply equal work.')
 fig,axs=plt.subplots(3,3,figsize=(18,14));metrics=[('crack:n_births','Discrete damage-network births'),('crack:n_merges','Discrete damage-network mergers'),('history_max_N20_200','History max 20–200-node count'),('final:N20_200','Final 20–200-node count'),('final:largest_fraction','Final largest / initial nodes (%)'),('final:debris12_detached_fraction','Final debris / detached mass (%)'),('withdrawal_medium_cohort_surviving','Medium cohort connected through end'),('secondary_after_withdrawal','Post-withdrawal secondary splits'),('final:W_J','Actual final input work (nJ)')]
 for ax,(key,label) in zip(axs.flat,metrics):
  vals=[]
  for k in names:
   v=cases[k]['crack_networks'][key.split(':')[1]] if key.startswith('crack:') else cases[k]['final'][key.split(':')[1]] if ':' in key else cases[k][key];vals.append(v)
  missing=[v is None for v in vals];vals=[0 if v is None else v for v in vals]
  scale=1e6 if key=='peak_reaction_N' else 100 if 'fraction' in key else 1e9 if key.endswith('W_J') else 1;vals=np.asarray(vals)*scale;ax.bar(['Uniform','Correlated'],vals,color=colors);ax.set_title(label)
  for i,v in enumerate(vals):ax.text(i,v,'UNDEFINED' if missing[i] else f'{v:.4g}',ha='center',va='bottom')
  ax.set_ylim(0,max(1,float(max(vals))*1.2))
 fig.suptitle('Same kinematic loading · one frozen heterogeneity realization · finite observation');fig.tight_layout(rect=(0,.035,1,1));fig.text(.5,.025,'Final detached / initial material: '+ ' · '.join(f"{t} {cases[k]['final']['detached_nodes']/n*100:.3f}%" for k,t in zip(names,['uniform','correlated'])),ha='center');save(fig,'Fig08_comparison','Exploratory comparison of one uniform and one fixed correlated-field case. Undefined debris ratios are labeled UNDEFINED; the denominator is detached mass, not the entire specimen. No claim of clot microstructure calibration, mesh convergence, ultrasound feasibility or V3 transfer.')
 base='../../../verification/'+TASK+'/'+RUN+'/'
 html='<!doctype html><meta charset="utf-8"><title>Clot multifragment crushing</title><style>body{font:17px sans-serif;max-width:1600px;margin:30px auto;background:#f6f7f8;color:#243447}section{background:white;padding:24px;margin:24px 0}img{width:100%}p{line-height:1.6}</style><h1>Local dynamic indentation and correlated material heterogeneity</h1><p>KINEMATIC_IMPACT_SURROGATE · EXPLORATORY_UNCALIBRATED · actual native accepted states</p>'
 html+='<p>'+ ' · '.join(f'<a href="{base+file}">{label}</a>' for file,label in [('REVIEW_CLOT_MULTIFRAGMENT_CRUSHING_MECHANISM_V1.md','中文科研报告'),('MULTIFRAGMENT_SIZE_DISTRIBUTION.csv','Fragment CSV'),('FINAL_STATUS.json','Final status')])+'</p>'
 for name,caption in figures:html+=f'<section><h2>{name}</h2><a href="{name}.svg"><img src="{name}.png"></a><p>{caption}</p></section>'
 (F/'OPEN_RESULTS.html').write_text(html);(O/'OPEN_RESULTS.html').write_text(f'<!doctype html><meta http-equiv="refresh" content="0;url=../../../results/{TASK}/{RUN}/OPEN_RESULTS.html">');dump(F/'FIGURE_MANIFEST.json',dict(figures=[dict(name=n,caption=c) for n,c in figures],coordinates='actual accepted FP64, displacement x1',camera=dict(elev=16,azim=-63,physical_ROI_m=ROI)));print('FIGURES_DONE',len(figures))
if __name__=='__main__':main()
