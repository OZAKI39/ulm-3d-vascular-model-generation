"""Scientific figures from actual accepted coordinates and original topology."""
from fragment_topology_analysis import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
plt.rcParams.update({'font.size':15,'axes.titlesize':17,'axes.labelsize':16,'legend.fontsize':12,'figure.facecolor':'white','savefig.facecolor':'white','svg.fonttype':'none','axes.grid':True,'grid.color':'#65717f','grid.alpha':.6})
ref=load_npz(O/'provenance/reference.npz');X=ref['X'];b=ref['bonds'];n=int(ref['clot_count']);names=['simultaneous','staggered'];titles=['SIMULTANEOUS','STAGGERED'];colors=['#2a679b','#bc5738'];results=read(O/'ANALYSIS.json');summaries={k:read(O/'cases'/k/'SUMMARY.json') for k in names};rows={k:histories(O/'cases'/k) for k in names};tools={k:DistributedLoading(ref,k) for k in names};figures=[]
ROI=np.array([[-800,-450,-650],[800,450,650]])*1e-6
registry={k:{int(cid):v for cid,v in read(O/'cases'/k/'COMPONENT_REGISTRY.json').items()} for k in names}
member_ids={k:{tuple(v['atom_ids']):cid for cid,v in registry[k].items()} for k in names}
def save(fig,name,caption):
    fig.text(.99,.004,'MULTISITE_MECHANICAL_SURROGATE · exploratory discrete specimen',ha='right',fontsize=10,color='#53606c')
    fig.savefig(F/(name+'.png'),dpi=230,bbox_inches='tight');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig);figures.append(dict(name=name,caption=caption))
def axis3(ax,roi=ROI):
    ax.view_init(elev=16,azim=-63);ax.set(xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)',xlim=roi[:,0]*1e6,ylim=roi[:,1]*1e6,zlim=roi[:,2]*1e6);ax.set_box_aspect(np.ptp(roi,axis=0));ax.tick_params(labelsize=11,pad=1)
    for index,axis in enumerate([ax.xaxis,ax.yaxis,ax.zaxis]):
        ticks=matplotlib.ticker.MaxNLocator(nbins=3 if index==1 else 4).tick_values(*(roi[:,index]*1e6))
        axis.set_ticks(ticks[(ticks>=roi[0,index]*1e6)&(ticks<=roi[1,index]*1e6)])
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.set_facecolor((.95,.96,.98,.35));axis._axinfo['grid']['color']=(.3,.36,.43,.55)
def sphere(ax,site,t,color):
    c=site.center(t);u=np.linspace(0,2*np.pi,17);v=np.linspace(0,np.pi,11)
    xx=c[0]+site.radius*np.outer(np.cos(u),np.sin(v));yy=c[1]+site.radius*np.outer(np.sin(u),np.sin(v));zz=c[2]+site.radius*np.outer(np.ones_like(u),np.cos(v))
    ax.plot_wireframe(xx*1e6,yy*1e6,zz*1e6,color=color,linewidth=.55,alpha=.45)
def point_colors(name,s):
    lab,sizes,att=graph.topology(ref,s['d']);c=np.tile(matplotlib.colors.to_rgba('#adb3ba'),(n,1));ids=np.zeros(n,int)
    for label in np.flatnonzero(~att):
        nodes=np.flatnonzero(lab==label);cid=member_ids[name][tuple((nodes+1).tolist())];c[nodes]=plt.cm.turbo(cid/max(len(registry[name]),1));ids[nodes]=cid
    return c,ids,lab,sizes,att
def state3(ax,name,s,title,show_tools=True):
    x=s['x'];c,ids,lab,sizes,att=point_colors(name,s);inside=np.all((x[:n]>=ROI[0])&(x[:n]<=ROI[1]),axis=1)
    ax.scatter(*(x[:n][inside]*1e6).T,s=8,c=c[inside],depthshade=False);ax.scatter(*(x[n:]*1e6).T,s=5,c='#5e6670',depthshade=False)
    live=np.flatnonzero((s['d']<1)&(ref['bond_types']==1));sample=live[::max(1,int(np.ceil(len(live)/500)))];edges=x[b[sample]];visible=np.all((edges>=ROI[0])&(edges<=ROI[1]),axis=(1,2));ax.add_collection3d(Line3DCollection(edges[visible]*1e6,colors='#606b77',linewidths=.3,alpha=.3))
    # All surviving internal bonds of >=14-node free components, no shell.
    for label in np.flatnonzero(~att & (sizes>=14)):
        nodes=np.flatnonzero(lab==label);_,internal=graph.boundary(ref,nodes);internal=internal[s['d'][internal]<1]
        ax.add_collection3d(Line3DCollection(x[b[internal]]*1e6,colors=[c[nodes[0]]],linewidths=.7,alpha=.65))
        center=x[nodes].mean(0)
        if np.all(center>=ROI[0]) and np.all(center<=ROI[1]):ax.text(*(center*1e6),f'ID {ids[nodes[0]]}',fontsize=11)
    if show_tools:
        for site,color in zip(tools[name].sites,colors):sphere(ax,site,float(s['time']),color)
    axis3(ax);ax.set_title(f'{title}\n{float(s["time"])*1e3:.4f} ms',pad=5)
    ax.text2D(.01,.92,f'Free {int(sizes[~att].sum())} nodes / {sum(~att)} components\nOutside ROI: {sum(~inside)}',transform=ax.transAxes,fontsize=11,bbox=dict(facecolor='white',edgecolor='none',alpha=.85))
def series(ax,key,scale=1,ylabel=None):
    for name,color,title in zip(names,colors,titles):ax.plot([q['time_s']*1e3 for q in rows[name]],[q[key]*scale if q[key] is not None else np.nan for q in rows[name]],color=color,label=title,lw=1.6)
    ax.axvline(2.8,color='#606060',ls='--',lw=1);ax.set(xlabel='Time (ms)',ylabel=ylabel or key)
    if key.startswith('N'):
        ymax=max(q[key]*scale for name in names for q in rows[name]);ax.set_ylim(0,max(1,ymax*1.15));ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True,nbins=4))
    elif 'fraction' in key:
        values=[q[key]*scale for name in names for q in rows[name] if q[key] is not None];ax.set_ylim(0,max(1,max(values,default=0)*1.08))
def main():
    F.mkdir(parents=True,exist_ok=True)
    fig=plt.figure(figsize=(18,9))
    for j,(key,title) in enumerate([('initial','External sites and paths'),('site0_peak','Actual contact at maximum indentation')]):
        ax=fig.add_subplot(1,2,j+1,projection='3d');s=load_npz(O/'cases/simultaneous/keys'/(key+'.npz'));state3(ax,'simultaneous',s,title)
        if j==0:
            cut=np.flatnonzero(ref['initial_cut'])[::18];ax.add_collection3d(Line3DCollection(X[b[cut]]*1e6,colors='#24262b',linewidths=.6,linestyles='--'))
            for site,color in zip(tools['simultaneous'].sites,colors):
                c=site.c0*1e6;v=site.direction*200;ax.quiver(*c,*v,color=color,arrow_length_ratio=.15,linewidth=2)
        else:
            for source,color in enumerate(colors):
                forces=s['site_forces'][source];sel=np.flatnonzero(np.linalg.norm(forces,axis=1)>0);sel=sel[::max(1,int(np.ceil(len(sel)/18)))];vec=forces[sel]/2e-6*100
                ax.quiver(*(s['x'][sel]*1e6).T,*vec.T,color=color,arrow_length_ratio=.22,linewidth=1.1)
            ax.text2D(.02,.03,'Force scale: 2 µN = 100 µm arrow',transform=ax.transAxes,fontsize=11)
    fig.suptitle('Two externally accessible, offset, frictionless normal contacts',y=.99);fig.subplots_adjust(wspace=.03)
    save(fig,'Fig01_sites_and_actual_contact','Left: original specimen, anchors, sampled original precuts and imposed path arrows. Right: actual accepted state and per-node contact forces (at most 18 arrows per source, physical scale shown). No tangential friction law; offset normal loading may produce mixed deformation. Attached points gray; free points colored by exact component ID.')
    oldrows=histories(HISTORY/'cases/uniform');fig,axs=plt.subplots(2,2,figsize=(16,10))
    axs[0,0].plot([q['time_s']*1e3 for q in oldrows],[q['indent_displacement_m']*1e6 for q in oldrows],color='#555555',ls=':',lw=2,label='Historical single: inward Z')
    for name,color in zip(names,colors):
        for j,ls in enumerate(['-','--']):axs[0,0].plot([q['time_s']*1e3 for q in rows[name]],[q[f'site{j}_displacement_m']*1e6 for q in rows[name]],color=color,ls=ls,label=f'{name}, source {j}')
    axs[0,0].set(xlabel='Time (ms)',ylabel='Inward travel incl. gap (µm)',title='Prescribed timing; different historical geometry');axs[0,0].legend(fontsize=11)
    for j,name in enumerate(names):
        for source,color in enumerate(colors):axs[0,1].plot([q['time_s']*1e3 for q in rows[name]],[q[f'site{source}_contact_nodes'] for q in rows[name]],color=color,ls=['-','--'][j],label=f'{name}, source {source}')
    axs[0,1].set(xlabel='Time (ms)',ylabel='Actual contact nodes');axs[0,1].legend(fontsize=11)
    series(axs[1,0],'N20_200',ylabel='True 20–200-node components');series(axs[1,1],'W_J',1e9,'Actual total external work (nJ)');axs[1,0].legend();fig.tight_layout()
    save(fig,'Fig02_historical_and_multisite_timing','A/B share sites and trajectories; B source 1 starts 0.4 ms later. The dotted historical single-source path is context, with a different contact position/direction. Equal travel does not imply equal force or work. Dashed vertical line marks the common 2.8 ms cohort, not guaranteed last contact.')
    historical=module('historical_geometry',HA/'indenter.py');historical.O=HISTORY;oldtool=historical.Indenter(ref);fig=plt.figure(figsize=(21,8));geometry_roi=ROI.copy();geometry_roi[1,2]=850e-6
    for j,label in enumerate(['Historical single','Simultaneous sites','Staggered sites']):
        ax=fig.add_subplot(1,3,j+1,projection='3d');ax.scatter(*(X[:n]*1e6).T,s=6,c='#adb3ba');ax.scatter(*(X[n:]*1e6).T,s=5,c='#5e6670')
        actors=[oldtool] if j==0 else tools[names[j-1]].sites
        for k,site in enumerate(actors):
            sphere(ax,site,0,colors[k]);start=site.center(0)*1e6;end=site.center(.0012 if j<2 or k==0 else .0016)*1e6;ax.quiver(*start,*(end-start),color=colors[k],arrow_length_ratio=.15,linewidth=2)
        axis3(ax,geometry_roi);ax.set_title(label+'\n'+(['Single inward Z','Both start at 0 ms','Source 1 delayed 0.4 ms'][j]))
    save(fig,'Fig02b_actual_loading_locations','Actual historical and frozen new external sphere centers and complete inward travel arrows on the same original geometry, physical ROI and camera. A/B sites are identical. Historical geometry is read from its accepted indenter implementation; it is not a same-geometry source-number control.')
    after=read(O/'CRACK_BRIDGE_AFTER_LOADING.json')['cases'];before=read(O/'CRACK_BRIDGE_AUDIT.json')['cases'];fig=plt.figure(figsize=(21,9));bond_samples={}
    for j,name in enumerate(['historical_single_uniform']+names):
        path=HISTORY/'cases/uniform/keys/final.npz' if j==0 else O/'cases'/name/'keys/final.npz';s=load_npz(path);q=(before['uniform'] if j==0 else after[name])['final']['regions'][0]
        ax=fig.add_subplot(1,3,j+1,projection='3d');ax.scatter(*(X[:n]*1e6).T,s=3,c='#c8cbd0',alpha=.22);ax.scatter(*(X[n:]*1e6).T,s=3,c='#66717d',alpha=.35)
        failed=np.flatnonzero((s['d']==1)&~ref['initial_cut']&(ref['bond_types']==1));failed=failed[::max(1,int(np.ceil(len(failed)/1800)))];live=np.asarray(q['original_surviving_bond_ids'])-1;live=live[::max(1,int(np.ceil(len(live)/350)))];bond_samples[name]=dict(failed_original_bond_ids=failed+1,bridge_original_bond_ids=live+1)
        ax.add_collection3d(Line3DCollection(X[b[failed]]*1e6,colors='#c44938',linewidths=.5,alpha=.5));ax.add_collection3d(Line3DCollection(X[b[live]]*1e6,colors='#2d719d',linewidths=.7,alpha=.7))
        interface=np.flatnonzero((s['d']==1)&~ref['initial_cut']&(ref['bond_types']==2));bond_samples[name]['failed_interface_original_bond_ids']=interface+1
        if len(interface):ax.add_collection3d(Line3DCollection(X[b[interface]]*1e6,colors='#9c3977',linewidths=1.2,alpha=.9))
        path=q['off_midplane_anchor_path']
        if path['found']:ax.plot(*(np.asarray(path['reference_coordinates_m'])*1e6).T,color='#16242c',linewidth=3,label='Live path to anchor')
        axis3(ax,np.array([X.min(0)-DX,X.max(0)+DX]));ax.set_title(name.replace('_',' ')+'\n'+f"Upper-plane bridges: {q['surviving_crossing_bonds']}")
    fig.legend(handles=[Line2D([0],[0],color='#c44938',label='Naturally failed bulk bonds (sample)'),Line2D([0],[0],color='#9c3977',label='Naturally failed interface bonds (all)'),Line2D([0],[0],color='#2d719d',label='Surviving cross-region bonds (sample)'),Line2D([0],[0],color='#16242c',lw=3,label='Full off-midplane live anchor path')],loc='lower center',bbox_to_anchor=(.5,.035),ncol=2);fig.suptitle('Final graph drawn at REFERENCE coordinates · no physical bond deletion');fig.subplots_adjust(wspace=.07,bottom=.16)
    save(fig,'Fig03_failed_network_and_surviving_bridges','Final actual damage/topology on reference geometry for comparable bridge locations. Original bulk/interface identity is retained; all interface failures are displayed separately. The original notch is not counted as naturally failed. Displayed bulk/bridge samples only reduce clutter; all original IDs and full counts remain in the evidence. Black path stays outside |Y| ≤ 50 µm and proves surviving three-dimensional connectivity to anchors. Pinned graph cuts are diagnostics, not fracture energy.')
    dump(F/'PLOTTED_BOND_SAMPLES.json',bond_samples)
    fig=plt.figure(figsize=(16,8))
    for j,name in enumerate(names):
        ax=fig.add_subplot(1,2,j+1,projection='3d');ax.scatter(*(X[:n]*1e6).T,s=3,c='#c8cbd0',alpha=.25);births=read(O/'cases'/name/'CRACK_NETWORK_AUDIT.json')['births']
        if births:
            p=np.array([q['reference_centroid_m'] for q in births])*1e6;sc=ax.scatter(*p.T,c=[q['time_s']*1e3 for q in births],cmap='viridis',s=35,vmin=0,vmax=3.6,depthshade=False)
        axis3(ax,np.array([X.min(0)-DX,X.max(0)+DX]));ax.set_title(f'{titles[j]}: {len(births)} failed-network births')
    if births:fig.colorbar(sc,ax=fig.axes,shrink=.55,pad=.03,label='Birth time (ms)')
    save(fig,'Fig03b_damage_nucleation_locations','Birth centroids use actual original bonds of new fully failed bulk networks, excluding precuts. A birth is a discrete failed-edge component event, not a resolved crack surface or a free fragment. Reference geometry is explicitly used.')
    keys=[('initial','Initial'),('first_softening','First softening'),('first_failure','First full failure'),('first_multinode','First ≥2-node free component'),('max_multifragment','Highest frozen objective tuple'),('final','Final accepted state')]
    for name,title in zip(names,titles):
        fig=plt.figure(figsize=(21,15))
        for j,(key,label) in enumerate(keys):
            ax=fig.add_subplot(2,3,j+1,projection='3d');p=O/'cases'/name/'keys'/(key+'.npz')
            if p.exists():
                s=load_npz(p);state3(ax,name,s,label)
                if key in ['first_softening','first_failure']:
                    ids=np.asarray(summaries[name]['events'][key]['original_bond_ids'])-1;points=s['x'][b[ids]].mean(1)*1e6;ax.scatter(*points.T,marker='*',s=100,c='black',depthshade=False)
            else:ax.set_axis_off();ax.text2D(.1,.5,label+'\nNOT OBSERVED',transform=ax.transAxes,fontsize=18)
        fig.suptitle(title+' · actual accepted displacement ×1 · common ROI and camera');fig.subplots_adjust(hspace=.28,wspace=.06)
        save(fig,'Fig04_'+name+'_evolution','Actual accepted FP64 states, displacement ×1, gray attached body and component-ID colors for detached points. Missing events are explicitly absent. All particles inside the fixed ROI are plotted; out-of-ROI counts are shown. A ≥2-node event can be unresolved debris. All internal live bonds are shown for ≥14-node free components, with no invented shell.')
    fig,axs=plt.subplots(2,2,figsize=(17,11));labels=['1–2','3–13','14–33','34–113','114–200','≥201']
    for j,(name,color) in enumerate(zip(names,colors)):
        sizes=np.array([q['nodes'] for q in results[name]['final_components']]);count=[];mass=[]
        for lo,hi in BINS:
            v=sizes[(sizes>=lo)&((sizes<=hi) if hi else True)];count.append(len(v));mass.append(v.sum()*RHO*V*1e9)
        axs[0,0].bar(np.arange(6)+(j-.5)*.36,count,width=.36,color=color,label=titles[j]);axs[0,1].bar(np.arange(6)+(j-.5)*.36,mass,width=.36,color=color)
    for ax in axs[0]:ax.set_xticks(range(6),labels);ax.set_xlabel('Nodes per actual free component')
    axs[0,0].set_ylabel('Final component count');axs[0,0].legend();axs[0,1].set_ylabel('Final detached mass (µg)');series(axs[1,0],'N20_200',ylabel='20–200-node component count');axs[1,0].axhline(5,c='black',ls=':',label='Exploratory target ≥5');axs[1,0].set_ylim(0,max(5.5,axs[1,0].get_ylim()[1]));axs[1,0].legend();series(axs[1,1],'mass_fraction_in_target_band',100,'Initial mass in target band (%)');fig.tight_layout()
    save(fig,'Fig05_size_distributions','Six frozen size bins from actual free connected components. Volumes use node count × reference node volume, and equivalent diameters are in CSV/JSON. Bins are exploratory discretization categories, not medical particle classes. Zero bars are retained.')
    fig,axs=plt.subplots(2,3,figsize=(20,10))
    for ax,(key,scale,label,gate) in zip(axs.flat,[('N_ge20',1,'Free components ≥20 nodes',None),('N_ge34',1,'Free components ≥34 nodes',None),('N_ge114',1,'Free components ≥114 nodes',None),('total_detached_mass_fraction',100,'Detached / initial mass (%)',10),('largest_free_fragment_fraction_detached',100,'Largest / detached mass (%)',50),('fine_debris_mass_fraction_detached',100,'1–2-node / detached mass (%)',20)]):
        series(ax,key,scale,label)
        if gate is not None:ax.axhline(gate,c='black',ls=':',lw=1);ax.set_ylim(0,max(ax.get_ylim()[1],gate*1.1))
    axs[0,0].legend();fig.tight_layout();save(fig,'Fig05b_independent_objective_constraints','Independent size, mass, dominance and fine-debris metrics; no weighted total masks failure. Undefined fractions before detachment are gaps. Dotted thresholds are frozen exploratory constraints; common cohort time 2.8 ms is dashed.')
    fig=plt.figure(figsize=(18,8))
    for j,name in enumerate(names):
        ax=fig.add_subplot(1,2,j+1,projection='3d');s=load_npz(O/'cases'/name/'keys/final.npz');state3(ax,name,s,titles[j]+' · final component identities',False)
    save(fig,'Fig06_fragment_ID_spatial_view','Exact member IDs define free components, colored consistently within each case across the evolution panels. Attached material remains neutral gray. Component IDs are labeled for ≥14-node groups; all component memberships are available in registry and history. Camera clipping is counted, not hidden.')
    fig,axs=plt.subplots(1,2,figsize=(17,8))
    for ax,name,title in zip(axs,names,titles):
        s=load_npz(O/'cases'/name/'keys/final.npz');c,_,_,_,_=point_colors(name,s);ax.scatter(s['x'][:n,0]*1e3,s['x'][:n,2]*1e3,s=8,c=c);ax.scatter(s['x'][n:,0]*1e3,s['x'][n:,2]*1e3,s=5,c='#5e6670');ax.set(xlabel='X (mm)',ylabel='Z (mm)',title=title+' · all final points');ax.set_aspect('equal',adjustable='box')
    fig.tight_layout();save(fig,'Fig06b_full_extent','All final points, without the 3D ROI clipping, in X–Z projection. Colors use exact surviving-component identity; overlap in projection does not merge fragments.')
    fig,axs=plt.subplots(2,2,figsize=(18,11))
    for j,name in enumerate(names):
        fragments=read(O/'cases'/name/'FRAGMENT_HISTORY.json')['fragments'];trace=read(O/'cases'/name/'FRAGMENT_TRACE.json');reg=registry[name];selected=[q for q in fragments if q['nodes']>=2];lifetimes={q['fragment_id']:q for q in fragments}
        for q in selected:
            axs[0,j].plot([q['birth_time_s']*1e3,q['last_observed_s']*1e3],[q['nodes']]*2,c=colors[j],lw=1.2);axs[0,j].scatter(q['birth_time_s']*1e3,q['nodes'],s=22,c=colors[j])
        for e in summaries[name]['secondary_events']:
            parent=reg[e['parent_component_id']]
            for child in e['children_component_ids']:axs[0,j].plot([lifetimes[e['parent_component_id']]['last_observed_s']*1e3,e['time_s']*1e3],[parent['nodes'],reg[child]['nodes']],c='#555555',lw=.8)
        axs[0,j].set(title=titles[j],xlabel='Birth / observation time (ms)',ylabel='Nodes per component',xlim=(0,3.6))
        axs[0,j].yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True,nbins=5))
        if not selected:axs[0,j].text(.1,.5,'No ≥2-node free component observed',transform=axs[0,j].transAxes)
        for cid in sorted(trace,key=lambda v:reg[int(v)]['nodes'],reverse=True)[:5]:
            ts=trace[cid];axs[1,j].plot([q['time_s']*1e3 for q in ts],[q['internal_survival'] for q in ts],label=f'ID {cid}, n={reg[int(cid)]["nodes"]}')
        if trace:axs[1,j].legend()
        else:axs[1,j].text(.1,.5,'No multinode cohesion trace',transform=axs[1,j].transAxes)
        axs[1,j].set(xlabel='Time (ms)',ylabel='Internal original bond survival',ylim=(-.04,1.04),xlim=(0,3.6))
        for ax in axs[:,j]:ax.axvline(2.8,c='#666666',ls='--')
    fig.tight_layout();save(fig,'Fig07_fragment_birth_split_and_cohesion','Exact original memberships define births and split ancestry. Lifelines are shown for all ≥2-node components; up to five largest have internal-cohesion traces. A two-node pair remains unresolved debris, even with 100% survival. Lifetimes are right-censored by the actual final time; the common 2.8 ms cohort is dashed.')
    fig,axs=plt.subplots(4,2,figsize=(18,19))
    for j,name in enumerate(names):
        t=[q['time_s']*1e3 for q in rows[name]]
        for source,color in enumerate(colors):
            axs[0,j].plot(t,[q[f'site{source}_reaction_norm_N']*1e6 for q in rows[name]],c=color,label=f'Source {source}');axs[1,j].plot(t,[q[f'site{source}_W_J']*1e9 for q in rows[name]],c=color,label=f'Source {source}')
        for key,label in [('W_J','W0+W1'),('U_J','Stored U'),('K_J','Kinetic K'),('D_J','Damage D'),('Uc_J','Particle contact'),('U_indenter_J','Tool contact')]:axs[2,j].plot(t,[q[key]*1e9 for q in rows[name]],label=label)
        for key,label in [('R_J','Signed residual R'),('Jc_J','Contact activation Jc')]:axs[3,j].plot(t,[q[key]*1e9 for q in rows[name]],label=label)
        for ax,ylabel in zip(axs[:,j],['Tool reaction norm (µN)','Per-source work (nJ)','Energy (nJ)','Ledger terms (nJ)']):ax.set(xlabel='Time (ms)',ylabel=ylabel);ax.legend(ncol=2);ax.axvline(tools[name].clear_time*1e3,c='#666666',ls='--')
        axs[0,j].set_title(titles[j])
    fig.tight_layout();save(fig,'Fig08_per_source_force_work_energy','Actual independent native source reactions and independently integrated work. W0+W1 = ΔU + ΔK + D + ΔUc − Jc + ΔUi + R; no damping. R is numerical accounting error, not dissipation. Dashed lines mark each protocol’s scheduled return; recontact remains possible.')
    base='../../../verification/'+TASK+'/'+RUN+'/'
    html='<!doctype html><meta charset="utf-8"><title>Clot multifragment breakup</title><style>body{font:18px sans-serif;max-width:1700px;margin:30px auto;background:#f3f5f7;color:#263445}section{background:white;padding:24px;margin:24px 0}img{width:100%}p{line-height:1.6}</style><h1>多区域加载与真实多块碎裂检验</h1><p>MULTISITE_MECHANICAL_SURROGATE · actual native accepted states · displacement ×1</p>'
    html+='<p>'+' · '.join(f'<a href="{base+file}">{label}</a>' for file,label in [('REVIEW_CLOT_MULTIFRAGMENT_BREAKUP_V1.md','中文科研报告'),('FINAL_STATUS.json','最终状态'),('MULTIFRAGMENT_RESULTS.csv','多块指标'),('FRAGMENT_SIZE_DISTRIBUTIONS.csv','尺寸分布'),('SURVIVING_BRIDGE_BONDS_AFTER_LOADING.csv','实际残余桥接键')])+'</p>'
    for q in figures:html+=f'<section><h2>{q["name"]}</h2><a href="{q["name"]}.svg"><img src="{q["name"]}.png"></a><p>{q["caption"]}</p></section>'
    (F/'OPEN_RESULTS.html').write_text(html);(O/'OPEN_RESULTS.html').write_text(f'<!doctype html><meta http-equiv="refresh" content="0;url=../../../results/{TASK}/{RUN}/OPEN_RESULTS.html">')
    dump(F/'FIGURE_MANIFEST.json',dict(figures=figures,actual_coordinate_scale=1,camera=dict(elev=16,azim=-63,ROI_m=ROI),reference_coordinate_exceptions=['Fig03_failed_network_and_surviving_bridges','Fig03b_damage_nucleation_locations']))
    print('FIGURES_DONE',len(figures))
if __name__=='__main__':main()
