"""Publication figures from frozen fields and accepted FP64 states only."""
from analyze import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.mplot3d.art3d import Poly3DCollection,Line3DCollection
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':14,'axes.titlesize':15,'axes.labelsize':13,'xtick.labelsize':11,'ytick.labelsize':11,'legend.fontsize':11,'svg.fonttype':'none'})
MANIFEST=[];COLORS=['#307d9a','#bf6846'];TITLES={'uniform':'A · Unidirectional','counter_rotating':'B · Counter-rotating'}
def export(fig,name,caption):
    F.mkdir(parents=True,exist_ok=True)
    for ext in ['png','svg']:fig.savefig(F/(name+'.'+ext),dpi=200,bbox_inches='tight',facecolor='white')
    plt.close(fig);MANIFEST.append(dict(name=name,caption=caption))
def axes3(ax,title,points):
    ax.set(xlabel='X (µm)',ylabel='Y (µm)',zlabel='Z (µm)',title=title);ax.view_init(15,-65)
    lo=points.min(0);hi=points.max(0);span=np.maximum(hi-lo,50);pad=span*.08
    ax.set_xlim(lo[0]-pad[0],hi[0]+pad[0]);ax.set_ylim(lo[1]-pad[1],hi[1]+pad[1]);ax.set_zlim(lo[2]-pad[2],hi[2]+pad[2]);ax.set_box_aspect(span)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.set_major_locator(MaxNLocator(3));axis._axinfo['grid'].update(color='#929ba6',linewidth=.7)
    ax.set_xlabel('X (µm)',labelpad=9);ax.set_ylabel('Y (µm)',labelpad=9);ax.set_zlabel('Z (µm)',labelpad=10)
def body(ax,ref,s):
    n=int(ref['clot_count']);x=s['x']*1e6
    ax.scatter(*x[:n].T,c='#b5bbc2',s=5,alpha=.24,depthshade=False)
    ax.scatter(*x[n:].T,c='#414f61',s=5,alpha=.3,depthshade=False)
def surface_panel(ax,ref,s,tool,title,zoom=False):
    ax.computed_zorder=False
    z=tool.evaluate(s['x'],s['velocity'],s['d'],float(s['time']));ids=np.flatnonzero(z['active']);norm=Normalize(0,1000)
    if not zoom:body(ax,ref,s)
    if len(ids):
        ax.add_collection3d(Poly3DCollection(s['x'][tool.tri[ids]]*1e6,facecolors=plt.cm.viridis(norm(np.linalg.norm(z['traction'][ids],axis=1))),edgecolors='#4b5e6b',linewidths=.35,alpha=.65,zorder=3))
        sample=ids[::max(1,len(ids)//20)];cent=z['centroid'][sample]*1e6
        ax.quiver(*cent.T,*(z['traction'][sample]*.08).T,color='#ed4b29',arrow_length_ratio=.35,linewidth=1.7,zorder=10)
        ax.quiver(*cent.T,*(z['normal'][sample]*35).T,color='#2e3e50',arrow_length_ratio=.28,linewidth=.8,zorder=9)
    points=s['x']*1e6
    if zoom and len(ids):points=np.vstack([s['x'][tool.tri[ids]].reshape(-1,3)*1e6,cent+z['normal'][sample]*35,cent+z['traction'][sample]*.08])
    axes3(ax,title,points)
    if zoom:
        ax.view_init(30,-65)
        if np.ptp(points[:,2])<100:ax.set_zticks([round(float(points[:,2].mean()))])
    ax.text2D(.01,.92,f't={float(s["time"])*1e3:.3f} ms\nretained faces {len(ids)}/182',transform=ax.transAxes,fontsize=11)
def pick(case,key):
    p=case/'keys'/(key+'.npz')
    return load_npz(p) if p.exists() else None
def nearest(case,t):
    states=[load_npz(p) for p in (case/'keys').glob('frame_*.npz')]
    return min(states,key=lambda s:abs(float(s['time'])-t))
def main(inputs_only=False):
    ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);tools={k:CurrentSurface(ref,k) for k in CASES};cases={k:O/'cases'/k for k in CASES}
    fig,axs=plt.subplots(1,3,figsize=(19,5.3),layout='constrained')
    for j,name in enumerate(CASES):
        tool=tools[name];s=tool.evaluate(ref['X'],np.zeros_like(ref['X']),ref['initial_cut'].astype(float),RISE);xy=s['centroid'][:,:2]*1e6;mag=np.linalg.norm(s['traction'],axis=1)
        axs[j].tripcolor(ref['X'][:,0]*1e6,ref['X'][:,1]*1e6,triangles=tool.tri,facecolors=mag,cmap='viridis',vmin=0,vmax=1000,edgecolors='#79919b',linewidth=.2)
        ids=np.arange(0,len(xy),3);axs[j].quiver(*xy[ids].T,*s['traction'][ids,:2].T,angles='xy',scale_units='xy',scale=14,color='#c74934',width=.004)
        axs[j].set(xlabel='X (µm)',ylabel='Y (µm)',title=TITLES[name]+'\nInitial flat surface at peak envelope');axs[j].set_aspect('equal')
        axs[j].text(.02,.02,f'actual max {mag.max():.1f} Pa',transform=axs[j].transAxes,fontsize=11,bbox=dict(fc='white',alpha=.8,ec='none'))
    fig.colorbar(ScalarMappable(norm=Normalize(0,1000),cmap='viridis'),ax=list(axs[:2]),shrink=.7,label='Stationary driving shear (Pa)')
    tt=np.linspace(0,END,701);env=[tools['uniform'].field.envelope(t) for t in tt]
    axs[2].plot(tt*1e3,np.array(env)*1000,c=COLORS[0]);axs[2].axvline(OFF*1e3,ls='--',c='.5');axs[2].set(xlabel='Time (ms)',ylabel='Stationary drive cap (Pa)',title='Shared envelope / lab-fixed region')
    axs[2].text(.04,.42,'µ = 1 mPa·s; h = 25 µm\nU cap = 25 m/s; Reₕ cap = 625\nAfter 2.8 ms: passive drag only\n\nMechanical requirement proxy\nNo bubble-streaming solve',transform=axs[2].transAxes,fontsize=12,bbox=dict(fc='white',alpha=.85,ec='none'))
    export(fig,'Fig1_two_shear_fields','Defined mechanical fields on actual original surface at rest. Arrows: 1000 Pa = 71.4 µm. Counter-rotating planar Gaussian streamfunction is inspired by vortex morphology, not a reproduction of a bubble solution.')
    if inputs_only:dump(O/'INPUT_FIGURE_MANIFEST.json',dict(figures=MANIFEST));return
    analysis=read(O/'ANALYSIS.json');hist={k:histories(v) for k,v in cases.items()};bridge=read(O/'CRACK_BRIDGE_AUDIT.json');networks=read(O/'DAMAGE_NETWORK_HISTORY.json')
    for name,case in cases.items():
        fig=plt.figure(figsize=(22,11));states=[pick(case,'initial'),nearest(case,.00075),pick(case,'first_softening'),pick(case,'final')]
        for j,(title,s) in enumerate(zip(['Initial','During loading','First softening','Final'],states),1):
            ax=fig.add_subplot(2,4,j,projection='3d')
            if s is None:ax.text2D(.2,.5,'NOT OBSERVED',transform=ax.transAxes);continue
            surface_panel(ax,ref,s,tools[name],title)
            ax=fig.add_subplot(2,4,j+4,projection='3d');surface_panel(ax,ref,s,tools[name],'Retained surface detail',zoom=True)
        fig.suptitle(TITLES[name]+' | current original surface, physical coordinates ×1',y=1.01)
        fig.text(.5,.01,'Red: actual slip shear, 1000 Pa → 80 µm | dark: unit normal → 35 µm | color: 0–1000 Pa\nNEW_CRACK_FACE_WETTING_UNRESOLVED; absent facets are not fabricated crack surfaces',ha='center',fontsize=12)
        fig.subplots_adjust(bottom=.15,wspace=.18,hspace=.42);export(fig,'Fig2_'+name+'_surface_following','Actual original-ID facets at each accepted state, area/normal updated. Top: full body; bottom: retained-facet detail with actual coordinates, no deformation magnification. Broken-edge unsupported triangles are omitted from loads, not removed from the material graph. Axes fit actual positions per panel.')
    fig=plt.figure(figsize=(18,10))
    for j,name in enumerate(CASES):
        case=cases[name];s=pick(case,'first_failure')
        if s is None:s=pick(case,'peak_load')
        z=tools[name].evaluate(s['x'],s['velocity'],s['d'],float(s['time']));x=s['x']*1e6
        ax=fig.add_subplot(2,2,2*j+1,projection='3d');surface_panel(ax,ref,s,tools[name],TITLES[name]+' · facet shear',zoom=True)
        ax=fig.add_subplot(2,2,2*j+2,projection='3d');ax.computed_zorder=False;f=s['external_force'];mag=np.linalg.norm(f,axis=1);ids=np.flatnonzero(mag>1e-18)
        sc=ax.scatter(*x[ids].T,c=mag[ids]*1e6,cmap='magma',s=24,depthshade=False,zorder=2);sample=ids[::2];ax.quiver(*x[sample].T,*(f[sample]*5e7).T,color='#c63e29',linewidth=1.4,arrow_length_ratio=.35,zorder=10);axes3(ax,TITLES[name]+f' · native node force\nt={float(s["time"])*1e3:.3f} ms | 1 µN → 50 µm',np.vstack([x[ids],x[sample]+f[sample]*5e7]));ax.view_init(30,-65);ax.set_zticks([round(float(x[ids,2].mean()))])
        fig.colorbar(sc,ax=ax,orientation='horizontal',shrink=.65,pad=.12,label='|Node force| (µN)')
    fig.subplots_adjust(wspace=.12,hspace=.42);export(fig,'Fig3_current_native_force_mapping','Actual native external-force arrays at each mode’s first failure, surface detail with physical coordinates. Positive 1/3 barycentric map uses current triangles; no signed affine V3 weights. Facet arrows follow Figure 2 scale; native force arrows 1 µN = 50 µm.')
    fig=plt.figure(figsize=(20,11))
    for row,name in enumerate(CASES):
        for col,key in enumerate(['first_softening','first_failure','final']):
            ax=fig.add_subplot(2,3,row*3+col+1,projection='3d');s=pick(cases[name],key)
            if s is None:ax.text2D(.2,.5,'NOT OBSERVED',transform=ax.transAxes);continue
            body(ax,ref,s);b=ref['bonds'];natural=~ref['initial_cut']&np.isin(ref['bond_types'],[1,2]);soft=np.flatnonzero(natural&(s['d']>0)&(s['d']<1));failed=np.flatnonzero(natural&(s['d']==1))
            for ids,color in [(soft,'#d99837'),(failed,'#b73431')]:
                ends=np.unique(b[ids]);ax.scatter(*(s['x'][ends]*1e6).T,c=color,s=8,alpha=.8,depthshade=False)
            if key=='final':
                ids=np.array(bridge['cases'][name]['final']['regions'][0]['original_surviving_bond_ids'])-1
                ax.add_collection3d(Line3DCollection(s['x'][b[ids]]*1e6,colors='#287c89',linewidths=.3,alpha=.13))
            axes3(ax,TITLES[name]+' · '+key.replace('_',' '),s['x']*1e6);ax.text2D(.01,.94,f'{float(s["time"])*1e3:.4f} ms | soft {len(soft)} | failed {len(failed)}',transform=ax.transAxes,fontsize=11)
    fig.suptitle('Actual damage endpoints and surviving 3D bridges | coordinates ×1',y=1.01);fig.text(.5,.01,'Orange: partial damage endpoints | red: naturally failed endpoints | teal: live original bridge bonds\nInitial notch excluded from new failures. Failed-edge connectivity is not a closed fracture surface.',ha='center',fontsize=12);fig.subplots_adjust(bottom=.12,hspace=.40);export(fig,'Fig4_crack_network_evolution','Current coordinates of actual damaged original bonds, anchors retained. A colored endpoint is not an inferred wetted crack face. Final bridges cross the historical diagnostic upper volume; no bonds are deleted by analysis.')
    fig,axs=plt.subplots(1,2,figsize=(15,5),layout='constrained')
    for name,color in zip(CASES,COLORS):
        q=networks[name]['all_material_including_interface'];hh=q['history'];axs[0].step([r['time_s']*1e3 for r in hh],[r['separate_failed_material_networks'] for r in hh],where='post',label=TITLES[name],color=color)
        rr=hist[name];axs[1].plot([r['time_s']*1e3 for r in rr],[r['N20_200'] for r in rr],label=TITLES[name],color=color)
    axs[0].set(xlabel='Time (ms)',ylabel='Separate failed-edge networks',title='Natural bulk + interface failure graph');axs[1].set(xlabel='Time (ms)',ylabel='True free 20–200 node components',title='Surviving material graph + anchor reachability');axs[1].set_ylim(-.1,max(1.1,max(q['history_max_N20_200'] for q in analysis.values())+1));axs[1].yaxis.set_major_locator(MaxNLocator(integer=True));axs[0].legend();axs[1].legend();export(fig,'Fig4b_failed_networks_vs_true_fragments','Independent graphs: failed bonds define damage clusters; surviving bonds define real material components. Initial cuts and contact do not count as natural fractures or cohesion.')
    states={k:pick(cases[k],'final') for k in CASES};allpoints=np.vstack([s['x'] for s in states.values()])*1e6
    fig=plt.figure(figsize=(17,7))
    for j,name in enumerate(CASES,1):
        ax=fig.add_subplot(1,2,j,projection='3d');s=states[name];x=s['x']*1e6;lab,sizes,att=graph.topology(ref,s['d']);attached=att[lab]
        ax.scatter(*x[:n][attached].T,c='#b7bec5',s=6,alpha=.28,depthshade=False);ax.scatter(*x[n:].T,c='#3d4d60',s=6,alpha=.3,depthshade=False)
        for i,label in enumerate(np.flatnonzero(~att)):
            ids=np.flatnonzero(lab==label);color='#b63b31' if len(ids)<=2 else plt.cm.tab20(i%20);ax.scatter(*x[ids].T,c=[color],s=24 if len(ids)>=3 else 16,depthshade=False)
        axes3(ax,TITLES[name]+' · final true components',allpoints)
        ax.text2D(.01,.94,f'free nodes={int(sizes[~att].sum())}; components={int((~att).sum())}\n20–200={int(np.sum((sizes[~att]>=20)&(sizes[~att]<=200)))}',transform=ax.transAxes,fontsize=12)
    fig.text(.5,.01,'Gray: attached body | dark: fixed anchors | red: 1–2 node fines | colors: ≥3 node components\nAll final nodes shown at actual positions ×1; same axes, no invented surfaces or display separation',ha='center',fontsize=12);fig.subplots_adjust(bottom=.17);export(fig,'Fig5_true_detached_fragments','Full actual final extents. Original surviving material bonds and anchor reachability alone identify components; contact/proximity never merges fragments. Colors may repeat.')
    table=list(csv.DictReader((O/'FRAGMENT_SIZE_DISTRIBUTION.csv').open()));fig,axs=plt.subplots(1,2,figsize=(16,5),layout='constrained')
    for j,name in enumerate(CASES):
        q=[r for r in table if r['case']==name and r['state']=='final'];xx=np.arange(6)+(j-.5)*.35
        axs[0].bar(xx,[float(r['component_count']) for r in q],.35,label=TITLES[name],color=COLORS[j]);axs[1].bar(xx,[float(r['mass_fraction_initial'])*100 for r in q],.35,label=TITLES[name],color=COLORS[j])
    for ax in axs:ax.set_xticks(np.arange(6),['1–2','3–13','14–33','34–113','114–200','>200']);ax.set_xlabel('Original material nodes per free component');ax.legend()
    axs[0].set(ylabel='Component count',title='Final size distribution');axs[1].set(ylabel='Initial material mass (%)',title='Actual reference-volume mass');export(fig,'Fig6_fragment_size_distribution','Zero bins remain zero. Equal original material volumes determine mass; anchors excluded. Requested bins cross the 20-node target threshold; exact 20–200 is reported separately.')
    fig,axs=plt.subplots(2,3,figsize=(20,10),layout='constrained')
    for j,name in enumerate(CASES):
        rr=hist[name];t=np.array([r['time_s'] for r in rr])*1e3
        for key,label in [('W_J','Net surface work'),('U_J','Stored LPS'),('K_J','Kinetic'),('D_J','Damage dissipation'),('Uc_J','Contact energy')]:axs[j,0].plot(t,[r[key]*1e9 for r in rr],label=label)
        for key,label in [('W_drive_J','Outer-flow drive'),('W_drag_J','Passive slip drag'),('W_J','Net surface work')]:axs[j,1].plot(t,[r[key]*1e9 for r in rr],label=label)
        axs[j,1].axvline(OFF*1e3,c='.5',ls='--',label='Driver off')
        for key,label in [('R_J','Numerical balance R'),('Jc_J','Contact activation J'),('work_path_vs_endpoint_J','Path − endpoint integral')]:axs[j,2].plot(t,[r[key]*1e9 for r in rr],label=label)
        for col,title in enumerate(['Energy partition','Drive / drag work','Independent accounting terms']):axs[j,col].set(xlabel='Time (ms)',ylabel='Energy (nJ)',title=TITLES[name]+' · '+title);axs[j,col].legend(fontsize=10)
    export(fig,'Fig7_actual_energy_ledger','Work uses current force along actual native drift. Passive drag is physical-proxy dissipation, separate from irreversible material damage. Numerical residual and contact activation remain explicit; no fitted energy correction.')
    fig,axs=plt.subplots(2,4,figsize=(20,10),layout='constrained')
    metricspec=[('history_max_N20_200','Maximum target fragment count','count'),('mass','Final target-band mass','% initial mass'),('largest','Largest free component, history','nodes'),('fine','Final 1–2 node fines','% detached mass'),('work','Actual net external work','nJ'),('damage','Material damage dissipation','nJ'),('residual','Maximum numerical residual','%'),('surface','Final retained original area','% reference area')]
    for ax,(key,title,unit) in zip(axs.ravel(),metricspec):
        values=[]
        for name in CASES:
            q=analysis[name];e=q['final'];v={'history_max_N20_200':q['history_max_N20_200'],'mass':e['mass_fraction20_200']*100,'largest':q['history_max_free_nodes'],'fine':None if e['fine_debris_mass_fraction_detached'] is None else e['fine_debris_mass_fraction_detached']*100,'work':e['W_J']*1e9,'damage':e['D_J']*1e9,'residual':q['max_relative_energy_residual']*100,'surface':e['active_reference_area_fraction']*100}[key];values.append(v)
        ax.bar(['A','B'],[v if v is not None else 0 for v in values],color=COLORS)
        for j,v in enumerate(values):ax.text(j,0 if v is None else v,'undefined' if v is None else f'{v:.4g}',ha='center',va='bottom',fontsize=12)
        ax.set(title=title,ylabel=unit);ax.margins(y=.2)
        if max((v or 0) for v in values)==0:ax.set_ylim(0,1);ax.set_yticks([0,1])
    fig.suptitle('Same frozen material / envelope / nominal drive cap; actual work differs',fontsize=17);export(fig,'Fig8_final_controlled_comparison','A versus B isolates the spatial mode within the new follower protocol. It does not isolate follower geometry relative to older dead-load experiments; actual slip peaks and work are not matched.')
    coverage=read(O/'FRAGMENT_SURFACE_COVERAGE_AUDIT.json')['cases']
    fig,axs=plt.subplots(2,3,figsize=(19,10),layout='constrained')
    for j,name in enumerate(CASES):
        rr=hist[name];t=[r['time_s']*1e3 for r in rr]
        axs[j,0].plot(t,[r['actual_peak_slip_shear_Pa'] for r in rr],label='Relative slip shear');axs[j,0].plot(t,[r['actual_peak_drive_shear_Pa'] for r in rr],label='Driving shear');axs[j,0].set(ylabel='Peak shear (Pa)',title=TITLES[name]+' · actual surface stress')
        axs[j,1].plot(t,[100*r['active_reference_area_fraction'] for r in rr],label='Retained reference area');axs[j,1].set(ylabel='Original loaded area (%)',title='Surface support retention')
        cc=sorted(coverage[name],key=lambda v:v['time_s']);ct=[v['time_s']*1e3 for v in cc]
        axs[j,2].plot(ct,[r['free_components_with_retained_faces'] for r in cc],label='With retained facets');axs[j,2].plot(ct,[r['free_components_without_retained_faces'] for r in cc],label='Without represented facets');axs[j,2].set(ylabel='Detached components',title='Fragment surface representation')
        for ax in axs[j]:ax.set_xlabel('Time (ms)');ax.axvline(OFF*1e3,c='.5',ls='--');ax.legend(fontsize=10)
    export(fig,'Fig9_surface_coverage_and_force','Retained original exterior faces only. Loss of all represented faces on a fragment means unresolved hydrodynamics, not a claim that its real fluid force is zero. New crack wetting remains unresolved.')
    dump(O/'FIGURE_MANIFEST.json',dict(figures=MANIFEST,actual_positions_scale=1,source_sha256=sha(A/'figures.py')));print('FIGURES_DONE',len(MANIFEST))
if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('--inputs-only',action='store_true');args=p.parse_args();main(args.inputs_only)
