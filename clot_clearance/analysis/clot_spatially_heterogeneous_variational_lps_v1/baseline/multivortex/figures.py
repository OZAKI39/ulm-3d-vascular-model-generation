from common import *
from analyze import read_history
from current_surface import CurrentSurface
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.lines import Line2D
import time,resource

plt.rcParams.update({'font.family':'DejaVu Sans','font.size':12,'axes.titlesize':14,'axes.labelsize':12,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'figure.facecolor':'white','savefig.facecolor':'white','axes.grid':True,'grid.color':'#737d88','grid.alpha':.65,'grid.linewidth':.65,'svg.fonttype':'none'})
DAMAGE_CMAP=LinearSegmentedColormap.from_list('supported_damage',['#306b98','#d9b76f','#b83342','#4f2046'])
NAMES={'counter_rotating':'B · Counter-rotating pair','checkerboard':'C · Checkerboard vortices'}
graph=module('figure_graph_v2',R/'analysis/v3_chunk_fragmentation_mechanism_v1/graph_ops.py')
def save(fig,name):
    fig.savefig(F/(name+'.png'),dpi=240,bbox_inches='tight');fig.savefig(F/(name+'.svg'),bbox_inches='tight');plt.close(fig)
def axes3(fig,pos,title,xyz,limits=None):
    ax=fig.add_subplot(pos,projection='3d',computed_zorder=False);ax.view_init(elev=16,azim=-64);ax.set_title(title,pad=15)
    lo,hi=(np.min(xyz,0),np.max(xyz,0)) if limits is None else limits
    margin=np.maximum((hi-lo)*.06,20.)
    for j,setlim in enumerate([ax.set_xlim,ax.set_ylim,ax.set_zlim]):setlim(lo[j]-margin[j],hi[j]+margin[j])
    ax.set_box_aspect(np.maximum(hi-lo,50));ax.set_xlabel('X (µm)',labelpad=8);ax.set_ylabel('Y (µm)',labelpad=8);ax.set_zlabel('Z (µm)',labelpad=8)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.pane.fill=False;axis.set_major_locator(MaxNLocator(4));axis._axinfo['grid'].update(color='#7b8792',linewidth=.7)
    return ax
def cloud(ax,s,ref,kind='damage',size=7):
    n=int(ref['clot_count']);x=s['x']*1e6
    ax.scatter(*x[n:].T,c='#7c858e',s=4,alpha=.25,depthshade=False)
    if kind=='components':
        lab,sz,att=graph.topology(ref,s['d']);colors=np.tile(matplotlib.colors.to_rgba('#afbdc6'),(n,1));palette=plt.get_cmap('tab20')
        for k,label in enumerate(np.flatnonzero(~att)):colors[lab==label]=np.r_[matplotlib.colors.hsv_to_rgb([(k*.61803398875)%1,.72,.78]),1.]
        ax.scatter(*x[:n].T,c=colors,s=size,depthshade=False,edgecolors='none')
    else:
        b=ref['bonds'];mask=~ref['initial_cut'];den=np.bincount(b[mask].ravel(),minlength=len(x));num=np.bincount(b[mask].ravel(),weights=np.repeat(s['d'][mask],2),minlength=len(x));value=num[:n]/np.maximum(den[:n],1)
        ax.scatter(*x[:n].T,c=value,s=size,cmap=DAMAGE_CMAP,vmin=0,vmax=1.,depthshade=False,edgecolors='none')
def detail_figure(ref,finals):
    # Honest magnification of the actual largest 3/2-node free components.
    fig=plt.figure(figsize=(14,7))
    for j,m in enumerate(CASES):
        s=finals[m];lab,sz,att=graph.topology(ref,s['d']);free_labels=np.flatnonzero(~att)
        if not len(free_labels):continue
        label=free_labels[np.argmax(sz[free_labels])];ids=np.flatnonzero(lab==label);xyz=s['x'][ids]*1e6;center=xyz.mean(0)
        ax=axes3(fig,121+j,NAMES[m]+f'\nLargest free component: {len(ids)} nodes',xyz,(center-125,center+125))
        ax.scatter(*xyz.T,s=110,color='#b14a57',edgecolors='#243541',depthshade=False)
        bonds=ref['bonds'];inside=np.isin(bonds,ids).all(1)&(s['d']<1)&(ref['bond_types']==1)
        for a,b in bonds[inside]:ax.plot(*s['x'][[a,b]].T*1e6,color='#294b69',lw=2)
        for point,node in zip(xyz,ids):ax.text(*(point+[2,2,3]),f'ID {node+1}',fontsize=10)
        ax.text2D(.03,.92,f'{inside.sum()} live internal bonds; no fitted shell',transform=ax.transAxes)
    fig.suptitle('Actual-coordinate detail · same 250 µm view width · displacement ×1',fontsize=15)
    save(fig,'Fig5b_Largest_component_detail')

def main():
    start=time.perf_counter();cpu=time.process_time();F.mkdir(parents=True,exist_ok=True)
    ref=load_npz(O/'provenance/reference.npz');n=int(ref['clot_count']);data={};hist={};states={};surfs={}
    for mode,file in zip(CASES,['SINGLE_VORTEX_NATIVE_RESULT.json','CHECKERBOARD_NATIVE_RESULT.json']):
        data[mode]=json.loads((O/file).read_text());hist[mode]=read_history(O/'cases'/mode/'history.csv');states[mode]={p.stem:load_npz(p) for p in (O/'cases'/mode/'keys').glob('*.npz') if p.stem!='failed_unaccepted'};surfs[mode]=CurrentSurface(ref,mode)
    # 1. Identical initial top footprint, unit-scale spatial templates.
    xx,yy=np.meshgrid(np.linspace(-325,325,23),np.linspace(-175,175,15));pos=np.c_[xx.ravel(),yy.ravel(),np.full(xx.size,375)]*1e-6;normal=np.tile([0,0,1.],(len(pos),1))
    fig,axs=plt.subplots(1,3,figsize=(17,4.4),layout='constrained')
    for ax,mode in zip(axs,['uniform']+CASES):
        flow=np.tile([1.,0,0],(len(pos),1)) if mode=='uniform' else surfs[mode].field.spatial(pos,normal)
        ax.quiver(pos[:,0]*1e6,pos[:,1]*1e6,flow[:,0],flow[:,1],np.linalg.norm(flow,axis=1),cmap='viridis',clim=(0,1),scale=20,width=.005)
        ax.set_aspect('equal');ax.set_xlabel('X (µm)');ax.set_ylabel('Y (µm)');ax.set_xlim(-350,350);ax.set_ylim(-200,200)
        ax.set_title('Historical A · Uniform\nOLD surface model; budget only' if mode=='uniform' else NAMES[mode]+'\nNew exposed-surface model')
    fig.suptitle('Spatial patterns on the actual 650 × 350 µm top · normalized vectors, not equal peak loads',fontsize=15)
    save(fig,'Fig1_Load_patterns')
    # 2. Categories use one unique physical triangle each. Detached is an overlay category.
    colors=['#3182bd','#f6ae2d','#2a9d66','#9844b7','#777777'];labels=['Original exposed','Open, not accessible','Accessible crack','Detached exposed','Closed / unresolved']
    for mode in CASES:
        snapshots=data[mode]['surface']['snapshots'];best=max(snapshots,key=lambda q:q['accessible_new_faces']);s=states[mode][best['key']];surf=surfs[mode];amp=float(s['field_amplitude_Pa']);surf.field.segment=(0,END,amp,amp);z=surf.evaluate(s['x'],s['velocity'],s['d'],float(s['time']))
        lab,sz,att=graph.topology(ref,s['d']);nodefree=np.r_[~att[lab],np.zeros(len(ref['X'])-n,bool)];detached=z['active']&nodefree[surf.tri[:,0]]
        category=np.full(len(surf.tri),-1);category[z['unresolved']]=4;category[z['new_open']]=1;category[z['original']]=0;category[z['accessible']]=2;category[detached]=3
        fig=plt.figure(figsize=(11,9.5));fig.subplots_adjust(bottom=.24,top=.86);ax=axes3(fig,111,NAMES[mode]+f" · t = {float(s['time'])*1e3:.3f} ms\nSupported current faces · quasistatic fluid-accessibility proxy",s['x']*1e6)
        shown=category>=0;collection=Poly3DCollection(s['x'][surf.tri[shown]]*1e6,facecolors=[colors[q] for q in category[shown]],edgecolors=(.1,.15,.2,.15),linewidths=.15,alpha=.7);ax.add_collection3d(collection)
        ax.scatter(*s['x'][:n].T*1e6,s=2,c='#203747',alpha=.25,depthshade=False)
        fig.legend([Line2D([0],[0],marker='s',ls='',color=c,markersize=9) for c in colors],labels,loc='upper center',ncol=3,bbox_to_anchor=(.5,.14))
        fig.text(.02,.025,'Actual node coordinates ×1. Unresolved faces are not loaded. No reconstructed fine-particle shells.',fontsize=10)
        save(fig,'Fig2_Current_surfaces_'+mode)
    # 3. Actual shear vectors, slip speed and nodal forces at nearest saved peak state.
    fig=plt.figure(figsize=(18,12))
    for row,mode in enumerate(CASES):
        target=max(hist[mode],key=lambda q:q['shear_peak_Pa'])['time_s'];s=min(states[mode].values(),key=lambda q:abs(float(q['time'])-target));surf=surfs[mode];active=s['surface_active'];ids=np.flatnonzero(active);ids=ids[::max(1,len(ids)//90)];pos=s['surface_centroid']*1e6;tau=np.linalg.norm(s['surface_traction'],axis=1);nv=s['surface_normal'];vf=s['surface_velocity'];vslip=np.linalg.norm(s['surface_outer_velocity']-vf+nv*np.sum(vf*nv,axis=1)[:,None],axis=1)
        xyz=s['x']*1e6
        for col in range(3):
            ax=axes3(fig,231+row*3+col,NAMES[mode]+f" · {float(s['time'])*1e3:.3f} ms\n"+['Facet shear (Pa)','Relative tangential speed (m/s)','Original-ID node force (N)'][col],xyz)
            ax.scatter(*xyz[:n].T,s=2,c='#b7c5ce',alpha=.15,depthshade=False)
            if col<2:
                val=tau if col==0 else vslip;sc=ax.scatter(*pos[active].T,c=val[active],s=7,cmap='viridis',depthshade=False)
                if col==0:
                    ids=np.flatnonzero(active&(tau>tau.max()*.05));ids=ids[::max(1,len(ids)//70)]
                    vec=s['surface_traction'][ids];vec=vec/np.maximum(np.linalg.norm(vec,axis=1)[:,None],1e-30)*60;ax.quiver(*pos[ids].T,*vec.T,color='#152f3c',linewidth=.9,arrow_length_ratio=.3)
            else:
                forces=np.linalg.norm(s['external_force'][:n],axis=1);sel=forces>forces.max()*.025;sc=ax.scatter(*xyz[:n][sel].T,c=forces[sel],s=10,cmap='plasma',depthshade=False);pick=np.flatnonzero(sel)[::max(1,int(sel.sum())//50)];vec=s['external_force'][pick]*50/2e-6;ax.quiver(*xyz[pick].T,*vec.T,color='#772b35',linewidth=.9,arrow_length_ratio=.3)
            fig.colorbar(sc,ax=ax,orientation='horizontal',fraction=.025,pad=.22,shrink=.7)
    fig.subplots_adjust(wspace=.14,hspace=.36,top=.86,bottom=.08);fig.suptitle('Actual accepted-state traction · current surface mapping · node-force arrows: 50 µm per 2 µN',fontsize=16,y=.98);save(fig,'Fig3_Fluid_traction')
    # 4. Six event frames per case; preserve absent event labels.
    events=['initial','first_softening','first_failure','first_multinode','max_multifragment','final'];titles=['Initial','First softening','First failure','First free block ≥2','Maximum target count','Final']
    for mode in CASES:
        fig=plt.figure(figsize=(18,10));limits=(np.min(np.concatenate([s['x'] for k,s in states[mode].items() if k in events]),0)*1e6,np.max(np.concatenate([s['x'] for k,s in states[mode].items() if k in events]),0)*1e6)
        for j,(key,title) in enumerate(zip(events,titles)):
            if key=='final' and key not in states[mode]:key='last_accepted'
            s=states[mode].get(key)
            if key=='max_multifragment' and data[mode]['native']['max_medium_components']==0:s=None
            ax=axes3(fig,231+j,title+(' · NOT_OBSERVED' if s is None else f"\n{float(s['time'])*1e3:.3f} ms"),ref['X']*1e6,limits)
            if s is not None:
                cloud(ax,s,ref,size=5)
                if key in ['first_softening','first_failure']:
                    event=data[mode]['native']['events'].get(key,{})
                    marked=np.unique(np.array(event.get('atom_ids',[]),int))-1
                    if len(marked):ax.scatter(*s['x'][marked].T*1e6,s=65,facecolors='none',edgecolors='#e94e3d',linewidths=1.5,depthshade=False)
        fig.suptitle(NAMES[mode]+' · actual damage evolution · displacement ×1',fontsize=17);fig.subplots_adjust(wspace=.18,hspace=.25,bottom=.13);cax=fig.add_axes([.35,.055,.3,.018]);fig.colorbar(ScalarMappable(norm=Normalize(0,1.),cmap=DAMAGE_CMAP),cax=cax,orientation='horizontal',label='Mean incident native damage (initial precut excluded)');save(fig,'Fig4_Damage_'+mode)
    # 5. Common camera, identical global extents, true component colors.
    finals={m:states[m].get('final',states[m].get('last_accepted')) for m in CASES};xyz=np.concatenate([s['x'] for s in finals.values()])*1e6;limits=(xyz.min(0),xyz.max(0))
    fig=plt.figure(figsize=(16,8))
    for j,mode in enumerate(CASES):
        q=data[mode];ax=axes3(fig,121+j,NAMES[mode]+f"\nFree sizes: {q['final_free_sizes'][:8]}"+(' …' if len(q['final_free_sizes'])>8 else ''),xyz,limits);cloud(ax,finals[mode],ref,'components',9)
    fig.suptitle('Actual material connected components · common camera and bounds · displacement ×1',fontsize=16);fig.text(.1,.03,'Gray: attached material / fixed anchors. Colored: true free components. No visual separation or proximity grouping.');save(fig,'Fig5_True_fragments')
    detail_figure(ref,finals)
    # 6. Six required bins, count and reference-mass fractions.
    import csv
    bins=list(csv.DictReader((O/'FRAGMENT_SIZE_DISTRIBUTION.csv').open()));fig,axs=plt.subplots(2,2,figsize=(14,8),layout='constrained')
    for j,m in enumerate(CASES):
        q=[r for r in bins if r['case']==m];labels=[r['bin'] for r in q];count=[int(r['count']) for r in q];mass=[100*float(r['fraction_initial']) for r in q]
        for ax,val,title in [(axs[0,j],count,'Free component count'),(axs[1,j],mass,'Mass / initial clot (%)')]:
            ax.bar(labels,val,color=['#a7b1ba','#779cb4','#65aaab','#2b9a75','#166f60','#a85468']);ax.set_title(NAMES[m]+'\n'+title);ax.set_xlabel('Material nodes per component');ax.set_ylim(0,max(1,max(val)*1.2));ax.grid(axis='x',visible=False)
            if title=='Free component count':ax.yaxis.set_major_locator(MaxNLocator(integer=True))
            for k,v in enumerate(val):ax.text(k,v,str(round(v,3)),ha='center',va='bottom',fontsize=10)
    save(fig,'Fig6_Fragment_sizes')
    # 7. Work channels and full energy accounting including signed residual.
    fig,axs=plt.subplots(2,2,figsize=(15,9),layout='constrained')
    for j,m in enumerate(CASES):
        h=hist[m];t=np.array([r['time_s'] for r in h])*1e3;w=data[m]['work']
        for key,label,color in [('W_drive_J','Drive','#14669a'),('W_drag_J','Passive drag','#a34352'),('W_J','Net','#161b22')]:axs[0,j].plot(t,[r[key]*1e9 for r in h],label=label,color=color,lw=1.8)
        axs[0,j].axhline(W_TARGET*1e9,color='#d88425',ls='--',label='Target 25.698197 nJ');axs[0,j].set_title(NAMES[m]+f"\nW actual = {w['actual_net_J']*1e9:.4f} nJ; error = {w['relative_error']*100:+.2f}%");axs[0,j].legend(loc='best');axs[0,j].set_ylabel('Actual work (nJ)')
        for key,label in [('U_J','Stored U'),('K_J','Kinetic K'),('D_J','Damage D'),('Uc_J','Contact U'),('R_J','Signed residual R')]:axs[1,j].plot(t,[r[key]*1e9 for r in h],label=label,lw=1.5)
        axs[1,j].legend(ncol=2);axs[1,j].set_ylabel('Energy (nJ)')
        for ax in axs[:,j]:ax.axvline(OFF*1e3,color='#777',ls=':',lw=1);ax.set_xlabel('Time (ms)')
    save(fig,'Fig7_Work_and_energy')
    # 8. Direct comparison with explicit percent definitions.
    fig,axs=plt.subplots(3,3,figsize=(15,11),layout='constrained')
    measures=[('Target count (20–200)',lambda q:q['objective']['n_fragments_20_200']),('Target mass / initial (%)',lambda q:100*q['objective']['mass_fraction_target_initial']),('Largest free block (nodes)',lambda q:max(q['final_free_sizes'],default=0)),('Fines / detached (%)',lambda q:None if q['objective']['fine_debris_mass_fraction_detached'] is None else 100*q['objective']['fine_debris_mass_fraction_detached']),('Target cohort survivors',lambda q:q['objective']['cohort_target_survivors_to_end']),('New accessible faces (final)',lambda q:q['surface']['final']['accessible_new_faces']),('Net work (nJ)',lambda q:q['work']['actual_net_J']*1e9),('Peak actual shear (Pa)',lambda q:q['work']['actual_peak_shear_Pa']),('Max |R| / scale (%)',lambda q:100*q['native']['max_relative_energy_residual'])]
    for index,(ax,(title,fun)) in enumerate(zip(axs.ravel(),measures)):
        values=[fun(data[m]) for m in CASES];v=[0 if a is None else a for a in values];ax.bar(['B pair','C checker'],v,color=['#447caa','#3c9d8e']);ax.set_title(title);ax.set_ylim(0,max(1,max(v)*1.25));ax.grid(axis='x',visible=False)
        if index in [0,2,4,5]:ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        if index==3:ax.set_ylim(0,110)
        for j,a in enumerate(values):ax.text(j,v[j],'UNDEFINED' if a is None else f'{a:.4g}',ha='center',va='bottom')
    save(fig,'Fig8_Comparison')
    costs=dict(wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024)
    previous=json.loads((O/'FIGURE_RESOURCES.json').read_text()) if (O/'FIGURE_RESOURCES.json').exists() else None
    runs=(previous.get('runs',[previous]) if previous else [])+[costs]
    dump(O/'FIGURE_RESOURCES.json',dict(runs=runs,wall_s_sum=sum(q['wall_s'] for q in runs),cpu_s_sum=sum(q['cpu_s'] for q in runs),peak_RSS_bytes=max(q['peak_RSS_bytes'] for q in runs)));print('FIGURES_COMPLETE',costs,flush=True)

if __name__=='__main__':
    import sys
    if '--detail-only' in sys.argv:
        start=time.perf_counter();cpu=time.process_time();ref=load_npz(O/'provenance/reference.npz')
        finals={m:load_npz(O/'cases'/m/'keys/final.npz') for m in CASES}
        detail_figure(ref,finals)
        dump(O/'FIGURE_DETAIL_FIX_RESOURCES.json',dict(wall_s=time.perf_counter()-start,cpu_s=time.process_time()-cpu,peak_RSS_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss*1024))
    else:main()
