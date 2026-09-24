"""Evidence-backed scientific figures; physical inlet coordinates and SI data."""
from pathlib import Path
import argparse,json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Circle,FancyArrowPatch
from .particle82a_admission import context
from .particle82a_pipeline import BASINS
from .particle82_provenance import atomic_json,sha256

BG='#090f19';INK='#e5edf7';MUTED='#a0acbd'
COLORS=['#e5a16d','#52c8da','#b59ae6','#637083']
METHOD_COLORS=['#55bed0','#dfaf66','#a896e1'];FAIL='#cd8586'
NAMES=['00_stage_overview','01_inlet_point_tracer_basins','02_inlet_max_passable_diameter_map',
 '03_passable_size_vs_sonovue','04_failure_reason_map_current_method','05_method_A_position_bias',
 '06_method_B_size_bias','07_method_C_inward_distance','08_A_B_C_admission_comparison',
 '09_A_B_C_full_trajectory_outlets','10_birth_plane_artifact_diagnosis','11_scientific_limitations']


def style():
    plt.rcParams.update({'figure.facecolor':BG,'axes.facecolor':BG,'savefig.facecolor':BG,
        'text.color':INK,'axes.labelcolor':INK,'xtick.color':MUTED,'ytick.color':MUTED,
        'axes.edgecolor':'#526075','font.family':'DejaVu Sans','font.size':12,
        'axes.titleweight':'bold','axes.grid':False,'legend.facecolor':BG,'legend.edgecolor':'#526075',
        'figure.dpi':120,'savefig.dpi':200})


def canvas(title,subtitle,cols=1):
    fig,axes=plt.subplots(1,cols,figsize=(14,8),squeeze=False)
    fig.subplots_adjust(left=.075,right=.92,bottom=.15,top=.76,wspace=.38)
    fig.text(.06,.94,title,fontsize=23,weight='bold')
    fig.text(.06,.87,subtitle,fontsize=12,color=MUTED)
    fig.text(.06,.055,'Particle-8.2A  |  Frozen FEM + original SonoVue input  |  diagnostic comparison',color=MUTED,fontsize=10)
    return fig,axes[0]


def tidy(ax):
    ax.spines[['top','right']].set_visible(False)


def outline(ax,geometry):
    for edge in geometry.perimeter:
        p=geometry.coordinates(edge)*1e6;ax.plot(p[:,0],p[:,1],color=INK,lw=.7,alpha=.55)
    ax.set_aspect('equal');ax.set_xlabel('Inlet u (µm)');ax.set_ylabel('Inlet v (µm)')


def save(fig,output,index,sources,extra=None):
    output=Path(output);(output/'figures').mkdir(parents=True,exist_ok=True);(output/'figure_sources').mkdir(exist_ok=True)
    name=NAMES[index];fig.savefig(output/'figures'/(name+'.png'));fig.savefig(output/'figures'/(name+'.pdf'));plt.close(fig)
    atomic_json(output/'figure_sources'/(name+'.json'),dict(source_files={str(p):sha256(p) for p in sources},
        geometry_coordinates_rotated=False,display_coordinates='ORTHONORMAL_INLET_UV_MICROMETERS',
        renderer_source_sha256=sha256(Path(__file__)),**(extra or {})))


def plot_all(admission,report,output,allow_partial=False):
    style();admission=Path(admission);report=Path(report);output=Path(output)
    c=context();g=c.geometry
    summary=json.loads((report/'ADMISSION_COMPARISON.json').read_text());source=report/'ADMISSION_COMPARISON.json'
    maps=np.load(admission/'map_n24.npz')['rows'];uv=g.coordinates(maps[:,:3])*1e6
    data=np.load(report/'admission_plot_data.npz');methods=summary['methods'];n=summary['common_event_count']
    # 00 — show algorithmic operations without suggesting forces or trajectories.
    fig,axes=canvas('Three entry strategies, two geometry questions',
        'A changes trial position   ·   B changes trial size   ·   C searches an inward representation',3)
    for i,(ax,title) in enumerate(zip(axes,['A  Fixed size','B  Fixed anchor','C  Fixed anchor + size'])):
        ax.set_xlim(-1,5);ax.set_ylim(-1,5);ax.set_aspect('equal');ax.axis('off');ax.set_title(title,color=METHOD_COLORS[i],pad=18)
        ax.plot([0,0],[0,4],color=MUTED,lw=2);ax.plot([4,4],[0,4],color=MUTED,lw=2)
        ax.plot([0,4],[0,0],color=INK,ls='--',lw=1);ax.text(2,-.6,'open CFD inlet',ha='center',fontsize=10)
        if i==0:
            for x,y,alpha in [(0.4,0,.3),(3.7,0,.4),(2,0,1)]:ax.add_patch(Circle((x,y),.65,fc=METHOD_COLORS[i],alpha=alpha))
            ax.text(2,3,'Resample position\nSame diameter',ha='center',va='center',linespacing=1.8)
        elif i==1:
            for r,alpha in [(.85,.2),(.6,.35),(.3,.9)]:ax.add_patch(Circle((.5,0),r,fc=METHOD_COLORS[i],alpha=alpha))
            ax.text(2,3,'Resample size\nSame anchor',ha='center',va='center',linespacing=1.8)
        else:
            ax.plot([1.8,1.9,2.2,2.1],[0,.5,1.2,2.3],color=METHOD_COLORS[i],lw=2)
            for x,y,alpha in [(1.8,0,.25),(2.2,1.2,1)]:ax.add_patch(Circle((x,y),.65,fc=METHOD_COLORS[i],alpha=alpha))
            ax.text(2,3.5,'First legal full center\nNo resampling',ha='center',va='center',linespacing=1.8)
    fig.text(.08,.15,'Question 1: Is the aperture locally too small?     Question 2: Does the current checker misuse the open cap?',fontsize=12)
    save(fig,output,0,[source],dict(schematic_not_physical_coordinates=True))
    # 01 — all high-resolution quadrature nodes, no assigned outlet quotas.
    fig,(ax,)=canvas('01  Inlet point-tracer outlet basins','Official inlet surface; unresolved paths remain grey')
    for k,label in enumerate(BASINS):
        mask=maps[:,5]==k;ax.scatter(uv[mask,0],uv[mask,1],s=.8,c=COLORS[k],rasterized=True,label=label.replace('_POINT_PATH',''))
    outline(ax,g);ax.legend(loc='upper left',bbox_to_anchor=(1.02,1),markerscale=5,fontsize=10)
    save(fig,output,1,[admission/'map_n24.npz'])
    fig,(ax,)=canvas('02  Maximum passable diameter at the open inlet','Dmax = 2 × minimum distance to true WALL / aperture perimeter; cap is not solid')
    sc=ax.scatter(uv[:,0],uv[:,1],c=maps[:,8],s=.8,cmap='viridis',vmin=0,vmax=maps[:,8].max(),rasterized=True)
    outline(ax,g);fig.colorbar(sc,ax=ax,pad=.05,label='Dmax (µm)')
    save(fig,output,2,[admission/'map_n24.npz'])
    fig,axes=canvas('03  Local clearance versus original SonoVue sizes','Flux-weighted Dmax distribution; shaded curve is the unmodified number-weighted SonoVue PDF',3)
    x=np.linspace(.001,5.3,1200)
    for k,ax in enumerate(axes):
        mask=maps[:,5]==k;ax.hist(maps[mask,8],bins=np.linspace(0,5.3,70),weights=maps[mask,4],density=True,
            histtype='step',lw=2,color=COLORS[k],label='Local Dmax')
        ax.fill_between(x,c.distribution.pdf(x),color=MUTED,alpha=.3,label='SonoVue diameter')
        ax.axvline(summary['official_sonovue']['quantiles_um']['D50'],ls=':',color=INK,lw=1)
        ax.set_title(BASINS[k],color=COLORS[k]);ax.set_xlabel('Diameter (µm)');tidy(ax)
        if k==0:ax.set_ylabel('Probability density (1/µm)');ax.legend(fontsize=9)
    save(fig,output,3,[admission/'map_n24.npz',source])
    # The retained original A candidates supply actual rejection map points.
    candidate_sample=[]
    from collections import Counter
    counts=Counter()
    for p in sorted((admission/'events').glob('events_*.npz')):
        a=np.load(p)['attempts'];a=a[(a[:,0]==0)&(a[:,8]!=0)]
        counts.update(a[:,8].astype(int).tolist())
        # Deterministic display thinning only, with complete counts preserved.
        candidate_sample.extend(a[::max(1,len(a)//8)][:8].tolist())
    a=np.asarray(candidate_sample);puv=g.coordinates(a[:,4:7])*1e6
    fig,(ax,)=canvas('04  Why current-method trial positions are rejected','Every A retry is retained in data; deterministic spatial display sample shown')
    cause_colors=['#cd8586','#e8b47d','#9d9ac8','#6daeb2','#9d9ac8','#9d9ac8','#9d9ac8','#9d9ac8','#9d9ac8']
    from .particle82a_pipeline import CAUSES
    for idx,k in enumerate(sorted(counts)):
        mask=a[:,8]==k;ax.scatter(puv[mask,0],puv[mask,1],s=3,c=cause_colors[idx],alpha=.5,label=CAUSES[k].replace('_',' '))
    outline(ax,g);ax.legend(loc='upper left',bbox_to_anchor=(1.02,1),fontsize=8,markerscale=3)
    save(fig,output,4,[source],dict(all_A_rejection_counts={CAUSES[k]:v for k,v in counts.items()},display_sample_count=len(a)))
    original=data['uv_um'];accepted=data['A_accepted'];final=g.coordinates(data['A_birth_m'][accepted])*1e6
    edges=[np.linspace(original[:,i].min(),original[:,i].max(),45) for i in range(2)]
    h0=np.histogram2d(*original.T,bins=edges,density=True)[0];h1=np.histogram2d(*final.T,bins=edges,density=True)[0]
    fig,axes=canvas('05  Method A changes the admitted position distribution','Each event keeps its first diameter; rejected positions are redrawn',3)
    for ax,h,title in zip(axes,[h0,h1,h1-h0],['Original anchors','Admitted positions','Admitted − original']):
        delta='−' in title;lim=np.max(np.abs(h)) if delta else max(h0.max(),h1.max())
        im=ax.pcolormesh(edges[0],edges[1],h.T,cmap='coolwarm' if delta else 'magma',vmin=-lim if delta else 0,vmax=lim,rasterized=True)
        outline(ax,g);ax.set_title(title,fontsize=13);fig.colorbar(im,ax=ax,shrink=.75,pad=.04,label='Density (1/µm²)')
    save(fig,output,5,[source,report/'admission_plot_data.npz'])
    fig,axes=canvas('06  Method B conditions sizes on anchor passability','Fixed anchor; accepted diameter distribution is not the original SonoVue distribution',3)
    for k,ax in enumerate(axes):
        mask=(data['anchor_basin']==k)&data['B_accepted'];diam=data['B_diameter_um'][mask]
        ax.plot(x,c.distribution.cdf(x),color=MUTED,lw=2,label='Original SonoVue CDF')
        if len(diam):
            ordered=np.sort(diam);ax.step(ordered,np.arange(1,len(diam)+1)/len(diam),where='post',color=COLORS[k],lw=2,label='B accepted CDF')
        ax.set_title(BASINS[k],color=COLORS[k]);ax.set_xlabel('Diameter (µm)');ax.set_ylim(0,1);tidy(ax)
        if k==0:ax.set_ylabel('Cumulative fraction');ax.legend(fontsize=9,loc='lower right')
    save(fig,output,6,[source,report/'admission_plot_data.npz'])
    fig,axes=canvas('07  Method C searches for the first full center','Fixed anchor and first diameter; inward distance follows the Frozen velocity streamline',2)
    for k in range(3):
        mask=(data['anchor_basin']==k)&data['C_accepted'];distance=data['C_s_birth_m'][mask]*1e6;radius=data['first_diameter_um'][mask]/2
        for ax,values in zip(axes,[distance,distance/radius]):
            if len(values):
                ordered=np.sort(values);ax.plot(ordered,np.arange(1,len(values)+1)/len(values),color=COLORS[k],lw=2,label=BASINS[k])
    axes[0].set_xlabel('First admissible inward distance (µm)');axes[1].set_xlabel('Inward distance / radius')
    for ax in axes:ax.set_ylabel('Cumulative admitted fraction');ax.set_ylim(0,1);ax.legend(fontsize=10);tidy(ax)
    save(fig,output,7,[source,report/'admission_plot_data.npz'])
    fig,axes=canvas('08  Paired A / B / C admission comparison',f'{n:,} identical initial events; accepted populations remain strategy-conditioned',3)
    axes[0].bar(list('ABC'),[methods[k]['acceptance_fraction']*100 for k in 'ABC'],color=METHOD_COLORS);axes[0].set_ylabel('Accepted events (%)');axes[0].set_ylim(0,100)
    bottom=np.zeros(3)
    for j,b in enumerate(BASINS):
        vals=np.array([methods[k]['accepted_birth_point_basin_fraction'][b] for k in 'ABC'])*100
        axes[1].bar(list('ABC'),vals,bottom=bottom,color=COLORS[j],label=b.replace('_POINT_PATH',''));bottom+=vals
    axes[1].set_ylabel('Admitted birth-point basin (%)');axes[1].legend(fontsize=8,loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2)
    rejection_keys=['INLET_PERIMETER_CLEARANCE_FAIL','WALL_INTERSECTION','NEARFIELD_HANDOFF_VIOLATION',
                    'NO_ADMISSIBLE_INWARD_LOCATION','NUMERICAL_GEOMETRY_UNRESOLVED']
    rejection_colors=['#d5a16f','#cf8587','#729fa7','#9285bd','#637083'];bottom=np.zeros(3)
    short_names=['Perimeter','Wall','Handoff','No legal center','Unresolved']
    for key,color,label in zip(rejection_keys,rejection_colors,short_names):
        values=[]
        for m in 'ABC':
            counts=methods[m]['rejection_counts'] if m!='C' else {k:v for k,v in methods[m]['status_counts'].items() if k!='ACCEPTED'}
            values.append(counts.get(key,0)/max(1,sum(counts.values()))*100)
        axes[2].bar(list('ABC'),values,bottom=bottom,color=color,label=label)
        bottom+=values
    axes[2].set_ylabel('Failure mix (%)');axes[2].set_ylim(0,100)
    axes[2].text(.5,-.16,'A/B: rejected trials; C: failed searches',transform=axes[2].transAxes,ha='center',fontsize=8,color=MUTED)
    axes[2].legend(fontsize=8,loc='lower center',bbox_to_anchor=(.5,1.01),ncol=2)
    for ax in axes:tidy(ax)
    save(fig,output,8,[source])
    trajectory_file=report/'FULL_TRAJECTORY_COMPARISON.json'
    if trajectory_file.exists():
        t=json.loads(trajectory_file.read_text());fig,axes=canvas('09  Full-vessel trajectory outcomes','Same P6.5 sphere dynamics and integration settings; numerical stops have unknown later fate',2)
        xpos=np.arange(3)
        for k,b in enumerate(BASINS[:3]):axes[0].bar(xpos+(k-1)*.23,[t[m]['outlet_counts'][b] for m in 'ABC'],.23,color=COLORS[k],label=b)
        axes[0].set_xticks(xpos,list('ABC'));axes[0].set_ylabel('Completed trajectories');axes[0].legend(fontsize=9)
        axes[1].bar(list('ABC'),[t[m]['stop_rate']*100 for m in 'ABC'],color=FAIL);axes[1].set_ylabel('Numerical safety stops (%)');axes[1].set_ylim(0,100)
        for ax in axes:tidy(ax)
        save(fig,output,9,[trajectory_file])
    elif not allow_partial:raise ValueError('Full trajectory comparison missing')
    fig,axes=canvas('10  Aperture geometry versus birth-plane representation','Current checker permits upstream overlap through open cap; full-center condition is separate',3)
    entry_file=report/'METHOD_C_ENTRY_SWEEP_AUDIT.json'
    entry=json.loads(entry_file.read_text())['summary'] if entry_file.exists() else None
    for k,ax in enumerate(axes):
        b=BASINS[k];base=summary['geometry_by_basin'][b]
        values=[base['open_aperture_passable_fraction']*100,base['center_on_plane_acceptance']*100,methods['C']['by_anchor_basin'][b]['acceptance_fraction']*100]
        ax.bar([0,1],values[:2],color=[MUTED,METHOD_COLORS[0]])
        if entry:
            clear=entry[b]['full_sweep_handoff_certified']/base['scheduled']*100
            ax.bar(2,clear,color=METHOD_COLORS[2],label='C: entry WALL certified')
            ax.bar(2,values[2]-clear,bottom=clear,color=FAIL,hatch='///',label='C: entry conflict / unverified')
        else:ax.bar(2,values[2],color=METHOD_COLORS[2])
        ax.set_xticks([0,1,2],['Open\naperture','Current\non-plane','C full\ncenter'])
        for i,v in enumerate(values):ax.text(i,v+.25,f'{v:.2f}%',ha='center',fontsize=10)
        ax.set_ylim(0,max(5,max(values)*1.2));ax.set_title(b,color=COLORS[k]);ax.set_ylabel('Paired first-draw events (%)');tidy(ax)
    if entry:axes[1].legend(loc='lower center',bbox_to_anchor=(.5,1.10),fontsize=8)
    fig.text(.06,.095,'C: a legal birth center alone does not establish a geometrically passable entry path.',fontsize=10,color=MUTED)
    save(fig,output,10,[source]+([entry_file] if entry else []))
    fig,(ax,)=canvas('11  Scientific interpretation and limits','Evidence dimensions are reported separately; no subjective overall ranking')
    ax.axis('off')
    rows=[['A','Fixed first size','Redraws position','Guard also filters admitted sizes'],
          ['B','Fixed anchor','Redraws size','Accepted sizes are conditional'],
          ['C','Fixed anchor + first size','Changes represented birth center','Legal center does not prove entry passage'],
          ['All','Same Frozen FEM + P6.5','Independent single-MB trajectories','No RBC coupling / no physical trapping claim']]
    table=ax.table(cellText=rows,colLabels=['Method','Preserved per event','Changed','Limit'],loc='center',cellLoc='left',colWidths=[.09,.25,.29,.37])
    table.auto_set_font_size(False);table.set_fontsize(12);table.scale(1,3.2)
    for (r,col),cell in table.get_celld().items():
        cell.set_facecolor('#131e2c' if r else '#203247');cell.set_edgecolor('#425066');cell.get_text().set_color(INK)
    fig.text(.09,.19,'Frozen geometry has finite resolution. Point-basin labels include unresolved paths.\nAdmitted-cohort fidelity differs from preserving each event’s original inputs.',fontsize=12,color=MUTED,linespacing=1.8)
    save(fig,output,11,[source])


def main():
    p=argparse.ArgumentParser();p.add_argument('--admission',required=True);p.add_argument('--report',required=True);p.add_argument('--output',required=True)
    p.add_argument('--allow-partial',action='store_true');a=p.parse_args();plot_all(a.admission,a.report,a.output,a.allow_partial)


if __name__=='__main__':main()
