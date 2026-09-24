"""Measured comparison figures. Missing P2 results are labelled, never imputed."""
from common import *
from p2_audit import weighted_quantiles
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

C1,C2='#225ea8','#d95f0e'
LABELS=['P1/P1 + VMS','P2/P1 Taylor–Hood']

def save(fig,name):
    fig.savefig(REPORT/'figures'/name,dpi=240,bbox_inches='tight');plt.close(fig)

def readcsv(name):return list(csv.DictReader((REPORT/'data'/name).open()))

def missing(ax):ax.text(.5,.5,'P2 steady solution unavailable\nNo effect claim can be made',ha='center',va='center',transform=ax.transAxes,color='#92400e',bbox=dict(facecolor='#fffbeb',edgecolor='#fcd34d',pad=8))

def main():
    plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
    p1=json.loads((REPORT/'data/p1_network_bc_local_conservation_baseline.json').read_text())
    p2path=REPORT/'data/p2_local_conservation_audit.json';p2=json.loads(p2path.read_text()) if p2path.exists() else None
    # Figure 02: same outlet definitions, no forced equality.
    fig,ax=plt.subplots(figsize=(8,4.5),layout='constrained');x=np.arange(3);names=['OUTLET_01','OUTLET_02','OUTLET_03']
    vals=[p1['outlet_fractions'][k]*100 for k in names];bars=ax.bar(x-.18,vals,.36,label=LABELS[0],color=C1);ax.bar_label(bars,fmt='%.4f%%',padding=4,fontsize=9)
    if p2:
        bars=ax.bar(x+.18,[p2['global_flow']['outlet_fractions'][k]*100 for k in names],.36,label=LABELS[1],color=C2);ax.bar_label(bars,fmt='%.4f%%',padding=4,fontsize=9)
    else:ax.text(.99,.97,'P2: unavailable',ha='right',va='top',transform=ax.transAxes,color=C2)
    ax.set(xticks=x,xticklabels=['O1','O2','O3'],ylabel='Outlet flow / inlet flow (%)',ylim=(0,62),title='Same network-derived outlet pressures');ax.legend(loc='upper left');save(fig,'02_global_flow_split.png')
    # Figure 03 includes all 612 candidates in the inset; only the 599 pre-existing
    # geometrically eligible sections define the root error metrics.
    one=readcsv('p1_root_section_flux.csv');two=readcsv('p2_root_section_flux.csv') if p2 else None
    fig,ax=plt.subplots(figsize=(11,5.5),layout='constrained');inset=ax.inset_axes([.59,.13,.38,.36])
    for rows,label,color in [(one,LABELS[0],C1),(two,LABELS[1],C2)]:
        if rows is None:continue
        rows=sorted(rows,key=lambda r:float(r['section_s_m']));s=np.array([float(r['section_s_m'])*1e6 for r in rows]);y=np.array([float(r['Q_over_Qin']) for r in rows]);ok=np.array([r['original_geometry_pass']=='True' for r in rows])
        ax.plot(s[ok],y[ok],color=color,lw=1.7,label=label);inset.plot(s,y,color=color,lw=1)
        inset.scatter(s[~ok],y[~ok],s=9,marker='x',color='#6b7280')
    ax.axhline(1,color='#374151',ls='--',lw=1,label='Q / Qin = 1');inset.axhline(1,color='#374151',ls='--',lw=.7)
    ax.set(xlabel='Root arclength from inlet (µm)',ylabel='Signed section flow / inlet flow',title='Internal section flux: identical physical sections')
    ax.grid(alpha=.18);ax.legend(loc='upper center',ncol=3,fontsize=10);inset.set_title('All 612 candidates (×: geometry excluded)',fontsize=8);inset.tick_params(labelsize=8)
    if p2 is None:ax.text(.02,.08,'P2 unavailable: baseline curve only',transform=ax.transAxes,color=C2)
    fig.text(.01,-.035,'Main curve: 599 sections with a wall-only boundary. Inset also shows 13 open-cap/intersecting cuts; no internal flux tolerance is applied.',fontsize=9)
    save(fig,'03_root_section_flux_comparison.png')
    # Figure 04: the same element-defined root bins, volume-weighted samples.
    geom=np.load(REPORT/'data/p1_geometry_diagnostics.npz');root=json.loads((OLD/'data/root_topology.json').read_text());edges=np.linspace(0,root['arclength_m'][-1],21)
    rows=[]
    for left,right in zip(edges[:-1],edges[1:]):
        mask=geom['root_mask']&(geom['root_arclength']>=left)&(geom['root_arclength']<right)
        if not mask.any():continue
        vals=geom['divergence'][mask];weights=geom['volumes'][mask];quant=weighted_quantiles(vals,weights,[.01,.1,.9,.99]);rows.append(dict(field='P1',s_left_m=left,s_right_m=right,mean=float(vals@weights/weights.sum()),**dict(zip(['P01','P10','P90','P99'],map(float,quant)))))
    if p2:rows=readcsv('root_divergence_bins.csv')
    else:csvwrite(REPORT/'data/root_divergence_bins_P1_only.csv',rows)
    fig,axes=plt.subplots(2,1,figsize=(10,6),sharex=True,sharey=True,layout='constrained')
    for ax,field,color,label in zip(axes,['P1','P2'],[C1,C2],LABELS):
        rr=[r for r in rows if r['field']==field]
        if rr:
            xx=np.array([(float(r['s_left_m'])+float(r['s_right_m']))*.5e6 for r in rr]);get=lambda k:np.array([float(r[k]) for r in rr])
            ax.fill_between(xx,get('P01'),get('P99'),color=color,alpha=.12,label='P01–P99');ax.fill_between(xx,get('P10'),get('P90'),color=color,alpha=.28,label='P10–P90');ax.plot(xx,get('mean'),color=color,label='Volume-weighted mean')
        else:missing(ax)
        ax.axhline(0,color='#6b7280',lw=.6);ax.set(ylabel='div(u) (s⁻¹)',title=label);ax.grid(alpha=.15)
    axes[0].legend(loc='upper right',ncol=3,fontsize=9);axes[-1].set_xlabel('Root arclength bin (µm)');save(fig,'04_root_divergence_comparison.png')
    fig,axes=plt.subplots(1,2,figsize=(11,4.3),sharey=True,layout='constrained')
    for ax,rows,label in zip(axes,[p1['six_slabs'],p2['six_slabs'] if p2 else None],LABELS):
        if rows:
            s=[v['section_s_m']*1e6 for v in rows];ax.plot(s,[v['loss_minus_wall_m3_s']*1e15 for v in rows],'-o',color=C1,label='Qin − Qsection − Qwall')
            ax.scatter(s,[v['negative_integral_div_m3_s']*1e15 for v in rows],marker='x',s=65,color=C2,label='−∫ div(u) dV',zorder=4)
        else:missing(ax)
        ax.set(xlabel='Slab terminal arclength (µm)',title=label);ax.grid(alpha=.2)
    axes[0].set_ylabel('Flux deficit (10⁻¹⁵ m³/s)');axes[0].legend(fontsize=9);save(fig,'05_gauss_slab_balance.png')
    metrics=[('Volume-weighted RMS divergence','s⁻¹',p1['divergence']['domain']['volume_weighted_RMS'],p2['divergence']['domain']['volume_weighted_RMS'] if p2 else None),
        ('Max eligible root section error','%',p1['root_section_max_error']*100,p2['root_section_max_error']*100 if p2 else None),
        ('RMS eligible root section error','%',p1['root_section_RMS_error']*100,p2['root_section_RMS_error']*100 if p2 else None),
        ('Global boundary mass error','relative',p1['global_mass_error'],p2['global_flow']['epsilon_mass'] if p2 else None)]
    fig,axes=plt.subplots(2,2,figsize=(9,6.5),layout='constrained')
    for ax,(title,unit,one,two) in zip(axes.ravel(),metrics):
        ax.bar([0],[one],color=C1,width=.55);ax.annotate(f'{one:.5g}',(0,one),xytext=(0,4),textcoords='offset points',ha='center',fontsize=9)
        if two is not None:ax.bar([1],[two],color=C2,width=.55);ax.annotate(f'{two:.5g}',(1,two),xytext=(0,4),textcoords='offset points',ha='center',fontsize=9)
        else:ax.text(1,one*.4,'N/A',ha='center',color=C2)
        ax.set(xticks=[0,1],xticklabels=['P1','P2'],ylabel=unit,title=title,xlim=(-.6,1.6));ax.margins(y=.2)
    save(fig,'06_local_conservation_summary.png')
    cost1=json.loads((REPORT/'data/p1_resource_metrics.json').read_text());est=json.loads((REPORT/'data/taylor_hood_resource_estimate.json').read_text())
    costpath=REPORT/'data/p2_resource_metrics.json';cost2=json.loads(costpath.read_text()) if costpath.exists() else None
    a=np.array([est['original_vertices'],est['P1_total_dof'],cost1['wall_time_s'],cost1['peak_tree_RSS_MiB_sampled'],cost1['peak_GPU_memory_MiB']])
    b=[est['tet10_nodes'],est['P2P1_total_dof'],cost2['wall_time_s'] if cost2 else None,cost2['peak_tree_RSS_MiB_sampled'] if cost2 else None,cost2['peak_GPU_memory_MiB'] if cost2 else None]
    fig,ax=plt.subplots(figsize=(10,4.8),layout='constrained');xx=np.arange(5)
    ax.bar(xx-.18,np.ones(5),.36,color=C1,label='P1 = 1');valid=[i for i,v in enumerate(b) if v is not None]
    ax.bar(xx[valid]+.18,[b[i]/a[i] for i in valid],.36,color=C2,label='P2 / P1')
    for i in range(5):
        if b[i] is None:ax.text(i+.18,.15,'N/A',ha='center',color=C2)
        else:ax.text(i+.18,b[i]/a[i]+.12,f'{b[i]/a[i]:.2f}×',ha='center',fontsize=10)
    ax.set(xticks=xx,xticklabels=['Nodes','Active DOF','Full-run time','Peak host RSS','Peak GPU memory'],ylabel='Ratio to P1 baseline',title='Cost of order elevation');ax.legend();ax.margins(y=.17)
    if cost2 is None:fig.text(.01,-.03,'Full P2 solve not completed. Smoke time and memory are documented separately and are not substituted for full-run costs.',fontsize=9)
    save(fig,'07_resource_cost_comparison.png')

if __name__=='__main__':main()
