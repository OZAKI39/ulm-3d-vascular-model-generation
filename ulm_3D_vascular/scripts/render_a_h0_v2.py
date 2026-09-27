"""Reproducible standalone scientific figures from saved H0/FEM evidence."""
from pathlib import Path
import json,csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize,LogNorm
from matplotlib.cm import ScalarMappable
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from matplotlib.patches import FancyBboxPatch,FancyArrowPatch

ROOT=Path(__file__).resolve().parents[1]
REPORT=ROOT/'reports/a_network_1d0d_boundary_v2_idealized'
DATA=REPORT/'data';FIG=REPORT/'figures'
COLORS=['#277da8','#dd7533','#48946a']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                     'axes.spines.right':False,'savefig.facecolor':'white','pdf.fonttype':42})

def read(name):return json.loads((DATA/name).read_text())
def save(fig,name):
    for suffix in ['png','pdf']:fig.savefig(FIG/(name+'.'+suffix),dpi=230,bbox_inches='tight')
    plt.close(fig)

def bars(groups,labels,name,title,footnote):
    fig,ax=plt.subplots(figsize=(11,6.2));x=np.arange(3);width=.76/len(groups)
    palette=['#9ba9b8','#398bb2','#d97837'][:len(groups)]
    for j,(g,label,color) in enumerate(zip(groups,labels,palette)):
        rects=ax.bar(x+(j-(len(groups)-1)/2)*width,100*np.array(g),width,label=label,color=color)
        ax.bar_label(rects,fmt='%.2f%%',padding=4,fontsize=11)
    ax.set_xticks(x,['O1','O2','O3']);ax.set_ylabel('Signed outlet flow / ROI inlet flow (%)')
    low=min(0,min(float(np.min(g))*100 for g in groups));high=max(100,max(float(np.max(g))*100 for g in groups)+15)
    ax.set_ylim(low-5 if low<0 else 0,high+5);ax.axhline(0,color='#465466',lw=.8)
    ax.grid(axis='y',alpha=.16);ax.set_axisbelow(True);ax.legend(loc='upper center',ncol=len(groups),fontsize=10,frameon=False)
    ax.set_title(title,pad=22,weight='bold');fig.text(.08,.015,footnote,fontsize=10,color='#425268')
    fig.subplots_adjust(bottom=.16,top=.86);save(fig,name)

def network_map(graph,sol,mode):
    xyz=graph['xyz_m']*1e6;e=graph['edges'];roi=graph['roi_internal_edge_mask'];ports=graph['port_indices']
    if mode=='pressure':
        values=sol['pressure_Pa'][e].mean(1)/1000;norm=Normalize(0,float(sol['pressure_Pa'].max()/1000));cmap='viridis';label='Gauge pressure (kPa)'
    else:
        values=np.abs(sol['edge_Q_m3s'])*1e15
        norm=LogNorm(max(1e-6,float(values[values>0].min())),float(values.max()));cmap='plasma';label='Edge flow magnitude (pL/s, logarithmic scale)'
    fig=plt.figure(figsize=(14,7.5));ax=fig.add_subplot(121,projection='3d');zoom=fig.add_subplot(122,projection='3d')
    for a,mask in [(ax,np.ones(len(e),bool)),(zoom,roi)]:
        lc=Line3DCollection(xyz[e[mask]],array=values[mask],norm=norm,cmap=cmap,linewidths=.85 if a is ax else 2.6,alpha=.95)
        a.add_collection3d(lc);coords=xyz if a is ax else xyz[np.unique(e[roi])]
        margin=np.ptp(coords,axis=0)*.1
        for setter,lo,hi,m in zip([a.set_xlim,a.set_ylim,a.set_zlim],coords.min(0),coords.max(0),margin):setter(lo-m,hi+m)
        a.set_box_aspect(np.maximum(np.ptp(coords,axis=0),1));a.view_init(elev=23,azim=-61)
        a.set_xlabel('x (um)',labelpad=8);a.set_ylabel('y (um)',labelpad=8);a.set_zlabel('z (um)',labelpad=8)
        a.tick_params(labelsize=8);a.grid(alpha=.15)
    ax.add_collection3d(Line3DCollection(xyz[e[roi]],colors='#dc407d',linewidths=2.2,alpha=.85))
    t=graph['terminal_indices'];ax.scatter(*xyz[t].T,c='#526273',s=8,alpha=.65,label='122 p_ref terminals',depthshade=False)
    source=int(graph['source_index']);ax.scatter(*xyz[source],marker='*',s=190,c='black',label='Idealized source 2410',depthshade=False)
    ax.text(*xyz[source], '  2410',fontsize=10,weight='bold');ax.legend(loc='upper left',frameon=False,fontsize=9)
    for j,name in enumerate(['IN','O1','O2','O3']):
        pos=xyz[ports[j]];zoom.scatter(*pos,s=65,c='black' if j==0 else COLORS[j-1],depthshade=False)
        zoom.text(*(pos+np.array([1.2,1.2,1.2])),name,fontsize=12,weight='bold')
    if mode=='flow':
        roi_edges=np.flatnonzero(roi)[::14]
        for idx in roi_edges:
            a,b=xyz[e[idx]];direction=(b-a)*np.sign(sol['edge_Q_m3s'][idx]);direction/=np.linalg.norm(direction)
            zoom.quiver(*((a+b)/2),*direction,length=3.5,color='#273444',arrow_length_ratio=.5,linewidth=1.2)
    ax.set_title('Analysis A only: 7,419 original nodes',fontsize=12)
    zoom.set_title('Retained ROI and four real ports',fontsize=12)
    fig.suptitle('Idealized H0 model | '+('pressure field' if mode=='pressure' else 'flow magnitude and ROI direction'),fontsize=17,weight='bold',y=.97)
    cb=fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),ax=[ax,zoom],orientation='horizontal',fraction=.04,pad=.12,shrink=.65)
    cb.set_label(label)
    fig.text(.05,.012,'Structural source 2410; all 122 other leaves at p_ref = 0. ROI inlet matched to 15.5136 pL/s.\nMagenta overview corridor = ROI. Three virtual cut nodes are added for the solve. Other 42 components are excluded.',fontsize=10)
    save(fig,'02_analysis_A_pressure' if mode=='pressure' else '03_analysis_A_flow')

def run():
    FIG.mkdir(exist_ok=True);s=read('final_summary.json');mass=s['network_mass_residual']
    graph=np.load(DATA/'analysis_A_H0_graph_si.npz');sol=np.load(DATA/'H0_solution_si.npz')
    residual=abs(sol['node_outflow_m3s'][sol['internal_node_indices']])/s['A_total_source_inflow']
    fig,axes=plt.subplots(1,2,figsize=(12,5))
    nonzero=residual[residual>0];axes[0].hist(np.log10(nonzero),bins=32,color='#398bb2',edgecolor='white')
    axes[0].set_xlabel('log10(|internal node imbalance| / A source flow)');axes[0].set_ylabel('Number of internal nodes')
    axes[0].set_title('Sparse network continuity residuals')
    axes[0].text(.04,.96,f'{len(residual)} internal nodes\n{np.count_nonzero(residual==0)} exact floating-point zeros\nNonzero residuals retained',transform=axes[0].transAxes,va='top',fontsize=10)
    vals=[mass['max_relative_nodal_imbalance_to_A_source'],mass['global_relative_residual'],mass['ROI_relative_residual']]
    axes[1].bar(['Max node','Full A','ROI'],vals,color=['#398bb2','#48946a','#dd7533'])
    axes[1].set_yscale('log');axes[1].set_ylim(1e-16,1e-8);axes[1].axhline(1e-9,ls='--',color='#aa4455',label='Network gate: 1e-9')
    for i,v in enumerate(vals):axes[1].text(i,v*1.35,f'{v:.2e}',ha='center',fontsize=10)
    axes[1].set_title('Relative continuity checks');axes[1].set_ylabel('Relative residual');axes[1].legend(frameon=False,fontsize=9)
    fig.suptitle('Idealized H0 model | mass balance',fontsize=16,weight='bold')
    fig.text(.075,.015,f'A source = {s["A_total_source_inflow"]*1e15:.4f} pL/s; ROI inlet = {s["ROI_inlet_Q"]*1e15:.4f} pL/s. These are distinct flows.',fontsize=10)
    fig.tight_layout(rect=[0,.07,1,.94]);save(fig,'01_network_mass_balance')
    network_map(graph,sol,'pressure');network_map(graph,sol,'flow')
    bars([s['current_3D_split'],s['H0_split']],['Current ROI 3D: all outlet p = 0','Full analysis-A: idealized H0'],
         '04_roi_split_current_vs_H0','The surrounding network changes the predicted flow split',
         'Different model classes: this comparison alone does not isolate the 3D boundary effect. All H0 ports flow outward.')
    loo=list(csv.DictReader((DATA/'terminal_leave_one_out_sensitivity.csv').open()));valid=[r for r in loo if r['status']=='PASS']
    fig,axes=plt.subplots(1,2,figsize=(13,5.8),gridspec_kw={'width_ratios':[1.5,1]})
    for j,name in enumerate(['O1','O2','O3']):
        delta=np.array([float(r[name+'_difference_pp']) for r in valid]);axes[0].scatter(range(1,len(valid)+1),delta,s=21,label=name,c=COLORS[j],alpha=.8)
    axes[0].set_xlabel('Removed terminal index (all 122 tested, in saved CSV order)');axes[0].set_ylabel('Change from H0 (percentage points)');axes[0].axhline(0,c='#586779',lw=.7)
    axes[0].set_title('Most terminals have almost no ROI effect');axes[0].legend(ncol=3,frameon=False)
    important=sorted(valid,key=lambda r:max(abs(float(r[n+'_difference_pp'])) for n in ['O1','O2','O3']),reverse=True)[:3]
    x=np.arange(3)
    for j,name in enumerate(['O1','O2','O3']):
        vals=[float(r[name+'_fraction'])*100 for r in important];rects=axes[1].bar(x+(j-1)*.25,vals,.25,color=COLORS[j],label=name)
        for rect,value in zip(rects,vals):
            label='~0' if abs(value)<1e-8 else f'{value:.1f}'
            axes[1].text(rect.get_x()+rect.get_width()/2,max(value,0)+1.5,label,ha='center',fontsize=8)
    axes[1].set_xticks(x,['Remove\n'+r['removed_terminal_id'] for r in important]);axes[1].set_ylabel('Outlet / ROI inlet (%)');axes[1].set_ylim(0,100)
    axes[1].set_title('A few local sinks strongly control the split')
    fig.suptitle('Terminal sensitivity | one reference-pressure sink removed at a time',fontsize=15,weight='bold')
    fig.text(.07,.018,'Each variant is solved again and matched to the same ROI inlet flow. 122 valid variants; no 3D sensitivity runs.\nTiny signed roundoff values are retained in CSV; they are not interpreted as resolved flow reversal.',fontsize=10)
    fig.tight_layout(rect=[0,.09,1,.93]);save(fig,'05_terminal_sensitivity')
    if s['new_3D_case_status']=='PASS':
        bars([s['current_3D_split'],s['H0_split'],s['new_3D_split']],
             ['Current 3D\nall outlets p = 0','Full-A H0\n1D/0D network','New 3D\nA-derived outlet pressures'],
             '06_three_way_flow_split','Outlet pressure changes the actual 3D flow split',
             'Current and new 3D use the same mesh, inlet, fluid, wall and numerical settings. Only outlet pressure values change.\nThe 1D/0D prediction is an idealized approximation; exact agreement with 3D is not an acceptance gate.')
    transfer=read('network_to_fem_pressure_transfer.json')['ports']
    fig,ax=plt.subplots(figsize=(14,7));ax.axis('off');ax.set_xlim(0,1);ax.set_ylim(0,1)
    ax.text(.5,.97,'Why are the new outlet pressures different?',ha='center',va='top',fontsize=21,weight='bold')
    boxes=[(.03,.67,.26,.18,'1  Analysis A','Keep the full saved network\nSource 2410; 122 reference terminals'),
           (.37,.67,.26,.18,'2  H0 pressure field','Solve the connected pipe network\nMatch the ROI inlet to 15.5136 pL/s'),
           (.71,.67,.26,.18,'3  Pressures at real cuts','O1  342.34 Pa\nO2  3065.61 Pa     O3  0 Pa'),
           (.53,.27,.40,.23,'4  Correct for added FEM tubes','Subtract the pressure lost along extensions\nO1  52.63 Pa   |   O2  429.17 Pa   |   O3  295.57 Pa\nRaw cap pressures: 289.71 / 2636.44 / -295.57 Pa'),
           (.04,.27,.40,.23,'5  Set the three FEM cap pressures','Add the same 295.57 Pa to all three\nO1  585.29 Pa   |   O2  2932.02 Pa   |   O3  0 Pa\nAll pressure differences stay unchanged')]
    for x,y,w,h,title,body in boxes:
        ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=.014',facecolor='#edf4f8',edgecolor='#6c96ad',lw=1.5))
        ax.text(x+w/2,y+h-.035,title,ha='center',va='top',fontsize=12,weight='bold')
        ax.text(x+w/2,y+h-.085,body,ha='center',va='top',fontsize=10.5,linespacing=1.65)
    for a,b in [((.30,.76),(.35,.76)),((.64,.76),(.69,.76)),((.84,.65),(.77,.52)),((.51,.385),(.46,.385))]:
        ax.add_patch(FancyArrowPatch(a,b,arrowstyle='-|>',mutation_scale=22,color='#376680',lw=2))
    ax.text(.5,.11,'Different saved downstream paths produce different pressures at the ROI cuts.\nO3 is itself a model terminal: its real-cut pressure is defined as zero; no downstream geometric R3 is invented.',ha='center',fontsize=12,linespacing=1.7)
    ax.text(.5,.025,'All pressures are model gauge pressures, not measured physiological mouse pressures.',ha='center',fontsize=11,color='#8d482d')
    save(fig,'07_pressure_boundary_explanation')
    print('Rendered',len(list(FIG.glob('*.png'))),'PNG figures and matching PDFs')

if __name__=='__main__':run()
