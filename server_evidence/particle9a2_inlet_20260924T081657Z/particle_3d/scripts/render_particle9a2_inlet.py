#!/usr/bin/env python3
"""Academic figures plus plain CSV source tables for the P9-A.2 inlet audit."""
from pathlib import Path
import csv,json,hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from particle_3d.injection_method_c import TruncatedSonoVue
from particle_3d.particle7_cases import SONOVUE

ROOT=Path(__file__).resolve().parents[2];R=ROOT/'particle_3d/reports/particle9a2_inlet_sampling';D=R/'data';F=R/'figures'
ROLES=['OUTLET_01','OUTLET_02','OUTLET_03','NO_EXIT'];COLORS=['#4477AA','#228833','#CC6677','#888888']
def read(p):return json.loads(Path(p).read_text())
def write(p,v):p.write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def csvwrite(p,rows):
    with p.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def style():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,'axes.labelsize':11,
        'figure.dpi':120,'savefig.dpi':300,'figure.facecolor':'white','axes.facecolor':'white',
        'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'ps.fonttype':42})

def main():
    style();F.mkdir(exist_ok=True);manifest=[]
    def save(fig,name,source):
        fig.savefig(F/(name+'.png'),bbox_inches='tight');fig.savefig(F/(name+'.pdf'),bbox_inches='tight');plt.close(fig)
        manifest.append(dict(name=name,sources=source,png_sha256=sha(F/(name+'.png')),pdf_sha256=sha(F/(name+'.pdf'))))
    a=read(R/'audit2000/admission/birth_ledger.json')['events'];b=read(R/'reference/method_b_500_birth_ledger.json')['events']
    da=np.array([e['diameter_um'] for e in a]);db=np.array([e['diameter_um'] for e in b]);d=TruncatedSonoVue(SONOVUE);old=d.original
    cap=read(D/'inlet_size_capacity.json')['D_geometry_max_m']*1e6
    edges=np.unique(np.r_[old.low,old.high[-1],4.,cap]);source=np.diff(d.cdf(edges*1e-6));original=np.diff(old.cdf(edges));emp=np.histogram(da,bins=edges)[0]/len(da);centers=(edges[1:]+edges[:-1])/2
    rows=[dict(low_um=float(x),high_um=float(y),original_mass=float(z),truncated_source_mass=float(s),entering_empirical_mass=float(e)) for x,y,z,s,e in zip(edges[:-1],edges[1:],original,source,emp)]
    csvwrite(D/'figure01_diameter_distribution.csv',rows)
    fig,axs=plt.subplots(1,2,figsize=(10.5,3.8),gridspec_kw={'width_ratios':[2.5,1]})
    for ax in axs:
        ax.stairs(original/np.diff(edges),edges,color='#4477AA',label='Original SonoVue',lw=1.8)
        ax.stairs(source/np.diff(edges),edges,color='#EE7733',label='Source: D ≤ 4 µm',lw=1.7,linestyle='--')
        ax.stairs(emp/np.diff(edges),edges,color='#228833',label='Entering: Method C (n = 2000)',fill=True,alpha=.23,lw=1.4)
        ax.set(xlabel='Diameter (µm)',ylabel='Probability density (µm⁻¹)');ax.grid(axis='y',alpha=.15)
    axs[0].axvline(cap,color='#555555',ls=':',lw=1,label=f'Inlet limit ≈ {cap:.4f} µm');axs[0].set_title('(a) Source and entering distributions');axs[0].legend(frameon=False,fontsize=8.5);axs[0].set_xlim(.7,5.3)
    axs[1].set(xlim=(3.65,4.35),ylim=(0,.09),title='(b) No clipping peak at 4 µm');axs[1].axvline(4.,color='#666666',ls=':',lw=1)
    fig.tight_layout();save(fig,'01_diameter_distributions',['figure01_diameter_distribution.csv','source_probability.json'])
    geom=np.load(D/'inlet_geometry.npz');tri=geom['triangles_m'];pts=tri.reshape(-1,3);origin=pts.mean(axis=0);basis=np.linalg.svd(pts-origin,full_matrices=False)[2][:2]
    projection=(tri-origin)@basis.T*1e6
    write(D/'inlet_projection.json',dict(origin_m=origin.tolist(),basis=basis.tolist(),triangles_um=projection.tolist()))
    fig,axs=plt.subplots(1,2,figsize=(10,4.6),sharex=True,sharey=True,layout='constrained');scatter=[]
    plotrows=[]
    for ax,es,ds,label in zip(axs,[b,a],[db,da],['(a) Legacy Method B (n = 500)','(b) Method C (n = 2000)']):
        ax.add_collection(PolyCollection(projection,facecolors='#fafafa',edgecolors='#cccccc',linewidths=.35))
        p=(np.array([e['birth_center_m'] for e in es])-origin)@basis.T*1e6
        sc=ax.scatter(p[:,0],p[:,1],c=ds,cmap='viridis',vmin=.75,vmax=cap,s=7 if len(es)>500 else 12,alpha=.8,linewidths=0,rasterized=True)
        scatter.append(sc);ax.set(aspect='equal',xlabel='Inlet coordinate 1 (µm)',title=label);ax.autoscale_view();ax.margins(.06)
        for e,xy,diam in zip(es,p,ds):plotrows.append(dict(method='B' if es is b else 'C',particle_id=e['particle_id'],x_um=xy[0],y_um=xy[1],diameter_um=diam))
    axs[0].set_ylabel('Inlet coordinate 2 (µm)');fig.colorbar(scatter[-1],ax=axs,label='Diameter (µm)',shrink=.8)
    csvwrite(D/'figure02_positions.csv',plotrows);save(fig,'02_inlet_positions',['figure02_positions.csv','inlet_projection.json'])
    ap=read(D/'audit2000_point.json');old_point=list(csv.DictReader((R/'reference/method_b_500_point.csv').open()));bc=np.array([sum(r['point_outlet']==k for r in old_point) for k in ROLES]);ac=np.array([ap['counts'][k] for k in ROLES])
    fig,axs=plt.subplots(1,2,figsize=(10,3.8),gridspec_kw={'width_ratios':[1,1.4]});x=np.arange(4);w=.36
    axs[0].bar(x-w/2,bc/5,w,color='#4477AA',label='Method B, n = 500');axs[0].bar(x+w/2,ac/20,w,color='#EE7733',label='Method C, n = 2000')
    axs[0].set(xticks=x,xticklabels=['O1','O2','O3','No exit'],ylabel='Birth-center point tracers (%)',title='(a) Natural point-basin distribution');axs[0].legend(frameon=False,fontsize=8.5);axs[0].set_ylim(0,112)
    cell=[[str(v) for v in bc],[str(v) for v in ac]];axs[1].axis('off');table=axs[1].table(cellText=cell,rowLabels=['Method B (500)','Method C (2000)'],colLabels=['O1','O2','O3','No exit'],loc='center',cellLoc='center',bbox=[.24,.3,.76,.4]);table.auto_set_font_size(False);table.set_fontsize(10)
    axs[1].set_title('(b) Actual counts; no outlet quotas');fig.tight_layout();save(fig,'03_point_basin_comparison',['audit2000_point.json','../reference/method_b_500_point.csv'])
    proposals=np.array([e['position_proposal_count'] for e in a]);low=np.array([e['flux_bounds']['feasible_flux_fraction_lower'] for e in a]);hi=np.array([e['flux_bounds']['feasible_flux_fraction_upper'] for e in a]);ix=np.argsort(da)
    csvwrite(D/'figure04_position_proposals.csv',[dict(particle_id=e['particle_id'],diameter_um=e['diameter_um'],proposals=e['position_proposal_count'],guard=e['flux_bounds']['position_guard'],feasible_flux_lower=e['flux_bounds']['feasible_flux_fraction_lower'],feasible_flux_upper=e['flux_bounds']['feasible_flux_fraction_upper']) for e in a])
    fig,axs=plt.subplots(1,2,figsize=(10,3.8));axs[0].scatter(da,proposals,s=9,c='#4477AA',alpha=.35,linewidths=0,rasterized=True);axs[0].set(xlabel='Fixed diameter (µm)',ylabel='Position proposals until acceptance',title='(a) Actual retries with exact acceleration',ylim=(.5,proposals.max()+1))
    axs[1].fill_between(da[ix],low[ix],hi[ix],color='#228833',alpha=.3,label='Certified lower / upper bounds');axs[1].set(xlabel='Fixed diameter (µm)',ylabel='Feasible fraction of original inlet flux',yscale='log',title='(b) Larger spheres have less feasible flux');axs[1].legend(frameon=False,fontsize=8.5)
    for ax in axs:ax.grid(axis='y',alpha=.15)
    fig.tight_layout();save(fig,'04_position_proposals',['figure04_position_proposals.csv'])
    draws=[v for e in a for v in e['source_diameter_draws']];rejected=np.array([v['diameter_m']*1e6 for v in draws if v['status']=='NO_FEASIBLE_INLET_POSITION_FOR_SIZE']);pr=read(D/'source_probability.json')
    csvwrite(D/'figure05_source_draws.csv',[dict(particle_id=e['particle_id'],attempt=v['diameter_draw_id'][-1],diameter_um=v['diameter_m']*1e6,status=v['status']) for e in a for v in e['source_diameter_draws']])
    fig,axs=plt.subplots(1,2,figsize=(10,3.8));axs[0].hist([da,rejected],bins=np.linspace(.75,4.,34),stacked=True,color=['#228833','#CC6677'],label=['Globally feasible → position sampling','No feasible inlet position for size']);axs[0].axvline(cap,color='#333333',ls=':',lw=1);axs[0].set(xlabel='Source draw diameter (µm)',ylabel='Number of source draws',title='(a) Global size rejection is explicit');axs[0].legend(frameon=False,fontsize=8)
    frac=len(rejected)/len(draws)*100;theory=np.mean(pr['global_impossible_source_probability_bracket'])*100
    axs[1].bar(['Observed draws','Source probability'],[frac,theory],color=['#CC6677','#999999'],width=.55)
    axs[1].set(ylabel='Globally impossible source sizes (%)',title='(b) Source and finite inlet capacity',ylim=(0,max(frac,theory)*1.3))
    axs[1].text(0,frac+.25,f'{len(rejected)} / {len(draws)}\n{frac:.2f}%',ha='center',fontsize=10);axs[1].text(1,theory+.25,f'{theory:.3f}%',ha='center',fontsize=10)
    fig.tight_layout();save(fig,'05_global_source_size_rejection',['figure05_source_draws.csv','inlet_size_capacity.json','source_probability.json'])
    if (D/'formal500_gate.json').exists():
        fp=read(D/'formal500_point.json');fg=read(D/'formal500_gate.json');flow=read(R/'reference/flow_split.json')['outlet_fractions']
        values=np.array([[flow[k]*100 for k in ROLES[:3]]+[0.],list(bc/5),[fp['counts'][k]/5 for k in ROLES],[fg['outlet_counts'][k]/5 for k in ROLES]])
        labels=['Fluid volume flow','Method B accepted-500 point','Method C accepted-500 point','Method C P9-A.2 final']
        csvwrite(D/'figure06_routing.csv',[dict(stage=label,**{k:float(v) for k,v in zip(ROLES,row)}) for label,row in zip(labels,values)])
        fig,(ax,tab)=plt.subplots(1,2,figsize=(13,4.3),gridspec_kw={'width_ratios':[2.3,1]},sharey=True);left=np.zeros(4)
        for j,(name,color) in enumerate(zip(['O1','O2','O3','No exit'],COLORS)):
            ax.barh(np.arange(4),values[:,j],left=left,color=color,label=name,height=.65)
            for i,v in enumerate(values[:,j]):
                if v>=8:ax.text(left[i]+v/2,i,f'{v:.1f}%',ha='center',va='center',color='white',fontsize=10)
            left+=values[:,j]
        ax.set(yticks=np.arange(4),yticklabels=labels,xlabel='Fraction (%)',xlim=(0,100));ax.invert_yaxis();ax.legend(ncol=4,frameon=False,loc='lower center',bbox_to_anchor=(.5,1.01))
        tab.set_xlim(-.55,3.55);tab.axis('off');tab.set_title('All fractions (%)',pad=17,fontsize=10)
        for j,name in enumerate(['O1','O2','O3','No exit']):
            tab.text(j,-.52,name,ha='center',va='center',color=COLORS[j],fontsize=10,weight='bold')
            for i in range(4):tab.text(j,i,f'{values[i,j]:.2f}',ha='center',va='center',fontsize=10)
        fig.suptitle('Outlet audit | inlet sampling and downstream routing',fontsize=12,y=1.01)
        fig.tight_layout();save(fig,'06_routing_stages',['figure06_routing.csv','formal500_gate.json','formal500_point.json'])
    if (D/'independent_flux_reference.json').exists():
        ref=read(D/'independent_flux_reference.json')['rows'];fig,axs=plt.subplots(1,2,figsize=(10.5,4));x=np.array([r['diameter_um'] for r in ref]);bottom=np.zeros(len(ref))
        for k,color in zip(ROLES,COLORS):
            y=np.array([r['conditional_point_counts'][k]/5 for r in ref]);axs[0].bar(x,y,width=.14,bottom=bottom,color=color,label=k.replace('OUTLET_0','O').replace('NO_EXIT','No exit'));bottom+=y
        axs[0].set(xlabel='Fixed diameter (µm)',ylabel='Point-basin fraction (%)',title='(a) Original flux sampler + direct rejection');axs[0].set_title('(a) Original flux sampler + direct rejection',pad=35);axs[0].legend(frameon=False,fontsize=8,ncol=4,loc='lower center',bbox_to_anchor=(.5,1.01))
        axs[1].hist(db,bins=np.linspace(.75,3.15,25),density=True,histtype='step',lw=1.8,color='#4477AA',label='Legacy B accepted 500');axs[1].hist(da,bins=np.linspace(.75,3.15,25),density=True,histtype='step',lw=1.8,color='#EE7733',label='Method C entering 2000');axs[1].set(xlabel='Diameter (µm)',ylabel='Probability density (µm⁻¹)',title='(b) Size conditioning changes the population');axs[1].legend(frameon=False,fontsize=8);axs[1].set_title('(b) Size conditioning changes the population',pad=35)
        fig.tight_layout();save(fig,'07_finite_size_basin_explanation',['independent_flux_reference.json','../reference/method_b_500_birth_ledger.json','../audit2000/admission/birth_ledger.json'])
    write(D/'figure_manifest.json',manifest)
    print('saved',len(manifest),'PNG + PDF figures')
if __name__=='__main__':main()
