"""All figures use saved full-population diagnostics; no new scientific runs."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

OUT=Path(__file__).resolve().parent;ROOT=OUT.parents[1]
FIG=ROOT/'results/higher_order_baqp_rk_oracle_v1/figures'
BLUE,ORANGE,GREEN,RED,GRAY='#176b9b','#c46522','#388368','#a73c50','#596777'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':11,'axes.spines.top':False,'axes.spines.right':False,
    'axes.grid':True,'grid.alpha':.17,'pdf.fonttype':42,'ps.fonttype':42,'legend.frameon':False})
QA={}
def read(n):return json.loads((OUT/n).read_text())
def rows(n):
    with (OUT/n).open() as f:return list(csv.DictReader(f))
def subset(r,**kw):return [v for v in r if all(v[k]==x for k,x in kw.items())]
def col(r,k):return np.array([float(v[k]) for v in r])
def save(fig,name,caption):
    fig.text(.08,.028,caption,fontsize=8.7,color=GRAY,va='bottom',linespacing=1.5)
    fig.canvas.draw();renderer=fig.canvas.get_renderer();w,h=fig.canvas.get_width_height();texts=list(fig.texts)
    for legend in fig.legends:texts+=legend.get_texts()
    for a in fig.axes:
        texts += [a.title,a.xaxis.label,a.yaxis.label]+list(a.texts)
        if a.get_legend():texts+=a.get_legend().get_texts()
        for axis,lim in [(a.xaxis,a.get_xlim()),(a.yaxis,a.get_ylim())]:
            texts.append(axis.get_offset_text())
            for t in axis.get_major_ticks()+axis.get_minor_ticks():
                if min(lim)<=t.get_loc()<=max(lim):texts.extend([t.label1,t.label2])
    outside=[]
    for t in texts:
        if not t.get_visible() or not t.get_text():continue
        b=t.get_window_extent(renderer)
        if b.x0 < -1 or b.y0 < -1 or b.x1>w+1 or b.y1>h+1:outside.append(t.get_text())
    QA[name]={'text_outside_canvas':outside,'all_values_from_saved_data':True,'random_or_synthetic_samples':False}
    if outside:raise RuntimeError((name,outside))
    fig.savefig(FIG/(name+'.png'),dpi=600,facecolor='white');fig.savefig(FIG/(name+'.pdf'),facecolor='white');plt.close(fig)
    print('FIGURE',name,flush=True)

def boundary():
    r=rows('F_DEPTH_PROFILES.csv');fig,ax=plt.subplots(1,3,figsize=(12.5,5.5))
    lowest=min(float(v['max_scaled']) for v in r if v['sample']=='Fbar' and float(v['max_scaled'])>0)/2
    fig.subplots_adjust(left=.08,right=.98,bottom=.25,top=.81,wspace=.28)
    for k,field in enumerate(['axial','shear','cross']):
        for model,c,m,ls in [('C1',BLUE,'o','-'),('PD2',ORANGE,'s','-'),('RK2_DIAG',GREEN,'^','--')]:
            q=subset(r,field=field,model=model,sample='Fbar');x=col(q,'depth_over_delta');y=col(q,'max_scaled')
            ax[k].semilogy(x,np.where(y>0,y,np.nan),ls,marker=m,mfc='none',ms=4,c=c,label=model)
        ax[k].axhline(1e-10,c=RED,ls=':',lw=1)
        ax[k].set(title=['(a) Quadratic axial','(b) Quadratic shear','(c) XZ cross term'][k],xlabel=r'Nearest-face depth $d/\delta$',ylim=(lowest,.3))
    ax[0].set_ylabel(r'Maximum $\|\bar F-F_{exact}\|_F/0.004$')
    fig.legend(*ax[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.53,.91),ncol=3)
    fig.suptitle('Full quadratic RK removes the truncated-family gradient bias',y=.975,fontsize=14)
    save(fig,'c1_vs_pd2_boundary_bias',
         'All 6,859 nodes; fixed coordinates, neighbors, horizon, quadratic kernel and physical volumes. Dotted line: 1e-10 gate.\n'
         'RK2_DIAG contains only three linear and three squared terms; it does not contain XY, XZ or YZ.')

def polynomial():
    with np.load(OUT/'REFERENCE_GEOMETRY.npz') as a:classes=a['geometry_class']
    with np.load(OUT/'POLYNOMIAL_REPRODUCTION.npz') as a:r=a['normalized_residual']
    cats=['interior-like','face-near','edge-near','corner-near'];values=np.stack([r[classes==c].max(axis=0) for c in cats],axis=1)
    fig,ax=plt.subplots(figsize=(8.5,7.4));fig.subplots_adjust(left=.12,right=.85,bottom=.19,top=.88)
    im=ax.pcolormesh(np.arange(5)-.5,np.arange(10)-.5,values,norm=LogNorm(vmin=values.min()*.8,vmax=values.max()*1.25),cmap='viridis',shading='flat',rasterized=False)
    ax.invert_yaxis()
    ax.set_xticks(range(4),['Interior','Face','Edge','Corner']);ax.set_yticks(range(9),['X','Y','Z',r'$X^2$',r'$Y^2$',r'$Z^2$','XY','XZ','YZ']);ax.grid(False)
    for i in range(9):
        for j in range(4):ax.text(j,i,f'{values[i,j]:.2e}',ha='center',va='center',color='white' if values[i,j]<np.sqrt(values.min()*values.max()) else '#162632',fontsize=10)
    cb=fig.colorbar(im,ax=ax,fraction=.045,pad=.05);cb.set_label('Maximum normalized residual');cb.solids.set_rasterized(False)
    fig.suptitle('All nine physical polynomial gradients are reproduced',y=.96,fontsize=13)
    ax.set_xlabel('Pure geometric class; every node in each class is included')
    save(fig,'pd2_polynomial_reproduction',
         'Normalization: 1 for linear derivatives; L = 1.125 mm for quadratic derivatives.\n'
         'All class maxima are below 1e-10. Global volume-weighted RMS also passes 1e-11.')

def conditioning():
    with np.load(OUT/'PD2_REFERENCE_CACHE.npz') as a:c=a['condition'];sigma=a['singular_values'][:,-1];rank=a['rank']
    with np.load(OUT/'REFERENCE_GEOMETRY.npz') as a:d=a['depth_over_delta'];classes=a['geometry_class']
    bins=np.rint(d*4.02).astype(int);levels=np.unique(bins);x=levels/4.02
    fig,ax=plt.subplots(2,2,figsize=(11,8.5));fig.subplots_adjust(left=.10,right=.97,bottom=.16,top=.87,hspace=.46,wspace=.33)
    for q,ls,c0,label in [(50,'-',BLUE,'Median'),(95,'--',GREEN,'P95'),(100,'-',ORANGE,'Maximum')]:
        ax[0,0].plot(x,[np.percentile(c[bins==z],q) for z in levels],ls,marker='o',ms=3,c=c0,label=label)
    ax[0,0].set(title='(a) Moment-matrix condition number',xlabel=r'Nearest-face depth $d/\delta$',ylabel=r'$\kappa_2(M_2)$');ax[0,0].legend(fontsize=9)
    ax[0,1].semilogy(x,[sigma[bins==z].min() for z in levels],'o-',c=BLUE,ms=4)
    ax[0,1].set(title='(b) Smallest singular value in each layer',xlabel=r'Nearest-face depth $d/\delta$',ylabel=r'Minimum $\sigma_{min}(M_2)$ [m$^3$]')
    cats=['interior-like','face-near','edge-near','corner-near'];labels=['Interior','Face','Edge','Corner'];counts=[int(((classes==g)&(rank<9)).sum()) for g in cats]
    ax[1,0].bar(range(4),counts,color=BLUE);ax[1,0].set_xticks(range(4),labels);ax[1,0].set_ylim(0,1)
    for k,v in enumerate(counts):ax[1,0].text(k,.05,str(v),ha='center',fontsize=13)
    ax[1,0].set(title='(c) Rank-deficient physical nodes',ylabel='Count (rank < 9)')
    ax[1,1].boxplot([c[classes==g] for g in cats],tick_labels=labels,showfliers=True,flierprops={'markersize':2,'alpha':.25},medianprops={'color':ORANGE})
    ax[1,1].set(title='(d) All-node condition distributions',ylabel=r'$\kappa_2(M_2)$')
    fig.suptitle('PD2 local moment systems remain full rank near faces, edges and corners',y=.96,fontsize=13)
    save(fig,'pd2_conditioning',
         'SVD rank tolerance: 9 eps sigma_max. No inverse, pseudo-inverse, regularization or horizon adjustment.\n'
         'All nodes are WELL_CONDITIONED under the predeclared report-only flag (kappa <= 1e6).')

def bond_energy():
    f=rows('F_DEPTH_PROFILES.csv');p=rows('PD2_PSI_PLUS_DIAGNOSTIC.csv')
    fig,ax=plt.subplots(1,2,figsize=(12,5.9));fig.subplots_adjust(left=.085,right=.97,bottom=.32,top=.76,wspace=.28)
    for field,c in [('axial',BLUE),('shear',ORANGE),('cross',GREEN)]:
        for model,ls,marker in [('C1','-','o'),('PD2','--','s')]:
            for a,rows0 in [(ax[0],f),(ax[1],p)]:
                q=subset(rows0,field=field,model=model,sample='Ftilde');y=col(q,'max_scaled')
                a.semilogy(col(q,'depth_over_delta'),np.where(y>0,y,np.nan),ls,c=c,marker=marker,ms=3,mfc='none',label=field.capitalize()+' / '+model)
    ax[0].axhline(1e-10,c=GRAY,ls=':',lw=1);ax[1].axhline(1e-8,c=GRAY,ls=':',lw=1)
    ax[0].set(title='(a) Bond-associated gradient',ylabel=r'Maximum $\|\widetilde F-F_{exact,mid}\|_F/0.004$')
    ax[1].set(title='(b) Bond positive-energy density',ylabel=r'Maximum $|\Delta\Psi^+|/\Psi_{ref}$')
    for a, records in zip(ax, [f, p]):
        positive=[float(v['max_scaled']) for v in records if v['sample']=='Ftilde' and float(v['max_scaled'])>0]
        a.set(xlabel=r'Bond-midpoint depth $d/\delta$',ylim=(min(positive)/2,.4))
    fig.suptitle('KINEMATIC / ENERGY-DRIVER DIAGNOSTIC\nNO HIGHER-ORDER FORCE USED',y=.965,fontsize=13,linespacing=1.6)
    fig.legend(*ax[0].get_legend_handles_labels(),loc='lower center',bbox_to_anchor=(.53,.115),ncol=3,fontsize=9)
    save(fig,'pd2_ftilde_and_energy_bias',
         'All 690,020 bonds. Solid: C1; dashed: PD2. Dotted lines: maximum-error gates.\n'
         'The same spectral split and one analytic physical-domain Psi reference per field are used.')

def julia():
    q=read('PERIDYNAMICSJL_CROSSCHECK.json')
    if q['status']!='PASS':return
    fig,ax=plt.subplots(1,3,figsize=(12,5.1));fig.subplots_adjust(left=.08,right=.98,bottom=.27,top=.81,wspace=.30)
    colors=[BLUE,ORANGE,GREEN,RED];names=['Interior','Face','Edge','Corner']
    for k,field in enumerate(['axial','shear','cross']):
        r=[v for v in q['records'] if v['field']==field]
        for n,v in enumerate(r):ax[k].scatter(np.asarray(v['python_F']).ravel(),np.asarray(v['julia_F']).ravel(),s=25,facecolors='none',edgecolors=colors[n],label=names[n])
        ax[k].plot([-.01,1.01],[-.01,1.01],'--',c=GRAY,lw=1)
        ax[k].set(title=f'({"abc"[k]}) {field.capitalize()}',xlabel='Python F component',ylabel='Peridynamics.jl F component')
        ax[k].text(.05,.83,'Max difference / 0.004\n'+f'{q["case_results"][field]["max_F_difference_scaled"]:.2e}',transform=ax[k].transAxes,fontsize=9)
    fig.suptitle('Independent full-PD2 gradient comparison with the installed package',y=.97,fontsize=13)
    fig.legend(*ax[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.54,.91),ncol=4,fontsize=9)
    save(fig,'python_vs_peridynamicsjl_pd2',
         'All nine tensor components at each of four representative nodes; many identity/zero components overlap.\n'
         'Same kernel, horizon and physical volumes. Package source unchanged; only gradient weights and F evaluated.')

if __name__=='__main__':
    for f in [boundary,polynomial,conditioning,bond_energy,julia]:f()
    (OUT/'FIGURE_LAYOUT_QA.json').write_text(json.dumps(QA,indent=2)+'\n')
