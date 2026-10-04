"""Scientific plots of measured fields; unrun candidates are explicitly N/A."""
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from scripts.activation_root_cause_analysis import read,save,arrays,ROOT,OUT,RESULTS,FULL,AXIAL,BoundaryProblem

FIG=RESULTS/'figures'
COLORS={5:'#1967a3',7:'#c15b22'}
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'savefig.dpi':180,'pdf.fonttype':42})
records=[]


def data(n,mode=FULL,e=.003):
    return read(OUT/'analysis'/f'C{n}_{mode}_{e:.4f}.json')


def table(n,mode=FULL):
    return list(csv.DictReader((RESULTS/'tables'/f'C{n}_{mode}_0.0030_AXIAL_BONDS.csv').open()))


def finish(fig,name,description):
    FIG.mkdir(parents=True,exist_ok=True)
    for ext in ('png','pdf'):fig.savefig(FIG/(name+'.'+ext),bbox_inches='tight')
    plt.close(fig);records.append(dict(name=name,formats=['PNG','PDF'],description=description))


def historical_activation():
    rows=[]
    for n,task,file in [(5,'pfpd_at1_c5_activation_stationarity_audit_v1','POINT.json'),(7,'pfpd_at1_c7_early_branch_v1','WEAK.json')]:
        base=ROOT/'verification'/task;r=read(base/'activation'/file)
        a=arrays(base/'intact'/r['label']/'STATE.npz');p=BoundaryProblem(n,r['strain'],FULL)
        st=p.state(a['x'],a['s']);mid=p.family.X[p.family.pairs].mean(axis=1)
        b=int(np.argmax(st.W0));d=.000125-abs(mid[b,2]);bracket=read(base/'activation/BRACKET.json')['interval']
        rows.append(dict(level=n,strain=r['strain'],bracket=bracket,reaction_stress_Pa=float(-st.force[p.top,2].sum()/.00025**2),
                         d_clamp_m=float(d),d_over_h=float(d/p.cloud.spacing),d_over_delta=float(d/p.family.horizon),
                         source=str((base/'intact'/r['label']/'STATE.npz').relative_to(ROOT)),new_equilibrium_solve=False))
    save(OUT/'BASELINE_ACTIVATION.json',rows);return rows


def main():
    decision=read(OUT/'PHASE_A_DECISION.json');scale=read(OUT/'PROFILE_SCALING.json');baseline=historical_activation()
    fig,axes=plt.subplots(1,2,figsize=(11,4),constrained_layout=True)
    for ax,n in zip(axes,(5,7)):
        d=data(n);raw=table(n);p=d['profile'];x=np.array([r['d_clamp_m'] for r in p])*1e6
        ax.scatter([float(r['d_clamp_m'])*1e6 for r in raw],[float(r['Psi_normalized']) for r in raw],s=8,alpha=.2,color=COLORS[n],label='Equilibrium: raw axial bonds')
        ax.plot(x,[r['mean_Psi_normalized'] for r in p],'o-',color=COLORS[n],label='Equilibrium: layer mean')
        ax.plot(x,[r['equilibrium_Psi_stats']['max']/read(OUT/'CONFIG.json')['Yc_Pa'] for r in p],'--',color=COLORS[n],label='Equilibrium: layer maximum')
        ax.axhline(d['affine']['Psi_plus_axial']['mean']/read(OUT/'CONFIG.json')['Yc_Pa'],color='#333333',ls=':',lw=2,label='Exact affine oracle')
        ax.set(title=f'C{n} | axial bonds | strain = 0.003',xlabel='Distance to nearest clamp (µm)',ylabel=r'$\Psi^+ / Y_c$')
        ax.grid(alpha=.2);ax.legend(fontsize=8)
    finish(fig,'baseline_affine_vs_equilibrium','Affine is uniform; hotspot develops in the equilibrium field.')

    fig,axes=plt.subplots(1,2,figsize=(11,4.5),constrained_layout=True)
    for ax,coord,label in zip(axes,['d_over_h','d_over_delta'],[r'$d_{clamp}/h$',r'$d_{clamp}/\delta$']):
        for n in (5,7):
            for e,ls in [(.0025,'--'),(.003,'-')]:
                p=data(n,e=e)['profile'];ax.plot([r[coord] for r in p],[r['mean_A_psi'] for r in p],ls,marker='o',ms=4,color=COLORS[n],label=f'C{n}, ε={e:.4f}')
        key='h' if coord=='d_over_h' else 'delta'
        vals=[r['metrics']['mean_A_psi']['RMSE_'+key] for r in scale['rows']]
        ax.text(.98,.95,f'RMSE: {vals[0]:.4f} / {vals[1]:.4f}',transform=ax.transAxes,ha='right',va='top')
        ax.set(xlabel=label,ylabel=r'Layer mean $A_\psi$',title='Raw profiles; linear interpolation on overlap')
        ax.grid(alpha=.2);ax.legend(fontsize=8,loc='upper center')
    domain=scale['rows'][0]['metrics']['mean_A_psi']['domain_delta']
    axes[1].axvspan(domain[0],domain[1],alpha=.07,color='grey')
    axes[1].annotate('Nearest C7 layer is outside\nthe common δ interval',xy=(.0829187,7.76),xytext=(.23,5.6),
                     arrowprops={'arrowstyle':'->','color':'#555555'},fontsize=8)
    fig.suptitle('δ scaling has lower RMSE; this does not establish peak or width convergence',fontsize=11)
    finish(fig,'baseline_profile_scaling','Both strains satisfy RMSE_delta < 0.5 RMSE_h; comparison-domain caveat shown.')

    fig,axes=plt.subplots(2,2,figsize=(10,7),constrained_layout=True)
    for rr,n in enumerate((5,7)):
        raw=table(n);d=data(n)
        for cc,(key,name,rhokey) in enumerate([('endpoint_support_mean','Mean endpoint weighted support','weighted_support'),('endpoint_shape_condition_max','Maximum endpoint shape condition','shape_condition')]):
            ax=axes[rr,cc];xx=np.array([float(r[key]) for r in raw]);yy=[float(r['A_psi']) for r in raw]
            if cc==0:xx/=max(xx);name+=' / max'
            ax.scatter(xx,yy,s=9,alpha=.2,color=COLORS[n]);sel=np.array([r['historical_activation']=='True' for r in raw])
            ax.scatter(xx[sel],np.array(yy)[sel],marker='*',s=130,color='#b32130',edgecolors='black',lw=.4,label='Historical activation bonds')
            ax.set(xlabel=name,ylabel=r'$A_\psi$',title=f'C{n}: Spearman ρ = {d["correlations"][rhokey]["rho"]:.3f}')
            ax.grid(alpha=.2);ax.legend(fontsize=8)
    finish(fig,'support_and_driver','Correlations are descriptive; highlighted hotspots do not have the weakest same-layer support.')

    fig,axes=plt.subplots(1,3,figsize=(12,4),constrained_layout=True)
    for ax,metric,title in zip(axes,['first','amp','reaction'],[r'First-layer max $\Psi^+/Y_c$',r'Maximum $A_\psi$','Reaction stress (Pa)']):
        for offset,mode,label,color in [(-.18,FULL,'Full clamp','#4b657e'),(.18,AXIAL,'Axial-only control','#e1a457')]:
            values=[]
            for n in (5,7):
                d=data(n,mode);values.append(d['first_layer_Psi_normalized']['max'] if metric=='first' else d['A_psi_all']['max'] if metric=='amp' else d['equilibrium']['reaction_stress_Pa'])
            ax.bar(np.arange(2)+offset,values,.35,label=label,color=color)
        ax.set_xticks([0,1],['C5','C7']);ax.set_title(title);ax.grid(axis='y',alpha=.2)
    axes[0].legend(fontsize=8);fig.suptitle('Strain = 0.003: releasing lateral constraints does not reduce peak drivers',fontsize=12)
    finish(fig,'full_clamp_vs_axial_only','Common-strain control results, not a validated candidate fix.')

    fig,axes=plt.subplots(2,2,figsize=(9,8),constrained_layout=True);norm=Normalize(0,1)
    for rr,n in enumerate((5,7)):
        for cc,mode in enumerate((FULL,AXIAL)):
            ax=axes[rr,cc];raw=table(n,mode);sel=[r for r in raw if abs(float(r['midpoint_y_m']))<1e-14]
            im=ax.scatter([float(r['midpoint_x_m'])*1e6 for r in sel],[float(r['midpoint_z_m'])*1e6 for r in sel],
                c=[float(r['Psi_normalized']) for r in sel],cmap='viridis',norm=norm,s=25)
            hot=[r for r in sel if r['historical_activation']=='True']
            ax.scatter([float(r['midpoint_x_m'])*1e6 for r in hot],[float(r['midpoint_z_m'])*1e6 for r in hot],facecolors='none',edgecolors='red',s=110,lw=1.3)
            for z in (-125,125):ax.axhline(z,color='#444444',lw=2)
            ax.set(xlim=(-135,135),ylim=(-140,140),aspect='equal',xlabel='x (µm)',ylabel='z (µm)',
                   title=f'C{n}: '+('full clamp' if mode==FULL else 'axial-only control'))
    fig.colorbar(im,ax=axes.ravel().tolist(),label=r'$\Psi^+/Y_c$ at ε=0.003',shrink=.8)
    fig.suptitle('Axial-bond midpoints in y=0 plane; red rings: historical first-layer hotspots',fontsize=10)
    finish(fig,'driver_spatial_map','Shared physical coordinates and color limits; loaded end planes marked.')

    fig,axes=plt.subplots(1,3,figsize=(11,4),constrained_layout=True)
    for ax,key,mult,title in zip(axes,['strain','reaction_stress_Pa','d_over_h'],[100,1,1],['Activation strain (%)','Activation reaction stress (Pa)',r'Activation $d/h$']):
        values=[r[key]*mult for r in baseline];ax.bar([0,1],values,color=[COLORS[5],COLORS[7]],width=.55)
        ax.set_xticks([0,1],['C5 baseline','C7 baseline']);ax.set_title(title)
        ax.text(.5,.9,'Candidate: NOT RUN',transform=ax.transAxes,ha='center',color='#7b2424',bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
        ax.grid(axis='y',alpha=.2)
    fig.suptitle('Historical activation retained; Phase A does not justify promoting the BC control',fontsize=11)
    finish(fig,'activation_before_after','Before values shown; after values unavailable because no candidate was selected.')

    fig,ax=plt.subplots(figsize=(7,4),constrained_layout=True)
    value=100*decision['C5_C7_ACTIVATION_DISCREPANCY_BASELINE'];ax.bar([0],[value],width=.5,color='#4b657e')
    ax.text(0,value+1,f'{value:.2f}%',ha='center');ax.text(1,15,'NOT RUN\nNo justified candidate',ha='center',color='#7b2424')
    ax.set_xticks([0,1],['Baseline','Candidate']);ax.set(xlim=(-.6,1.6),ylim=(0,45),ylabel='Symmetric relative discrepancy (%)',title='C5–C7 activation sensitivity')
    ax.text(.5,.97,'Denominator = mean(C5, C7)\nEquivalent signed change relative to C5: −30.235%',transform=ax.transAxes,ha='center',va='top',fontsize=9)
    ax.grid(axis='y',alpha=.2)
    finish(fig,'activation_discrepancy','Uses the requested mean denominator; candidate is not plotted as zero.')

    # Additional physical mechanism view using only already accepted states.
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for ax,n in zip(axes,(5,7)):
        d=data(n);g=d['centerline_gaps'];z=[.5*(r['z_left_m']+r['z_right_m'])*1e6 for r in g]
        ax.plot(z,[r['engineering_strain']/.003 for r in g],'o-',color=COLORS[n],label='Equilibrium gap strain / global strain')
        ax.axhline(1,color='black',ls=':',label='Affine reference');ax.axhline(0,color='grey',lw=.6)
        ax.set(title=f'C{n}: centerline nearest-neighbor intervals',xlabel='Interval midpoint z (µm)',ylabel='Local / global engineering strain')
        ax.legend(fontsize=8);ax.grid(alpha=.2)
    finish(fig,'centerline_displacement_gaps','Boundary displacement gaps persist as h shrinks; some interior intervals compress.')
    save(OUT/'FIGURES.json',dict(figures=records,early_phase_support='NOT_GENERATED: no new positive-phase states',
                               activation_after='NOT_RUN, explicitly marked in figures 6 and 7'))


if __name__=='__main__':main()
