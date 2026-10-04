"""Eight black-background review figures; read saved scientific arrays only."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pd_clot.pfpd_candidate.coupling_stability_audit import AuditProblem,physics_config,read_json,read_npz
from pd_clot.pfpd_candidate.branch_stability_tracking import OUT
from scripts.analyze_branch_stability_tracking import all_rows,mode_record

font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/arial.ttf')
plt.rcParams.update({'font.family':'Arial','font.size':12,'font.weight':'normal','axes.titleweight':'normal',
    'figure.facecolor':'black','axes.facecolor':'black','savefig.facecolor':'black','text.color':'white',
    'axes.labelcolor':'white','axes.edgecolor':'#999999','xtick.color':'white','ytick.color':'white',
    'grid.color':'#555555','grid.alpha':.45,'legend.facecolor':'black','legend.edgecolor':'none','pdf.fonttype':42,
    'axes.spines.top':False,'axes.spines.right':False})
BLUE='#49b8ff';ORANGE='#ffb35b';PINK='#ec81b3';GREEN='#75dda5';GREY='#b5bec9'
F=OUT/'figures';F.mkdir(exist_ok=True)
A,B=all_rows('A'),all_rows('B');brackets={k:read_json(OUT/f'critical_bracket/SCHEDULE_{k}.json') for k in ['A','B']}
lo,hi=brackets['B']['final_interval'];damage=read_json(OUT/'comparisons/DAMAGE_ONSET.json')


def clean(ax):ax.grid(True,linewidth=.5);ax.tick_params(labelsize=11)
def shade(ax):ax.axvspan(lo*100,hi*100,color=PINK,alpha=.25,label='Stability-loss interval')
def finish(fig,name,note):
    fig.text(.5,.015,note,ha='center',va='bottom',fontsize=10.5,color=GREY)
    fig.savefig(F/(name+'.png'),dpi=180,bbox_inches='tight',pad_inches=.18)
    fig.savefig(F/(name+'.pdf'),bbox_inches='tight',pad_inches=.18);plt.close(fig)


def stability():
    fig,ax=plt.subplots(figsize=(11,6.5));fig.subplots_adjust(bottom=.2)
    for rows,marker,label in [(A,'o','Schedule A'),(B,'x','Schedule B')]:
        for state,color in [('STABLE_SAMPLED_STATE',GREEN),('UNSTABLE_SAMPLED_STATE',ORANGE),('UNRESOLVED',GREY)]:
            r=[v for v in rows if v['classification']==state and v['indicator_kind']=='reduced_eigenvector']
            if r:ax.scatter([v['strain']*100 for v in r],[v['stability_indicator'] for v in r],color=color,marker=marker,s=44,label=label+': '+state.split('_')[0].lower())
    upper=brackets['B']['upper']
    ax.scatter([100*hi],[upper['stability_indicator']],marker='v',s=100,color=PINK,label='Feasible cone negative witness')
    ax.scatter([100*hi],[upper['reduced_eigenvalues_min'][1]],facecolors='none',edgecolors=PINK,s=100,label='Reduced Hessian excludes weak bounds')
    ax.axhline(0,color='white',ls='--',lw=1);shade(ax);clean(ax)
    ax.set_xlabel('Engineering strain (%)');ax.set_ylabel('Scaled minimum eigenvalue / cone curvature witness')
    ax.set_title('Stability changes when weak-bound directions become admissible');ax.legend(fontsize=9.5,loc='center right')
    finish(fig,'stability_eigenvalue_vs_strain','Triangle: a feasible negative Rayleigh quotient, not the exact minimum over the cone\nNo smooth zero crossing is inferred across the active-set switch')


def bracket_zoom():
    fig,ax=plt.subplots(figsize=(11,6));fig.subplots_adjust(bottom=.2)
    for k,rows,color,marker in [('A',A,BLUE,'s'),('B',B,ORANGE,'x')]:
        r=[v for v in rows if lo-3*(hi-lo)<=v['strain']<=hi]
        ax.scatter([v['strain']*100 for v in r],[v['stability_indicator'] for v in r],color=color,marker=marker,s=70,label='Schedule '+k)
    shade(ax);ax.axhline(0,color=GREY,ls='--');ax.set_xlim((lo-2*(hi-lo))*100,(hi+.6*(hi-lo))*100)
    ax.set_xlabel('Engineering strain (%)');ax.set_ylabel('Scaled eigenvalue / cone curvature witness');clean(ax)
    ax.set_title('Both schedules give the same nonsmooth stability interval');ax.legend(loc='center left')
    ax.text(.03,.94,f'Interval: [{lo:.10f}, {hi:.10f}]\nWidth: {hi-lo:.7g} engineering strain',transform=ax.transAxes,va='top',fontsize=11)
    finish(fig,'critical_bracket_zoom','Markers are actual solved states; no spline or continuous-root interpolation is used')


def force_phase():
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True);fig.subplots_adjust(hspace=.28,bottom=.15)
    for rows,color,marker,name in [(A,BLUE,'o','A'),(B,ORANGE,'x','B')]:
        x=[r['strain']*100 for r in rows]
        axs[0].plot(x,[r['reaction_N']*1e6 for r in rows],color=color,marker=marker,ms=4,label='Schedule '+name)
        axs[1].plot(x,[r['metrics']['maximum_s'] for r in rows],color=color,marker=marker,ms=4,label=name+': maximum phase')
        axs[1].plot(x,[r['metrics']['mean_s'] for r in rows],color=color,ls='--',label=name+': weighted mean phase')
    for ax in axs:shade(ax);clean(ax)
    axs[0].set_title('Stationary reaction and phase along the prescribed branch');axs[0].set_ylabel('Tensile reaction (uN)');axs[0].legend(fontsize=10)
    d=damage['B']['damage_onset_interval'];axs[1].axvspan(d[0]*100,d[1]*100,color=GREEN,alpha=.15,label='Damage onset: sampling interval')
    axs[1].set_ylabel('Phase');axs[1].set_xlabel('Engineering strain (%)');axs[1].legend(fontsize=9.5)
    finish(fig,'branch_force_phase','First sampled phase > 1e-8: A at 1.05%, B at 1.025%\nThe analytic local AT1 threshold is 1.00%; sampling alone does not prove a causal ordering')


def active():
    fig,axs=plt.subplots(1,2,figsize=(13,5.8));fig.subplots_adjust(bottom=.23,wspace=.25)
    for ax in axs:
        for key,color,label,marker in [('LOWER_STRONG_ACTIVE',GREEN,'Strong lower','o'),('LOWER_WEAK_NEAR_SWITCH',PINK,'Weak lower / near switch','v'),('INTERIOR',BLUE,'Interior','s')]:
            ax.scatter([r['strain']*100 for r in B],[r['active']['counts'][key] for r in B],color=color,label=label,marker=marker,s=40)
        shade(ax);ax.set_ylim(-8,200);ax.set_xlabel('Engineering strain (%)');clean(ax)
    axs[0].set_ylabel('Number of phase bonds');axs[0].set_title('Sampled active sets: Schedule B');axs[0].legend(fontsize=9)
    axs[1].set_xlim(.9975,1.003);axs[1].set_title('All 185 bonds become weak at the switch')
    finish(fig,'active_set_evolution','Strong lower bounds prevent phase opening in the critical cone\nWeak lower bounds permit one-sided phase opening; sampled interior states occur above 1.00%')


def landscape():
    fig,axs=plt.subplots(1,3,figsize=(15,6));fig.subplots_adjust(bottom=.27,wspace=.34,top=.82)
    labels=[('last_stable',1e7,1e26,'Last stable state'),('onset_cone_direction',1e5,1e23,'Onset: one-sided cone direction'),('nearest_interior_unstable',1e5,1e22,'Nearest interior unstable state')]
    for ax,(label,xscale,escale,title) in zip(axs,labels):
        r=read_json(OUT/f'energy_landscape/{label}.json');valid=[p for p in r['points'] if p['feasible']]
        alpha=np.array([p['alpha'] for p in valid]);energy=np.array([p['delta_energy_J'] for p in valid])
        ax.scatter(alpha*xscale,energy*escale,color=BLUE,s=40,label='Energy evaluations')
        a=np.linspace(min(alpha),max(alpha),200)
        ax.plot(a*xscale,r['energy_scale_J']*(r['first_variation']*a+.5*r['curvature']*a*a)*escale,color=ORANGE,ls='--',label='Hessian quadratic reference')
        ax.axhline(0,color=GREY,lw=.7);ax.set_xlabel(f'Mode amplitude (1e-{int(np.log10(xscale))})')
        ax.set_ylabel(f'Energy change (1e-{int(np.log10(escale))} J)');clean(ax)
        ax.set_title(title+f"\nStrain {r['strain']*100:.7g}%",fontsize=12)
    axs[0].legend(fontsize=9);fig.suptitle('Admissible energy probes change from upward to downward curvature',fontsize=15)
    finish(fig,'curvature_along_critical_mode','Left: |alpha| >= 1e-6 exits the original pure-tension domain; only smaller valid probes are plotted\nMiddle: negative alpha violates the phase lower bound. Right: both signs are admissible\nTiny energy differences are sensitive to float64 cancellation; dashed curves are references, not replacement data')


def critical_map():
    row=brackets['B']['upper'];d=mode_record(row);p=AuditProblem(physics_config(),row['strain'],d['lower'])
    v=d['critical_mode'];phase=np.zeros(len(p.family.X));weight=np.zeros_like(phase)
    for j,(i,k) in enumerate(p.family.pairs):
        phase[[i,k]]+=p.weights[j]*v[p.nx+j];weight[[i,k]]+=p.weights[j]
    phase/=weight;phase/=max(abs(phase));u=np.zeros_like(d['x']);u[p.free]=v[:p.nx].reshape(-1,3)*p.length_scale
    magnitude=np.linalg.norm(u,axis=1);magnitude/=magnitude.max();points=(d['x']-p.family.X.min(axis=0))*1e6
    fig=plt.figure(figsize=(13,7));fig.subplots_adjust(left=.02,right=.96,bottom=.15,top=.77,wspace=.16)
    for i,(values,cmap,title) in enumerate([(phase,'viridis','Phase opening: weighted particle aggregation'),(magnitude,'viridis','Displacement-mode magnitude')]):
        ax=fig.add_subplot(1,2,i+1,projection='3d');ax.set_proj_type('ortho')
        sc=ax.scatter(*points.T,c=values,vmin=0,vmax=1,cmap=cmap,s=190,depthshade=False,edgecolors='#bbbbbb',linewidths=.4)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:axis.set_pane_color((0,0,0,1));axis._axinfo['grid']['color']=(.35,.35,.35,.5)
        ax.set_xlabel('x (um)',labelpad=9);ax.set_ylabel('y (um)',labelpad=9);ax.set_zlabel('z (um)',labelpad=8)
        ax.set_box_aspect([1,1,1]);ax.view_init(elev=23,azim=-53);ax.set_title(title,fontsize=12,pad=14)
        fig.colorbar(sc,ax=ax,shrink=.6,pad=.12,label='Normalized mode amplitude')
    fig.suptitle('NUMERICAL / VARIATIONAL CRITICAL MODE',y=.96,fontsize=16)
    fig.text(.5,.88,'NOT PHYSICAL CRACK PREDICTION',ha='center',fontsize=14,color=ORANGE)
    finish(fig,'critical_mode_map','At 1.00% strain: a feasible one-sided negative-curvature witness, not a unique crack shape\nThe saved coordinates are unchanged. Both components are normalized independently; no force arrows')


def quality():
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True);fig.subplots_adjust(hspace=.3,bottom=.15)
    for rows,color,marker,name in [(A,BLUE,'o','A'),(B,ORANGE,'x','B')]:
        x=[r['strain']*100 for r in rows]
        axs[0].semilogy(x,[r['metrics']['force_residual'] for r in rows],marker,color=color,label='Schedule '+name)
        axs[1].plot(x,[r['metrics']['projected_phase_KKT'] for r in rows],marker,color=color,label='Schedule '+name)
    axs[0].axhline(1e-9,color=PINK,ls='--',label='Original force tolerance');axs[0].set_ylim(1e-17,3e-9)
    axs[1].set_yscale('symlog',linthresh=1e-16);axs[1].axhline(1e-8,color=PINK,ls='--',label='Phase KKT tolerance');axs[1].set_ylim(-1e-16,3e-8)
    axs[0].set_title('Unstable states still satisfy stationary equilibrium gates')
    axs[0].set_ylabel('Scaled force residual');axs[1].set_ylabel('Projected phase KKT residual');axs[1].set_xlabel('Engineering strain (%)')
    for ax in axs:shade(ax);clean(ax);ax.legend(fontsize=10)
    finish(fig,'stationarity_quality','Force, phase bounds, irreversibility, fixed boundaries and the original model-domain gates pass\nStationary equilibrium is checked independently of curvature classification')


def consistency():
    fig,axs=plt.subplots(1,3,figsize=(14,6));fig.subplots_adjust(bottom=.24,wspace=.4,top=.82)
    for i,name in enumerate(['A','B']):
        b=brackets[name];L,U=b['final_interval'];mid=(L+U)/2
        axs[0].errorbar(i,mid*100,yerr=(U-L)*50,fmt='o',color=BLUE,capsize=9)
        for state,offset,col,marker in [('lower',-.08,GREEN,'o'),('upper',.08,ORANGE,'v')]:
            r=b[state];axs[1].scatter(i+offset,r['metrics']['maximum_s'],color=col,marker=marker,s=60,label=state.capitalize()+' bracket endpoint' if i==0 else None)
            axs[2].scatter(i+offset,r['reaction_N']*1e6,color=col,marker=marker,s=60)
    axs[0].set_ylabel('Critical interval: engineering strain (%)');axs[0].set_title('Same interval')
    axs[1].set_yscale('symlog',linthresh=1e-14);axs[1].axhline(1e-8,color=PINK,ls='--',label='Damage detection tolerance')
    axs[1].set_ylim(-1e-14,3e-8);axs[1].set_ylabel('Maximum phase');axs[1].set_title('Phase remains below detection');axs[1].legend(fontsize=8.5,loc='center')
    axs[2].set_ylabel('Reaction (uN)');axs[2].set_title('Same reaction at matched endpoints')
    for ax in axs:ax.set_xticks([0,1],['Schedule A','Schedule B']);ax.set_xlim(-.4,1.4);clean(ax)
    fig.suptitle('Tracking resolution changes neither the bracket nor the matched states',fontsize=15)
    finish(fig,'schedule_consistency','A: strain step 5e-4 | B: strain step 2.5e-4 | Both refined to width 7.8125e-6\nMaximum phase at the upper endpoint is about 4.53e-14, far below the 1e-8 damage detection tolerance')


if __name__=='__main__':
    for plot in [stability,bracket_zoom,force_phase,active,landscape,critical_map,quality,consistency]:plot();print(plot.__name__,flush=True)
