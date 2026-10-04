"""Eight static review figures from saved diagnostics; no new model calculations."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from pd_clot.pfpd_candidate.coupling_stability_audit import OUT,read_json,read_npz,load_audit_state

FONT=Path('/mnt/c/Windows/Fonts/arial.ttf')
if not FONT.exists():raise RuntimeError('Arial font required for the existing visual style')
font_manager.fontManager.addfont(str(FONT))
plt.rcParams.update({'font.family':'Arial','font.size':12,'font.weight':'normal',
    'axes.titleweight':'normal','axes.labelweight':'normal','figure.facecolor':'black',
    'axes.facecolor':'black','savefig.facecolor':'black','text.color':'white','axes.labelcolor':'white',
    'axes.edgecolor':'#999999','xtick.color':'white','ytick.color':'white','grid.color':'#555555',
    'grid.alpha':.45,'legend.facecolor':'black','legend.edgecolor':'none','pdf.fonttype':42,
    'axes.spines.top':False,'axes.spines.right':False})
BLUE='#49b8ff';ORANGE='#ffb35b';PINK='#ec81b3';GREEN='#75dda5';GREY='#b5bec9'
F=OUT/'figures';F.mkdir(exist_ok=True)
paths=read_json(OUT/'energy_path/SAVED_PATHS.json')
failure=[r for r in paths['loading_unloading'] if r['load_step']==10]
phase=[r for r in failure if r['kind']=='phase_after']
extra=[r for r in paths['loading_unloading_extra_iterations'] if r['load_step']==8]
stable=[r for r in paths['loading_unloading'] if r['load_step']==4]
spectra=read_json(OUT/'spectrum/EPSILON_STABILITY.json')
alignment=read_json(OUT/'spectrum/MODE_ALIGNMENT.json')['rows']
comparison=read_json(OUT/'monolithic/ACCEPTANCE_AND_COMPARISON.json')
reference=read_json(OUT/'energy_path/REFERENCE_STATES.json')


def finish(fig,name,note=None):
    if note:fig.text(.5,.02,note,ha='center',va='bottom',fontsize=11,color=GREY)
    fig.savefig(F/(name+'.png'),dpi=180,bbox_inches='tight',pad_inches=.18)
    fig.savefig(F/(name+'.pdf'),bbox_inches='tight',pad_inches=.18)
    plt.close(fig)


def clean(ax):ax.grid(True,linewidth=.5);ax.tick_params(labelsize=11)


def amplification():
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True);fig.subplots_adjust(hspace=.28,bottom=.13)
    for rows,color,label in [(phase,BLUE,'Failed step: standard path'),([r for r in extra if r['kind']=='phase_after'],ORANGE,'Failed step: extra-iteration path'),(stable,GREEN,'Stable control A')]:
        x=[r['iteration'] for r in rows]
        axs[0].plot(x,[r['delta_s_L2'] for r in rows],'-o',ms=3,color=color,label=label+' (L2)')
        axs[0].plot(x,[r['phase_increment'] for r in rows],'--',color=color,alpha=.65)
        axs[1].plot(x,[np.nan if r['L2_ratio'] is None else r['L2_ratio'] for r in rows],'-o',ms=3,color=color,label=label)
        axs[1].plot(x,[np.nan if r['Linf_ratio'] is None else r['Linf_ratio'] for r in rows],'--',color=color)
    axs[0].set_yscale('symlog',linthresh=1e-10);axs[0].set_ylabel('Phase increment norm')
    axs[0].set_title('Phase increments grow after the first correction');axs[0].legend(fontsize=10)
    axs[1].axhline(1,color=GREY,ls=':');axs[1].set_ylim(-.08,2.15);axs[1].set_ylabel('Successive norm ratio');axs[1].set_xlabel('Stagger iteration')
    for ax in axs:clean(ax)
    finish(fig,'iteration_amplification','Solid: L2 | Dashed: maximum norm | Stable increments are exactly zero; their ratio is undefined')


def spectrum():
    fig,axs=plt.subplots(1,3,figsize=(14,5.8));fig.subplots_adjust(bottom=.22,wspace=.3)
    t=np.linspace(0,2*np.pi,400)
    for ax,key,color in zip(axs,['A','B','C'],[GREEN,ORANGE,BLUE]):
        b=read_npz(OUT/f'spectrum/{key}/eps_3e-06/MATRICES.npz');v=b['eigenvalues']
        ax.plot(np.cos(t),np.sin(t),color=GREY,ls='--',lw=1,label='Unit circle')
        ax.scatter(v.real,v.imag,s=25,color=color,alpha=.7,label='Eigenvalues')
        ax.set_aspect('equal');ax.set_xlim(-1.15,2.12);ax.set_ylim(-1.15,1.15);ax.set_xlabel('Real part');clean(ax)
        ax.set_title(f"{key}: "+('Feasible one-sided map' if key=='A' else 'Conditional frozen branch')+f"\nLargest magnitude: {abs(v[0]):.6f}",fontsize=12)
    axs[0].set_ylabel('Imaginary part')
    finish(fig,'fixed_point_spectrum','A: marginal (identity phase transport) | B and C: ordinary map spectrum UNRESOLVED\nBranch switching prevents interpreting these conditional eigenvalues as an ordinary Jacobian spectrum')


def observed_mode():
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True);fig.subplots_adjust(hspace=.28,bottom=.14)
    rows=alignment[2:];x=[r['iteration'] for r in rows];lam=rows[0]['conditional_lambda_magnitude']
    axs[0].plot(x,[r['L2_ratio'] for r in rows],color=BLUE,marker='o',ms=3,label='Observed L2 amplification')
    axs[0].plot(x,[r['Linf_ratio'] for r in rows],color=ORANGE,ls='--',label='Observed maximum-norm amplification')
    axs[0].axhline(lam,color=GREEN,ls=':',label=f'Conditional branch |lambda| = {lam:.6f}')
    axs[0].set_ylabel('Amplification');axs[0].set_title('Observed growth aligns with one selected branch');axs[0].legend(fontsize=10)
    # Plot the small angular deficit so near-unity alignment remains visible.
    axs[1].semilogy(x,[max(1-r['mode_cosine'],1e-15) for r in rows],'-o',ms=3,color=BLUE,label='1 - cosine with conditional mode')
    axs[1].semilogy(x,[max(1-r['successive_cosine'],1e-15) for r in rows],'--',color=ORANGE,label='1 - cosine between successive updates')
    axs[1].set_ylabel('Cosine deficit (smaller = closer)');axs[1].set_xlabel('Stagger iteration');axs[1].legend(fontsize=10)
    for ax in axs:clean(ax)
    finish(fig,'observed_vs_eigenmode','Directional agreement supports this branch interpretation; opposite-direction probes disagree\nAssociation is not a proof of a differentiable fixed-point map')


def energy():
    fig,axs=plt.subplots(2,2,figsize=(13,9));fig.subplots_adjust(hspace=.43,wspace=.3,bottom=.15)
    base=reference['previous_accepted_at_own_load']['total_J'];first=phase[0]['total_J'];x=np.arange(len(phase))
    ax=axs[0,0]
    ax.scatter([-2],[0],color=GREEN,label='Previous accepted, strain 0.013')
    ax.scatter([-1],[(comparison['previous_affine_lift_energy_J']-base)*1e12],color=GREY,label='Affine lift to strain 0.014')
    ax.plot(x,[(r['total_J']-base)*1e12 for r in phase],color=BLUE,label='Stagger phase-after, strain 0.014')
    ax.set_title('Energy relative to previous accepted state');ax.set_ylabel('Energy difference (pJ)');ax.set_xlabel('Stagger iteration (-2: previous, -1: lift)');ax.legend(fontsize=9)
    ax=axs[0,1]
    ax.plot(x,[(r['total_J']-first)*1e18 for r in phase],color=BLUE,label='Phase-after energy')
    ax.plot(x[1:],[(r['mechanical_end_energy_J']-first)*1e18 for r in phase[1:]],color=ORANGE,ls='--',label='Mechanical-end energy (iteration >= 1)')
    for method,col,marker in [('L-BFGS-B',GREEN,'o'),('SLSQP',PINK,'s')]:
        r=next(r for r in comparison['rows'] if r['method']==method and r['first_order_passed'] and r['strain']==.014)
        ax.scatter([25 if method=='L-BFGS-B' else 26],[r['relative_to_failed_first_phase_J']*1e18],color=col,marker=marker,label=method+' stationary saddle')
    ax.set_title('Same-load detail: energy falls without convergence');ax.set_ylabel('Energy relative to first phase state (aJ)');ax.set_xlabel('Iteration; optimizer endpoints at right');ax.legend(fontsize=9)
    ax=axs[1,0]
    ax.plot(x,[r['mechanical_energy_change_J']*1e18 for r in phase],color=ORANGE,label='Mechanical substep')
    ax.plot(x,[r['phase_energy_change_J']*1e18 for r in phase],color=BLUE,label='Phase substep')
    ax.set_yscale('symlog',linthresh=1e-7);ax.axhline(0,color=GREY,ls=':');ax.set_title('Substep energy changes');ax.set_ylabel('Energy change (aJ)');ax.set_xlabel('Stagger iteration');ax.legend(fontsize=10)
    ax.set_yticks([-1e5,-1e3,-10,-.1,-.001,-1e-5,0]);ax.minorticks_off()
    ax=axs[1,1];rows=[r for r in comparison['rows'] if r['strain']==.014]
    for i,r in enumerate(rows):
        ax.scatter(i,r['relative_to_failed_first_phase_J']*1e18,marker='x' if r['domain_exit'] else 'o',color=ORANGE if r['domain_exit'] else GREEN,s=55)
    ax.axhline(0,color=GREY,ls=':');ax.set_yscale('symlog',linthresh=.05);ax.set_xticks(range(8),['L1','L2','L3','L4','S1','S2','S3','S4'])
    ax.set_yticks([-10,-1,0,1,10,100,1000]);ax.minorticks_off()
    ax.set_title('All joint optimizer attempts');ax.set_ylabel('Energy relative to first phase state (aJ)');ax.set_xlabel('L: L-BFGS-B | S: SLSQP | 1-4: initial states')
    for ax in axs.flat:clean(ax)
    finish(fig,'monolithic_vs_staggered_energy','Circles: stationary saddles | Crosses: last valid evaluations before domain exit, not solutions\nThe old-load subtraction is a plotting baseline; stability comparisons use the same target load')


def curvature():
    fig,ax=plt.subplots(figsize=(11,6));fig.subplots_adjust(bottom=.18)
    curv=read_json(OUT/'hessian/EPSILON_STABILITY.json')
    for key,strain,label,color,marker in [('A',.009,'A: stable control',GREEN,'o'),('B',.013,'B: last accepted',ORANGE,'o'),('C',.014,'C: early failed iterate',BLUE,'o'),('L-BFGS-B_failed_iteration_000_phase_after',.014,'L-BFGS-B stationary saddle',PINK,'x'),('SLSQP_failed_iteration_000_phase_after',.014,'SLSQP stationary saddle',GREY,'+')]:
        vals=np.array(curv[key]['minimum_eigenvalues']);med=np.median(vals)
        ax.errorbar(strain*100,med,yerr=[[med-min(vals)],[max(vals)-med]],fmt=marker,color=color,ms=10,capsize=6,label=label)
    ax.axhline(0,color='white',ls='--');ax.set_xlim(.85,1.46);ax.set_ylim(-.055,.075);ax.set_xlabel('Engineering strain (%)');ax.set_ylabel('Minimum reduced Hessian eigenvalue (scaled)')
    ax.set_title('Negative coupled curvature is already present at the last accepted state');ax.legend(fontsize=11);clean(ax)
    finish(fig,'smallest_hessian_eigenvalue','Three finite-difference scales agree; error bars are smaller than markers\nOnly sampled loads are shown. The exact onset of instability has not been located')


def mode_map():
    p,d=load_audit_state('B');b=read_npz(OUT/'hessian/B/eps_3e-06/BLOCKS.npz');v=b['minimum_mode_scaled']
    phase_mode=v[p.nx:];aggregate=np.zeros(len(p.family.X));den=np.zeros_like(aggregate)
    for j,(i,k) in enumerate(p.family.pairs):
        aggregate[i]+=p.weights[j]*phase_mode[j];aggregate[k]+=p.weights[j]*phase_mode[j];den[[i,k]]+=p.weights[j]
    aggregate/=den;aggregate/=max(abs(aggregate))
    displacement=np.zeros_like(p.family.X);displacement[p.free]=v[:p.nx].reshape(-1,3)*p.length_scale
    magnitude=np.linalg.norm(displacement,axis=1);magnitude/=magnitude.max();points=(d['x']-p.family.X.min(axis=0))*1e6
    fig=plt.figure(figsize=(13,7));fig.subplots_adjust(left=.02,right=.96,bottom=.14,top=.78,wspace=.16)
    for i,(values,cmap,title,lim) in enumerate([(aggregate,'coolwarm','Phase mode: weighted bond-to-particle average',(-1,1)),(magnitude,'viridis','Displacement mode magnitude',(0,1))]):
        ax=fig.add_subplot(1,2,i+1,projection='3d');ax.set_facecolor('black');ax.set_proj_type('ortho')
        sc=ax.scatter(*points.T,c=values,cmap=cmap,vmin=lim[0],vmax=lim[1],s=190,edgecolors='#bbbbbb',linewidths=.4,depthshade=False)
        for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
            axis.set_pane_color((0,0,0,1));axis._axinfo['grid']['color']=(.35,.35,.35,.5)
        ax.set_xlabel('x (um)',labelpad=9);ax.set_ylabel('y (um)',labelpad=9);ax.set_zlabel('z (um)',labelpad=8)
        ax.set_box_aspect([1,1,1]);ax.view_init(elev=23,azim=-53);ax.set_title(title,fontsize=12,pad=14)
        fig.colorbar(sc,ax=ax,shrink=.6,pad=.12,label='Normalized mode amplitude')
    fig.suptitle('Numerical stability mode at B: negative Hessian eigenvector',y=.96,fontsize=16)
    fig.text(.5,.88,'NUMERICAL STABILITY MODE  |  NOT PHYSICAL CRACK PREDICTION',ha='center',fontsize=13,color=ORANGE)
    finish(fig,'unstable_mode_map','Both panels show the same coupled mode, independently normalized for legibility\nCoordinates are the saved state; no displaced geometry or force arrows are shown')


def subdivision():
    fig,axs=plt.subplots(1,3,figsize=(14,6));fig.subplots_adjust(bottom=.24,wspace=.33,top=.78)
    mono=next(r for r in comparison['rows'] if r['strain']==.014 and r['method']=='L-BFGS-B' and r['first_order_passed'])
    metrics=[('Endpoint maximum phase',mono['metrics']['maximum_s']),('Endpoint energy (pJ)',mono['metrics']['total_J']*1e12),('State distance to reference\nscaled L2 norm',0)]
    for ax,(label,value) in zip(axs,metrics):
        ax.scatter([4],[value],color=ORANGE,s=70);ax.set_xticks(range(5),['1 step','2 steps','4 steps','8 steps','Joint\nreference']);ax.set_ylabel(label)
        low,high=ax.get_ylim();span=max(abs(value)*.2,.1);ax.set_ylim(value-span,value+span)
        for i in range(4):ax.text(i,.5,'N/A',transform=ax.get_xaxis_transform(),ha='center',color=GREY)
        ax.set_xlim(-.5,4.6);clean(ax)
    fig.suptitle('DIAGNOSTIC ONLY: no staggered path reaches the same endpoint',fontsize=16)
    fig.text(.5,.82,'Target strain 0.014 | All subdivisions stop before accepting their first substep',ha='center',fontsize=12)
    finish(fig,'substep_diagnostic','N/A: endpoint state does not exist; partial states at different loads are not substituted\nThe joint reference is a stationary saddle, not a stable minimum. Its self-distance is zero\n2 / 4 / 8 subdivisions attempt strain 0.0135 / 0.01325 / 0.013125; all hit the unchanged pure-tension gate')


def ordering():
    fig,axs=plt.subplots(2,1,figsize=(11,8),sharex=True);fig.subplots_adjust(hspace=.28,bottom=.14)
    orders=read_json(OUT/'energy_path/TENSION_GATE_ORDERING.json')
    for rows,key,color,label in [(failure,'loading_unloading',BLUE,'Standard path'),(extra,'loading_unloading_extra_iterations',ORANGE,'Extra-iteration path')]:
        ph=[r for r in rows if r['kind']=='phase_after']
        axs[0].semilogy([r['iteration'] for r in ph],[r['phase_increment'] for r in ph],'-o',ms=3,color=color,label=label)
        axs[1].plot([r['iteration'] for r in rows],[r['minimum_Green_strain'] for r in rows],'-o',ms=3,color=color,label=label)
        o=orders[key];axs[0].axvline(o['first_clearly_amplifying_iteration'],color=color,ls=':',alpha=.8)
        axs[1].axvline(o['first_domain_exit_iteration'],color=color,ls=':',alpha=.8)
    axs[0].set_ylabel('Maximum phase increment');axs[0].set_title('Phase amplification appears before the domain exit');axs[0].legend()
    axs[1].set_yscale('symlog',linthresh=1e-17);axs[1].axhline(-1e-9,color=PINK,ls='--',label='Pure-tension threshold: -1e-9')
    axs[1].set_ylim(-5e-9,1e-16)
    axs[1].set_ylabel('Minimum Green strain');axs[1].set_xlabel('Stagger iteration');axs[1].legend(fontsize=10)
    for ax in axs:clean(ax)
    finish(fig,'tension_gate_ordering','Clear amplification: iterations 4 / 6 | First domain exit: iterations 24 / 27\nTemporal ordering supports a later domain violation; it does not establish causality')


if __name__=='__main__':
    for plot in [amplification,spectrum,observed_mode,energy,curvature,mode_map,subdivision,ordering]:
        plot();print(plot.__name__,flush=True)
