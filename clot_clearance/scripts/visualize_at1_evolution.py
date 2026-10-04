"""Seven review figures from saved data, including the failed coupled gates."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from scripts.verify_at1_evolution import ROOT, OUT, config, write_json

FIG = OUT/'figures'
RESULTS = ROOT/'results/pfpd_at1_evolution_v1'
BLUE, ORANGE, GREEN, RED = '#58c8ff','#ffbb65','#76dfae','#ff777f'


def read(path): return json.loads(Path(path).read_text())


def style():
    arial = Path('/mnt/c/Windows/Fonts/arial.ttf')
    if not arial.exists(): raise FileNotFoundError('Requested Arial font is unavailable')
    font_manager.fontManager.addfont(str(arial))
    plt.rcParams.update({'font.family':'Arial','font.size':13,'font.weight':'normal',
        'axes.titleweight':'normal','axes.labelweight':'normal','axes.titlesize':16,
        'figure.facecolor':'black','axes.facecolor':'black','savefig.facecolor':'black',
        'text.color':'white','axes.labelcolor':'white','axes.edgecolor':'#999999',
        'xtick.color':'white','ytick.color':'white','grid.color':'#555555','grid.alpha':.35,
        'axes.spines.top':False,'axes.spines.right':False,'legend.facecolor':'black',
        'legend.edgecolor':'none','pdf.fonttype':42,'figure.constrained_layout.use':True})
    FIG.mkdir(exist_ok=True)
    write_json(FIG/'STYLE.json',dict(background='black',font='Arial',weight='regular',
        arrows=False,displacement_scale=1,scientific_arrays_modified=False,
        phase_particle_aggregation='BAQP pair-energy weights to both endpoints; visualization only'))


def finish(fig,name):
    for ax in fig.axes:
        if ax.name != '3d': ax.grid(True,alpha=.25)
    fig.savefig(FIG/(name+'.png'),dpi=200)
    fig.savefig(FIG/(name+'.pdf'))
    plt.close(fig)


def main():
    style(); c = config()
    fig,ax = plt.subplots(figsize=(10,6))
    s = np.linspace(0,1,601)
    for ratio,color in zip([.5,1.,1.5,3.],[BLUE,ORANGE,GREEN,RED]):
        ax.plot(s,(1-s)**2*ratio+2*s,color=color,label=f'Y/Yc = {ratio:g}')
        minimum = 0. if ratio <= 1 else 1-1/ratio
        ax.scatter([minimum],[(1-minimum)**2*ratio+2*minimum],color=color,s=65,zorder=5)
    ax.scatter([],[],color='white',s=55,label='Constrained minimum')
    ax.set(xlabel='Fracture phase s',ylabel='Candidate energy density / Yc',
           title='AT1 energy landscape',xlim=(-.02,1.02))
    ax.legend(ncol=2,fontsize=11)
    finish(fig,'at1_energy_landscape')

    rows = [r for r in read(OUT/'minimization/CASES.json') if r['s_old']==0]
    fig,axes = plt.subplots(1,2,figsize=(13,5.5))
    for ax,maximum in zip(axes,[3,100]):
        r = np.linspace(0,maximum,3001)
        analytic = np.where(r>1,1-1/np.maximum(r,1),0.)
        ax.plot(r,analytic,color=BLUE,label='AT1 analytic')
        rr = [v for v in rows if v['ratio']<=maximum]
        ax.scatter([v['ratio'] for v in rr],[v['selected_s'] for v in rr],
                   facecolor='none',edgecolor=ORANGE,s=58,label='SciPy bounded minimizer',zorder=5)
        ax.axvline(1,color=GREEN,ls='--',label='Threshold = 1')
        ax.set(xlabel='Undamaged energy / Yc',ylabel='Virgin fracture phase s',ylim=(-.035,1.02))
        if maximum>3:ax.set_xscale('symlog',linthresh=1);ax.set_xlim(0,100)
    axes[0].legend(fontsize=10,loc='lower right')
    axes[0].set_title('Damage starts above the threshold')
    axes[1].set_title('Full scalar range')
    finish(fig,'at1_analytic_vs_numeric')

    bond = read(OUT/'analytic/SINGLE_BOND_PATH.json'); steps = [r['step'] for r in bond]
    yc = c['phase']['Yc_design_Pa']
    fig,axes = plt.subplots(3,1,figsize=(12,10),sharex=True)
    axes[0].plot(steps,[r['ratio'] for r in bond],'-o',color='white')
    axes[0].axhline(1,color=GREEN,ls='--');axes[0].set_ylabel('Driver Y/Yc')
    axes[1].plot(steps,[r['s_new'] for r in bond],'-o',color=BLUE,label='Irreversible phase')
    axes[1].set_ylabel('Fracture phase s');axes[1].legend(loc='upper left')
    for key,label,color in [('elastic_after','Elastic',BLUE),('crack_after','Crack',ORANGE),('total_after','Total',GREEN)]:
        axes[2].plot(steps,[r[key]/yc for r in bond],'-o',color=color,label=label)
    axes[2].set(xlabel='Single-bond load step',ylabel='Energy density / Yc');axes[2].legend(ncol=3)
    for ax in axes: ax.axvspan(8,12,color='#cccccc',alpha=.1)
    fig.suptitle('Single-bond loading, unloading and reloading\nCoupled cube unloading was not reached',fontsize=16)
    finish(fig,'loading_unloading_irreversibility')

    fixed = read(OUT/'fixed_x/SUMMARY.json'); e = fixed['energies_J']
    fig,axes = plt.subplots(1,2,figsize=(12,5.5))
    changes = [(e[k+'_after']-e[k+'_before'])*1e12 for k in ['elastic','crack','total']]
    axes[0].bar(['Elastic','Crack','Total'],changes,color=[BLUE,ORANGE,GREEN])
    axes[0].axhline(0,color='white',lw=.7);axes[0].set_ylabel('Energy change (pJ)')
    axes[0].set_title('Fixed-position phase update')
    names=['K','Fbar','Ftilde'];values=[fixed['maximum_changes'][k] for k in names]
    axes[1].scatter(names,values,color=BLUE,s=70)
    axes[1].set(ylim=(-.1,1.),ylabel='Maximum absolute change',title='Kinematics are bitwise identical')
    for i in range(3):axes[1].text(i,.05,'0',ha='center',fontsize=18)
    axes[1].text(1,.55,'Shape derivatives also unchanged',ha='center',fontsize=12)
    finish(fig,'fixed_x_energy_audit')

    fig,axes = plt.subplots(2,1,figsize=(12,8),sharex=True)
    for name,label,color in [('tensile_cube','Monotonic cube: last accepted step',GREEN),
                              ('loading_unloading','Loading path: failed step',RED),
                              ('loading_unloading_extra_iterations','Extra iterations: failed step',ORANGE)]:
        history = read(RESULTS/name/'STAGGERED_HISTORY.json'); last = history[-1]['load_step']
        subset = [r for r in history if r['load_step']==last]
        ix = [r['iteration'] for r in subset]
        axes[0].semilogy(ix,np.maximum([r['relative_force_residual'] for r in subset],1e-18),'-o',
                         ms=3,color=color,label=label)
        axes[1].semilogy(ix,np.maximum([r['maximum_delta_s'] for r in subset],1e-18),'-o',ms=3,color=color)
    axes[0].axhline(c['mechanical_equilibrium']['relative_force_residual_tolerance'],color='white',ls='--',label='Acceptance tolerance')
    axes[1].axhline(c['tolerances']['phase_step'],color='white',ls='--')
    axes[0].set(ylabel='Post-phase force residual',title='Staggered convergence and the failed robustness check')
    axes[0].legend(fontsize=10,loc='upper left')
    axes[1].set(xlabel='Staggered iteration within the load step',ylabel='Maximum phase increment')
    finish(fig,'staggered_convergence')

    mono = read(RESULTS/'tensile_cube/LOAD_PATH.json'); lur = read(RESULTS/'loading_unloading/LOAD_PATH.json')
    fig,axes = plt.subplots(2,1,figsize=(11,9),sharex=True)
    disp = np.array([r['displacement_m'] for r in mono])*1e6
    axes[0].plot(disp,[r['reaction_N']*1e6 for r in mono],'-o',color=BLUE,label='Monotonic cube: accepted states')
    axes[0].plot([r['displacement_m']*1e6 for r in lur],[r['reaction_N']*1e6 for r in lur],':s',
                 color=ORANGE,label='Loading path: accepted prefix only',ms=4)
    onset = next(r for r in mono if r['maximum_s']>c['tolerances']['phase_step'])
    axes[0].scatter([onset['displacement_m']*1e6],[onset['reaction_N']*1e6],s=150,
                    facecolor='none',edgecolor=GREEN,label='First sampled damage above tolerance',zorder=5)
    axes[0].set(ylabel='Grip reaction (µN)',title='Tensile response with demonstration fracture energy')
    axes[0].legend(fontsize=10)
    axes[1].plot(disp,[r['maximum_s'] for r in mono],'-o',color=BLUE,label='Maximum bond phase')
    axes[1].plot(disp,[r['mean_s'] for r in mono],'--s',color=ORANGE,label='Weighted mean bond phase',ms=4)
    axes[1].axhline(c['phase']['sc'],color=RED,ls='--',label='High-damage gate')
    axes[1].set(xlabel='Prescribed displacement (µm)',ylabel='Fracture phase s',ylim=(-.04,1.03))
    threshold_disp=c['phase']['threshold_engineering_strain']*c['geometry']['grip_center_separation_m']*1e6
    for ax in axes:ax.axvline(threshold_disp,color=GREEN,ls=':',lw=1)
    axes[1].legend(fontsize=10,loc='lower right')
    finish(fig,'tensile_force_displacement')

    last = mono[-1]['load_step']
    with np.load(RESULTS/'tensile_cube'/f'load_{last:03d}/accepted_state.npz') as state:
        x,phase = state['x'],state['s_pf_bond']
    with np.load(RESULTS/'tensile_cube/reference.npz') as ref:
        pairs,weights,X = ref['pairs'],ref['weights'],ref['X']
    total,denom = np.zeros(len(x)),np.zeros(len(x))
    for end in [0,1]:
        np.add.at(total,pairs[:,end],weights*phase);np.add.at(denom,pairs[:,end],weights)
    point_phase=total/denom
    np.savez_compressed(FIG/'phase_map_data.npz',x_actual=x,X=X,s_pf_bond=phase,
                        D_pf_particle=point_phase,weights=weights)
    fig=plt.figure(figsize=(10,8));ax=fig.add_subplot(projection='3d')
    scatter=ax.scatter(x[:,0]*1e3,x[:,1]*1e3,x[:,2]*1e3,c=point_phase,cmap='inferno',
                       vmin=0,vmax=c['phase']['sc'],s=180,depthshade=False,edgecolors='#bbbbbb',linewidths=.3)
    ax.set(xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)')
    ax.set_box_aspect(np.ptp(x,axis=0));ax.view_init(elev=22,azim=-58)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
        axis.set_pane_color((0,0,0,0));axis._axinfo['grid'].update(color=(.3,.3,.3,.5),linewidth=.5)
    fig.colorbar(scatter,ax=ax,shrink=.7,pad=.1,label='Weighted particle phase (visualization only)')
    fig.suptitle('Last accepted state before the high-damage gate\nActual displacement scale = 1×\nDemonstration fracture parameter · Not clot calibration',fontsize=15)
    finish(fig,'phase_map_before_high_damage_gate')
    print('Saved seven PNG/PDF pairs; failed coupled gates are shown explicitly.')


if __name__ == '__main__': main()
