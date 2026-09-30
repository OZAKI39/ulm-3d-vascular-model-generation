"""Scientific study plots. Values are never altered to imply fragmentation."""
import argparse,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.patches import Rectangle,Circle
from .analysis import collect

ROOT=Path(__file__).resolve().parents[2]
COLORS=['#66e0c2','#e9b44c','#c69bff','#ff879c']


def style():
    font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/arial.ttf')
    plt.rcParams.update({'font.family':'Arial','font.weight':'normal','font.size':11,
        'figure.facecolor':'black','axes.facecolor':'black','savefig.facecolor':'black',
        'text.color':'white','axes.labelcolor':'white','axes.edgecolor':'#647487',
        'xtick.color':'white','ytick.color':'white','grid.color':'#263545','axes.titleweight':'normal',
        'axes.labelweight':'normal','legend.frameon':False,'legend.labelcolor':'white'})


def create(root,output):
    root=Path(root);out=Path(output);out.mkdir(parents=True,exist_ok=True);style()
    data=collect(root,out/'data');rows={r['run']:r for r in data['runs']}
    history={n:json.loads((root/n/'history.json').read_text()) for n in rows}

    def save(fig,name):
        fig.savefig(out/f'{name}.png',dpi=180,bbox_inches='tight')
        fig.savefig(out/f'{name}.pdf',bbox_inches='tight');plt.close(fig)

    mesh=[n for n in ['coarse','medium','fine'] if n in rows]
    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for k,name in enumerate(mesh):
        cal=json.loads((root/name/'CALIBRATION.json').read_text());h=cal['history'];color=COLORS[k]
        axes[0].plot(cal['spacing_m']*1e3,cal['Gc_measured_J_m2'],'o',color=color,label=name.capitalize())
        N=[a['coupon_cycle'] for a in h]
        axes[1].plot(N,[a['damage_dissipation_J']*1e9 for a in h],color=color,label=f'{name.capitalize()} dissipation')
        axes[1].plot(N,[a['external_work_J']*1e9 for a in h],color=color,ls='--',alpha=.6)
    axes[0].axhline(.01,color='white',ls='--',label='Demo target')
    axes[0].set(xlabel='Particle spacing (mm)',ylabel='Measured Gc (J/m²)',title='Fracture energy calibration',ylim=(.0095,.0105))
    axes[1].set(xlabel='Coupon cycle',ylabel='Energy (nJ)',title='Dissipation and force-integrated work')
    for ax in axes:ax.grid(alpha=.5);ax.legend(fontsize=9)
    save(fig,'fracture_energy_calibration')

    metrics=[('attached_volume_fraction','Attached volume fraction',1),('detached_volume_fraction','Detached volume fraction',1),
             ('mean_damage','Mean damage',1),('damage_dissipation_estimate_J','Damage dissipation (J)',1),
             ('fraction_of_damage_inside_influence_region','Damage inside influence region',1),('maximum_displacement_m','Maximum displacement (µm)',1e6)]
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout='constrained')
    for ax,(key,label,factor) in zip(axes.flat,metrics):
        ax.plot([rows[n]['spacing_mm'] for n in mesh],[rows[n][key]*factor for n in mesh],'o-',color=COLORS[0]);ax.set(xlabel='Particle spacing (mm)',ylabel=label);ax.grid(alpha=.5)
    fig.suptitle('Particle-spacing sensitivity — failure regime not established',fontsize=15)
    save(fig,'mesh_convergence')

    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for k,name in enumerate(['coarse','cycle_500','cycle_250']):
        if name not in history:continue
        h=history[name];N=[r['represented_cycles'] for r in h]
        for ax,(key,label,_) in zip(axes.flat,[metrics[2],metrics[3],metrics[0],metrics[1]]):
            ax.plot(N,[r[key] for r in h],color=COLORS[k],label=f'ΔN = {rows[name]["DeltaN"]}');ax.set(xlabel='Represented cycles',ylabel=label);ax.grid(alpha=.5)
    axes[0,0].legend();fig.suptitle('Cycle-jump sensitivity',fontsize=15);save(fig,'cycle_jump_convergence')

    fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
    for k,name in enumerate(['coarse','dt_half']):
        if name not in history:continue
        h=history[name];N=[r['represented_cycles'] for r in h]
        for ax,key,label,factor in zip(axes,['mean_damage','maximum_displacement_m','relative_numerical_energy_residual'],['Mean damage','Maximum displacement (µm)','Signed energy residual (%)'],[1,1e6,100]):
            ax.plot(N,[r[key]*factor for r in h],color=COLORS[k],label=f'dt = {rows[name]["dt_s"]*1e6:g} µs');ax.set(xlabel='Represented cycles',ylabel=label);ax.grid(alpha=.5)
    axes[0].legend();fig.suptitle('Time-step sensitivity — no failure event in this exposure',fontsize=15);save(fig,'time_step_sensitivity')

    h=history['coarse'];N=[r['represented_cycles'] for r in h]
    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for key,label,color in [('elastic_energy_J','Elastic',COLORS[0]),('stabilization_energy_J','Stabilization',COLORS[1]),('kinetic_energy_J','Kinetic',COLORS[2])]:
        axes[0,0].plot(N,[r[key] for r in h],color=color,label=label)
    for key,label,color in [('external_work_J','Surface work',COLORS[0]),('damping_dissipation_J','Damping',COLORS[1]),('transport_work_J','Transport work',COLORS[2])]:
        axes[0,1].plot(N,[r[key] for r in h],color=color,label=label)
    axes[1,0].plot(N,[r['damage_dissipation_estimate_J'] for r in h],color=COLORS[0]);axes[1,0].set_title('Accumulated damage dissipation')
    axes[1,1].plot(N,[r['numerical_energy_residual_J'] for r in h],color=COLORS[1]);axes[1,1].set_title('Signed numerical energy residual')
    for ax in axes.flat:ax.set(xlabel='Represented cycles',ylabel='Energy (J)');ax.grid(alpha=.5)
    axes[0,0].legend();axes[0,1].legend();fig.suptitle('Energy on the simulated representative mechanical path',fontsize=15);save(fig,'energy_history')

    fig,axes=plt.subplots(1,3,figsize=(15,4.5),layout='constrained')
    for k,name in enumerate(mesh):
        h=history[name];N=[r['represented_cycles'] for r in h]
        axes[0].plot(N,[r['fraction_of_damage_inside_influence_region'] for r in h],color=COLORS[k],label=name.capitalize())
        axes[1].plot(N,[None if r['damage_weighted_distance_from_source_m'] is None else r['damage_weighted_distance_from_source_m']*1e3 for r in h],color=COLORS[k])
    axes[0].set(xlabel='Represented cycles',ylabel='Fraction of damage inside region',ylim=(0,1));axes[0].legend()
    axes[1].set(xlabel='Represented cycles',ylabel='Damage-weighted distance (mm)')
    c=json.loads((root/'coarse/CONFIG.json').read_text());z=np.load(root/'coarse/states.npz');D=z['damage'][-1].reshape(c['clot']['cells']).mean(axis=1).T
    origin=np.array(c['clot']['origin_m'])*1e3;size=np.array(c['clot']['cells'])*c['clot']['particle_spacing_m']*1e3;center=np.array(c['streaming']['bubble_center_m'])*1e3
    im=axes[2].imshow(D,origin='lower',extent=[origin[0],origin[0]+size[0],origin[2],origin[2]+size[2]],cmap='inferno',vmin=0,aspect='equal')
    axes[2].plot(center[0],center[2],'D',mfc='none',mec='white');axes[2].add_patch(Circle(center[[0,2]],c['localization']['influence_radius_m']*1e3,fill=False,color='white',alpha=.4))
    axes[2].set(xlabel='X (mm)',ylabel='Z (mm)',title='Y-averaged damage, coarse',xlim=(-1.4,.8),ylim=(-1.4,-.05));fig.colorbar(im,ax=axes[2],shrink=.7,label='Actual damage')
    save(fig,'localization_metrics')

    fig,axes=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    latest=json.loads((root/'coarse/components.json').read_text())[-1]['components'];free=[r for r in latest if not r['attached_to_base']]
    if free:
        axes[0].hist([r['particle_count'] for r in free]);axes[0].set(xlabel='Particles per component',ylabel='Component count')
    else:
        axes[0].text(.5,.5,'No detached components\nwithin 25,000 represented cycles\n\nSize distribution: not defined',ha='center',va='center',transform=axes[0].transAxes);axes[0].set_xticks([]);axes[0].set_yticks([])
    axes[0].set_title('New regularized case')
    legacy=json.loads((ROOT/'verification/regularization/LEGACY_RESOLUTION_DIAGNOSTIC.json').read_text())
    old=[r for r in legacy['components'] if not r['attached_to_base']];counts=np.array([r['particle_count'] for r in old]);unique,number=np.unique(counts,return_counts=True)
    axes[1].bar(np.arange(len(unique)),number,color=COLORS[3]);axes[1].set_xticks(np.arange(len(unique)),[str(v) for v in unique]);axes[1].set(xlabel='Particles per numerical component',ylabel='Count',title='Legacy forced-failure case: different loading')
    for k,n in enumerate(number):axes[1].text(k,n+2,str(n),ha='center')
    save(fig,'fragment_size_distribution')
    distributions={name:json.loads((root/name/'components.json').read_text())[-1] for name in mesh}
    (out/'data/fragment_size_distributions.json').write_text(json.dumps(distributions,indent=2)+'\n')

    fig,axes=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    for ax,key,label in zip(axes.flat,['Gc_demo_J_m2','m_E','C_E','streaming_amplitude_m_s'],['Gc demo (J/m²)','Damage exponent','Cyclic coefficient','Streaming amplitude (m/s)']):
        prefix='streaming_velocity_scale_m_s' if key=='streaming_amplitude_m_s' else key
        values=[rows['coarse']]+[r for n,r in rows.items() if n.startswith('sweep_'+prefix+'_')]
        values.sort(key=lambda r:r[key]);ax.plot([r[key] for r in values],[r['mean_damage'] for r in values],'o-',color=COLORS[0]);ax.set(xlabel=label,ylabel='Final mean damage');ax.grid(alpha=.5)
    fig.suptitle('Parameter sensitivity: all sampled outcomes retained',fontsize=15);save(fig,'parameter_sensitivity')
    return data


def main():
    p=argparse.ArgumentParser();p.add_argument('--runs',type=Path,default=ROOT/'results/streaming_regularized_demo');p.add_argument('--output',type=Path,default=ROOT/'visualization/streaming_regularized_demo/figures');a=p.parse_args();create(a.runs,a.output)


if __name__=='__main__':main()
