"""Render signed audit data without rescaling displacement or changing mechanics."""
from pathlib import Path
import json
import hashlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, colors

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'verification/damage_energy_audit_v1'


def read(p):
    return json.loads(p.read_text())


def write(p, obj):
    with p.open('x') as f:
        json.dump(obj, f, indent=2, allow_nan=False)


def main():
    s = read(OUT / 'comparisons/MECHANISM_SUMMARY.json')
    plan = read(OUT / 'provenance/PREDECLARED_PLAN.json')
    # Preserve the first analysis, including its over-strict near-zero summary flag.
    rows = [read(p) for p in (OUT / 'fixtures').glob('*/uniform_*/SUMMARY.json')
            if p.parent.parent.name != 'pattern_study']
    zero = [r for r in rows if r['deformation'] == 'undeformed']
    positive = [r for r in rows if r['deformation'] != 'undeformed']
    zero_tol = plan['undeformed_absolute_tolerance_J']
    correction = dict(
        original_uniform_monotonic_flag=s['fixtures']['uniform_monotonic_with_roundoff'],
        explanation='Original summary applied a relative positive-energy criterion to the zero-energy state. Three undeformed signed changes are negative roundoff. Raw data and original summary are retained; all original declared tolerances are unchanged.',
        zero_state_declared_absolute_tolerance_J=zero_tol,
        zero_state_stays_zero=all(max(abs(r['U_before_J']), abs(r['U_after_J'])) < zero_tol for r in zero),
        zero_state_negative_changes_J=[r['delta_U_total_J'] for r in zero if r['delta_U_total_J'] < 0],
        positive_energy_fixture_count=len(positive),
        positive_energy_uniform_monotonic=all(r['delta_U_total_J'] >= 0 for r in positive),
        positive_energy_maximum_relative_scaling_residual=max(abs(r['uniform_scaling_residual_J']) / abs(r['U_before_J']) for r in positive),
        positive_energy_scaling_tolerance=plan['uniform_energy_relative_tolerance'],
        conclusion='Uniform reduction is dissipative for positive-energy full-rank fixtures; undeformed energy remains zero within the predeclared zero-state tolerance.')
    interpretation = OUT / 'comparisons/UNIFORM_SCALING_INTERPRETATION.json'
    if interpretation.exists():
        assert read(interpretation) == correction
    else:
        write(interpretation, correction)

    font = Path('/mnt/c/Windows/Fonts/arial.ttf')
    font_manager.fontManager.addfont(str(font))
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 13, 'font.weight': 'normal',
        'axes.titleweight': 'normal', 'axes.labelweight': 'normal', 'figure.facecolor': '#000000',
        'axes.facecolor': '#000000', 'savefig.facecolor': '#000000', 'text.color': 'white',
        'axes.labelcolor': 'white', 'xtick.color': '#dddddd', 'ytick.color': '#dddddd',
        'axes.edgecolor': '#888888', 'grid.color': '#353535', 'grid.alpha': .65,
        'pdf.fonttype': 42, 'ps.fonttype': 42, 'legend.frameon': False})
    manifest = []

    def setup(title, ylabel, xlabel, size=(10.5, 6.5)):
        fig, ax = plt.subplots(figsize=size, layout='constrained')
        ax.set_title(title, pad=18)
        ax.set(xlabel=xlabel, ylabel=ylabel)
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(True)
        return fig, ax

    def save(fig, name, source, note):
        fig.canvas.draw()
        for ext in ['png', 'pdf']:
            p = OUT / 'figures' / f'{name}.{ext}'
            if p.exists():
                raise FileExistsError(p)
            fig.savefig(p, dpi=180)
        manifest.append(dict(name=name, source=source, note=note))
        plt.close(fig)

    actual = read(OUT / 'actual_run_replay/damage_jump_audit.json')
    fig, axes = plt.subplots(1, 3, figsize=(16, 6.2), layout='constrained', sharey=True)
    series = [('softening', 'Softening', '#d7d7d7'), ('reconstruction', 'Reconstruction', '#e79d57'),
              ('stabilization', 'Stabilization', '#699cbc'), ('total', 'Total', '#ffffff')]
    for ax, case, dt in zip(axes, ['dt_erosion_0p02', 'dt_erosion_0p01', 'dt_erosion_0p005'], ['0.02', '0.01', '0.005']):
        data = sorted([r for r in actual if r['case'] == case], key=lambda r: r['erosion_time'])
        for key, label, color in series:
            ax.plot([r['erosion_time'] for r in data], [r[f'delta_U_{key}_J'] / 1e-18 for r in data],
                    marker='o' if key == 'total' else '.', markersize=4, color=color, label=label,
                    linestyle='--' if key == 'total' else '-', linewidth=1.8)
        ax.axhline(0, color='#999999', linewidth=.7)
        ax.set_title(f'Slow step {dt} s')
        ax.set_xlabel('Erosion exposure time (s)')
        ax.grid(True)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel('Signed energy change per jump (10^-18 J)')
    axes[1].legend(loc='upper right', fontsize=11)
    fig.suptitle('Reconstruction outweighs softening and stabilization\nNegative values mean stored energy increases', fontsize=16)
    save(fig, 'energy_decomposition_actual_jumps', 'actual_run_replay/damage_jump_audit.json', '35 jumps; each protocol is separate, not cumulative across protocols.')

    fig, ax = setup('The actual damage pattern raises energy at every tested amplitude',
                    'Fixed-position signed energy change (J)', 'Maximum bond damage increment')
    for key, label, color, marker in [('uniform','Uniform','#dddddd','o'), ('one_bond','Single bond','#aebca6','s'),
            ('coherent_patch','Local patch','#699cbc','^'), ('actual_pattern','Actual pattern','#e79d57','D')]:
        data = sorted([r for r in s['patterns'] if r['pattern'] == key], key=lambda r:r['damage_amplitude'])
        ax.plot([r['damage_amplitude'] for r in data], [r['delta_U_total_J'] for r in data],
                color=color, marker=marker, label=label)
    ax.set_xscale('log'); ax.set_yscale('symlog', linthresh=1e-21)
    ax.axhline(0, color='white', linewidth=.6)
    ax.legend(loc='upper left', fontsize=11)
    ax.text(.98, .43, 'Same fixed equilibrium coordinates\nEqual maximum damage; total damage differs',
            transform=ax.transAxes, ha='right', fontsize=11)
    save(fig, 'uniform_vs_local_damage', 'fixtures/pattern_study/*/SUMMARY.json', 'Signed symlog, no absolute-value transformation. Same-sum uniform comparisons are in raw summaries.')

    z = np.load(OUT / 'actual_run_replay/dt_erosion_0p01/jump_001/particle_audit.npz')
    categories = [('INTERIOR','Interior','#d7d7d7','o'), ('FREE_SURFACE','Free surface','#699cbc','s'),
                  ('FIXED_BASE_NEAR','Fixed base near','#9f9a80','^'), ('SOURCE_NEAR','Source near','#e79d57','D')]
    fig, ax = setup('Larger reconstruction changes are associated with energy increases',
                    'Particle signed energy change (10^-20 J)', 'Relative change in F (10^-9)')
    for key, label, color, marker in categories:
        m = z['exclusive_category'] == key
        ax.scatter(z['relative_F_change'][m]/1e-9, z['delta_U_total'][m]/1e-20,
                   s=28, marker=marker, color=color, alpha=.8, label=label)
    ax.axhline(0, color='white', linewidth=.7); ax.legend(fontsize=11)
    ax.text(.02,.04,'First baseline jump · Spearman ρ = −0.766\nAssociation does not establish causation',transform=ax.transAxes,fontsize=11)
    save(fig, 'F_reconstruction_vs_energy', 'actual_run_replay/dt_erosion_0p01/jump_001/particle_audit.npz', 'All 360 particles; exclusive geometric categories.')

    fig, ax = setup('Condition changes have a weaker association with reconstruction energy',
                    'Reconstruction contribution (10^-20 J)', 'Change in K condition number (10^-4)')
    for key, label, color, marker in categories:
        m=z['exclusive_category']==key
        ax.scatter(z['K_condition_change'][m]/1e-4,z['delta_U_reconstruction'][m]/1e-20,
                   s=28,marker=marker,color=color,alpha=.8,label=label)
    ax.axhline(0,color='white',linewidth=.7); ax.legend(fontsize=11)
    ax.text(.02,.035,'First baseline jump · Spearman ρ = −0.386\nAll audited supports remain full rank',transform=ax.transAxes,fontsize=11)
    save(fig, 'K_condition_vs_energy', 'actual_run_replay/dt_erosion_0p01/jump_001/particle_audit.npz', 'No claim of singular K or causal attribution from correlation.')

    f=np.load(OUT/'bond_sensitivity/baseline_first/finite_difference.npz')
    values=f['estimate_J'][f['sign_resolved']]/1e-16
    fig, ax=setup('Some sampled bonds have a resolved negative damage force',
                   'Number of resolved sampled bonds', 'Y = −∂U/∂D at fixed x (10^-16 J per unit D)')
    ax.hist(values,bins=35,color='#699cbc',edgecolor='#000000',linewidth=.5)
    ax.axvline(0,color='#e79d57',linewidth=2,label='Y = 0')
    ax.legend(loc='upper right')
    ax.text(.02,.95,'Negative: 300 / 951 = 31.55% of sample\nResolved: 906; unresolved: 45\nSelection favors damage and energy hotspots\nSample fraction is not a population estimate',
            transform=ax.transAxes,va='top',fontsize=12)
    save(fig,'thermodynamic_force_distribution','bond_sensitivity/baseline_first/finite_difference.npz','Only resolved estimates enter histogram; all sampled bonds are in the denominator.')

    fig=plt.figure(figsize=(11.8,7.6),layout='constrained')
    ax=fig.add_subplot(projection='3d')
    xyz=z['x']*1e3
    val=z['delta_U_reconstruction']/1e-20
    v=float(np.max(np.abs(val)))
    cmap=colors.LinearSegmentedColormap.from_list('signed_audit',['#e9a265','#747474','#6da9cf'])
    pc=ax.scatter(*xyz.T,c=val,cmap=cmap,vmin=-v,vmax=v,s=25,depthshade=False)
    ax.scatter(*xyz[z['fixed']].T,s=43,facecolors='none',edgecolors='#eeeeee',linewidths=.7,label='Fixed particles',depthshade=False)
    cfg=read(ROOT/'configs/streaming_erosion_equilibrium_guarded_v1.json')
    source=np.array(cfg['streaming']['bubble_center_m'])*1e3
    ax.scatter(*source,s=130,marker='*',color='white',label='Streaming source',depthshade=False)
    ax.set(xlabel='X (mm)',ylabel='Y (mm)',zlabel='Z (mm)')
    ranges=np.ptp(np.vstack([xyz,source]),axis=0)
    ax.set_box_aspect(ranges)
    ax.view_init(elev=22,azim=-63)
    for axis in [ax.xaxis,ax.yaxis,ax.zaxis]:
        axis.set_pane_color((0,0,0,1))
        axis._axinfo['grid']['color']=(.25,.25,.25,.65)
    ax.legend(loc='upper left',fontsize=12)
    cbar=fig.colorbar(pc,ax=ax,shrink=.64,pad=.03)
    cbar.set_label('Reconstruction contribution (10^-20 J)\nNegative: energy increase')
    ax.set_title('Reconstruction energy change in the clot\nActual displacement scale 1× · First baseline jump',pad=22)
    save(fig,'spatial_energy_anomaly_map','actual_run_replay/dt_erosion_0p01/jump_001/particle_audit.npz','True saved x; source is upstream left; fixed-particle rings; no arrows or microbubble spheres.')

    rows=sorted([r for r in s['stabilization_sensitivity'] if r['state']=='actual'],key=lambda r:r['factor'])
    fig,ax=setup('Stabilization offsets part of the increase; it is not its source',
                  'Signed energy change (10^-18 J)', 'Stabilization coefficient', size=(10.5,6.5))
    ax.bar(np.arange(4),[r['delta_U_total_J']/1e-18 for r in rows],color='#699cbc',width=.6)
    ax.set_xticks(np.arange(4),['0','0.5× (α = 1.5)','1× (α = 3)','2× (α = 6)'])
    ax.set_ylim(-3.2, .1)
    ax.axhline(0,color='white',linewidth=.7)
    ax.text(.02,.08,'DIAGNOSTIC ONLY\nNOT A CANDIDATE MODEL / NOT USED IN PRODUCTION\nFixed coordinates, fixed old/new damage · No re-equilibration',
            transform=ax.transAxes,fontsize=11)
    save(fig,'stabilization_sensitivity','counterfactuals/stabilization_sensitivity/actual_*/SUMMARY.json','Only audit coefficient varies. All four signed changes remain negative.')

    neg=int(np.flatnonzero(f['negative_resolved'])[np.argmin(f['estimate_J'][f['negative_resolved']])])
    pos=int(np.nanargmax(np.where(f['sign_resolved'],f['estimate_J'],np.nan)))
    unres=int(np.flatnonzero(~f['sign_resolved'])[0])
    fig,axes=plt.subplots(1,3,figsize=(15,5.5),layout='constrained')
    for ax,k,label,col in zip(axes,[neg,pos,unres],['Resolved negative','Resolved positive','Numerically unresolved'],['#e79d57','#699cbc','#bbbbbb']):
        ax.plot(f['epsilons'],f['Y_J'][k]/1e-16,'o-',color=col)
        ax.set_xscale('log'); ax.axhline(0,color='white',linewidth=.6)
        ax.set_title(f'{label}\nBond {f["bond_ids"][k]}')
        ax.set_xlabel('Finite difference epsilon');ax.grid(True)
        ax.spines[['top','right']].set_visible(False)
    axes[0].set_ylabel('Raw Y estimate (10^-16 J per unit D)')
    fig.suptitle('Multiple step sizes expose unresolved derivative estimates',fontsize=16)
    save(fig,'finite_difference_stability','bond_sensitivity/baseline_first/finite_difference.npz','Representative curves, no smoothing. Predeclared stable-triple rule retains unresolved cases.')
    write(OUT/'figures/FIGURE_MANIFEST.json',dict(figures=manifest,font_file=str(font),font_sha256=hashlib.sha256(font.read_bytes()).hexdigest(),
        background='#000000',font='Arial',font_weight='normal',arrow_count=0,displacement_scale=1,production_changes=False))
    print(json.dumps(dict(figures=len(manifest),uniform_interpretation=correction)))


if __name__=='__main__':
    main()
