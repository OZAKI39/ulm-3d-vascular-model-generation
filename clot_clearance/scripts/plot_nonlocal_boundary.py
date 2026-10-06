"""Publication figures from saved NONLOCAL_BOUNDARY_FORMULATION_V1 data only."""
import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT/'verification/nonlocal_boundary_formulation_v1'
OUT = ROOT/'results/nonlocal_boundary_formulation_v1/figures'
COLORS = {'B0': '#bb5538', 'B1': '#087f8c', 'B2': '#7661a5'}


def load(name): return json.loads((DATA/name).read_text())


def rows(name):
    with (DATA/name).open() as f:
        return [{k: float(v) if v else np.nan for k, v in r.items()} for r in csv.DictReader(f)]


def finish(fig, name):
    fig.savefig(OUT/f'{name}.png', dpi=600, facecolor='white')
    fig.savefig(OUT/f'{name}.pdf', facecolor='white')
    plt.close(fig)


def schematic():
    cfg = load('PREDECLARED_CONFIG.json')
    a = np.load(DATA/'NORMAL_AFFINE.npz'); X = a['X']*1e3
    cut = abs(X[:, 1]) < 1e-9; xx = X[cut]; direct = a['direct_region'][cut]
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.3), layout='constrained', gridspec_kw={'width_ratios': [1, 1, 1.12]})
    ax = axes[0]
    for mask, color, label in [(direct, '#ddc799', r'Direct region: $d<\delta$'), (~direct, '#456a84', r'Geometric core: $d\geq\delta$')]:
        ax.scatter(xx[mask, 0], xx[mask, 2], s=10, color=color, label=label)
    bound = cfg['L_m']*1e3/2
    for z in [-bound, bound]: ax.plot([-bound, bound], [z, z], color='#333333', lw=1.6)
    ax.text(0, bound+.042, r'$\Gamma_D$: prescribed displacement', ha='center', fontsize=8)
    ax.text(-bound-.035, 0, r'$\Gamma_N$: traction', ha='right', va='center', rotation=90, fontsize=8)
    ax.set(title=r'(a) Physical domain $\Omega$: central slice', xlabel='Reference x [mm]', ylabel='Reference z [mm]', aspect='equal', xlim=(-bound-.08,bound+.025), ylim=(-bound-.025,bound+.095))
    ax.legend(loc='upper center', bbox_to_anchor=(.5, -.18), frameon=False, fontsize=8)
    ax = axes[1]; force = a['load'][cut][:, [0, 2]]
    mag = np.linalg.norm(force, axis=1); scale = mag.max()/.043
    ax.scatter(xx[:, 0], xx[:, 2], s=7, color='.8')
    ax.quiver(xx[:, 0], xx[:, 2], force[:, 0]/scale, force[:, 1]/scale,
              color=COLORS['B1'], angles='xy', scale_units='xy', scale=1., width=.004)
    ax.set(title='(b) Frozen equivalent volume force', xlabel='Reference x [mm]', ylabel='Reference z [mm]', aspect='equal')
    ax.text(.5, -.2, 'Arrows use one visual scale\nDependency support: all 6,859 particles', ha='center', va='top', transform=ax.transAxes, fontsize=8)
    ax = axes[2]; ax.axis('off'); ax.set_title('(c) Boundary-data conversion')
    labels = [r'Classical data: $u_D$, $t=P N$', r'Analytic local solution $u_L$', r'Unchanged BAQP: $r_h(u_L)$', r'Frozen force $f^{BC}=-r_h(u_L)$', r'Solve $r_h(u)+f^{BC}=0$']
    for k, label in enumerate(labels):
        y = .91-k*.18
        box = FancyBboxPatch((.04, y-.105), .92, .105, boxstyle='round,pad=.012', transform=ax.transAxes, facecolor='#edf5f5', edgecolor='#41818a')
        ax.add_patch(box); ax.text(.5, y-.0525, label, ha='center', va='center', fontsize=10, transform=ax.transAxes)
        if k < 4: ax.annotate('', xy=(.5, y-.16), xytext=(.5, y-.12), xycoords='axes fraction', arrowprops={'arrowstyle': '->', 'color': '.3'})
    ax.text(.5, -.09, 'Independent gates: force, moment and work\nIntact affine mechanics; no phase evolution', ha='center', va='top', fontsize=9, transform=ax.transAxes)
    fig.suptitle(r'Nonlocal boundary formulation | $N=19$, $h=62.5\,\mu$m, $\delta=251.25\,\mu$m', fontsize=13)
    finish(fig, 'boundary_formulation_schematic')


def normal_comparison():
    cases = {c: load(f'normal_{c}.json') for c in COLORS}
    fig, axes = plt.subplots(2, 2, figsize=(9.4, 6.6), layout='constrained')
    for ax, section, metric, title, gate in zip(axes.ravel(),
        ['load_checks', 'load_checks', 'fields', 'fields'],
        ['resultant_error', 'virtual_work_error', 'F_error_core_RMS', 'F_error_boundary_RMS'],
        ['(a) Resultant force', '(b) Affine virtual work', '(c) Geometric-core gradient RMS', '(d) Boundary-region gradient RMS'],
        [1e-6, 1e-6, 1e-5, 1e-5]):
        values = [cases[c][section][metric] for c in COLORS]
        ax.bar(list(COLORS), np.maximum(values, 1e-16), color=list(COLORS.values()), width=.56)
        ax.axhline(gate, color='.25', ls='--', lw=1, label=f'Gate: {gate:.0e}')
        ax.set(yscale='log', ylim=(1e-16, max(2., max(values)*30)), title=title, ylabel='Normalized error')
        for i, v in enumerate(values): ax.text(i, max(v, 1e-16)*1.6, f'{v:.2e}', ha='center', fontsize=8)
        ax.legend(frameon=False, fontsize=8, loc='upper right'); ax.grid(axis='y', alpha=.15)
    fig.suptitle('Normal affine benchmark | total external load includes constraint reactions', fontsize=12)
    fig.supxlabel(r'B0: surface baseline     B1: volume force     B2: locality-strong reference'+'\n'+r'Gradient error normalized by $\|F_* - I\|_F$; core means $d\geq\delta$, not a closed bulk stencil.', fontsize=9)
    finish(fig, 'normal_traction_candidate_comparison')


def support():
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7), layout='constrained')
    for col, kind in enumerate(['normal', 'shear']):
        data = rows(f'{kind.upper()}_BOUNDARY_FORCE_SUPPORT.csv'); ref = load(f'{kind.upper()}_REFERENCE.json')
        d = np.array([r['depth_over_delta'] for r in data]); dh = np.array([r['depth_h'] for r in data])
        ax = axes[0, col]
        for key, label, style in [('load_mean_N', 'Mean per particle', '-o'), ('load_max_N', 'Maximum particle', '--s')]:
            ax.semilogy(d, [max(r[key], 1e-25) for r in data], style, ms=4, label=label)
        ax.axvline(1, color='.5', ls=':', lw=1); ax.set(title=f'({"ab"[col]}) {kind.capitalize()} equivalent load', xlabel=r'Nearest-boundary depth $d/\delta$', ylabel=r'$\|f_i^{BC}\|$ [N]')
        ax.legend(frameon=False, fontsize=8); ax.grid(alpha=.15)
        ax = axes[1, col]; ax.step(dh, [r['cumulative_L1_fraction'] for r in data], where='post', color=COLORS['B1'], lw=2)
        ax.scatter(dh, [r['cumulative_L1_fraction'] for r in data], s=15, color=COLORS['B1'])
        q = ref['B1_support_quantiles']; note = '\n'.join(f'{name}: {v["depth_h"]} h = {v["depth_m"]*1e6:.1f} µm' for name, v in q.items())
        ax.text(.97, .07, note, transform=ax.transAxes, va='bottom', ha='right', fontsize=9, bbox={'facecolor':'white','alpha':.9,'edgecolor':'none'})
        ax.set(title=f'({"cd"[col]}) Cumulative load magnitude', xlabel=r'Nearest-boundary depth $d/h$', ylabel=r'$\sum_{d_i\leq d}\|f_i\|\,/\,\sum_i\|f_i\|$', ylim=(0, 1.04)); ax.grid(alpha=.15)
    fig.suptitle('B1 load support | distribution of the conversion, not a physical boundary-layer thickness', fontsize=12)
    fig.supxlabel('The last depth bin contains one center particle; its near-zero force reflects symmetry, not full stencil closure.',fontsize=8)
    finish(fig, 'equivalent_boundary_force_profile')


def normal_profiles():
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7), layout='constrained')
    for c, color in COLORS.items():
        data = rows(f'normal_{c}_FIELD_PROFILE.csv'); d = [r['depth_over_delta'] for r in data]
        axes[0, 0].plot(d, [r['Fzz_mean'] for r in data], '-o', ms=3, color=color, label=c)
        axes[0, 1].plot(d, [r['Psi_plus_mean_Pa'] for r in data], '-o', ms=3, color=color, label=c)
        a = np.load(DATA/f'normal_{c}_STATE.npz')
        plane = np.rint((.001125/2-abs(a['X'][:, 2]))/62.5e-6).astype(int)
        residual = a['force_internal']+a['external']; scale = load('NORMAL_REFERENCE.json')['B1_conversion']['force_scale_N']
        active_layers, rms = [], []
        for k in np.unique(plane):
            active=(plane==k)&~a['fixed']
            if active.any():
                active_layers.append(k/4.02); rms.append(np.linalg.norm(residual[active])/np.sqrt(active.sum())/scale)
        axes[1, 0].semilogy(active_layers, np.maximum(rms,1e-18), '-o', ms=3, color=color, label=c)
        exact_F=np.diag([1.,1.,1.003])
        error = np.linalg.norm(a['Fbar']-exact_F, axis=(1,2))/.003
        depth = np.rint(np.min(.001125/2-abs(a['X']),axis=1)/62.5e-6).astype(int)
        values = [np.sqrt(np.mean(error[depth == k]**2)) for k in np.unique(depth)]
        axes[1, 1].semilogy(np.unique(depth)/4.02, np.maximum(values, 1e-16), '-o', ms=3, color=color, label=c)
    axes[0, 0].axhline(1.003, color='k', ls='--', lw=1, label='Affine target')
    axes[0, 1].axhline(load('NORMAL_REFERENCE.json')['response']['W'], color='k', ls='--', lw=1, label='Affine target')
    for ax, title, ylabel in zip(axes.ravel(), ['(a) Mean reconstructed axial gradient', '(b) Mean nodal positive energy', '(c) Free-node residual RMS', '(d) Gradient error by boundary distance'],
        [r'$\overline{F}_{zz}$', r'$\Psi^+$ [Pa]', r'Force residual / $(\|P\|_F L^2)$', r'RMS $\|\overline{F}-F_*\|_F\,/\,\|F_*-I\|_F$']):
        ax.set(title=title, xlabel=r'Depth from nearest z-grip plane / $\delta$', ylabel=ylabel); ax.grid(alpha=.15); ax.legend(frameon=False, fontsize=8)
    axes[1,1].set_xlabel(r'Nearest-boundary depth $d/\delta$ (all six faces)')
    fig.suptitle('Normal equilibrium fields | B2 constrains a physical one-horizon collar', fontsize=12)
    fig.supxlabel('Residual panel excludes layers with no free nodes. Gradient errors use the imposed-gradient scale.',fontsize=8)
    finish(fig, 'normal_affine_field_profiles')


def shear_consistency():
    ref = load('SHEAR_REFERENCE.json'); cases = {c: load(f'shear_{c}.json') for c in COLORS}
    fig, axes = plt.subplots(2, 2, figsize=(11.5, 7.2), layout='constrained')
    ax = axes[0, 0]; faces = ref['face_targets']; ix = np.arange(6)
    for offset, key, label, color in [(-.18, 'normal_force', 'Normal magnitude', '#db9f43'), (.18, 'tangential_force', 'Tangential magnitude', '#087f8c')]:
        ax.bar(ix+offset, [np.linalg.norm(r[key])*1e6 for r in faces], .34, label=label, color=color)
    ax.set(xticks=ix, xticklabels=[r['face'] for r in faces], title='(a) Exact face loads: every face included', ylabel='Force magnitude [µN]'); ax.legend(frameon=False, fontsize=8)
    inset=ax.inset_axes([.365,.55,.285,.30])
    inset.bar(ix,[np.linalg.norm(r['normal_force'])*1e9 for r in faces],color='#db9f43',width=.65)
    inset.set(xticks=ix,xticklabels=[r['face'] for r in faces],ylim=(0,13),yticks=[0,6,12])
    inset.tick_params(labelsize=6.5,pad=1);inset.set_title('Normal component [nN]',fontsize=8,pad=3)
    ax = axes[0, 1]; channels = [('normal', ref['traction_channels']['normal']), ('tangent', ref['traction_channels']['tangential']), ('combined', ref['B1_conversion'])]
    for j, (key, label) in enumerate([('resultant_error','Force'), ('moment_error','Moment'), ('virtual_work_error','Work')]):
        ax.bar(np.arange(3)+(j-1)*.24, [max(c[key], 1e-18) for _, c in channels], .23, label=label)
    ax.axhline(1e-6, ls='--', color='.2', lw=1); ax.set(yscale='log', ylim=(1e-18, 1e-4), xticks=range(3), xticklabels=[c for c, _ in channels], title='(b) B1 independent traction channels', ylabel='Normalized error'); ax.legend(frameon=False, fontsize=8)
    ax = axes[1, 0]; variants = ['A_xx', 'A_xz', 'A_zx', 'A_zz']
    for j, (c, color) in enumerate(COLORS.items()):
        work = {w['variation']: w for w in cases[c]['load_checks']['work']}
        ax.bar(np.arange(4)+(j-1)*.24, [max(work[v]['normalized_error'], 1e-16) for v in variants], .23, label=c, color=color)
    ax.axhline(1e-6, color='.2', ls='--', lw=1); ax.set(yscale='log', ylim=(1e-16, 10), xticks=range(4), xticklabels=[r'$A_{xx}$',r'$A_{xz}$',r'$A_{zx}$',r'$A_{zz}$'], title='(c) Nonzero affine virtual-work variations', ylabel='Relative work error'); ax.legend(frameon=False, fontsize=8)
    ax = axes[1, 1]; target = ref['B1_conversion']['target_moment_Nm'][1]
    for j, (c, color) in enumerate(COLORS.items()):
        value = cases[c]['load_checks']['discrete_moment_Nm'][1]
        ax.bar(j, value*1e15, color=color, width=.55)
        ax.text(j, value*1e15, f'{value:.3e}', ha='center', va='bottom', fontsize=8)
    ax.axhline(target*1e15, ls='--', color='k', label='Exact reference moment')
    ax.set(xticks=range(3), xticklabels=list(COLORS), title='(d) Moment about reference centroid', ylabel=r'$M_y$ [$10^{-15}$ N m]'); ax.margins(y=.3); ax.legend(frameon=False, fontsize=8)
    for ax in axes.ravel(): ax.grid(axis='y', alpha=.15); ax.set_axisbelow(True)
    fig.suptitle('Affine shear consistency | reference-coordinate moment is small and nonzero', fontsize=12)
    finish(fig, 'shear_traction_consistency')


def provenance():
    (OUT/'FIGURE_PROVENANCE.json').write_text(json.dumps({
        'generator':'scripts/plot_nonlocal_boundary.py', 'dpi':600,
        'formats':['PNG','PDF'], 'input_directory':str(DATA),
        'data_type':'Computed intact affine numerical benchmarks with analytic continuum targets; not experimental clot data',
        'schematic_arrow_scale':'Single visual gain on saved normal force vectors',
        'log_plot_floors':'1e-16/1e-18 as labeled numerical-display floors; source metrics unchanged'}, indent=2)+'\n')


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family':'DejaVu Sans', 'font.size':9, 'axes.titlesize':10,
                         'axes.spines.top':False, 'axes.spines.right':False, 'pdf.fonttype':42})
    for function in [schematic, normal_comparison, support, normal_profiles, shear_consistency]: function()
    provenance()


if __name__ == '__main__': main()
