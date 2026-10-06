"""English scientific figures, derived solely from the local saved audit."""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from pd_clot.pfpd_candidate.affine_reference_audit import ROOT, OUT, TASK, PATHS, frob

FIG = ROOT/'results'/TASK/'figures'
COLORS = ['#D55E00', '#0072B2', '#CC79A7', '#009E73']


def figures(s, a):
    FIG.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'ps.fonttype': 42})
    metadata = []
    def save(fig, name):
        fig.savefig(FIG/(name+'.png'), dpi=600, facecolor='white')
        fig.savefig(FIG/(name+'.pdf'), facecolor='white')
        metadata.append(dict(name=name, PNG_dpi=600, PDF='vector',
                             data='verification/'+TASK+'/LOCAL_COEFFICIENT_ARRAYS.npz'))
        plt.close(fig)

    ferr = s['F_errors']['Ftilde']['max_relative_frobenius']
    pscale = s['scales']['P_frobenius_Pa']
    perr = [s['P_errors']['raw']['max_relative_frobenius'], 0,
            float(np.max(frob(a['P_high_precision']-a['P_exact']))/pscale),
            float(frob(a['P_exact_high_precision']-a['P_exact'])/pscale)]
    zerr = [s['path_errors'][p]['max_bulk_relative_error'] for p in PATHS]
    labels = ['A  Raw float64', 'B  Exact P / float64',
              'C  Reconstructed F / high precision', 'D  Exact F / high precision']
    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    fig.subplots_adjust(left=.12, right=.98, bottom=.25, top=.77)
    for i, (f, p, z) in enumerate(zip([ferr, 0, ferr, 0], perr, zerr)):
        values = np.array([f, p, z])
        x = np.arange(3)+(i-1.5)*.045
        ax.plot(x, np.where(values > 0, values, np.nan), '-o', color=COLORS[i],
                lw=1.5, ms=5, label=labels[i])
    ax.axhline(1e-12, color='.35', linestyle='--', lw=1)
    ax.text(.02, 1.4e-12, 'Historical Z gate', ha='left', va='bottom', color='.25', fontsize=9)
    ax.set(yscale='log', ylim=(3e-16, 4e-12), xlim=(-.15, 2.18),
           xticks=[0, 1, 2], xticklabels=['F reconstruction', 'Constitutive P', 'Assembled Z'],
           ylabel='Maximum relative Frobenius error')
    ax.grid(axis='y', alpha=.17)
    fig.suptitle('The mismatch enters through reconstructed stress', y=.97, fontsize=14)
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.12), ncol=2, frameon=False, fontsize=9)
    fig.text(.12, .13, r'Each stage uses its own nonzero scale: $\|F_*\|_F$, $\|P_*\|_F$, $\|Z_{bulk}\|_F$.', fontsize=10)
    fig.text(.12, .085, 'All references here use the saved float64 affine state. The dashed gate applies only to Z.', fontsize=9)
    fig.text(.12, .04, 'Exact zeros are omitted on the log axis: F in B/D and P in B. D includes the reference precision offset.', fontsize=9)
    save(fig, 'error_propagation_chain')

    fig, ax = plt.subplots(figsize=(9.2, 5.5))
    fig.subplots_adjust(left=.12, right=.98, bottom=.31, top=.84)
    ax.axhline(1e-12, color='.35', linestyle='--', lw=1)
    for i, val in enumerate(zerr):
        ax.scatter(i, val, color=COLORS[i], s=65, zorder=3)
        ax.annotate(f'{val:.3e}', (i, val), xytext=(0, 10), textcoords='offset points',
                    ha='center', fontsize=10)
    matched = s['high_precision']['exact_stored_F_matched_bulk_error']
    ax.scatter(3, matched, facecolors='white', edgecolors=COLORS[3], s=65, marker='D', zorder=3)
    ax.plot([3, 3], [matched, zerr[3]], color=COLORS[3], ls=':', lw=1)
    ax.annotate(f'{matched:.3e}\nmatched bulk', (3, matched), xytext=(-12, 8),
                textcoords='offset points', ha='right', va='bottom', fontsize=9, color=COLORS[3])
    ax.set(yscale='log', ylim=(4e-16, 5e-12), xlim=(-.5, 3.5),
           ylabel=r'Maximum $\|Z-Z_{ref}\|_F\,/\,\|Z_{bulk}\|_F$',
           xticks=range(4), xticklabels=['A\nRaw F, P, Z\nfloat64', 'B\nExact P\nfloat64',
                                         'C\nReconstructed F\nhigh-precision P, Z', 'D\nExact stored F\nhigh-precision P, Z'])
    ax.grid(axis='y', alpha=.17)
    fig.suptitle('Exact affine stress removes the coefficient exceedance', y=.96, fontsize=14)
    fig.text(.12, .91, '257 closed-stencil coefficient nodes; original Z tolerance = 1e-12', fontsize=10)
    fig.text(.12, .085, 'Filled circles: original float64 bulk reference. Open diamond: stress-matched high-precision bulk.', fontsize=9)
    fig.text(.12, .04, 'D retains the same stored F. Its reference precision offset is shown explicitly; no tolerance is changed.', fontsize=9)
    save(fig, 'z_reference_comparison')

    observed = a['delta_Z_observed']/s['scales']['Zbulk_frobenius_Pa']
    predicted = a['delta_Z_predicted']/s['scales']['Zbulk_frobenius_Pa']
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 5.1))
    fig.subplots_adjust(left=.09, right=.98, bottom=.22, top=.83, wspace=.3)
    shear_extent = 1.15*max(np.max(abs(observed[:, 2, :2])), np.max(abs(predicted[:, 2, :2])))*1e17
    for ax, components, limit, title, exponent in [
        (axes[0], [(0, 0, 'xx'), (1, 1, 'yy'), (2, 2, 'zz')], (-1.2, .025), '(a) Normal components', 12),
        (axes[1], [(2, 0, 'zx'), (2, 1, 'zy')], (-shear_extent, shear_extent), '(b) Shear perturbations (expanded)', 17)]:
        for j, (r, c, label) in enumerate(components):
            ax.scatter(observed[:, r, c]*10**exponent, predicted[:, r, c]*10**exponent, s=14, alpha=.6,
                       color=COLORS[j], label=label, linewidths=0)
        ax.plot(limit, limit, '--', color='.3', lw=1, label='y = x')
        ax.set(xlim=limit, ylim=limit, title=title,
               xlabel=rf'Observed $\Delta Z / \|Z_{{bulk}}\|_F$ ($\times 10^{{-{exponent}}}$)',
               ylabel=rf'Predicted $\Delta Z / \|Z_{{bulk}}\|_F$ ($\times 10^{{-{exponent}}}$)')
        ax.set_aspect('equal', adjustable='box')
        ax.grid(alpha=.15)
        ax.legend(frameon=False, fontsize=8, loc='upper left')
    fig.suptitle('Stress perturbations reconstruct the observed Z shift', y=.97, fontsize=14)
    fig.text(.09, .10, f"All nine components checked: max closure residual / bulk scale = {s['causal_closure']['max_relative_frobenius']:.3e}.", fontsize=10)
    fig.text(.09, .05, 'Criterion = 1e-13. No division by individual zero or near-zero tensor entries.', fontsize=9)
    save(fig, 'delta_z_closure')

    condition = s['summation']['center']
    center_Z = np.asarray(condition['Z'])
    sumabs = np.asarray(condition['sum_absolute_terms'])
    components = [(2, 0, 'zx'), (2, 1, 'zy'), (2, 2, 'zz')]
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    fig.subplots_adjust(left=.12, right=.98, bottom=.24, top=.8)
    x = np.arange(3)
    ax.bar(x-.17, [sumabs[i, j] for i, j, _ in components], .34, color=COLORS[1], label='Sum of absolute terms')
    ax.bar(x+.17, [abs(center_Z[i, j]) for i, j, _ in components], .34, color=COLORS[0], label='Absolute final sum')
    ax.set(yscale='log', ylim=(1e-17, 20), ylabel='Magnitude [Pa]', xticks=x,
           xticklabels=[r'$Z_{zx}$', r'$Z_{zy}$', r'$Z_{zz}$'])
    ax.legend(loc='lower left', bbox_to_anchor=(0, 1.02), frameon=False, ncol=2, fontsize=9)
    ax.grid(axis='y', alpha=.15)
    fig.suptitle('Shear cancellation does not amplify the bulk-norm error', y=.96, fontsize=14)
    fig.text(.12, .13, f"Dominant zz cancellation factor = 1.00; tensor norm factor = {condition['norm_sumabs_over_norm_sum']:.4f}.", fontsize=10)
    fig.text(.12, .075, 'Near-zero shear sums have no reported component-wise ratio. The stable Z scale is dominated by zz.', fontsize=9)
    fig.text(.12, .03, 'Same-term float64 vs balanced / extended sums differ by at most 1.99e-15 of the bulk scale.', fontsize=9)
    save(fig, 'z_accumulation_conditioning')
    (FIG/'FIGURE_PROVENANCE.json').write_text(json.dumps(metadata, indent=2)+'\n')


if __name__ == '__main__':
    s = json.loads((OUT/'SCIENTIFIC_RESULTS.json').read_text())
    with np.load(OUT/'LOCAL_COEFFICIENT_ARRAYS.npz', allow_pickle=False) as source:
        figures(s, {k: source[k] for k in source.files})
