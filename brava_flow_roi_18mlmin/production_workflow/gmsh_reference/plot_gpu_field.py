"""Scientific field/steady evidence only; no CPU/GPU performance comparison."""
from pathlib import Path
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]


def main():
    field = json.loads((ROOT / 'FIELD_VALIDATION.json').read_text())
    assert field.get('GPU_manufactured_status') == 'GPU_MANUFACTURED_STEADY_PASS'
    assert field['local_field_accuracy_assessed']
    execution = json.loads((ROOT / 'manufactured_pipe/R1_GMSH_V2/cfd_gpu/reports/execution.json').read_text())
    rows = [row for row in field['groups'] if row['role'] == 'PARTICLE_RELEVANT']
    distances = [row['wall_distance_m'] * 1e6 for row in rows]
    plt.rcParams.update({'font.size': 12})
    fig, axes = plt.subplots(2, 4, figsize=(20, 9), constrained_layout=True)
    for ax, key, threshold, label in [
        (axes[0, 0], 'E_u', 1e-5, 'Saved-state velocity change'),
        (axes[0, 1], 'E_Q', 1e-6, 'Saved-state flow change')]:
        intervals = execution['intervals']
        ax.semilogy([r['step'] for r in intervals],
                    [r[key] if r[key] > 0 else float('nan') for r in intervals], 'o-', color='#246ba0')
        zero_steps = [r['step'] for r in intervals if r[key] == 0]
        if zero_steps:
            ax.text(.04, .12, f'Exact zero at step(s) {zero_steps}\nNot drawn on logarithmic axis',
                    transform=ax.transAxes, fontsize=9)
        ax.axhline(threshold, ls='--', color='#ba2525', label='Frozen threshold')
        qualified = execution['stop']['step']
        ax.axvspan(intervals[-5]['previous_step'], qualified, alpha=.1, color='green', label='Five qualifying intervals')
        ax.set(title=label, xlabel='Native step', ylabel='Relative change')
        ax.legend(fontsize=9)
    for ax, key, limit, label in [
        (axes[0, 2], 'u', 5., 'Velocity u'),
        (axes[0, 3], 'G', 8., 'Full gradient G (Frobenius)'),
        (axes[1, 0], 's', 8., 'Tangential shear vector s'),
        (axes[1, 1], 'curl', 8., 'Vorticity curl(u)')]:
        values = [row['metrics'][key]['relative_L2'] * 100 for row in rows]
        ax.plot(distances, values, 'o-', color='#246ba0')
        ax.axhline(limit, ls='--', color='#ba2525', label=f'Hard gate {limit:g}%')
        for x, y in zip(distances, values):
            ax.annotate(f'{y:.2f}%', (x, y), xytext=(0, 8), textcoords='offset points', ha='center', fontsize=10)
        ax.set_xscale('log')
        ax.margins(x=.12)
        ax.set_xticks(distances, ['0.377', '1', '2', '5', '10'])
        ax.set_ylim(0, max(limit * 1.35, max(values) * 1.25))
        ax.set(title=label, xlabel='Exact-cylinder wall distance (µm)', ylabel='Relative L2 error (%)')
        ax.legend(fontsize=9)
    bulk = [row for row in field['groups'] if row['role'] == 'BULK']
    values = [row['metrics']['u']['relative_L2'] * 100 for row in bulk]
    ax = axes[1, 2]
    ax.bar(range(len(bulk)), values, color=['#ba2525' if v > 5 else '#246ba0' for v in values])
    ax.axhline(5, ls='--', color='#ba2525', label='Hard gate 5%')
    for x, y in enumerate(values):
        ax.text(x, y + .12, f'{y:.6f}%', ha='center', fontsize=9)
    ax.set_xticks(range(len(bulk)), [f"{row['wall_distance_m'] * 1e6:.2f}" for row in bulk])
    ax.set(title='Bulk velocity u', xlabel='Exact-cylinder wall distance (µm)',
           ylabel='Relative L2 error (%)', ylim=(0, max(values) * 1.3))
    ax.legend(fontsize=9)
    ax = axes[1, 3]
    integral = [max(s['Q_relative_error'] for s in field['internal_sections']) * 100,
                field['dp_dz_relative_error'] * 100]
    ax.bar([0, 1], integral, color='#246ba0', width=.5)
    for x, (value, limit) in enumerate(zip(integral, [5., 8.])):
        ax.hlines(limit, x - .35, x + .35, ls='--', color='#ba2525')
        ax.text(x, value + .2, f'{value:.3f}%', ha='center')
        ax.text(x, limit + .2, f'Gate {limit:g}%', ha='center', color='#ba2525')
    ax.set_xticks([0, 1], ['Max internal Q', 'Pressure gradient'])
    ax.set(title='Integral field gates', ylabel='Relative error (%)', ylim=(0, 10))
    for ax in axes.flat:
        ax.grid(alpha=.25)
    fig.suptitle('Frozen Gmsh V2: GPU steady field validation\n' + field['status'], fontsize=17)
    out = ROOT / 'gmsh_reference/figures'
    fig.savefig(out / 'gpu_v2_field_validation.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
