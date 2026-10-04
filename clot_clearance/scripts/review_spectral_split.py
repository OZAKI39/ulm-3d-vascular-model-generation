"""Render the six requested figure groups and the evidence-only Chinese review.

This script reads saved states. It never launches mechanics, phase or stability
solves, and never writes outside this audit's newly created artifact directories.
"""
from scripts.spectral_split_audit import ROOT, OUT, RESULTS, read, save, np
from pathlib import Path
import hashlib
import os
import subprocess
import time
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.colors import TwoSlopeNorm

COLORS = {'C3': '#64caff', 'C5': '#ffcc6b', 'C7': '#ef94d2'}
FIGURES = RESULTS / 'figures'
LABELS = {
    'SVK_SPECTRAL_SPLIT': 'VERIFIED',
    'OPEN_SOURCE_REFERENCE': 'USED',
    'INTACT_S0_EQUIVALENCE': 'VERIFIED',
    'PURE_TENSION_HARD_GATE': 'REMOVED_FOR_NEW_CANDIDATE_ONLY',
    'C5_INTACT_EQUILIBRIUM': 'VERIFIED',
    'C7_INTACT_EQUILIBRIUM': 'VERIFIED',
    'MIXED_STRAIN_STATES': 'ADMITTED',
    'AT1_DRIVER': 'PSI_PLUS',
    'C3_STABILITY_REGRESSION': 'PASSED',
    'C5_STABILITY': 'UNRESOLVED',
    'C7_STABILITY': 'NOT_RUN',
    'SPATIAL_DISCRETIZATION_AUDIT': 'RESUMED_BUT_INCOMPLETE',
    'PRODUCTION MODEL': 'UNCHANGED',
    'PHYSICAL_CLOT_CALIBRATION': 'NOT PERFORMED',
}


def parse_tests(path):
    tree = ET.parse(path)
    cases = tree.findall('.//testcase')
    failures = {}
    for case in cases:
        failure = case.find('failure')
        if failure is not None:
            failures[case.attrib['classname'] + '::' + case.attrib['name']] = {
                'message': failure.attrib.get('message'), 'traceback': failure.text}
    suites = tree.findall('.//testsuite')
    skipped = sum(int(s.attrib.get('skipped', 0)) for s in suites)
    errors = sum(int(s.attrib.get('errors', 0)) for s in suites)
    return dict(passed=len(cases)-len(failures)-skipped-errors, failed=len(failures),
                skipped=skipped, errors=errors, failures=failures,
                wall_seconds=sum(float(s.attrib.get('time', 0)) for s in suites))


def test_summary():
    baseline = parse_tests(OUT / 'tests/BASELINE.xml')
    final = parse_tests(OUT / 'tests/FINAL.xml')
    previous = parse_tests(ROOT / 'verification/pfpd_at1_spatial_discretization_audit_v2/tests/FINAL.xml')
    result = dict(baseline=baseline, final=final,
        final_known_failure_ID_message_traceback_identical=baseline['failures'] == final['failures'],
        previous_V2_known_failure_ID_message_traceback_identical=baseline['failures'] == previous['failures'],
        targeted={name: parse_tests(OUT / f'tests/{name}.xml') for name in
                  ['CONSTITUTIVE', 'DOMAIN', 'EQUILIBRIUM', 'BRANCH']})
    assert result['final_known_failure_ID_message_traceback_identical']
    assert result['previous_V2_known_failure_ID_message_traceback_identical']
    assert final['failed'] == 3 and final['errors'] == final['skipped'] == 0
    save(OUT / 'tests/TEST_SUMMARY.json', result)
    return result


def preservation():
    before = read(OUT / 'provenance/PRE_RUN_MANIFEST.json')
    hashes, links, changed, missing = {}, {}, [], []
    for name, expected in before['files'].items():
        path = ROOT / name
        if not path.is_file() or path.is_symlink():
            missing.append(name)
            continue
        digest = hashlib.sha256(path.read_bytes()).hexdigest()
        hashes[name] = digest
        if digest != expected:
            changed.append(name)
    for name, expected in before['symlinks'].items():
        path = ROOT / name
        if not path.is_symlink():
            missing.append(name)
            continue
        links[name] = os.readlink(path)
        if links[name] != expected:
            changed.append(name)
    post = dict(root=str(ROOT), files=hashes, symlinks=links, timestamp=time.time(),
                scope='Exactly the pre-run file/symlink set, including untracked historical artifacts')
    save(OUT / 'provenance/POST_RUN_MANIFEST.json', post)
    result = dict(checked_files=len(hashes), checked_symlinks=len(links), changed=changed,
                  missing=missing, passed=not changed and not missing,
                  production_files_modified=0 if not changed and not missing else None,
                  historical_files_modified=0 if not changed and not missing else None)
    save(OUT / 'provenance/PRESERVATION_CHECK.json', result)
    assert result['passed'], result
    return result


def summarize():
    intact = read(OUT / 'equilibrium/ALL.json')
    c3 = read(OUT / 'critical_brackets/C3_A.json')
    c5 = read(OUT / 'critical_brackets/C5_A.json')
    assert all(r['accepted'] for r in intact)
    assert read(OUT / 'constitutive/VERIFICATION.json')['gate1_passed']
    assert read(OUT / 'constitutive/VERIFICATION.json')['gate2_passed']
    assert read(OUT / 'comparisons/C3_STABILITY_REGRESSION.json')['passed']
    assert not c5['complete'] and c5['interval'] is None
    assert read(OUT / 'GATES.json')['C7_heavy'] == 'NOT_RUN'
    accepted = {'C3': [r for r in c3['states'] if r['stationary']],
                'C5': [r for r in c5['states'] if r['stationary']], 'C7': []}
    save(OUT / 'comparisons/ACCEPTED_BRANCH_STATES.json', accepted)
    save(OUT / 'critical_brackets/C5_ATTEMPT_03_AT1_PREDICTOR.json', c5)
    resources, levels = {}, []
    for level in COLORS:
        eq = [r for r in intact if r['level'] == level]
        branch = [read(p) for p in (OUT / 'branch_tracking').glob(f'{level}_*/STATIONARY.json')]
        init = [read(p) for p in (OUT / 'branch_tracking').glob(f'{level}_*/FROZEN_PHASE_INITIALIZER.json')]
        resources[level] = dict(intact_equilibrium_wall_seconds=sum(r['wall_seconds'] for r in eq),
            stationary_solve_wall_seconds=sum(r['wall_seconds'] for r in branch),
            frozen_phase_initialization_wall_seconds=sum(r['wall_seconds'] for r in init),
            intact_peak_RSS_bytes=max(r['peak_RSS_bytes'] for r in eq),
            accepted_intact_states=len(eq), accepted_coupled_states=len(accepted[level]),
            stationary_solves_including_rejected=len(branch))
        row = next(r for r in eq if r['strain'] == .01)
        levels.append(dict(level=level, strain=.01, intact_equilibrium='VERIFIED',
            min_principal_green_strain=row['min_principal_green_strain'],
            negative_principal_particle_fraction=row['negative_principal_particle_fraction'],
            negative_principal_bond_fraction=row['negative_principal_bond_fraction'],
            Psi_plus_fraction=row['Psi_plus_fraction'], Psi_minus_fraction=row['Psi_minus_fraction'],
            nominal_stress_Pa=row['nominal_stress_Pa'],
            critical_interval=c3['interval'] if level == 'C3' else None))
    save(OUT / 'comparisons/LEVEL_SUMMARY.json', levels)
    save(RESULTS / 'LEVEL_SUMMARY.json', levels)
    save(OUT / 'provenance/COMPUTE_SUMMARY.json', dict(levels=resources,
        timing_scope='Measured calls only. Stationary sums exclude subsequent C3 cone/Hessian audits; no total CPU-time or GPU speedup inferred.',
        GPU_AVAILABLE_BUT_NOT_BENEFICIAL='NOT_ESTABLISHED_NO_BENCHMARK',
        GPU_backend_used=False, GPU_speedup=None,
        backend_reason='No installed CuPy/JAX/Torch batch eigensolver. CUDA driver and numba.cuda available; no custom GPU eigensolver or large stack installed.'))
    # The predeclared file retains copied V2 descriptions. Record the actual
    # candidate overrides explicitly rather than silently rewriting that archive.
    save(OUT / 'provenance/EFFECTIVE_CONFIG_NOTES.json', dict(
        predeclared_file='PREDECLARED_CONFIG.json',
        inherited_descriptive_fields=['domain_policy', 'scope', 'initial_window', 'jacobian', 'created_utc'],
        effective_domain='Negative principal Green strain diagnostic only for the new candidate; finite positive determinant, phase bounds, fixed boundary and old force/KKT gates retained.',
        effective_scope='Four intact strains only; declared spectral branch window plus recorded downward onset extension. C7 heavy work stopped after C5 unresolved.',
        effective_jacobian='Existing dense branch Jacobian for C3; unmodified V2 matrix_free_residual_jacobian for C5, dogbox/lsmr as in V2.',
        constitutive_probe_s0_equivalence=True, numerical_tolerances_changed=False,
        clarification='This file clarifies inherited text after the run; it does not change the numerical configuration. The copied created_utc is the V2 declaration time, not this run start; this run pre-manifest timestamp is authoritative.'))
    for directory in [OUT, RESULTS]:
        save(directory / 'FINAL_STATUS.json', LABELS)
    return intact, accepted, resources, levels


def style():
    font_manager.fontManager.addfont('/mnt/c/Windows/Fonts/arial.ttf')
    plt.rcParams.update({'font.family': 'Arial', 'font.weight': 'normal',
        'axes.titleweight': 'normal', 'figure.facecolor': '#080b10',
        'axes.facecolor': '#080b10', 'savefig.facecolor': '#080b10',
        'text.color': '#eeeeee', 'axes.labelcolor': '#eeeeee', 'xtick.color': '#dddddd',
        'ytick.color': '#dddddd', 'axes.edgecolor': '#7c8796', 'grid.color': '#39414c',
        'font.size': 11, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42})


def finish(fig, name):
    FIGURES.mkdir(parents=True, exist_ok=True)
    for ext in ['png', 'pdf']:
        fig.savefig(FIGURES / f'{name}.{ext}', dpi=190, bbox_inches='tight')
    plt.close(fig)


def figures(intact, accepted):
    style()
    rows = read(OUT / 'constitutive/VERIFICATION.json')['tensors']
    names = ['pure_tension', 'pure_compression', 'mixed', 'shear']
    selected = [next(r for r in rows if r['case'] == name and abs(r['lambda_Pa']) < 1e-10) for name in names]
    fig, ax = plt.subplots(figsize=(8.4, 4.5), layout='constrained')
    plus = np.array([r['Psi_plus'] for r in selected])
    ax.bar(range(4), plus, width=.58, color=COLORS['C3'], label='Psi plus')
    ax.bar(range(4), [r['Psi_minus'] for r in selected], bottom=plus,
           width=.58, color=COLORS['C7'], label='Psi minus')
    ax.set(xticks=range(4), xticklabels=['Pure tension', 'Pure compression', 'Mixed', 'Shear type'],
           ylabel='Energy density (Pa)', title='Spectral split of prescribed Green strain tensors',
           ylim=(0, 1.2*max(r['Psi_plus']+r['Psi_minus'] for r in selected)))
    ax.legend(frameon=False, loc='upper right'); ax.grid(axis='y', alpha=.3)
    finish(fig, 'spectral_split_energy')

    fig, ax = plt.subplots(figsize=(8.4, 4.7), layout='constrained')
    for level, color in COLORS.items():
        rows = [r for r in intact if r['level'] == level]
        ax.plot([r['strain'] for r in rows], [r['min_principal_green_strain'] for r in rows],
                'o-', color=color, label=level)
    ax.axhline(-1e-9, color='#eeeeee', ls='--', lw=1, label='Old -1e-9 gate: diagnostic only')
    ax.set(xlabel='Engineering strain', ylabel='Minimum principal Green strain (BAQP)',
           title='Accepted intact mechanical equilibria, s = 0 fixed')
    ax.grid(alpha=.3); ax.legend(frameon=False, loc='lower left')
    finish(fig, 'mixed_strain_equilibrium')

    fig, axes = plt.subplots(1, 2, figsize=(10.4, 5.6), layout='constrained')
    data = {n: np.load(OUT / f'equilibrium/C{n}/0.00900000/STATE.npz') for n in [5, 7]}
    low = min(float(d['particle_min_green'].min()) for d in data.values()) * 1e4
    high = max(float(d['particle_min_green'].max()) for d in data.values()) * 1e4
    for ax, (n, d) in zip(axes, data.items()):
        X = d['particle_X'] * 1e6
        mask = abs(X[:, 1]) < 1e-7
        points = X[mask]; values = d['particle_min_green'][mask] * 1e4
        order = np.lexsort((points[:, 0], points[:, 2]))
        coordinates = np.linspace(-125, 125, n)
        edges = np.r_[-125, (coordinates[1:] + coordinates[:-1])/2, 125]
        plot = ax.pcolormesh(edges, edges, values[order].reshape(n, n), cmap='coolwarm',
                             norm=TwoSlopeNorm(vmin=low, vcenter=0, vmax=high), rasterized=True)
        ax.scatter(points[:, 0], points[:, 2], color='white', s=9, alpha=.7)
        ax.hlines([-125, 125], -125, 125, color='#f3d779', lw=3, clip_on=False)
        ax.vlines([-125, 125], -125, 125, color='white', lw=1, ls='--', clip_on=False)
        ax.text(.5, 1.02, 'Top clamp', transform=ax.transAxes, ha='center', fontsize=10)
        ax.text(.5, .015, 'Bottom clamp', transform=ax.transAxes, ha='center', fontsize=10,
                bbox=dict(facecolor='#080b10', edgecolor='none', alpha=.85, pad=3))
        ax.text(1.01, .5, 'Free side', transform=ax.transAxes, va='center', rotation=90, fontsize=10)
        ax.set(xlim=(-125, 125), ylim=(-125, 125), xlabel='x (micrometre)', ylabel='z (micrometre)', aspect='equal')
        ax.set_title(f'C{n}: particle Fbar, y = 0', pad=29)
    fig.colorbar(plot, ax=axes, shrink=.66, pad=.06, label='Minimum Green strain (x 1e-4)')
    fig.suptitle('Compression location at engineering strain 0.009\nParticle field; BAQP extrema are reported separately', fontsize=13)
    finish(fig, 'compression_location')

    fig, ax = plt.subplots(figsize=(8.4, 4.7), layout='constrained')
    for level, rows in accepted.items():
        if rows:
            ax.plot([r['strain'] for r in rows], [r['metrics']['Psi_plus_fraction'] for r in rows],
                    'o-', color=COLORS[level], label=f'{level}: accepted branch states')
    ax.text(.04, .55, 'C7 coupled branch: NOT RUN\nC5 stops before accepted phase activation', transform=ax.transAxes, color=COLORS['C7'])
    ax.set(xlabel='Engineering strain', ylabel='Psi plus / (Psi plus + Psi minus)', ylim=(.9965, 1.0003),
           title='Undegraded active energy fraction on accepted coupled states')
    ax.grid(alpha=.3); ax.legend(frameon=False, loc='lower right')
    finish(fig, 'active_energy_fraction')

    fig, ax = plt.subplots(figsize=(8.4, 4.7), layout='constrained')
    bracket = read(OUT / 'critical_brackets/C3_A.json')['interval']
    mid = sum(bracket)/2
    ax.errorbar(125, mid, yerr=np.array([[mid-bracket[0]], [bracket[1]-mid]]), fmt='o', color=COLORS['C3'], capsize=7)
    ax.text(121, .010012, 'C3 verified\n[0.0099921875, 0.0100000000]', ha='right', color=COLORS['C3'])
    for x, level, note in [(62.5, 'C5', 'UNRESOLVED'), (250/6, 'C7', 'NOT RUN')]:
        ax.text(x, .5, level + '\n' + note, transform=ax.get_xaxis_transform(), ha='center', color=COLORS[level], fontsize=10)
    ax.set(xlabel='Particle spacing h (micrometre)', ylabel='Critical engineering strain',
           xlim=(30, 140), ylim=(.00997, .010035), title='Critical interval after split: spatial comparison incomplete')
    ax.grid(alpha=.3)
    finish(fig, 'critical_interval_after_split')

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), layout='constrained', sharey=True)
    for level, color in COLORS.items():
        rows = [r for r in intact if r['level'] == level]
        axes[0].plot([r['strain'] for r in rows], [r['nominal_stress_Pa'] for r in rows], 'o-', color=color, label=level)
        rows = accepted[level]
        if rows:
            axes[1].plot([r['strain'] for r in rows], [r['nominal_stress_Pa'] for r in rows], 'o-', color=color, label=level)
    axes[0].set(title='Intact equilibrium: s = 0 fixed', ylabel='Nominal reaction stress (Pa)')
    axes[1].set(title='Accepted coupled stationary states')
    axes[1].text(.04, .75, 'C7 branch: NOT RUN\nRejected C5 states excluded', transform=axes[1].transAxes, fontsize=10)
    for ax in axes:
        ax.set_xlabel('Engineering strain'); ax.grid(alpha=.3); ax.legend(frameon=False, loc='lower right')
    finish(fig, 'reaction_stress_after_equilibrium')
    save(OUT / 'provenance/FIGURES.json', dict(directory=str(FIGURES), groups=6,
        files=sorted(p.name for p in FIGURES.iterdir()),
        source_policy='Only accepted fixed-phase mechanical equilibria and accepted coupled stationary states. Missing fine critical intervals never interpolated.',
        style='Dark background, Arial regular, English annotations, no arrows'))


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |', '| ' + ' | '.join(['---']*len(headers)) + ' |'] +
                     ['| ' + ' | '.join(str(v) for v in row) + ' |' for row in rows])


def report(intact, accepted, resources, levels, tests, protected):
    maxima = read(OUT / 'constitutive/VERIFICATION.json')['maxima']
    c3 = read(OUT / 'critical_brackets/C3_A.json')
    c5 = read(OUT / 'critical_brackets/C5_A.json')
    last = accepted['C5'][-1]
    attempts = [read(OUT / f'branch_tracking/{label}/STATIONARY.json') for label in ['C5_A_001', 'C5_W_005', 'C5_P_005']]
    baseline, final = tests['baseline'], tests['final']
    diag = read(OUT / 'comparisons/C5_ACTIVATION_LINEARIZATION_DIAGNOSTIC.json')
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip()
    top = subprocess.check_output(['git', 'rev-parse', '--show-toplevel'], cwd=ROOT, text=True).strip()
    eq_table = table(['层级', '应变', '力残差', '应力 Pa', '弹性能 J', '最小 BAQP Green', '负粒子/负键 %', 'Ψ+ 比例'], [
        [r['level'], f"{r['strain']:.5f}", f"{r['force_residual']:.3e}", f"{r['nominal_stress_Pa']:.7g}",
         f"{r['elastic_energy_J']:.6e}", f"{r['min_principal_green_strain']:.6e}",
         f"{100*r['negative_principal_particle_fraction']:.3f} / {100*r['negative_principal_bond_fraction']:.3f}",
         f"{r['Psi_plus_fraction']:.9f}"] for r in intact])
    motion_table = table(['层级', '应变', '最大位移 µm', '最大偏离 affine µm', '最小 J'], [
        [r['level'], f"{r['strain']:.5f}", f"{r['maximum_displacement_m']*1e6:.6g}",
         f"{r['maximum_deviation_from_affine_m']*1e6:.6g}", f"{r['minimum_J']:.9f}"] for r in intact])
    locations, groups, percentiles = [], [], []
    for n in [5, 7]:
        r = next(r for r in intact if r['level'] == f'C{n}' and r['strain'] == .009)
        for name, key, value in [('particle Fbar', 'most_compressive_particle', r['minimum_particle_green_strain']),
                                 ('BAQP Ftilde', 'most_compressive_BAQP', r['min_principal_green_strain'])]:
            loc = r[key]
            locations.append([f'C{n}', name, f'{value:.7e}', ', '.join(f'{x*1e6:.3f}' for x in loc['position_m']),
                              f"{loc['clamp_distance_m']*1e6:.3f}", f"{loc['side_distance_m']*1e6:.3f}"])
        for name, group in r['spatial_groups'].items():
            groups.append([f'C{n}', name, group['count'], f"{group['min']:.6e}", f"{100*group['negative_fraction']:.3f}"])
        for key in ['particle_percentiles', 'bond_percentiles']:
            a = r[key]
            percentiles.append([f'C{n}', key.replace('_percentiles', '')] + [f'{a[i]:.6e}' for i in [0, 2, 4, 6, 8]])
    attempt_table = table(['初值', '应变', 'force', 'KKT', '优化器停止', '物理验收'], [
        [name, f"{r['strain']:.11f}", f"{r['metrics']['force_residual']:.6e}",
         f"{r['metrics']['projected_phase_KKT']:.6e}", r['reason'], 'REJECTED']
        for name, r in zip(['原位移 predictor', '固定相场 ADR', 'ADR + 一次现有 AT1 predictor'], attempts)])
    res_table = table(['层级', 'intact 秒', 'stationary 秒', '初值 ADR 秒', 'intact 进程峰值 RSS MiB'], [
        [level, f"{r['intact_equilibrium_wall_seconds']:.3f}", f"{r['stationary_solve_wall_seconds']:.3f}",
         f"{r['frozen_phase_initialization_wall_seconds']:.3f}", f"{r['intact_peak_RSS_bytes']/2**20:.1f}"] for level, r in resources.items()])
    labels = '\n'.join(f'{k}: {v}' for k, v in LABELS.items())
    failures = '\n'.join('- `' + name + '`' for name in baseline['failures'])
    text = f'''# PFPD_AT1_SPECTRAL_SPLIT_V1 审核报告

结论：Q1 为是，C3/C5/C7 共 12 个固定相场完整态全部达到原力平衡要求。Q2 部分通过：C3 临界区间完全复现；C5 在相场激活前已有合法 stationary states，但首次激活附近的 coupled solve 未通过 force/KKT；C7 按 gate 停止重型分支和稳定性计算。空间离散稳定性审计仍未完成。

Repository：`{top}`；项目目录：`{ROOT}`；HEAD：`{head}`。本轮仅增加独立候选模块、审计脚本、测试和新输出目录，没有提交或替换 production。

## 1. 为什么 old pure-tension gate 不再合适

全局拉伸边界不保证每个 BAQP 局部 Green 张量三个主值均非负。旧候选把最小主值低于 −1e-9 的 trial 整体拒绝，阻止了细离散中的混合拉压完整态平衡。新候选按张量谱分解保留压缩弹性能，只退化拉伸能量，所以负主值可以作为合法诊断。旧 candidate 和旧测试的硬门槛保持原样，−1e-9 仍用于本报告负比例统计，不是新物理拒绝条件。

## 2. 使用了哪篇论文/开源实现

- [PhaseFieldX theory](https://phasefieldx.readthedocs.io/en/stable/theory/phase_field_fracture/main.html) 与 [energy-splitting example](https://phasefieldx.readthedocs.io/en/stable/auto_examples/Elasticity/plot_theory_anisotropic.html) 提供 Miehe spectral / Amor vol-dev 公式对照。实际源码固定 [a9714c122497f8829860efc68d776e4324d48eca](https://github.com/CastillonMiguel/phasefieldx/tree/a9714c122497f8829860efc68d776e4324d48eca)，MIT，`ORACLE_USED`。
- [Miehe et al. 2010](https://doi.org/10.1016/j.cma.2010.04.011)：仅采用拉压分裂概念；本轮未取到论文全文，公式以用户指定内容与 PhaseFieldX 可读源码核对。未引入 AT2、history field、残余刚度或谱平滑。
- [Rahimi & Moutsanidis 2022](https://doi.org/10.1016/j.cma.2022.115191)：全文未获得。只检查 [Stony Brook 文献元数据](https://researchconnect.stonybrook.edu/en/publications/a-smoothed-particle-hydrodynamics-approach-for-phase-field-modeli/)，不从摘要补写 Green/SVK 应力公式。
- [Nardi invariant split](https://github.com/deborahnardi/phase-field-invariant-split/tree/dadd9deae360d54807180d6ad6346e8bb9939636)：固定 `dadd9deae360d54807180d6ad6346e8bb9939636`，MIT。阅读 README、LICENSE 和 `src/getSpectralConstitutiveTensor.cpp`；例程依赖 Element/Solid/material/PETSc 及自定义张量工具，不是可直接编译的独立函数。本轮仅源码参考，`ORACLE_NOT_NEEDED`，未移植或声称数值对照通过。

PhaseFieldX oracle 只执行固定源码中的五个张量代数函数，以 NumPy 替代 UFL 标量/张量运算，以输入 E 替代 epsilon(u)；谱绝对值独立使用 scipy.linalg.sqrtm(E@E)。这不是完整 UFL 特征分解验证，更不是 finite-strain FEM solver equivalence。原源码、许可、SHA 和逐例结果见 `references/`；没有 vendor FEM solver。Amor 只保留 tensor 对照表，没有第二套候选 campaign。

## 3. spectral split 是什么

E = (FᵀF − I)/2 = Q diag(e) Qᵀ；E± = Q diag(max/min(e,0)) Qᵀ。

注意体积项用的是 max/min(tr(E),0)，不是 tr(E±)。

Ψ± = λ〈tr E〉±²/2 + μ tr(E±²)；S± = λ〈tr E〉±I + 2μE±。

Ψ(F,s) = (1−s)²Ψ+ + Ψ−；S(F,s) = (1−s)²S+ + S−；P = FS。

AT1 的 Y = Ψ+，q = −2(1−s)Ψ+ + 2Yc。仍采用原 AT1 irreversibility、sc=0.95、h=1、c0=0.25、delta=0.25125 mm、Gc_demo=1.26884390625e-5 J/m²、Yc=0.1010025 Pa。原材料 μ=1000 Pa、bulk=666.6666666666666 Pa、λ=0；用于构成验证的额外 bulk=2300 Pa 仅是非零 λ 测试输入，不用于 cube 计算。

新增 `pd_clot/pfpd_candidate/svk_spectral_split.py` 提供全部要求 API。审计通过串行、可恢复的 constitutive callback 复用原 BAQP 装配，外层 g=1，避免重复退化；原 state.W0 在该候选中只是 Ψ+ 的兼容别名。能量和力来自已正确分裂、退化的 Ψ/P。旧 BAQP/ADR/branch/KKT/HVP/cone/localization 源码未修改。该 scoped callback 不是线程安全生产接口，本轮始终串行调用。

## 4. 为什么 s=0 时 intact mechanics 不变

E=E++E−，谱正负部分正交，且 〈tr E〉+²+〈tr E〉−²=(tr E)²，因此 Ψ++Ψ−=Ψ0、S++S−=S0。g(0)=1，于是 P=F(S++S−)=P0。

128 个 deterministic、det(F)>0 的随机 F 同时比较 energy/S/P，最大归一化误差 {maxima['intact_equivalence']:.6e}，通过 2e-12 门槛。此后 s=0 ADR 复用原 SVK 力回调；最终状态另用新 split 重新评价并诊断，不改变完整态问题。

## 5. constitutive verification

8 类 Green 张量包括 pure tension/compression、mixed、shear、hydrostatic tension/compression、repeated eigenvalues、near-zero；分别用原材料和非零 λ 材料验证。

{table(['检查', '最大误差', '门槛'], [
    ['Ψ 分解', f"{maxima['Psi_decomposition']:.6e}", '2e-12'],
    ['S 分解', f"{maxima['stress_decomposition']:.6e}", '2e-12'],
    ['能量方向导数', f"{maxima['energy_derivative']:.6e}", '5e-6'],
    ['s=0 invariance', f"{maxima['intact_equivalence']:.6e}", '2e-12'],
    ['PhaseFieldX tensor oracle', f"{maxima['oracle']:.6e}", '2e-8']])}

离开 switching surface 的导数对照使用 1e-5、3e-6、1e-6 三档 epsilon，并覆盖 s=0 和 s=0.37。误差归一化分母及逐次值保存在 `constitutive/VERIFICATION.json`；导数误差以方向导数/应力范数尺度归一，oracle 对接近零分量用 max(norm,1)，不夸称全部是相对误差。

switching state 只记录单侧诊断，不把切线不光滑当成实现错误。重复主值的特征基旋转不改变重构张量。微型 BAQP pure compression fixture 的 Ψ+ 和相场增量为零；simple shear 的 Ψ+、Ψ− 均正。未另跑大型压缩、刚体全矩阵或旧校准 campaign。

![Spectral split energy](../../results/pfpd_at1_spectral_split_v1/figures/spectral_split_energy.png)

## 6. C3/C5/C7 intact mixed-strain equilibrium

L=0.25 mm boundary-node cube，V=1.5625e-11 m³，A=6.25e-8 m²，delta 固定。C3/C5/C7 分别为 27/125/343 粒子及 185/5392/46059 键。本轮只跑指定四个 intact 应变，s=0 固定，无 phase updates。全部 12 状态通过原 ADR criteria、force≤1e-9、固定边界严格相等、finite/positive J 与现有 safety 检查。

{eq_table}

{motion_table}

表中 force 是原归一化指标；详细反力、绝对力残差与能量见 `equilibrium/C*/.../SUMMARY.json`，每点完整状态和 ADR 历史同目录保存。

![Mixed strain equilibrium](../../results/pfpd_at1_spectral_split_v1/figures/mixed_strain_equilibrium.png)

## 7. C5/C7 原 blocker 是否解除

原 `PURE_TENSION_DOMAIN_VIOLATION` 已在新候选解除。C5/C7 完整态可以在负局部主应变下达到平衡，`MIXED_STRAIN_EQUILIBRIUM_CONFIRMED`。这不等于 coupled phase stationarity 或 mesh objectivity 已通过。旧候选仍拒绝同类 mixed trial；新增回归测试明确验证新旧行为差别。

## 8. 负 principal strain 在哪里、多大

以下为应变 0.009 的 accepted intact equilibrium；坐标和距离均为 reference configuration。Particle 使用 Fbar，BAQP 使用 bond-associated Ftilde，两者不能混为同一个应变场。键中点仅用于 BAQP 空间定位，不能把一个键的最小主值直接称为该中点处唯一连续应变值。

{table(['层级', '场', '最小主值', '位置 x,y,z µm', '到 clamp µm', '到侧面 µm'], locations)}

{table(['层级', '粒子组', '数量', '最小主值', '负比例 %'], groups)}

{table(['层级', '场', 'min', 'p5', 'p50', 'p95', 'max'], percentiles)}

clamp 和侧面组在棱边可重叠，比例分母是各组自己的粒子数。C5 的最负 particle 位于 clamp/侧面交界，但最负 BAQP 位于内部键中点；C7 同样需要区分两种场。负应变既不限于表面，也不是 damage 指标。各完整态均保存九档 percentiles 和 raw<0 比例；正文采用 −1e-9 排除浮点零附近符号噪声。

![Compression location](../../results/pfpd_at1_spectral_split_v1/figures/compression_location.png)

## 9. Psi+ / Psi- 如何分布

能量比例按原 BAQP 体积权重求和，即 ΣwΨ+ / Σw(Ψ++Ψ−)。这是未退化弹性成分比例，不以含 crack energy 的总势能作分母。每个 accepted coupled state 和 intact equilibrium 均保存 total_Psi_plus_J、total_Psi_minus_J 和两种比例。

应变 0.01：C5 负粒子比例 80%、负键比例约 74.50%，但 Ψ− 仅约 0.2422%；C7 两个负比例约 53.64% 和 63.89%，但 Ψ− 仅约 0.4959%。因此负主应变出现频率不能替代压缩能量权重。C3 的压缩分量仅浮点量级。

下图只用 accepted coupled states；C5 仅画到相场激活前，C7 不以完整态冒充 coupled branch。完整态的四个载荷能量比例见第 6 节。

![Active energy fraction](../../results/pfpd_at1_spectral_split_v1/figures/active_energy_fraction.png)

## 10. C3 old stability 是否复现

C3 完整复现 [0.0099921875, 0.01]，宽度 7.8125e-6。两端 energy 最大相对差 3.330669e-16，reaction 最大相对差 2.220446e-16；pure-tension driver 对照也通过 2e-12。Ψ− 在该分支为浮点零量级。

最后稳定点为 185 条 strong lower-active 键；0.01 端为 185 条 weak lower 键，力和 KKT 均合格。复用原 critical cone 得到可行负曲率 witness，三档 epsilon 曲率为 {c3['first_unstable']['negative_witness']['curvatures']}，物理可行探针也通过。没有用正的随机方向样本证明稳定。

复用原 V2 HVP 的 C3 dense 对照最大相对误差 1.9165873e-5，小于原 5e-5 门槛；本轮没有预生成 C5/C7 dense Hessian。

## 11. C5/C7 branch/stability 是否建立

C5 在 intact strain=0.005 时 max(Ψ+)/Yc 已约 1.237，因此相场激活明显低于原 0.009 窗口，按用户允许向低应变扩展。低载荷比例估计仅用于选采样点，不能替代临界区间。

复用原 branch solver：C3 使用原 dense Jacobian，C5 沿用 V2 matrix-free residual Jacobian 与 dogbox/lsmr。force=1e-9、KKT=1e-8、ftol/xtol/gtol=1e-13 以及曲率门槛均未调节。保留三次初值尝试的失败记录：

{attempt_table}

第二次以原 ADR 在前一载荷 phase 固定时生成力平衡初值；第三次额外调用一次现有 AT1 constrained update，仅作为合法初猜，前一 accepted phase 的不可逆下界保持不变。没有另写分支求解器、KKT、stagger iteration、phase history，未提交任何被拒绝状态的 damage。

C5 最后一个 accepted coupled stationary state 位于 ε={last['strain']:.12f}，force={last['metrics']['force_residual']:.6e}，KKT={last['metrics']['projected_phase_KKT']:.6e}，max(s)={last['metrics']['maximum_s']:.6e}；全为 strong lower active。它尚未接受曲率分类，不能写成 STABLE。

首次激活附近 ε={attempts[-1]['strain']:.12f} 的最后一次求解仍为 force={attempts[-1]['metrics']['force_residual']:.6e}、KKT={attempts[-1]['metrics']['projected_phase_KKT']:.6e}。优化器 ftol/xtol 终止不等于 stationary。该状态只归档诊断，不进入 accepted branch 图、临界锥或模态结论。

拒绝态的四个 residual-Jacobian 方向/伴随检查记录于 `comparisons/C5_ACTIVATION_LINEARIZATION_DIAGNOSTIC.json`（forward 最大相对误差约 2.59e-7，伴随归一化误差约 5.12e-10），没有证据把障碍直接称为 Jacobian bug；这些局部检查也不证明 solver 一定收敛。

C5 没有 accepted weak/interior threshold state，因此其 heavy stability 未运行，临界区间为 null，科学状态 UNRESOLVED。C7 只保留四个已通过的完整态诊断，coupled branch/HVP/cone/mode 全部 NOT_RUN。脚本内细级 operator/block adapter 尚未在 C5/C7 可接受激活态执行，不能标为已验证。

## 12. spatial consistency 是否终于可判断

目前不能判断稳定性区间、onset reaction 或临界模态的空间一致性。C3 有有效区间，C5/C7 没有；不插值、不填零、不把 onset 采样估计当临界应变。

![Critical interval](../../results/pfpd_at1_spectral_split_v1/figures/critical_interval_after_split.png)

可判断的是 intact equilibrium reaction 对离散仍明显敏感：ε=0.01 时 C3/C5/C7 为 20.301、12.6634、9.53305 Pa。它们均是真实固定相场平衡反力，不是 affine 非平衡反力；这也不是 fracture onset 一致性结论。三个级别不足以据此宣称收敛阶或 mesh objectivity。

![Reaction stress](../../results/pfpd_at1_spectral_split_v1/figures/reaction_stress_after_equilibrium.png)

## 13. CPU/GPU 使用情况和实际加速

CPU：Intel i7-13700HX，环境可见 12 physical / 24 logical cores；保留至少两核，NUMBA_NUM_THREADS=10，BLAS/OMP/MKL=1。CPU NumPy batched eigh 与现有 NumPy BAQP 是 canonical backend。当前 BAQP/spectral 路径没有 Numba parallel kernel，不把线程环境配置误称成已获得并行加速。Python 3.13.11，NumPy 2.5.3，SciPy 1.18.1，Numba 0.67.0。

GPU：RTX 4060 Laptop，8188 MiB；CUDA UMD 13.4，numba.cuda available，nvcc 与 CuPy/JAX/Torch 不在环境中。没有为本轮安装大 CUDA stack 或新写 GPU eigensolver，GPU backend 未使用，CPU/GPU benchmark 未做，加速比 N/A。没有数据声称 GPU 已测得不划算；准确原因是没有可直接复用并验证的批量线代 backend。

{res_table}

这些是实际调用 wall time；stationary 时间不含后续 C3 Hessian/cone 审计，RSS 为 intact 运行进程累计峰值，不是各状态独占内存，也不把 artifact 时间跨度当 CPU 时间。硬件原始输出保存在 `provenance/HARDWARE.json`。C7 stability 没有启动。

## 14. tests

开始 full suite 一次：{baseline['passed']} passed / {baseline['failed']} failed；结束 full suite 一次：{final['passed']} passed / {final['failed']} failed。新增 40 项测试全部通过；无 skip/xfail 或新 failure。开发中只执行 constitutive/domain/equilibrium/branch 定向组，没有重复运行旧完整审计 campaign。

3 个原有失败的 ID、message、traceback 与本轮 baseline、前一 V2 final 完全一致：

{failures}

`tests/BASELINE.log`、`FINAL.log`、JUnit XML、`TEST_SUMMARY.json` 保存全部证据。

SHA256 保全校验覆盖本轮开始前 {protected['checked_files']} 个文件和 {protected['checked_symlinks']} 个符号链接，包括既有未跟踪结果；changed=0、missing=0。Production files modified=0；historical files modified=0。旧配置、Gc、L、delta、sc、材料及旧 pure-tension gate 未改。新增文件清单和 SHA 见 `provenance/NEW_SOURCE_FILES.json`。

`PREDECLARED_CONFIG.json` 保留从 V2 复制的少量旧说明字段；本轮实际覆盖关系、C3/C5 Jacobian 后端及 domain policy 在 `EFFECTIVE_CONFIG_NOTES.json` 明确列出，原数值容差不变。

## 15. scientific status

```text
{labels}
```

图仅生成六组 PNG/PDF，目录 `results/pfpd_at1_spectral_split_v1/figures/`。新模块及已通过 gates 的数值验证已完成；空间离散稳定性结论仍缺 C5/C7 有效 coupled endpoints。

## 16. 唯一下一步

人工审核后，只开展 C5 首次 AT1 active-set transition 的 stationary termination 审计：在同一已归档 ε≈0.004514117556 邻域解释为何 ftol/xtol 终止而 force/KKT 不达标，并确认怎样在原物理参数及验收门槛下获得合法 stationary state。现有局部 Jacobian 检查通过，不能预设根因。未通过该 gate 前不启动 C7 stability，也不进入任何新的材料、耦合或 post-instability 分支任务。本轮到此停止，等待人工审核。
'''
    (OUT / 'REVIEW_SPECTRAL_SPLIT_V1.md').write_text(text)


def source_inventory():
    names = ['pd_clot/pfpd_candidate/svk_spectral_split.py', 'scripts/spectral_split_audit.py',
             'scripts/run_spectral_branch.py', 'scripts/review_spectral_split.py',
             'tests/test_svk_spectral_split.py', 'tests/test_pfpd_spectral_domain.py',
             'tests/test_pfpd_mixed_strain_equilibrium.py', 'tests/test_spectral_branch_gates.py']
    before = read(OUT / 'provenance/PRE_RUN_MANIFEST.json')
    assert all(name not in before['files'] for name in names)
    save(OUT / 'provenance/NEW_SOURCE_FILES.json',
         {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in names})


if __name__ == '__main__':
    tests = test_summary()
    intact, accepted, resources, levels = summarize()
    figures(intact, accepted)
    protected = preservation()
    source_inventory()
    report(intact, accepted, resources, levels, tests, protected)
    (OUT / 'provenance/STATUS_AFTER.txt').write_text(
        subprocess.check_output(['git', 'status', '--short'], cwd=ROOT, text=True))
    print('REPORT', OUT / 'REVIEW_SPECTRAL_SPLIT_V1.md')
    print('FIGURES', FIGURES)
    print('TESTS', tests['final']['passed'], 'passed', tests['final']['failed'], 'known failed')
    print('PRESERVATION', protected)
