"""Run/review only the frozen local coefficient audit; no campaign dispatch."""
import argparse
import csv
import hashlib
import os
import resource
import subprocess
import sys
import threading
import time
import xml.etree.ElementTree as ET

import numpy as np
import psutil

from pd_clot.pfpd_candidate.affine_reference_audit import (
    ROOT, OUT, OLD, TASK, PATHS, compute, dump, historical_integrity, load,
)
from scripts.plot_affine_coefficient_reference import FIG, figures


def tables(s, a):
    f = s['F_errors']['Ftilde']['max_relative_frobenius']
    p = s['P_errors']['raw']['max_relative_frobenius']
    z = s['causal_closure']['observed_delta_Z']['max_relative_frobenius']
    summ = max(s['summation'][k]['max_relative_frobenius'] for k in ['pairwise', 'longdouble_same_terms'])
    budget = [
        dict(stage='F reconstruction', F_error=f, P_error=None, Z_error=None,
             reference_scale='||F_exact||_F', amplification_factor=None,
             interpretation='Maximum Ftilde relative Frobenius error; all 50356 unique bonds. Fbar is separately reported.'),
        dict(stage='constitutive P', F_error=f, P_error=p, Z_error=None,
             reference_scale='||F_exact||_F; ||P_exact_float64||_F', amplification_factor=p/f,
             interpretation='Ratio of maximum normalized errors, not a paired local condition number; axial tangent condition is 335.834.'),
        dict(stage='Z assembly', F_error=None, P_error=p, Z_error=z,
             reference_scale='||P_exact_float64||_F; ||Z_bulk_float64||_F', amplification_factor=z/p,
             interpretation='Raw minus exact-P Z. Stable-norm propagation attenuates rather than amplifies the relative P error.'),
        dict(stage='float64 summation', F_error=None, P_error=0., Z_error=summ,
             reference_scale='||Z_bulk_float64||_F', amplification_factor=summ/z,
             interpretation='Maximum same-term balanced/longdouble shift; amplification column is fraction of raw-minus-exact-P Z shift.'),
        dict(stage='high-precision reference', F_error=0., P_error=None,
             Z_error=s['high_precision']['exact_stored_F_matched_bulk_error'],
             reference_scale='Stress-matched bulk in numerator; original ||Z_bulk_float64||_F denominator', amplification_factor=None,
             interpretation='Spatial/weight floor. D versus historical bulk is separately 1.35291e-14 because constitutive reference precision differs.'),
    ]
    dump(OUT/'ERROR_BUDGET.json', dict(rows=budget,
         convention='All numeric errors are dimensionless relative Frobenius norms; maxima need not occur on the same bond/node. Rows are diagnostic, not an additive budget.',
         physical_scales=s['scales']))
    def csv_file(path, rows):
        with path.open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
    csv_file(OUT/'ERROR_BUDGET.csv', budget)
    csv_file(OUT/'P_SPATIAL_ORIENTATION.csv', s['P_spatial_orientation_groups'])
    rows = []
    for kind, error in s['F_errors'].items():
        rows.append(dict(quantity=kind, **{k: error[k] for k in [
            'count', 'max_component_absolute', 'max_frobenius', 'mean_frobenius', 'rms_frobenius', 'rms_component']},
            **error['frobenius_percentiles']))
    error = s['P_errors']['raw']
    rows.append(dict(quantity='P_Pa', **{k: error[k] for k in [
        'count', 'max_component_absolute', 'max_frobenius', 'mean_frobenius', 'rms_frobenius', 'rms_component']},
        **error['frobenius_percentiles']))
    csv_file(OUT/'F_P_ERROR_DISTRIBUTION.csv', rows)
    zscale = s['scales']['Zbulk_frobenius_Pa']
    rows = []
    for i, pid in enumerate(a['Z_node_ids']):
        row = dict(ID=int(pid))
        for name in PATHS:
            row[name] = float(np.sqrt(np.sum((a[name][i].astype(np.longdouble)-a['Zbulk'])**2))/zscale)
        for label in ['delta_Z_observed', 'delta_Z_predicted']:
            row[label+'_norm_Pa'] = float(np.sqrt(np.sum(a[label][i]**2)))
        row['closure_residual_bulk_scaled'] = float(np.sqrt(np.sum((a['delta_Z_predicted'][i]-a['delta_Z_observed'][i])**2))/zscale)
        rows.append(row)
    csv_file(OUT/'Z_REFERENCE_COMPARISON.csv', rows)
    return budget


def report(s, budget, tests, resources):
    st = s['status']; f = s['F_errors']; p = s['P_errors']; c = s['causal_closure']; hp = s['high_precision']
    ref = load(OUT/'EXACT_AFFINE_REFERENCE.json')
    table = '\n'.join(f"| {name} | {s['path_errors'][name]['max_bulk_relative_error']:.10e} | {s['path_errors'][name]['worst_node_ID']} |" for name in PATHS)
    distributions = '\n'.join(f"| {name} | {d['max_component_absolute']:.9e} | {d['max_frobenius']:.9e} | {d['mean_frobenius']:.9e} | {d['rms_frobenius']:.9e} | {d['frobenius_percentiles']['p50']:.9e} | {d['frobenius_percentiles']['p95']:.9e} |" for name, d in [*f.items(), ('P [Pa]', p['raw'])])
    groups = '\n'.join(f"| {r['group']} | [{r['lower']:.2f}, {min(r['upper'], 1.) if r['group']=='abs_bond_cos_z' else r['upper']:.2f}] | {r['count']} | {r['mean_delta_Pzz_Pa']:.6e} | {r['rms_P_error_Pa']:.6e} |" for r in s['P_spatial_orientation_groups'])
    status_lines = '\n'.join(f'{k}: {v}' for k, v in st.items())
    text = f'''# PFPD_BAQP_AFFINE_COEFFICIENT_REFERENCE_AUDIT_V1

**Decision: {st['Z_MISMATCH_SOURCE']}.** The raw excess is explained by reconstructed-stress floating-point propagation within the explicitly measured Z assembly floor. Exact affine stress recovers bulk coefficients. The dominant Z component is well conditioned; cancellation-sensitive summation is not the main source.

**Previous raw reconstructed-Z gate remains FAIL under its original 1e-12 criterion.** This document is a new interpretation layer. No old result, review, tolerance, production formula or accepted state was changed.

## 1. Why the previous audit stopped at NOT_ESTABLISHED

The [previous audit](../{OLD.name}/REVIEW_FULLY_CLOSED_FORCE_STENCIL_V1.md) measured N27 center force 1.945942391685588e-22 N, local normalized 8.265455792242755e-15, passing 1e-10. Its conjunction also required raw coefficients below 1e-12. Raw Z = 1.1362855097181066e-12 exceeded that threshold. Its `N27_COEFFICIENT_GATE=FAIL`, `BULK_AFFINE_FORCE_PATCH_TEST=NOT_ESTABLISHED` and primary `INCONCLUSIVE` labels remain historical facts. They are not retroactively replaced by this audit.

## 2. What is unchanged from the previous audit

h=6.25e-5 m, delta=0.00025125 m, epsilon=0.003, mu=1000 Pa, lambda=0, bulk=666.6666666666666 Pa; quadratic kernel; phase disabled, s=0. The saved N27 center D1/D2/D3 closure is reused. Saved integer offsets, oracle geometry, weights, F/P and Z are checked bitwise against the unchanged local helper before reuse. Coordinates were not previously saved; the small N27 coordinate array is reconstructed identically. Only 257 Z nodes, 1797 Fbar nodes and {s['scope']['unique_physical_bonds']:,} unique incident physical bonds are evaluated (65,792 directed Z contributions). No full-domain bond state or new force/equilibrium/activation campaign is built.

Historical input/results protection covers {historical_integrity()['protected_files']} directly relevant files, including all previous audit files/figures and the imported production constitutive, BAQP and kernel sources. [HISTORICAL_INTEGRITY.json](HISTORICAL_INTEGRITY.json) records exact SHA256 equality. No historical cleanup manifests or full pytest were run. New computation provenance is in [SOURCE_PROVENANCE.json](SOURCE_PROVENANCE.json).

## 3. Exact affine F / P reference

Starting from the current SVK energy, with lambda=0,

$$E=\\tfrac12(F^T F-I),\\quad S=2\\mu E,\\quad P=FS,\\quad W=\\mu E:E.$$

For $F=\\operatorname{{diag}}(1,1,1+\\epsilon)$,

$$E_{{zz}}=\\epsilon+\\epsilon^2/2,\\quad S_{{zz}}=2\\mu E_{{zz}},\\quad P_{{zz}}=(1+\\epsilon)S_{{zz}},\\quad W=\\mu E_{{zz}}^2.$$

All other E/S/P components are zero; Fxx=Fyy=1. For intended decimal epsilon=0.003: Ezz=0.0030045, Szz=6.009 Pa, Pzz=6.027027 Pa and W=0.00902702025 Pa. This independent closed form is saved in [EXACT_AFFINE_REFERENCE.json](EXACT_AFFINE_REFERENCE.json).

Finite representation is explicitly separated: stored float64 Fzz = `{ref['stored_Fzz_exact_decimal']}`; current float64 constitutive gives Pzz=6.027026999999707 Pa. Path B uses exactly that saved P to reproduce the prior diagnostic. Path D promotes exactly the same stored F and evaluates the analytic stress in long double. The intended-decimal reference is an additional diagnostic, never silently substituted for the historical reference. Current P differs from the closed form at the same stored F by {ref['current_P_vs_stored_F_closed_form_relative']:.6e} of the P scale, and from the intended-decimal reference by {ref['current_P_vs_decimal_closed_form_relative']:.6e}.

Formula sources: unchanged [energy.py](../../pd_clot/pfpd_candidate/energy.py), [baqp.py](../../pd_clot/pfpd_candidate/baqp.py), and [force_stencil_audit.py](../../pd_clot/pfpd_candidate/force_stencil_audit.py). The reference is independently derived here from those frozen constitutive definitions; no new external physical model is introduced.

## 4. Magnitude of reconstructed F error

All 1797 required Fbar tensors and all {s['scope']['unique_physical_bonds']:,} distinct Ftilde tensors are included, with each physical bond counted once.

| Quantity | Max component | Max Frobenius | Mean Frobenius | RMS Frobenius | Median Frobenius | P95 Frobenius |
|---|---:|---:|---:|---:|---:|---:|
{distributions}

F is dimensionless; P errors are in Pa. Mean signed Fbar zz error = {f['Fbar']['mean_signed_tensor'][2][2]:.9e}; Ftilde zz = {f['Ftilde']['mean_signed_tensor'][2][2]:.9e}. The error is mainly a systematic small negative axial bias, not zero-mean scatter. Full component means, RMS and p0/p25/p50/p75/p90/p95/p99/p100 appear in [SCIENTIFIC_RESULTS.json](SCIENTIFIC_RESULTS.json) and [F_P_ERROR_DISTRIBUTION.csv](F_P_ERROR_DISTRIBUTION.csv). Per-node/per-bond tensors are preserved in the compact [LOCAL_COEFFICIENT_ARRAYS.npz](LOCAL_COEFFICIENT_ARRAYS.npz); long-double arrays retain their dtype there.

## 5. How F error changes P

Maximum relative P Frobenius error is {p['raw']['max_relative_frobenius']:.9e}; RMS is {p['raw']['rms_relative_frobenius']:.9e}, referenced to ||Pexact||F={s['scales']['P_frobenius_Pa']:.15g} Pa. Delta Pzz is negative on {p['delta_Pzz_negative']:,} bonds, positive on {p['delta_Pzz_positive']:,}, and exactly zero on {p['delta_Pzz_zero']:,}. Mean delta Pzz={p['raw']['mean_signed_tensor'][2][2]:.9e} Pa. Orientation and position bins are descriptive and are not independent causal experiments:

| Group | Interval | Bonds | Mean delta Pzz [Pa] | RMS P error [Pa] |
|---|---:|---:|---:|---:|
{groups}

The tangent at the stored exact F is

$$\\delta E=\\operatorname{{sym}}(F^T\\delta F),\\quad \\delta P_{{lin}}=\\delta F S+2\\mu F\\delta E.$$

All bonds are evaluated. Max ||delta Praw - delta Plin||F / ||Pexact||F = {s['tangent']['raw_delta_P_minus_tangent']['max_relative_frobenius']:.6e}; with high-precision stresses this becomes {s['tangent']['high_precision_delta_P_minus_tangent']['max_relative_frobenius']:.6e}. The raw RMS remainder is {s['tangent']['raw_remainder_over_delta_P_rms']:.3%} of the RMS perturbation. This supports first-order sensitivity plus arithmetic roundoff, with no resolvable nonlinear physical discrepancy at this scale.

For the axial mode $P_{{zz}}=\\mu(a^3-a)$, $dP_{{zz}}/da=\\mu(3a^2-1)$, so $|a P'_{{zz}}/P_{{zz}}|={s['tangent']['estimated_relative_condition_number_axial']:.6f}$. Stress measures a small departure from the identity. A roughly 1e-15 F error can therefore produce roughly 1e-12 **relative stress error**; the stress error itself is not merely 1e-15. Transverse bonds retain more axial reconstruction error than nearly axial bonds under the existing bond correction.

## 6. How delta P propagates into Z

At fixed geometry and weights, define $Q_{{ij}}=I-\\xi_{{ij}}\\otimes\\xi_{{ij}}/r_{{ij}}^2$ and

$$Z_i(P)=\\sum_j\\bar\\omega_{{ij}}P_{{ij}}Q_{{ij}}V_j,\\qquad
\\Delta Z_{{pred,i}}=\\sum_j\\bar\\omega_{{ij}}(P_{{raw,ij}}-P_*)Q_{{ij}}V_j.$$

The same float64 endpoint order, projections, neighbor volumes and weights are used for raw, exact-P, and delta-P assembly. No boundary normalization is changed. Every bond P is verified bitwise against the unchanged selective evaluator. The stable scale is ||Zbulk||F = {s['scales']['Zbulk_frobenius_Pa']:.15g} Pa. Ratios never divide by an individual near-zero tensor entry.

## 7. Does delta P reconstruct the observed delta Z?

**{st['DELTA_Z_FROM_DELTA_P_RECONSTRUCTION']}**, over all 257 tensors and all nine components. Observed delta Z = Zraw - ZexactP has max Frobenius norm {c['observed_delta_Z']['max_frobenius']:.9e} Pa and max component {c['observed_delta_Z']['max_component_absolute']:.9e} Pa. Max observed norm / bulk scale = {c['observed_delta_Z']['max_relative_frobenius']:.9e}.

The max causal residual ||delta Zpred - delta Zobs||F is {c['max_frobenius']:.9e} Pa; max component is {c['max_component_absolute']:.9e} Pa; relative to the stable scale it is **{c['max_relative_frobenius']:.9e} <= 1e-13**. RMS residual/bulk = {c['rms_relative_frobenius']:.9e}. Floating-point subtraction and separately rounded assemblies need not be exactly linear. The residual is reported, not forced to zero. Its maximum is {c['max_relative_frobenius']/c['observed_delta_Z']['max_relative_frobenius']:.3%} of the maximum observed shift (different maxima may occur at different nodes).

## 8. Role of floating-point summation

Holding the complete float64 terms fixed, sequential versus balanced summation differs by at most {s['summation']['pairwise']['max_relative_frobenius']:.9e} of the bulk scale ({s['summation']['pairwise']['max_frobenius']:.9e} Pa); sequential versus long-double accumulation differs by {s['summation']['longdouble_same_terms']['max_relative_frobenius']:.9e} ({s['summation']['longdouble_same_terms']['max_frobenius']:.9e} Pa). Classification: **{st['SUMMATION_ORDER_EFFECT']}** relative to the 1e-12 exceedance and predeclared 1e-13 diagnostic criterion.

At the center, zz sum(abs(terms))/abs(sum(terms)) is approximately 1.00. ||sum(abs(terms))||F/||sum(terms)||F={s['summation']['center']['norm_sumabs_over_norm_sum']:.9f}; sum(||term||F)/||sum(terms)||F={s['summation']['center']['sum_term_frobenius_over_norm_sum']:.9f}. zx/zy undergo symmetry cancellation, but their final values are approximately 1e-16 Pa and do not dominate the stable tensor norm. Ratios for unsafe denominators are null, with a declared safety threshold of {s['summation']['center']['safe_denominator_Pa']:.6e} Pa. Thus the proposed cancellation-amplification story is **not supported for the dominant Z norm**. The main amplification occurs in the relative F-to-P constitutive sensitivity.

## 9. High-precision reference result

This platform's long double has {hp['longdouble_mantissa_bits']} significand bits and epsilon={hp['longdouble_eps']:.9e}; its {hp['storage_bits']}-bit storage does not mean IEEE quad arithmetic. mpmath was not installed; no package was installed. Pure 3x3 algebra and summation use long double, while physical geometry, K/B, weights and projection remain the identical float64 inputs.

| Path | Max Z error / original bulk scale | Worst node ID |
|---|---:|---:|
{table}

All four primary rows compare to the **same saved float64 bulk tensor**. Path C improves constitutive and accumulation precision but retains the reconstructed F error, so its approximately 1e-12 discrepancy persists. A separate C diagnostic casts high-precision P back to float64 before historical Z assembly and gives {s['summation']['high_precision_P_with_float64_assembly_max_Z_error']:.9e}, confirming the same conclusion without changing accumulation precision.

For D, the high-precision bulk versus historical float64 bulk differs by {hp['stored_F_high_precision_bulk_vs_historical']:.9e} of the stable scale. Against the **matched-stress bulk**, D's spatial/weight discrepancy is only **{hp['exact_stored_F_matched_bulk_error']:.9e}**. Both references are reported; the numerical floor is not hidden by a precision-dependent reference offset. The optional intended-decimal F diagnostic has matched bulk error {hp['intended_decimal_F_matched_bulk_error']:.9e}, and historical-reference error {hp['intended_decimal_F_vs_historical_bulk_error']:.9e}. Analytic versus general long-double P at the same stored F differs by {hp['analytic_vs_general_high_precision_P_error']:.9e} Pa.

## 10. Is the 1.136e-12 mismatch physical or numerical?

**{st['Z_MISMATCH_SOURCE']}**. Case A is satisfied: raw and exact-P historical values reproduce, exact-P Z is below the original coefficient tolerance, high-precision matched reference is at the fixed-input arithmetic floor, and delta-P reconstructs the raw shift within the declared causal tolerance. No structural coefficient discrepancy is detected. This conclusion is limited to the frozen, closed, intact affine stencil.

[ERROR_BUDGET.json](ERROR_BUDGET.json) and [ERROR_BUDGET.csv](ERROR_BUDGET.csv) record reconstruction, constitutive, Z assembly, summation and reference contributions. Entries use explicitly named relative Frobenius scales. The rows are not an additive error budget. Maximum relative P/F error ratio={budget[1]['amplification_factor']:.6f}; Z-shift/P ratio={budget[2]['amplification_factor']:.6f}. This is consistent with constitutive sensitivity followed by a well-conditioned stable-norm Z assembly.

## 11. Updated scientific interpretation

```text
{status_lines}
```

These are this audit's labels. In particular `BULK_OPERATOR_CONSISTENCY=SUPPORTED` does not rewrite the prior raw-gate failure or make a universal operator theorem.

## 12. What this means for finite-domain boundary influence

The previous unchanged N27 noncentral one-hop/two-hop observations have local normalized force approximately 1.53e-2 / 3.53e-7, while the fully closed center passes force and now has an independently explained coefficient reference. Taken together, the evidence supports finite-domain multi-hop influence as the source of prior intact mechanical contamination. This is a scientific inference combining the saved stencil audit with the present reference audit; no new boundary force campaign is claimed. Center inversion symmetry remains a caveat for inferring minimal closure depth.

## 13. What remains unresolved for clot-wall modeling

The correct physical nonlocal loading/interaction-domain condition is unresolved. This audit does not establish that surface Dirichlet is mathematically invalid, that C5/C7 activation is mesh objective, that clot-wall fracture is validated, or that 3delta is the required physical collar thickness. Graph dependency reach is not a physical boundary prescription. Current C5/C7 states and results are unchanged; fracture, stability, phase and CFD were not run.

Six targeted tests only: **{tests['passed']} passed, {tests['failed']} failed, {tests['errors']} errors, {tests['skipped']} skipped**; [TARGETED_TESTS.log](TARGETED_TESTS.log). The tests cover analytic SVK, raw reproduction including saved representative tensors, exact-P reproduction, causal assembly, accumulation/reference precision, and unchanged historical bytes. No old test assertions were changed.

Measured numerical computation/export + final figure/test/report phase wall time: {resources['wall_seconds']:.3f} s. CPU workers=1, BLAS threads=1, peak process-tree RSS={resources['peak_RSS_GiB']:.3f} GiB; GPU=NOT_USED_BY_DESIGN. Resource scope and sampling method are in [RESOURCES.json](RESOURCES.json); authoring/inspection time is not computational runtime.

Scientific figures (English, PNG 600 dpi and vector PDF):

''' + '\n'.join(f'- [{name} — PNG](../../results/{TASK}/figures/{name}.png) · [PDF](../../results/{TASK}/figures/{name}.pdf)' for name in [
       'error_propagation_chain', 'z_reference_comparison', 'delta_z_closure', 'z_accumulation_conditioning']) + f'''

## 14. One next development task

`NONLOCAL_BOUNDARY_FORMULATION_V1` is recommended for **separate human review**, conditional on this Case A interpretation being accepted. Define the physical load and interaction domain, then compare its traction/virtual-work target without tuning collars to match activation thresholds. This task is **not started**. No ghosts, volume constraints, force-density collars, surface corrections or boundary fixes are implemented here.

Reproduce only this local audit from the project root:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 MPLCONFIGDIR=/tmp/clot_affine_reference_mpl /home/lzy/projects/computation_examples/_shared/python_environment/bin/python -B -m scripts.affine_coefficient_reference_audit --recompute
```

Omit `--recompute` to re-export figures/report and run the six focused tests from the saved local tensors. Historical SHA256 protection is checked before and after. Stop for human review.
'''
    (OUT/'REVIEW_AFFINE_COEFFICIENT_REFERENCE_V1.md').write_text(text)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--recompute', action='store_true')
    args = parser.parse_args()
    assert os.environ.get('OPENBLAS_NUM_THREADS') == '1'
    assert os.environ.get('OMP_NUM_THREADS') == '1'
    assert historical_integrity()['PASS']
    start = time.perf_counter()
    peak = [0]
    stop = threading.Event()
    process = psutil.Process()
    def monitor():
        while not stop.is_set():
            total = 0
            for p in [process]+process.children(recursive=True):
                try: total += p.memory_info().rss
                except (psutil.NoSuchProcess, psutil.AccessDenied): pass
            peak[0] = max(peak[0], total)
            stop.wait(.02)
    thread = threading.Thread(target=monitor, daemon=True); thread.start()
    if args.recompute or not (OUT/'SCIENTIFIC_RESULTS.json').exists():
        s, a, reference = compute()
        np.savez_compressed(OUT/'LOCAL_COEFFICIENT_ARRAYS.npz', **a)
        dump(OUT/'SCIENTIFIC_RESULTS.json', s)
        dump(OUT/'EXACT_AFFINE_REFERENCE.json', reference)
        compute_resources = dict(wall_seconds=0, peak_RSS_GiB=0,
                                 note='Compute included in this invocation, not added twice.')
    else:
        s = load(OUT/'SCIENTIFIC_RESULTS.json')
        with np.load(OUT/'LOCAL_COEFFICIENT_ARRAYS.npz', allow_pickle=False) as source:
            a = {k: source[k] for k in source.files}
        compute_resources = load(OUT/'COMPUTE_RESOURCES.json')
    budget = tables(s, a)
    figures(s, a)
    test_start = time.perf_counter()
    run = subprocess.run([sys.executable, '-B', '-m', 'pytest', '-q', '-p', 'no:cacheprovider',
                          'tests/test_affine_coefficient_reference_audit.py',
                          '--junitxml='+str(OUT/'TARGETED_TESTS.xml')], cwd=ROOT,
                         text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    (OUT/'TARGETED_TESTS.log').write_text(run.stdout)
    suites = ET.parse(OUT/'TARGETED_TESTS.xml').getroot()
    counts = {k: sum(int(x.attrib.get(k, 0)) for x in suites.iter('testsuite'))
              for k in ['tests', 'failures', 'errors', 'skipped']}
    tests = dict(passed=counts['tests']-counts['failures']-counts['errors']-counts['skipped'],
                 failed=counts['failures'], errors=counts['errors'], skipped=counts['skipped'],
                 exit_code=run.returncode, wall_seconds=time.perf_counter()-test_start)
    dump(OUT/'TARGETED_TESTS.json', tests)
    integrity = historical_integrity(); dump(OUT/'HISTORICAL_INTEGRITY.json', integrity)
    if run.returncode or not integrity['PASS']:
        print(run.stdout); raise SystemExit('Focused validation failed; inspect saved evidence.')
    dump(OUT/'STATUS.json', s['status'])
    sources = ['pd_clot/pfpd_candidate/affine_reference_audit.py',
               'scripts/affine_coefficient_reference_audit.py', 'scripts/plot_affine_coefficient_reference.py',
               'tests/test_affine_coefficient_reference_audit.py']
    dump(OUT/'SOURCE_PROVENANCE.json', {p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in sources})
    resources = dict(wall_seconds=compute_resources['wall_seconds']+time.perf_counter()-start,
                     compute_phase=compute_resources, final_phase_wall_seconds=time.perf_counter()-start,
                     CPU_workers=1, BLAS_threads=1, GPU='NOT_USED_BY_DESIGN',
                     peak_RSS_GiB=max(compute_resources['peak_RSS_GiB'], peak[0]/1024**3,
                                      resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2),
                     method='20 ms process-tree RSS sampling for final phase; OS ru_maxrss for compute. Separate phase maxima, not sums.',
                     scope='Numerical compute plus local array export, final plots, six targeted tests and report; excludes interactive authoring and preview.')
    report(s, budget, tests, resources)
    resources['wall_seconds'] = compute_resources['wall_seconds']+time.perf_counter()-start
    resources['final_phase_wall_seconds'] = time.perf_counter()-start
    dump(OUT/'RESOURCES.json', resources)
    stop.set(); thread.join()
    starting = load(OUT/'STARTING_GIT_STATE.json')
    st = s['status']
    lines = [f"branch: {starting['branch']}", f"HEAD: {starting['HEAD']}",
             'epsilon: 0.003', 'mu: 1000 Pa', 'lambda: 0 Pa', 'h: 6.25e-05 m', 'delta: 0.00025125 m']
    lines += [f'{k}: {v}' for k, v in st.items()]
    lines += [f"delta_Z reconstruction error / bulk: {s['causal_closure']['max_relative_frobenius']:.12e}",
              f"float64 vs pairwise / bulk: {s['summation']['pairwise']['max_relative_frobenius']:.12e}",
              f"float64 vs longdouble / bulk: {s['summation']['longdouble_same_terms']['max_relative_frobenius']:.12e}",
              f"targeted tests: {tests['passed']} passed / {tests['failed']} failed",
              f"wall time: {resources['wall_seconds']:.3f} s", 'CPU workers: 1; BLAS threads: 1',
              f"peak RAM: {resources['peak_RSS_GiB']:.3f} GiB", 'GPU: NOT_USED_BY_DESIGN',
              f'REVIEW: {OUT}/REVIEW_AFFINE_COEFFICIENT_REFERENCE_V1.md', f'figures: {FIG}']
    terminal = '\n'.join(lines)+'\n'
    (OUT/'TERMINAL_SUMMARY.txt').write_text(terminal)
    print(terminal, end='')


if __name__ == '__main__':
    main()
