"""Local, intact N27 coefficient reference diagnostics; never a solver or BC.

All historical kernels and the two-pass endpoint order are imported unchanged.
Long double is used only for algebra and sums; geometry and weights stay frozen.
"""
import hashlib
import json
from pathlib import Path

import numpy as np

from .energy import constitutive
from .force_stencil_audit import AffineStencil, H, DELTA, EPSILON, material

TASK = 'pfpd_baqp_affine_coefficient_reference_audit_v1'
PREVIOUS = 'pfpd_baqp_fully_closed_force_stencil_audit_v1'
ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'verification' / TASK
OLD = ROOT / 'verification' / PREVIOUS
PATHS = ('RAW_RECONSTRUCTED_FLOAT64', 'EXACT_AFFINE_P_FLOAT64',
         'RECONSTRUCTED_F_HIGH_PRECISION_P', 'EXACT_AFFINE_HIGH_PRECISION_REFERENCE')
LD = np.longdouble


def load(path):
    return json.loads(Path(path).read_text())


def serial(x):
    if isinstance(x, np.ndarray):
        # Long-double tensors are additionally retained losslessly in NPZ.
        return np.asarray(x, dtype=float).tolist()
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, np.floating):
        return float(x)
    if isinstance(x, np.bool_):
        return bool(x)
    raise TypeError(type(x))


def dump(path, data):
    Path(path).write_text(json.dumps(data, default=serial, indent=2, allow_nan=False)+'\n')


def frob(x):
    return np.sqrt(np.sum(np.asarray(x)**2, axis=(-2, -1)))


def pairwise(terms):
    """Deterministic adjacent balanced reduction, preserving every term."""
    a = np.array(terms, copy=True)
    while len(a) > 1:
        pairs = len(a)//2
        b = a[:2*pairs:2] + a[1:2*pairs:2]
        a = np.concatenate((b, a[-1:]), axis=0) if len(a) % 2 else b
    return a[0]


def svk_extended(F):
    """Independent lambda=0 algebra with no float64 conversion or linalg."""
    F = np.asarray(F, dtype=LD)
    E = (np.swapaxes(F, -1, -2)@F-np.eye(3, dtype=LD))/LD(2)
    S = LD(2000)*E
    return E, S, F@S, LD(1000)*np.sum(E*E, axis=(-2, -1))


def closed_form(epsilon, dtype=LD):
    e = dtype(epsilon)
    F = np.diag(np.array([1, 1, 1+e], dtype=dtype))
    E = np.zeros((3, 3), dtype=dtype)
    E[2, 2] = e+e*e/dtype(2)
    S = dtype(2000)*E
    P = F@S
    W = dtype(1000)*E[2, 2]**2
    return dict(F=F, E=E, S=S, P=P, W=W, Pzz=P[2, 2])


def stats(delta, reference_scale):
    a = np.asarray(delta, dtype=LD)
    norms = frob(a)
    percentile = dict(zip(['p0', 'p25', 'p50', 'p75', 'p90', 'p95', 'p99', 'p100'],
                          np.percentile(np.asarray(norms, float), [0, 25, 50, 75, 90, 95, 99, 100])))
    return dict(count=len(a), max_component_absolute=np.max(abs(a)),
                max_frobenius=np.max(norms), mean_frobenius=np.mean(norms),
                rms_frobenius=np.sqrt(np.mean(norms**2)),
                rms_component=np.sqrt(np.mean(a*a)), mean_signed_tensor=np.mean(a, axis=0),
                frobenius_percentiles=percentile, reference_scale=reference_scale,
                max_relative_frobenius=np.max(norms)/reference_scale,
                rms_relative_frobenius=np.sqrt(np.mean(norms**2))/reference_scale)


def historical_integrity():
    original = load(OUT/'HISTORICAL_INPUT_HASHES.json')
    differences = [p for p, sha in original.items()
                   if not (ROOT/p).is_file() or hashlib.sha256((ROOT/p).read_bytes()).hexdigest() != sha]
    return dict(protected_files=len(original), differences=differences,
                PASS=not differences)


def compute():
    cfg = load(OUT/'PREDECLARED_CONFIG.json')
    old_cfg = load(OLD/'PREDECLARED_CONFIG.json')
    for key, expected in [('h_m', H), ('delta_m', DELTA), ('engineering_strain', EPSILON),
                          ('mu_Pa', 1000.), ('lambda_Pa', 0.), ('kernel', 'quadratic')]:
        assert old_cfg[key] == expected
    assert material().lame == 0 and cfg['coefficient_tolerance'] == 1e-12
    oracle = load(OLD/'BULK_LATTICE_ORACLE.json')
    previous = load(OLD/'N27.json')
    saved_center = next(r for r in previous['observations'] if r['label'] == 'center')
    assert saved_center['closure']['FULL_FORCE_STENCIL_CLOSED']

    # Coordinates/Fbar/bond tensors were not saved: reconstruct only the local
    # coefficient objects, never the full N27 campaign, force or closure graph.
    grid = AffineStencil(27)
    for key in ['offsets', 'xi', 'omega', 'wbar', 'K_inf', 'F', 'P', 'Z_inf']:
        assert np.array_equal(grid.oracle[key], oracle[key]), key
    grid.offsets = np.asarray(oracle['offsets'], dtype=int)
    grid.oracle = {k: np.asarray(v) if isinstance(v, list) else v for k, v in oracle.items()}
    center = saved_center['ID']
    assert np.array_equal(grid.X[center], saved_center['X_m'])
    first, second = grid.coefficient_sets(center)
    assert (len(first), len(second)) == (257, 1797)
    Fexact = np.asarray(oracle['F'])
    Pexact = np.asarray(oracle['P'])
    Zbulk = np.asarray(oracle['Z_inf'])
    fscale, pscale, zscale = frob(Fexact), frob(Pexact), frob(Zbulk)
    Fbar = np.array([grid.Fbar(int(pid)) for pid in second])
    Fb_lookup = {int(pid): k for k, pid in enumerate(second)}

    rows = []
    pairs = []
    for pid in first:
        ids, xi, r, w = grid.row(int(pid))
        row_pairs = np.column_stack((np.minimum(pid, ids), np.maximum(pid, ids)))
        pairs.append(row_pairs)
        wb = grid.weight(int(pid), ids, w, 'original')
        projection = np.eye(3)-xi[:, :, None]*xi[:, None, :]/r[:, None, None]**2
        rows.append((ids, projection, wb))
    bonds, inverse = np.unique(np.concatenate(pairs), axis=0, return_inverse=True)
    row_bonds = inverse.reshape(len(first), 256)
    lo, hi = bonds.T
    av = .5*(Fbar[[Fb_lookup[int(p)] for p in lo]]+Fbar[[Fb_lookup[int(p)] for p in hi]])
    xi = grid.X[hi]-grid.X[lo]
    eta = grid.x[hi]-grid.x[lo]
    r2 = np.linalg.norm(xi, axis=1)**2
    residual = eta-np.einsum('bij,bj->bi', av, xi)
    Ft = av+residual[:, :, None]*xi[:, None, :]/r2[:, None, None]
    Praw = constitutive(Ft, material())[1]
    _, _, Php, _ = svk_extended(Ft)
    # Preserve the SAME stored affine F as paths A/B, promoted losslessly.
    # Also retain the intended decimal-strain analytic reference separately.
    stored_ref = closed_form(LD(Fexact[2, 2])-LD(1))
    decimal_ref = closed_form('0.003')
    Phpexact = stored_ref['P']
    Ehp, Shp, Pgeneralhp, Whp = svk_extended(Fexact)
    dF = Ft-Fexact
    dP = Praw-Pexact
    dElinear = (np.swapaxes(dF.astype(LD), -1, -2)@Fexact.astype(LD)
                +Fexact.astype(LD).T@dF.astype(LD))/LD(2)
    dPlinear = dF.astype(LD)@stored_ref['S']+Fexact.astype(LD)@(LD(2000)*dElinear)
    dPhp = Php-Phpexact

    paths = {p: [] for p in PATHS}
    extra = {p: [] for p in ['delta_Z_predicted', 'delta_Z_observed', 'Z_pairwise',
                             'Z_longdouble_same_terms', 'Z_exact_decimal', 'Z_Php_float64_assembly']}
    center_condition = None
    for k, (pid, (ids, projection, wb)) in enumerate(zip(first, rows)):
        b = row_bonds[k]
        volume = grid.V[ids, None, None]
        # Recheck P against the unchanged, previously validated incident routine.
        assert np.array_equal(Praw[b], grid.bond_stress(bonds[b, 0], bonds[b, 1]))
        def terms(P, high=False):
            dtype = LD if high else np.float64
            return (wb.astype(dtype)[:, None, None]*(np.asarray(P, dtype=dtype)
                    @projection.astype(dtype)))*volume.astype(dtype)
        traw = terms(Praw[b])
        raw = np.add.accumulate(traw, axis=0)[-1]
        exact = np.add.accumulate(terms(Pexact), axis=0)[-1]
        reconstructed_hp = np.add.accumulate(terms(Php[b], True), axis=0)[-1]
        exact_hp = np.add.accumulate(terms(Phpexact, True), axis=0)[-1]
        for name, Z in zip(PATHS, [raw, exact, reconstructed_hp, exact_hp]):
            paths[name].append(Z)
        extra['delta_Z_predicted'].append(np.add.accumulate(terms(dP[b]), axis=0)[-1])
        extra['delta_Z_observed'].append(raw-exact)
        extra['Z_pairwise'].append(pairwise(traw))
        extra['Z_longdouble_same_terms'].append(np.add.accumulate(traw.astype(LD), axis=0)[-1])
        extra['Z_exact_decimal'].append(np.add.accumulate(terms(decimal_ref['P'], True), axis=0)[-1])
        extra['Z_Php_float64_assembly'].append(np.add.accumulate(terms(Php[b]), axis=0)[-1])
        if int(pid) == center:
            sumabs = np.sum(abs(traw).astype(LD), axis=0)
            safe = LD(100)*np.finfo(float).eps*LD(zscale)
            ratios = [[float(sumabs[i, j]/abs(raw[i, j])) if abs(raw[i, j]) > safe else None
                       for j in range(3)] for i in range(3)]
            center_condition = dict(Z=raw, sum_absolute_terms=sumabs, safe_denominator_Pa=safe,
                                    component_cancellation_factor=ratios,
                                    unsafe_denominator_policy='null: component too small relative to stable bulk scale',
                                    norm_sumabs_over_norm_sum=frob(sumabs)/frob(raw),
                                    sum_term_frobenius_over_norm_sum=np.sum(frob(traw))/frob(raw))

    paths = {k: np.array(v) for k, v in paths.items()}
    extra = {k: np.array(v) for k, v in extra.items()}
    raw, exact, reconstructed_hp, exact_hp = [paths[p] for p in PATHS]
    # Fixed saved geometry in the infinite-lattice row for a matched-stress
    # high-precision reference; no change to K, B, weights or projection.
    xib = np.asarray(oracle['xi']); rb = np.linalg.norm(xib, axis=1)
    projb = np.eye(3)-xib[:, :, None]*xib[:, None, :]/rb[:, None, None]**2
    def matched_bulk(P):
        t = np.asarray(oracle['wbar'], dtype=LD)[:, None, None]*(P@projb.astype(LD))*LD(H**3)
        return np.add.accumulate(t, axis=0)[-1]
    Zbulkhp = matched_bulk(Phpexact)
    Zbulkdecimal = matched_bulk(decimal_ref['P'])

    path_errors = {}
    for name, Z in paths.items():
        error = frob(Z.astype(LD)-Zbulk.astype(LD))/LD(zscale)
        worst = int(np.argmax(error))
        path_errors[name] = dict(max_bulk_relative_error=error[worst], worst_node_ID=int(first[worst]),
                                 rms_bulk_relative_error=np.sqrt(np.mean(error**2)),
                                 center_bulk_relative_error=error[np.where(first == center)[0][0]])
    closure = extra['delta_Z_predicted']-extra['delta_Z_observed']
    closure_stats = stats(closure, zscale)
    closure_stats['observed_delta_Z'] = stats(extra['delta_Z_observed'], zscale)
    closure_stats['predicted_delta_Z'] = stats(extra['delta_Z_predicted'], zscale)
    closure_stats['PASS'] = bool(closure_stats['max_relative_frobenius'] <= cfg['causal_closure_tolerance_Zbulk_scaled'])
    summation = dict(pairwise=stats(extra['Z_pairwise']-raw, zscale),
                     longdouble_same_terms=stats(extra['Z_longdouble_same_terms']-raw.astype(LD), zscale),
                     high_precision_P_with_float64_assembly_max_Z_error=np.max(frob(extra['Z_Php_float64_assembly']-Zbulk))/zscale,
                     convention='Same float64 P, weights, projection, and precomputed terms; only reduction changes',
                     center=center_condition)
    fstats = dict(Fbar=stats(Fbar-Fexact, fscale), Ftilde=stats(dF, fscale))
    pstats = dict(raw=stats(dP, pscale), high_precision=stats(dPhp, pscale),
                  constitutive_float64_minus_high_precision=stats(Praw.astype(LD)-Php, pscale),
                  delta_Pzz_positive=int(np.sum(dP[:, 2, 2] > 0)),
                  delta_Pzz_negative=int(np.sum(dP[:, 2, 2] < 0)),
                  delta_Pzz_zero=int(np.sum(dP[:, 2, 2] == 0)))
    direction_cos = abs(xi[:, 2])/np.sqrt(r2)
    mid = .5*(grid.X[lo]+grid.X[hi])
    groups = []
    for label, values, edges in [('abs_bond_cos_z', direction_cos, [0, .25, .5, .75, 1.0000001]),
                                 ('midpoint_z_over_h', mid[:, 2]/H, [-8.01, -4, 0, 4, 8.01])]:
        for left, right in zip(edges[:-1], edges[1:]):
            mask = (values >= left) & (values < right)
            if np.any(mask):
                groups.append(dict(group=label, lower=left, upper=right, count=int(mask.sum()),
                                   mean_delta_Pzz_Pa=np.mean(dP[mask, 2, 2]),
                                   rms_P_error_Pa=np.sqrt(np.mean(frob(dP[mask])**2)),
                                   max_P_error_Pa=np.max(frob(dP[mask]))))
    tangent = dict(raw_delta_P_minus_tangent=stats(dP.astype(LD)-dPlinear, pscale),
                   high_precision_delta_P_minus_tangent=stats(dPhp-dPlinear, pscale),
                   estimated_relative_condition_number_axial=float(Fexact[2, 2]*(1000*(3*Fexact[2, 2]**2-1))/Pexact[2, 2]),
                   formula='dP = dF S + F [2 mu sym(F^T dF)] at lambda=0',
                   raw_remainder_over_delta_P_rms=float(np.sqrt(np.mean(frob(dP.astype(LD)-dPlinear)**2))/np.sqrt(np.mean(frob(dP)**2))))
    hp = dict(longdouble_eps=np.finfo(LD).eps, longdouble_mantissa_bits=np.finfo(LD).nmant+1,
              storage_bits=np.dtype(LD).itemsize*8, mpmath='NOT_INSTALLED_NOT_REQUIRED',
              exact_stored_F_matched_bulk_error=np.max(frob(exact_hp-Zbulkhp))/LD(zscale),
              intended_decimal_F_matched_bulk_error=np.max(frob(extra['Z_exact_decimal']-Zbulkdecimal))/LD(zscale),
              intended_decimal_F_vs_historical_bulk_error=np.max(frob(extra['Z_exact_decimal']-Zbulk.astype(LD)))/LD(zscale),
              stored_F_high_precision_bulk_vs_historical=frob(Zbulkhp-Zbulk.astype(LD))/LD(zscale),
              analytic_vs_general_high_precision_P_error=frob(Phpexact-Pgeneralhp),
              exact_stored_F=stored_ref, intended_decimal_F=decimal_ref,
              explanation='D preserves the stored float64 F exactly; matched bulk removes only the reference constitutive precision offset. Decimal epsilon is a separate analytic reference, not substituted into A/B.')
    historical_error = load(OLD/'STATUS.json')['Z_BULK_MATCH']
    saved_exact = previous['center_stress_roundoff_diagnostic']['constant_P_Z_max_bulk_relative_error']
    reproduction = dict(raw_saved_error=historical_error, exact_P_saved_error=saved_exact,
                         raw_difference=path_errors[PATHS[0]]['max_bulk_relative_error']-historical_error,
                         exact_P_difference=path_errors[PATHS[1]]['max_bulk_relative_error']-saved_exact)
    consistent = all(abs(reproduction[k]) <= cfg['saved_reproduction_absolute_tolerance_on_relative_error']
                     for k in ['raw_difference', 'exact_P_difference'])
    summation_max = max(summation[k]['max_relative_frobenius'] for k in ['pairwise', 'longdouble_same_terms'])
    case_a = (consistent and path_errors[PATHS[1]]['max_bulk_relative_error'] <= 1e-12
              and hp['exact_stored_F_matched_bulk_error'] <= cfg['high_precision_noise_criterion_Zbulk_scaled']
              and closure_stats['PASS'] and historical_integrity()['PASS'])
    case_b = path_errors[PATHS[1]]['max_bulk_relative_error'] > 1e-12
    status = dict(HISTORICAL_RAW_Z_GATE='FAIL_UNCHANGED', RAW_Z_ERROR=path_errors[PATHS[0]]['max_bulk_relative_error'],
                  EXACT_AFFINE_P_Z_ERROR=path_errors[PATHS[1]]['max_bulk_relative_error'],
                  HIGH_PRECISION_Z_ERROR=path_errors[PATHS[3]]['max_bulk_relative_error'],
                  HIGH_PRECISION_MATCHED_REFERENCE_Z_ERROR=hp['exact_stored_F_matched_bulk_error'],
                  MAX_F_RECONSTRUCTION_ERROR=max(fstats[k]['max_component_absolute'] for k in fstats),
                  MAX_P_ERROR=pstats['raw']['max_component_absolute'],
                  DELTA_Z_FROM_DELTA_P_RECONSTRUCTION='PASS' if closure_stats['PASS'] else 'FAIL',
                  SUMMATION_ORDER_EFFECT='NEGLIGIBLE' if summation_max <= cfg['summation_negligible_criterion_Zbulk_scaled'] else 'MATERIAL',
                  Z_MISMATCH_SOURCE='RECONSTRUCTED_STRESS_FLOATING_POINT_PROPAGATION' if case_a else ('STRUCTURAL_COEFFICIENT_MISMATCH' if case_b else 'INCONCLUSIVE'),
                  BULK_COEFFICIENT_PHYSICS='CONSISTENT_WITH_EXACT_AFFINE_REFERENCE' if case_a else ('INCONSISTENT' if case_b else 'INCONCLUSIVE'),
                  BAQP_BULK_OPERATOR_DEFECT='NOT_SUPPORTED' if case_a else ('SUPPORTED' if case_b else 'INCONCLUSIVE'),
                  BULK_OPERATOR_CONSISTENCY='SUPPORTED' if case_a else 'NOT_ESTABLISHED',
                  FINITE_DOMAIN_MULTI_HOP_BOUNDARY_INFLUENCE='SUPPORTED' if case_a else 'NOT_ESTABLISHED',
                  BOUNDARY_MODEL_DEVELOPMENT_READINESS='READY_FOR_SEPARATE_REVIEW' if case_a else 'NOT_READY',
                  CURRENT_C5_C7_ACTIVATION_RESULTS='UNCHANGED', BOUNDARY_FIX='NOT_IMPLEMENTED', FRACTURE='NOT_RUN')
    summary = dict(configuration=cfg, scope=dict(Z_nodes=len(first), Fbar_nodes=len(second),
                   unique_physical_bonds=len(bonds), oriented_Z_contributions=len(first)*256,
                   saved_closure=saved_center['closure'], full_domain_bond_state=False,
                   geometry_reconstruction='N27 coordinates only; all coefficient work restricted to the center stencil'),
                   scales=dict(F_frobenius=fscale, P_frobenius_Pa=pscale, Zbulk_frobenius_Pa=zscale),
                   path_errors=path_errors, F_errors=fstats, P_errors=pstats, P_spatial_orientation_groups=groups,
                   tangent=tangent, causal_closure=closure_stats, summation=summation, high_precision=hp,
                   saved_reproduction=reproduction, status=status)
    arrays = dict(Z_node_ids=first, Fbar_node_ids=second, bonds=bonds, bond_xi=xi, bond_midpoints=mid,
                  Fbar=Fbar, Ftilde=Ft, Praw=Praw, P_high_precision=Php, delta_P=dP,
                  delta_P_tangent=dPlinear, Zbulk=Zbulk, Zbulk_high_precision=Zbulkhp,
                  F_exact=Fexact, P_exact=Pexact, P_exact_high_precision=Phpexact,
                  **paths, **extra)
    reference = dict(F=decimal_ref['F'], E=decimal_ref['E'], S=decimal_ref['S'], P=decimal_ref['P'],
                     W=decimal_ref['W'], Pzz=decimal_ref['Pzz'],
                     decimal_exact_strings=dict(epsilon='0.003', Fzz='1.003', Ezz='0.0030045',
                                                Szz='6.009', Pzz='6.027027', W='0.00902702025'),
                     units=dict(F='dimensionless', E='dimensionless', S='Pa', P='Pa', W='Pa'),
                     stored_float64_F=Fexact, current_constitutive_P=Pexact,
                     stored_float64_F_high_precision={k: v for k, v in stored_ref.items()},
                     stored_Fzz_exact_decimal=format(float(Fexact[2, 2]), '.55f'),
                     current_P_vs_decimal_closed_form_relative=frob(Pexact.astype(LD)-decimal_ref['P'])/pscale,
                     current_P_vs_stored_F_closed_form_relative=frob(Pexact.astype(LD)-Phpexact)/pscale,
                     derivation='Ezz=epsilon+epsilon^2/2; Szz=2 mu Ezz; Pzz=(1+epsilon)Szz; W=mu Ezz^2; all other E,S,P entries zero')
    return summary, arrays, reference
