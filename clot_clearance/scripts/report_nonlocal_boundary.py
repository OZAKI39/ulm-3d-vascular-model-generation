"""Assemble an auditable report without rerunning any scientific solve."""
import hashlib
import json
import os
import platform
import sys
from importlib.metadata import version
from pathlib import Path
import numpy as np
from types import SimpleNamespace
from .nonlocal_boundary_audit import OUT, RESULTS, ROOT, load, save, csvsave, integrity, surrogate, field_profiles
from pd_clot.pfpd_candidate.nonlocal_boundary.geometry import BlockBoundary


def f(v): return 'N/A' if v is None else f'{v:.6e}'


def table(headers, rows):
    return '\n'.join(['| '+' | '.join(headers)+' |', '| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(str(c) for c in row)+' |' for row in rows])


def main():
    cfg = load(OUT/'PREDECLARED_CONFIG.json')
    save(OUT/'ENVIRONMENT.json',dict(python=sys.version,executable=sys.executable,
        platform=platform.system(),architecture=platform.machine(),
        packages={name:version(name) for name in ['numpy','scipy','matplotlib','pytest','psutil']},
        thread_environment={name:os.environ.get(name) for name in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']},
        observation='Postprocessing environment; no dependency installation during this task'))
    refs = {k: load(OUT/f'{k.upper()}_REFERENCE.json') for k in ['normal','shear']}
    cases = {k: {c: load(OUT/f'{k}_{c}.json') for c in ['B0','B1','B2']} for k in refs}
    # Derive profile columns from saved states only; no new force evaluation.
    # Retain the original plane-RMS column and add an explicit free-only RMS.
    g=BlockBoundary(np.full(3,-cfg['L_m']/2),np.full(3,cfg['L_m']/2))
    for kind in refs:
        for candidate in cases[kind]:
            with np.load(OUT/f'{kind}_{candidate}_STATE.npz') as a:
                family=SimpleNamespace(X=a['X'],horizon=cfg['delta_m'])
                state=SimpleNamespace(force=a['force_internal'].copy(),Fbar=a['Fbar'])
                profile,_=field_profiles(g,family,a['x'],state,a['external'],a['fixed'],surrogate(kind),a['Psi_plus_nodal'],cfg['h_m'])
                layer=np.rint(g.distance(a['X'])/cfg['h_m']).astype(int)
                error=np.linalg.norm(a['Fbar']-surrogate(kind).F,axis=(1,2))/.003
                distance_profile=[]
                for k in np.unique(layer):
                    mask=layer==k
                    distance_profile.append(dict(depth_h=int(k),depth_m=float(k*cfg['h_m']),depth_over_delta=float(k*cfg['h_m']/cfg['delta_m']),
                        nodes=int(mask.sum()),F_error_RMS=float(np.sqrt(np.mean(error[mask]**2))),F_error_max=float(error[mask].max()),
                        Fzz_mean=float(a['Fbar'][mask,2,2].mean()),Fxz_mean=float(a['Fbar'][mask,0,2].mean()),
                        Psi_plus_nodal_mean_Pa=float(a['Psi_plus_nodal'][mask].mean()),Psi_plus_nodal_max_Pa=float(a['Psi_plus_nodal'][mask].max())))
            csvsave(OUT/f'{kind}_{candidate}_FIELD_PROFILE.csv',profile)
            csvsave(OUT/f'{kind}_{candidate}_DISTANCE_PROFILE.csv',distance_profile)
    work, balance, profiles, targets = [], [], [], []
    for kind in refs:
        for c, row in cases[kind].items():
            q = row['load_checks']
            balance.append(dict(benchmark=kind, candidate=c, load_definition='Applied load plus all constraint reactions',
                resultant_absolute_error_N=q['resultant_absolute_error_N'], resultant_error=q['resultant_error'],
                moment_absolute_error_Nm=q['moment_absolute_error_Nm'], moment_error=q['moment_error'],
                virtual_work_error=q['virtual_work_error'], force_scale_N=q['force_scale_N'], moment_scale_Nm=q['moment_scale_Nm'],
                solver_converged=row['solve']['converged']))
            for key in ['target_resultant_N','discrete_resultant_N','target_moment_Nm','discrete_moment_Nm']:
                for axis,value in zip('xyz',q[key]): balance[-1][key+'_'+axis]=value
            balance[-1]['moment_error_relative_to_nonzero_target']=q['moment_absolute_error_Nm']/np.linalg.norm(q['target_moment_Nm']) if np.linalg.norm(q['target_moment_Nm'])>0 else None
            for r in q['work']: work.append(dict(benchmark=kind,candidate=c,channel='total_external',**r))
        for channel, q in [('combined',refs[kind]['B1_conversion'])]+list(refs[kind]['traction_channels'].items()):
            if 'work' not in q: continue
            for r in q['work']: work.append(dict(benchmark=kind,candidate='B1_frozen',channel=channel,**r))
        q=refs[kind]['B1_conversion']
        balance.append(dict(benchmark=kind,candidate='B1_frozen',load_definition='Frozen converted load without constraint reactions',
            resultant_absolute_error_N=q['resultant_absolute_error_N'],resultant_error=q['resultant_error'],
            moment_absolute_error_Nm=q['moment_absolute_error_Nm'],moment_error=q['moment_error'],
            virtual_work_error=q['virtual_work_error'],force_scale_N=q['force_scale_N'],moment_scale_Nm=q['moment_scale_Nm'],solver_converged='NOT_APPLICABLE'))
        for key in ['target_resultant_N','discrete_resultant_N','target_moment_Nm','discrete_moment_Nm']:
            for axis,value in zip('xyz',q[key]): balance[-1][key+'_'+axis]=value
        balance[-1]['moment_error_relative_to_nonzero_target']=q['moment_absolute_error_Nm']/np.linalg.norm(q['target_moment_Nm']) if np.linalg.norm(q['target_moment_Nm'])>0 else None
        import csv
        with (OUT/f'{kind.upper()}_BOUNDARY_FORCE_SUPPORT.csv').open() as stream:
            profiles.extend(dict(benchmark=kind,**r) for r in csv.DictReader(stream))
        for r in refs[kind]['face_targets']:
            flat = dict(benchmark=kind, face=r['face'], area_m2=r['area_m2'])
            for key in ['traction','traction_normal','traction_tangential','force','normal_force','tangential_force','moment']:
                for axis,value in zip('xyz',r[key]): flat[f'{key}_{axis}']=value
            targets.append(flat)
        save(OUT/f'{kind.upper()}_BENCHMARK.json',dict(reference=refs[kind], candidates=cases[kind],
            acceptance_gates=cfg, scope='Intact affine finite-block mechanics only',
            solver_diagnostic_note='Legacy surface_traction_mean/max fields are unpopulated placeholders for the nodal-force adapter. Use face_targets and saved external nodal loads; zeros in those legacy fields do not imply zero classical shear.'))
    csvsave(OUT/'VIRTUAL_WORK_RESULTS.csv',work); csvsave(OUT/'RESULTANT_BALANCE.csv',balance)
    csvsave(OUT/'BOUNDARY_FORCE_SUPPORT.csv',profiles); csvsave(OUT/'CLASSICAL_FACE_TARGETS.csv',targets)
    status = dict(BOUNDARY_API='IMPLEMENTED',DELIA_BAQP_VOLUME_FORCE='IMPLEMENTED')
    for kind in refs:
        q=refs[kind]['B1_conversion']; solved=cases[kind]['B1']
        for name,key in [('TRACTION_RESULTANT','resultant_PASS'),('TRACTION_MOMENT','moment_PASS'),('VIRTUAL_WORK','virtual_work_PASS')]:
            status[f'{kind.upper()}_{name}']='PASS' if q[key] else 'FAIL'
        pass_all = all(q[key] for key in ['resultant_PASS','moment_PASS','virtual_work_PASS'])
        recovered = (solved['solve']['converged'] and
                     solved['fields']['F_error_max'] <= cfg['field_accuracy_gate_relative_to_imposed_gradient'] and
                     all(solved['load_checks'][key] for key in ['resultant_PASS','moment_PASS','virtual_work_PASS']))
        status[f'{kind.upper()}_BOUNDARY_CANDIDATE']='VERIFIED' if pass_all and recovered and refs[kind]['permutation']['PASS'] else 'NOT_VERIFIED'
    both=all(status[f'{k.upper()}_BOUNDARY_CANDIDATE']=='VERIFIED' for k in refs)
    status.update(VOLUME_DIRICHLET_REFERENCE='IMPLEMENTED', SCABBIA_GHOST_CONTROL='DEFERRED',
        PERIFIT='COMPATIBILITY_REVIEW_ONLY',BAQP_SPECIFIC_SURFACE_CORRECTION='NOT_IMPLEMENTED',
        BOUNDARY_DATA_CONVERSION_STATUS='VERIFIED_FOR_INTACT_AFFINE_MECHANICS' if both else ('NORMAL_ONLY' if status['NORMAL_BOUNDARY_CANDIDATE']=='VERIFIED' else 'NOT_VERIFIED'),
        SURFACE_OPERATOR_CORRECTION_STATUS='SEPARATE_ISSUE',C5_C7_ACTIVATION='NOT_RUN',FRACTURE='NOT_RUN',CLOT_WALL_PROTOTYPE='NOT_RUN')
    save(OUT/'STATUS.json',status)
    resource_rows = {p.name:load(p) for p in sorted(OUT.glob('RESOURCES_*.json'))}
    resources=dict(stages=resource_rows, CPU_workers=1, BLAS_threads=1, GPU='GPU_NOT_USED_BY_DESIGN',
        maximum_concurrent_compute_processes=3,
        concurrency_scope='One full-size scientific solver worker; occasional small tangent timing, focused tests or rendering of saved data overlap, never another memory-heavy 3D solve',
        measured_stage_wall_seconds=sum(r['wall_seconds'] for r in resource_rows.values()),
        peak_RSS_bytes=max(r['peak_RSS_bytes'] for r in resource_rows.values()),
        minimum_available_RAM_fraction=min(r['minimum_available_RAM_fraction'] for r in resource_rows.values()),
        timing_scope='Summed measured compute stages including interrupted attempt; excludes literature reading and authoring')
    save(OUT/'RESOURCE_SUMMARY.json',resources)
    protected=integrity(); save(OUT/'HISTORICAL_INTEGRITY.json',protected); assert protected['PASS']
    initial_adapter=OUT/'references/MECHANICS_ADAPTER_INITIAL.py.txt'
    current_adapter=ROOT/'pd_clot/pfpd_candidate/nonlocal_boundary/mechanics_adapter.py'
    save(OUT/'SOLVER_PROVENANCE.json',dict(
        normal_cases=dict(assembly='NumPy indexed sums',adapter_snapshot=str(initial_adapter.relative_to(ROOT)),sha256=hashlib.sha256(initial_adapter.read_bytes()).hexdigest()),
        shear_cases=dict(assembly='CSR incidence sums',adapter_source=str(current_adapter.relative_to(ROOT)),sha256=hashlib.sha256(current_adapter.read_bytes()).hexdigest()),
        nonlinear_solver='Unchanged regularization.equilibrium._scipy code object; TRF/LSMR',
        gates_and_tolerances_changed=False,
        restart='Normal states saved completely; process interrupted during initial evaluation of shear B0, before any shear state was accepted; all shear cases then run with CSR assembly',
        raw_traceback_note='Intentional interrupt can display updated on-disk source lines against the old in-memory code; archived initial source supplies the original line mapping'))
    metrics = {}
    for kind in refs:
        metrics[kind]=table(['Candidate','Solver force gate','Force error','Moment error','Max work error','Core F RMS','Boundary F RMS','Max F error'],[
            [c,'PASS' if r['solve']['converged'] else 'INCOMPLETE',f(r['load_checks']['resultant_error']),f(r['load_checks']['moment_error']),
             f(r['load_checks']['virtual_work_error']),f(r['fields']['F_error_core_RMS']),f(r['fields']['F_error_boundary_RMS']),f(r['fields']['F_error_max'])]
            for c,r in cases[kind].items()])
    support_table=table(['Load','d50 [µm / h]','d90 [µm / h]','d95 [µm / h]','d99 [µm / h]'],[
        [kind]+[f"{refs[kind]['B1_support_quantiles'][q]['depth_m']*1e6:.2f} / {refs[kind]['B1_support_quantiles'][q]['depth_h']}" for q in ['d50','d90','d95','d99']] for kind in refs])
    field_table=table(['Load / method','Psi+ max [Pa]','Psi+ P95 [Pa]','Energy [J]','Top reaction x [N]','Top reaction z [N]','First interior jump norm / affine'],[
        [kind+' / '+c,f(r['fields']['Psi_plus_max_Pa']),f(r['fields']['Psi_plus_P95_Pa']),f(r['fields']['Psi_total_J']),f(r['reaction_top_N'][0]),f(r['reaction_top_N'][2]),
         f(float(np.linalg.norm(r['first_interior_jumps'][0]['displacement_jump_m'])/np.linalg.norm(r['first_interior_jumps'][0]['surrogate_jump_m'])))]
        for kind in refs for c,r in cases[kind].items()])
    channel_table=table(['Shear channel','Force error','Moment error','Max work error'],[
        [name,f(q['resultant_error']),f(q['moment_error']),f(q['virtual_work_error'])] for name,q in refs['shear']['traction_channels'].items()])
    refinement=[]
    if (OUT/'REFINEMENT.json').exists():
        refinement=load(OUT/'REFINEMENT.json')
    refinement_note='Both main B1 cases passed before these four N21 evaluations.' if refinement else 'Refinement NOT_RUN: main B1 acceptance was not established; only N19 rows are shown below.'
    refine_table=table(['Resolution / load','h [µm]','delta [µm]','Force error','Moment error','Work error','d99 / delta'],[
        [r['geometry']['name']+' / '+r['kind'], f(r['geometry']['h_m']*1e6), f(r['geometry']['delta_m']*1e6),
         f(r['B1_conversion']['resultant_error']),f(r['B1_conversion']['moment_error']),f(r['B1_conversion']['virtual_work_error']),
         f(r['B1_support_quantiles']['d99']['depth_over_delta'])] for r in list(refs.values())+refinement])
    run_table=table(['Case','BAQP evaluations','Tangent actions','Wall [s]','Termination'],[
        [kind+' / '+c,r['solve']['BAQP_evaluations'],r['solve']['tangent_actions'],f(r['solve']['wall_seconds']),r['solve']['status']]
        for kind in refs for c,r in cases[kind].items()])
    normal0, normal1 = cases['normal']['B0'], cases['normal']['B1']
    field_interpretation = ('Both B0 and B1 meet the mechanical force gate. The normal core gradient RMS changes from '+
        f(normal0['fields']['F_error_core_RMS'])+' to '+f(normal1['fields']['F_error_core_RMS'])+
        '; boundary RMS changes from '+f(normal0['fields']['F_error_boundary_RMS'])+' to '+
        f(normal1['fields']['F_error_boundary_RMS'])+'. The positive-energy maximum relative to the affine target changes from '+
        f(normal0['fields']['Psi_plus_max_Pa']/refs['normal']['response']['W'])+' to '+
        f(normal1['fields']['Psi_plus_max_Pa']/refs['normal']['response']['W'])+'.') if all(r['solve']['converged'] for r in [normal0,normal1]) else (
        'At least one normal comparison has not reached the force gate; its endpoint fields are diagnostic only, and a converged B0-to-B1 field-quality comparison is not established.')
    tests=(OUT/'TARGETED_TESTS.log').read_text().strip()
    start=load(OUT/'STARTING_GIT_STATE.json')
    text=f'''# NONLOCAL_BOUNDARY_FORMULATION_V1 — review

Date: 2026-10-05. Model: `PFPD_AT1_SPECTRAL_CANDIDATE_V1`. **Intact mechanics only.**

Decision: **{status['BOUNDARY_DATA_CONVERSION_STATUS']}**. Scope is the two analytic affine block benchmarks and the exact discretization documented below. This is not validation of general nonlinear clot, fracture, adhesion or CFD boundary conditions. The construction identity is not the acceptance test. B0/B2 virtual-work failures are retained below.

## 1. Why bulk BAQP development is no longer the current blocker

The preceding affine/coefficient audits support the existing variational BAQP operator and explain the small raw Z mismatch through reconstructed-stress floating-point propagation. Their original raw FAIL records remain unchanged. This round changes boundary data, not bulk mechanics. {protected['protected_files']} protected input/source files retain their SHA256 values. See [integrity record](HISTORICAL_INTEGRITY.json); historical C5/C7 scientific states are neither changed nor rerun.

## 2. What the four papers actually contribute

The complete source map, equation locations, access limitations and distinctions are in [LITERATURE_TO_BAQP_MAPPING.md](LITERATURE_TO_BAQP_MAPPING.md). D'Elia supplies the conversion architecture; Yu supplies a traction-consistency reference; Scabbia informs an independent ghost-control design; Shojaei/perifit concerns a separate operator correction. Public full texts were read before coding because the uploaded PDF bytes were unavailable locally. The local PDF provenance does not certify identity with those uploads.

## 3. What was deliberately not copied

No scalar diffusion theorem is asserted for BAQP. No Yu LPS traction formula, OSB dilatation/weighted-volume ghost formula, or perifit weight is inserted. No kernel, K/F/Z coefficient rule, material model, phase law or historical assertion is changed. No arbitrary 3-delta physical collar is introduced. [Compatibility notes](PERFIT_BAQP_COMPATIBILITY_NOTES.md) identify future energy/variation obligations.

## 4. Boundary API design

The independent namespace `pd_clot/pfpd_candidate/nonlocal_boundary/` provides Cartesian geometry, face partitions, normals/tangents, distance/projection, physical/interior/exterior interaction masks, analytic local solutions, frozen volume data and diagnostics. [BOUNDARY_SPECIFICATION.md](BOUNDARY_SPECIFICATION.md) gives equations and node roles. Main geometry: N19, 6,859 physical particles, 690,020 bonds, L=1.125 mm, h=62.5 µm, delta=251.25 µm, delta/h=4.02. SVK is explicit: mu=1000 Pa, lambda=0, K=666.6666667 Pa.

Pre-run peak estimates are 1.41 GiB for N19, 2.46 GiB for fixed-delta N21, and 1.85 GiB for fixed-ratio N21. The estimates and all gates preceded results in [PREDECLARED_CONFIG.json](PREDECLARED_CONFIG.json). The geometric core has 729 nodes with d>=delta; it is not a completely closed bulk stencil.

## 5. D'Elia-style BAQP volume-force conversion

B1 applies the existing internal-force evaluator to u_L and freezes its negative on a topology-derived support. Nodal N and volumetric N/m³ units are explicit. Load updates during equilibrium are forbidden. The complete support covers this N19 domain; omitted-support force is {f(refs['normal']['out_of_support_residual_N'])} N in normal loading and {f(refs['shear']['out_of_support_residual_N'])} N in shear. The exact cancellation identity is retained as a check, not scientific proof. Analytic continuum targets, all affine work modes, relabeling, perturbed-start recovery and limited refinement provide additional evidence.

This finite-domain in-domain force representation differs from the original exterior nonlocal-flux problem. It is a restricted surrogate-dependent construction, not a theorem that any prescribed traction admits a unique such representation.

## 6. Normal traction benchmark

F=diag(1,1,1.003); Ezz=0.0030045, Szz=6.009 Pa, Pzz=6.027027 Pa, W=0.00902702025 Pa. Every face uses PN. Each z-face resultant has magnitude {f(np.linalg.norm(refs['normal']['face_targets'][-1]['force']))} N; closed-surface target resultant is zero. E/S/P/W and each exact face load are saved in [NORMAL_BENCHMARK.json](NORMAL_BENCHMARK.json) and [CLASSICAL_FACE_TARGETS.csv](CLASSICAL_FACE_TARGETS.csv).

{metrics['normal']}

F errors are relative to the imposed gradient norm 0.003. Force and moment use stable stress-area and stress-volume scales; nonzero work uses its own target. All three load columns include applied forces plus constraint reactions. Frozen B1 load checks are reported separately in section 9.

## 7. Volume-Dirichlet reference

B2 is **LOCALITY-STRONG REFERENCE**: all physical nodes d<delta on all six faces are prescribed to the local surrogate; no external nodes, ghost energy or extra material are introduced. There are {cases['normal']['B2']['constraint_nodes']} constrained physical particles and {6859-cases['normal']['B2']['constraint_nodes']} free particles. Constraint matching is exact to {f(cases['normal']['B2']['constrained_displacement_error_m'])} m. This is an in-domain support adaptation, not a demonstrated equivalent exterior-domain boundary-value problem. Small constrained-region error alone would not justify selecting it. One-horizon constraints need not remove every multi-hop influence in the free core.

The normal and shear B2 work errors are {f(cases['normal']['B2']['load_checks']['virtual_work_error'])} and {f(cases['shear']['B2']['load_checks']['virtual_work_error'])}; both exceed 1e-6 despite mechanical convergence. The prescribed support was not widened after seeing these outcomes.

## 8. Scabbia control if implemented

`SCABBIA_GHOST_CONTROL=DEFERRED` (not implemented in V1). Taylor continuation defines affine ghost positions, but does not define ghost-side BAQP K/F/Z, the physical-only quadrature energy, or the pullback of reactions when fictitious positions depend on physical DOFs. Treating those nodes as ordinary extra material would change the problem. No ghost run or unverifiable independent-control claim is made.

## 9. Normal resultant / moment / virtual work

For the frozen B1 load: force error {f(refs['normal']['B1_conversion']['resultant_error'])}, moment error {f(refs['normal']['B1_conversion']['moment_error'])}, maximum work error {f(refs['normal']['B1_conversion']['virtual_work_error'])}. The predeclared gates are 1e-6. All nine affine A basis modes and three translations are saved in [VIRTUAL_WORK_RESULTS.csv](VIRTUAL_WORK_RESULTS.csv). Zero closed-surface force alone is not sufficient: a deliberately zero load fails nonzero-work unit tests.

Normal load relabeling changes the force by {f(refs['normal']['permutation']['relative_force_L2_error'])} relative L2, below 1e-10. The upper-half volume-load versus top-face target discrepancy is {f(refs['normal']['half_space_vs_face_diagnostic'][-1]['stable_scaled_difference'])} on the stress-area scale. That latter check is diagnostic; no unique face assignment is assumed.

## 10. Shear traction benchmark

Shear was run after normal force/work gates passed. F=I+0.003 ex tensor ez. Nonzero P components [Pa] are Pxx=Pzz=0.009, Pxz=3.000027, Pzx=3. The exact energy is 0.00450002025 Pa. Lateral x-face traction is included. All tensors, targets and control states are in [SHEAR_BENCHMARK.json](SHEAR_BENCHMARK.json).

{metrics['shear']}

## 11. Shear resultant / moment / virtual work

Frozen B1 errors: force {f(refs['shear']['B1_conversion']['resultant_error'])}, reference moment {f(refs['shear']['B1_conversion']['moment_error'])}, work {f(refs['shear']['B1_conversion']['virtual_work_error'])}. The reference-coordinate target My is {f(refs['shear']['B1_conversion']['target_moment_Nm'][1])} N m, not zero. Relabeling force error is {f(refs['shear']['permutation']['relative_force_L2_error'])}.

Moment gates use the stable stress-volume scale specified before the run. The balance CSV also gives error relative to the small nonzero moment itself: shear B0 differs by {cases['shear']['B0']['load_checks']['moment_absolute_error_Nm']/np.linalg.norm(refs['shear']['B1_conversion']['target_moment_Nm']):.2%} of that target, despite passing the stable moment gate. The separate virtual-work and known-field gates remain necessary.

{channel_table}

These channels use constant Piola-stress diagonal/off-diagonal components in the unchanged force assembly. They verify Cartesian normal/tangential accounting at frozen geometry; they are not a nonlinear constitutive superposition theorem. Face-wise normal/tangent vectors are retained separately in the target CSV.

## 12. Spatial support of the equivalent nonlocal load

{support_table}

The distribution is weighted by summed nodal force magnitudes. Direct-interaction nodes number {refs['normal']['direct_horizon_nodes']}; the graph support contains all 6,859 nodes. d99 near one horizon does not imply a compact one-horizon physical BC. A small interior tail is retained. [BOUNDARY_FORCE_SUPPORT.csv](BOUNDARY_FORCE_SUPPORT.csv) records both d/h and d/delta.

The final depth bin contains only the center particle. Its near-zero force is consistent with inversion symmetry and must not be used to assert three-hop stencil closure. The favorable N19 half-space diagnostics do not retroactively validate the much smaller historical C5/C7 specimens.

Small predeclared refinement: fixed L, N19 versus N21, with fixed delta or fixed delta/h. Only affine conversion/work/support is reevaluated at N21; these are not additional equilibrium convergence solves.

{refinement_note}

{refine_table}

Two resolutions cannot establish asymptotic compatibility, mesh objectivity, or fracture convergence. Grid-quantized d99 can change even while the force distribution is close.

## 13. Does the new method reduce the old artificial boundary field?

{field_table}

{field_interpretation}

The first-jump ratio uses the centerline first interior layer relative to its affine displacement increment. Top constraint reaction is separate from already-applied B1 volume force; a nearly zero B1 grip reaction does not mean absent physical loading. Full vectors, reactions and free-node residual profiles are retained in candidate JSON/CSV/NPZ files. These Psi+ changes are intact mechanical diagnostics, not fracture corrections.

Each `*_DISTANCE_PROFILE.csv` records gradient RMS/maximum and nodal energy against the nearest of all six physical faces. `*_FIELD_PROFILE.csv` separately records profiles from the z-grip planes; its original plane-RMS residual column is retained and an explicit free-only RMS/count is supplied. Layers without free nodes have an empty free-only RMS, not an artificially perfect equilibrium metric.

Known-solution recovery is accepted only when the existing solver reaches its 1e-8 force gate and the maximum reconstructed F error is <=1e-5 of the imposed gradient. The final statuses below enforce that additional criterion. All six equilibrium solves meet the force gate, supporting interpretation of the large B0/B2 work and field discrepancies as differences in their discretized boundary treatments at this grid.

{run_table}

The adapter runs the existing SciPy solver code, supplying an analytic BAQP tangent action; no new nonlinear iteration is written. Its derivative is checked against force finite differences and reciprocal bilinear products. Matrix-free condition number is explicitly unavailable. A short initial attempt was interrupted to fix the reporting adapter's float conversion; it accepted no state and its raw log/resources remain saved. After the normal cases, an equivalent CSR bond-to-node aggregation was enabled for subsequent cases. [TANGENT_SPARSE_ASSEMBLY_CHECK.json](TANGENT_SPARSE_ASSEMBLY_CHECK.json) compares its action to the archived initial adapter, with relative difference 3.29e-16 on N9. No tolerance, force law or nonlinear algorithm changed; the initial adapter source is retained under `references/`.

## 14. What remains surface/operator correction rather than BC conversion?

`SURFACE_OPERATOR_CORRECTION_STATUS=SEPARATE_ISSUE`. Affine conversion under the unchanged variational operator is the demonstrated object. It does not remove all effects of truncated horizons, prove a boundary law for arbitrary nonaffine fields, or establish any need for a weighted BAQP operator. A manufactured affine load can preserve global work even when local force support is distributed. Future evidence must separate external data conversion from operator changes.

## 15. perifit compatibility assessment

[PERFIT_BAQP_COMPATIBILITY_NOTES.md](PERFIT_BAQP_COMPATIBILITY_NOTES.md) describes candidate BAQP targets: shape invertibility/completeness, homogeneous energy, a common force first variation, and Z-mediated contributions. Symmetric scalar weights alone do not prove all tensor and energy conditions. No fit or kernel modification is performed.

## 16. Implication for future clot-wall mechanics

The API accepts geometric displacement regions and pressure-plus-tangential traction data. Clot-wall adhesion would require an independent interface model. The present Cartesian affine surrogate cannot represent a real clot's curved, heterogeneous or nonaffine solution. Mixed loading, physical clot-wall prototypes, real-clot material calibration, activation, stability, fracture, erosion, microbubble loading and CFD remain NOT_RUN.

## 17. One next development task

After human review only: **NONLOCAL_BOUNDARY_NONAFFINE_LOCALIZATION_V1** — test a known nonaffine elastic field with compatible body force and surface data, and determine which part of the finite-domain residual is genuinely boundary data. Require spatially resolved traction/work and interior-load accounting before claiming a transferable nonlocal Neumann law. Do not start this task automatically.

## Final status and audit evidence

{table(['Status','Value'],status.items())}

Focused tests only:

```text
{tests}
```

An initial test run had two new-test API-name typos; these were fixed without changing scientific assertions, and its log is preserved. No full pytest suite was run. Historical missing artifacts and old failures were not modified.

Resources: one CPU worker, one BLAS thread, GPU_NOT_USED_BY_DESIGN. Summed measured stage wall time {resources['measured_stage_wall_seconds']:.2f} s, including the interrupted attempt; peak recorded RSS {resources['peak_RSS_bytes']/2**30:.3f} GiB; minimum available system RAM {resources['minimum_available_RAM_fraction']:.1%}. Authoring/literature time is outside this compute timing. [RESOURCE_SUMMARY.json](RESOURCE_SUMMARY.json) preserves per-stage measurements. No memory-heavy solves overlap.

Figures: [results/nonlocal_boundary_formulation_v1/figures](../../results/nonlocal_boundary_formulation_v1/figures), English labels, 600 dpi PNG and vector PDF, generated from saved numerical data by `scripts/plot_nonlocal_boundary.py`.

The observed Python/package versions are in [ENVIRONMENT.json](ENVIRONMENT.json); no dependencies were installed. [SOLVER_PROVENANCE.json](SOLVER_PROVENANCE.json) identifies the adapter version used for each loading family.

**Stopped for human review. No subsequent scientific stage is authorized by this report.**
'''
    (OUT/'REVIEW_NONLOCAL_BOUNDARY_FORMULATION_V1.md').write_text(text)
    figure_names=['boundary_formulation_schematic','normal_traction_candidate_comparison','equivalent_boundary_force_profile','normal_affine_field_profiles','shear_traction_consistency']
    figure_index=table(['Figure','PNG 600 dpi','Vector PDF'],[[name,f'[PNG](figures/{name}.png)',f'[PDF](figures/{name}.pdf)'] for name in figure_names])
    (RESULTS/'README.md').write_text('# NONLOCAL_BOUNDARY_FORMULATION_V1\n\n'+
        'Review: [REVIEW_NONLOCAL_BOUNDARY_FORMULATION_V1.md](../../verification/nonlocal_boundary_formulation_v1/REVIEW_NONLOCAL_BOUNDARY_FORMULATION_V1.md).\n\n'+
        'Five figure pairs in `figures/` are generated from the saved benchmark JSON/NPZ/CSV files. Status: `'+status['BOUNDARY_DATA_CONVERSION_STATUS']+'`.\n\n'+figure_index+'\n\n'+
        'The schematic uses the actual central particle slice and one visual arrow scale. Force-support quantiles describe this conversion, not a physical collar thickness. The geometric core is not a closed bulk stencil. No fracture or CFD result is shown.\n')
    paths=list((ROOT/'pd_clot/pfpd_candidate/nonlocal_boundary').glob('*.py'))+[ROOT/'scripts/nonlocal_boundary_audit.py',ROOT/'scripts/plot_nonlocal_boundary.py',ROOT/'scripts/report_nonlocal_boundary.py',ROOT/'tests/test_nonlocal_boundary_formulation_v1.py']
    paths += [p for p in OUT.rglob('*') if p.is_file() and p.name != 'OUTPUT_MANIFEST.json']
    paths += [p for p in RESULTS.rglob('*') if p.is_file()]
    save(OUT/'OUTPUT_MANIFEST.json',{str(p.relative_to(ROOT)):dict(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sorted(paths)})
    print(json.dumps(status,indent=2))


if __name__ == '__main__': main()
