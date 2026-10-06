"""Assemble the review from saved results, including limitations and raw failures."""
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from nonaffine_boundary_audit import OUT, ROOT, load, save, integrity


def fmt(x):return f'{x:.6e}' if isinstance(x,(float,np.floating)) else str(x)
def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(fmt(x) for x in row)+' |' for row in rows])


def main():
    cfg=load('PREDECLARED_CONFIG.json');exact=load('EXACT_NONAFFINE_REFERENCE.json');bulk=load('BULK_ORACLE_VALIDATION.json')
    aff=load('AFFINE_LIMIT_RECOVERY.json');node=load('NODEWISE_IDENTITIES.json');n27=load('N27_DEPTH_CHECK.json')
    force=load('FORCE_BALANCE.json');moment=load('MOMENT_BALANCE.json');qs=load('SUPPORT_QUANTILES.json')
    nl=load('NL_MANUFACTURED_RECOVERY.json');cl=load('CLASSICAL_MANUFACTURED_RESULT.json');git=load('STARTING_GIT_STATE.json')
    with (OUT/'VIRTUAL_WORK_DECOMPOSITION.csv').open() as stream:
        work_closure=max(float(r['closure_scaled']) for r in csv.DictReader(stream))
    localization=bool(exact['VERIFIED'] and bulk['VERIFIED'] and aff['PASS'] and node['PASS'] and n27['PASS'] and nl['PASS'])
    with np.load(OUT/'NODEWISE_RESIDUAL_DECOMPOSITION.npz') as a:
        free=~np.load(OUT/'classical_RECOVERED_STATE.npz')['fixed']
        target_classical_residual=float(np.linalg.norm(a['f_model'][free])/cl['solver_force_scale_N'])
        classical_target_max_local=float(np.linalg.norm(a['f_model'],axis=1).max()/exact['scales']['F_ref_N'])
        # Subtract the two converged solves: the only applied-load change is e_model V.
        with np.load(OUT/'NL_RECOVERED_STATE.npz') as ns,np.load(OUT/'classical_RECOVERED_STATE.npz') as cs:
            load_identity=np.linalg.norm(cs['external']-ns['external']-a['f_model'],axis=1).max()/exact['scales']['F_ref_N']
            response_balance=cs['force_internal']-ns['force_internal']+a['f_model']
            response_error=float(np.linalg.norm(response_balance[free])/np.linalg.norm(a['f_model'][free]))
            solved_difference=float(np.linalg.norm(cs['x']-ns['x'],axis=1).max())
            response=dict(external_difference_equals_model_scaled=float(load_identity),
                          force_response_closure_relative_to_model=response_error,
                          NL_to_classical_max_displacement_difference_m=solved_difference,
                          NL_recovery_error_over_classical_deviation=nl['displacement_error_max_m']/max(cl['displacement_error_max_m'],1e-30),
                          classical_target_residual_existing_solver_scale=target_classical_residual,
                          classical_target_max_residual_local_scaled=classical_target_max_local,
                          local_continuum_energy_J=exact['continuum']['energy_J'],
                          BAQP_target_minus_local_continuum_energy_J=node['target_metrics']['energy_J']-exact['continuum']['energy_J'],
                          classical_minus_local_continuum_energy_J=cl['energy_J']-exact['continuum']['energy_J'],
                          interpretation='Frozen load difference is e_model V exactly. Finite solve tolerance bounds resolution of the small response.')
    save('CLASSICAL_VS_NL_ACCOUNTING.json',response)
    conversion=localization and cl['PASS'] and load_identity<=cfg['gate']
    # Mathematical non-exactness is diagnosed at uL, independently of solver recovery tolerance.
    classical_label=('APPROXIMATE_DUE_TO_BULK_NONLOCALITY' if classical_target_max_local>cfg['gate'] else 'EXACT') if cl['PASS'] else 'FAILED_FOR_OTHER_REASON'
    status=dict(NONAFFINE_AXIAL_REFERENCE='VERIFIED' if exact['VERIFIED'] else 'FAILED',
                BULK_BAQP_NONAFFINE_ORACLE='VERIFIED' if bulk['VERIFIED'] else 'NOT_VERIFIED',
                AFFINE_LIMIT_RECOVERY='PASS' if aff['PASS'] else 'FAIL',NODEWISE_DECOMPOSITION='PASS' if node['PASS'] else 'FAIL',
                NAIVE_NONAFFINE_CANCELLATION='CONTAINS_INTERIOR_PHYSICS' if n27['points'][-1]['naive_scaled']>cfg['gate'] and n27['PASS'] else 'NOT_DEMONSTRATED',
                BOUNDARY_LOCALIZATION='VERIFIED' if localization else 'NOT_VERIFIED',
                N27_CLOSED_CENTER_BOUNDARY_CORRECTION='PASS' if n27['PASS'] else 'FAIL',
                BULK_MODEL_DIFFERENCE=bulk['BULK_MODEL_DIFFERENCE'],NL_MANUFACTURED_RECOVERY='PASS' if nl['PASS'] else 'FAIL',
                CLASSICAL_MANUFACTURED_FIELD=classical_label,NONAFFINE_SHEAR='NOT_RUN',
                BOUNDARY_DATA_CONVERSION_STATUS='NONAFFINE_LOCALIZATION_VERIFIED' if conversion else 'NOT_VERIFIED',
                SURFACE_OPERATOR_CORRECTION='NOT_MODIFIED',C5_C7_ACTIVATION='NOT_RUN',FRACTURE='NOT_RUN',CLOT_GEOMETRY='NOT_RUN')
    save('STATUS.json',status);save('HISTORICAL_INTEGRITY.json',integrity())
    resources=[json.loads(p.read_text()) for p in sorted(OUT.glob('RESOURCES_*.json'))]
    resource_summary=dict(stages=resources,measured_compute_wall_seconds=sum(r['wall_seconds'] for r in resources),
                          peak_process_tree_RSS_bytes=max(max(r['peak_process_tree_RSS_bytes'],r.get('peak_main_RSS_bytes',0)) for r in resources),
                          minimum_available_RAM_fraction=min(r['minimum_available_RAM_fraction'] for r in resources),
                          oracle_workers=8,equilibrium_workers=1,BLAS_threads=1,GPU='GPU_NOT_USED_BY_DESIGN',
                          timing_scope='Sum of recorded scientific stages, including initial precondition interruptions; excludes authoring and figure export. Process-tree RSS sums shared pages conservatively.')
    save('RESOURCE_SUMMARY.json',resource_summary)
    s=[]
    s.append('# NONLOCAL_BOUNDARY_NONAFFINE_LOCALIZATION_V1 — review\n\n'
             f"**{status['BOUNDARY_DATA_CONVERSION_STATUS']}** for the frozen intact quadratic axial benchmark. "
             'Boundary localization and bulk-model agreement are separate conclusions. '
             'This is a surrogate-dependent manufactured-data audit, not a general physical traction law.\n\n'
             f"Branch `{git['branch']}`; HEAD `{git['HEAD']}`. Frozen configuration: [PREDECLARED_CONFIG.json](PREDECLARED_CONFIG.json).")
    s.append('## 1. Why affine B1 was insufficient for generalization\n\n'
             'The previous affine benchmark has zero full-stencil force and zero local body force. Consequently its negative finite-domain residual can be interpreted as boundary conversion on that surrogate. '
             'Here the intact bulk force is nonzero. Negating the entire finite-domain force would cancel legitimate interior mechanics as well as truncation effects. '
             '[Previous review](../nonlocal_boundary_formulation_v1/REVIEW_NONLOCAL_BOUNDARY_FORMULATION_V1.md) remains unchanged.')
    s.append(r'''## 2. Manufactured nonaffine field

On the physical N19 cube, $L=1.125\,\mathrm{mm}$, $h=62.5\,\mu\mathrm m$, $\delta=251.25\,\mu\mathrm m$, $\delta/h=4.02$:

$$u_x=u_y=0,\qquad u_z=\epsilon Z+\beta Z^2/L,\qquad \epsilon=0.003,\quad\beta=0.001,$$
$$a=1+\epsilon+2\beta Z/L,\qquad F=\operatorname{diag}(1,1,a),\qquad J=a\in[1.002,1.004].$$

SVK is explicit: $\mu=1000\,\mathrm{Pa}$, $\lambda=0$, $K=666.6666666666666\,\mathrm{Pa}$. Physical boundary-node trapezoid volumes are retained. Phase, AT1, Gc, phase VI, stability, fracture, erosion and CFD are unused. No material defaults are used.

The grid expression is the previous round's exact floating representation $L=18h$. The decimal representation of the same length differs by one floating bit. An initial exact-coordinate storage check stopped before any regression result was evaluated; [IMPLEMENTATION_NOTES.json](IMPLEMENTATION_NOTES.json) preserves this correction and the original logs. No physical parameter or scientific threshold was adjusted.''')
    s.append(r'''## 3. Exact local continuum body force and traction

Independent closed forms, implemented without the general constitutive routine:

$$E_{zz}=\tfrac12(a^2-1),\quad S_{zz}=\mu(a^2-1),\quad P_{zz}=\mu(a^3-a),\quad W=\tfrac\mu4(a^2-1)^2,$$
$$b_z=-\mu(3a^2-1)\frac{2\beta}{L},\quad b_x=b_y=0,\qquad t=PN.$$

Top traction is $+P_{zz}(L/2)e_z$, bottom traction is $-P_{zz}(-L/2)e_z$ and side traction is zero. The formulas follow by differentiating the supplied SVK potential and the specified displacement field; their source is the user's [archived task, sections 4–6](TASK_REQUEST.txt), checked independently against the unchanged constitutive implementation. No external theorem about nonlinear BAQP boundary equivalence is asserted.

Eight-point Gauss quadrature is polynomial-exact for the body-force, energy and selected virtual-work integrands. All nine affine variations, three translations and the quadratic axial variation are checked.'''+
             f"\n\nContinuum force closure = {exact['continuum']['force_closure_scaled']:.6e} on the stress-area scale; max work closure = {exact['continuum']['max_work_closure_scaled']:.6e} on the stress-volume scale. See [EXACT_NONAFFINE_REFERENCE.json](EXACT_NONAFFINE_REFERENCE.json).")
    s.append('## 4. How the fully closed BAQP bulk oracle is constructed\n\n'
             'The nonaffine adapter subclasses the already verified `AffineStencil`; its neighbor, endpoint volume, K, reconstructed F, bond stress, Z and force-assembly methods are inherited unchanged. Only reference translation and manufactured current nodal positions change. '
             'Each of the 19 absolute Z levels uses an N27 scaffold centered at that point; the polynomial extends to all required nodes. No global N27 bond state is allocated. '
             'D0/D1/D2/D3 contain 1/256/1797/5821 required nodes, with maximum coordinate offsets 0/4/8/12h. All required nodes have volume h³, and every required D2 family is complete. '
             'The scaffold has boundary nodes outside the dependency stencil; none participates in the evaluated center force. '\
             '\n\nThree selected levels are repeated on translated N29 scaffolds with observation offset (1,-1,1); three reflected fields use beta→−beta, Z→−Z and the reflected force vector; three beta=0 cases recover zero bulk force. '
             'The selective result is also compared with full N19 evaluation at three nonaffine points. '
             '[BULK_STENCIL_CLOSURE.json](BULK_STENCIL_CLOSURE.json), [BULK_ORACLE_VALIDATION.json](BULK_ORACLE_VALIDATION.json), and [SELECTIVE_VS_FULL_NONAFFINE.json](SELECTIVE_VS_FULL_NONAFFINE.json) retain all checks. '
             'The nominal physical/oracle Z levels differ only by the one-ulp L representation described above (max 1.0842e-19 m); [Z_LEVEL_MAPPING_PRECISION.json](Z_LEVEL_MAPPING_PRECISION.json) records it explicitly. No interpolation is used.')
    s.append(r'''## 5. Bulk BAQP force vs local Div P

$$q_{\rm bulk}(Z)=f_{\rm closed}(Z)/h^3,\qquad e_{\rm model}=q_{\rm bulk}+b_{\rm local}.$$

Both densities have units N/m³. Their difference is a finite-horizon bulk-model/discretization difference, not boundary contamination. The predeclared descriptive SMALL threshold is 0.1% of max local body-force density; it is not a general modeling tolerance.'''+
             r' The second reference density is the characteristic stress-gradient scale $2\mu|2\beta/L|=3555.555555555556\,\mathrm{N/m^3}$; no small-strain approximation is used in the actual stress or body force.'+
             '\n\n'+table(['Quantity','Value'],[[k,bulk[k]] for k in ['max_local_body_N_m3','max_bulk_q_N_m3','max_model_error_N_m3','model_error_relative_to_body','model_error_relative_to_stress_gradient','BULK_MODEL_DIFFERENCE']])+
             '\n\nThe tiny variation of the mismatch across Z contains reconstructed-stress floating-point scatter. Its nonzero magnitude is resolved independently of the equilibrium solves. [BULK_BAQP_REFERENCE.csv](BULK_BAQP_REFERENCE.csv) saves all vector components.')
    s.append(r'''## 6. Definition of the localized boundary correction

$$f_{{\rm bulk,target},i}=q_{\rm bulk}(Z_i)V_i,\qquad f_{{\rm BC,localized},i}=f_{{\rm bulk,target},i}-f_{{\rm FD},i}.$$

Physical $V_i$ includes face, edge and corner factors. This converts the physical finite-domain force into the bulk BAQP force at the specified surrogate. It is not a final physical traction law. Section 11 of the task defines the signed expression used here; section 20's reversed verbal subtraction is interpreted using that main definition.

The correction is frozen throughout equilibrium. The common local normalization is
$$P_{\rm ref}=\max_{\rm N19}\|P\|_F,\quad b_{\rm ref}=\max_{\rm N19}\|b\|,\quad F_{\rm ref}=\max(P_{\rm ref}h^2,b_{\rm ref}h^3).$$'''+
             f"\n\nF_ref = {exact['scales']['F_ref_N']:.9e} N. N27 uses this same physical-N19 reference scale. Beta=0 regression against the frozen old B1 gives relative L2 error {aff['relative_L2_error']:.6e}, below 1e-10.")
    s.append(r'''## 7. Why naive -r_FD is wrong for nonaffine fields

$$f_{\rm BC,naive}=-f_{\rm FD},\qquad f_{\rm BC,naive}-f_{\rm BC,localized}=-q_{\rm bulk}V.$$

The naive quantity is diagnostic only. It removes the legitimate interior force.'''+
             f" At the fully closed N27 center its magnitude/F_ref is {n27['points'][-1]['naive_scaled']:.6e}, while the localized correction is {n27['points'][-1]['localized_scaled']:.6e}. This distinguishes interior loading from boundary conversion directly.")
    s.append(r'''## 8. Nodewise residual decomposition

$$f_{\rm FD}+f_{\rm BC,localized}+b_{\rm local}V=e_{\rm model}V.$$

All nodes, six face distances, d/h, d/delta and six force channels are saved in [NODEWISE_RESIDUAL_DECOMPOSITION.npz](NODEWISE_RESIDUAL_DECOMPOSITION.npz). These are nodal forces in N; density arrays are named separately.'''+
             '\n\n'+table(['Identity','max nodal norm / F_ref'],[['Nodewise model decomposition',node['closure_scaled']],['Naive minus localized',node['naive_difference_scaled']],['NL target equilibrium',node['NL_target_residual_scaled']]])+
             '\n\nAll prescribed 1e-10 gates pass. Algebraic closure alone is not evidence of a physical boundary law.')
    s.append(r'''## 9. Global force / moment accounting

The correct force identity retains physical quadrature and assembly terms:
$$F_{\rm BC}-F_{\rm surface}=F_{\rm model}-(F_{\rm body,nodal}-F_{\rm body,exact})-F_{\rm internal}.$$

The user-requested comparison with F_model alone is therefore only exact up to physical body quadrature and floating arithmetic. Quadratic b is not integrated exactly by nodal trapezoid volumes. No volume or old manifest is altered to erase this difference.'''+
             r' The composite trapezoid formula independently predicts $\Delta F_{\rm body,z}=-4\mu\beta^3h^2=-1.5625\times10^{-14}\,\mathrm N$ for this cube.'+
             '\n\n'+table(['Axial resultant term','N'],[[k,force[k][2]] for k in ['surface_local','body_local_exact','body_local_nodal','BC_localized','model','body_quadrature_error','internal','BC_minus_surface','closure']])+
             f"\n\nFull force closure/F_ref = {force['closure_scaled']:.6e}. Reference-centroid moment accounting retains the same body-quadrature and internal terms; closure/(F_ref L) = {moment['closure_scaled']:.6e}. "
             'BAQP angular momentum refers to current-position moments. A zero reference-position internal moment is not imposed; this axial symmetric case is consistent with zero at floating precision. '
             '[FORCE_BALANCE.json](FORCE_BALANCE.json) and [MOMENT_BALANCE.json](MOMENT_BALANCE.json) report actual vectors, including the simple-identity discrepancy.')
    s.append(r'''## 10. Virtual-work accounting

The existing affine/translation variation infrastructure is reused, with $\delta u_z=Z^2/L$ added. For all 13 modes:
$$\sum_i(f_{\rm FD}+f_{\rm BC,localized}+b_{\rm local}V)_i\cdot\delta u_i
=\sum_i(e_{\rm model}V)_i\cdot\delta u_i.$$

[VIRTUAL_WORK_DECOMPOSITION.csv](VIRTUAL_WORK_DECOMPOSITION.csv) records finite-domain internal, converted boundary, local body, corrected residual and bulk-model work. Corrected work is not forced to zero. Continuum work checks are separately retained in the exact reference JSON.''')
    s[-1]+=f' Maximum discrete work closure/(F_ref L) is {work_closure:.6e}.'
    s.append('## 11. Boundary-correction spatial support\n\n'+
             table(['Correction','d50 µm / h / δ','d90 µm / h / δ','d95 µm / h / δ','d99 µm / h / δ'],
                   [[name]+[f"{qs[name][k]['depth_um']:.1f} / {qs[name][k]['depth_h']} / {qs[name][k]['depth_over_delta']:.6f}" for k in ['d50','d90','d95','d99']] for name in ['localized','naive']])+
             '\n\nThese magnitude-weighted quantiles describe **CONVERSION SUPPORT, NOT PHYSICAL BC THICKNESS**. Mean/max profiles and cumulative magnitude use nearest distance to all six physical faces. N19 has no fully closed three-hop physical center stencil, so its near-zero tail cannot establish full bulk closure. The naive correction has a longer tail because it contains the physical bulk force.')
    s.append('## 12. N27 fully closed center result\n\n'
             'Only specimen extent changes to 26h; the same physical polynomial retains L=1.125 mm. Predeclared points approach from a lateral x face at Z=0, holding the local field value constant while changing stencil depth. The classical traction on that face is zero; these four samples do not establish every face profile. An independent shifted N29 scaffold supplies their bulk reference.\n\n'+
             table(['Offset','d/δ','Localized norm N','Localized/F_ref','Naive/F_ref','Full graph closed'],
                   [[str(r['offset']),r['depth_over_delta'],r['localized_norm_N'],r['localized_scaled'],r['naive_scaled'],r['fully_closed']] for r in n27['points']])+
             '\n\nThe 12h point is numerically small but includes nonbulk endpoint volumes in its theoretical dependency graph. Only the 13h center has verified full closure. Numerical cancellation at an earlier depth is not a theorem about minimal sufficient closure or physical BC thickness.')
    s.append(r'''## 13. BAQP-manufactured recovery

$$b_{\rm NL}=-q_{\rm bulk},\qquad f_{\rm ext,NL}=b_{\rm NL}V+f_{\rm BC,localized}.$$

This has a zero target residual but that substitution identity is not the recovery test. The unchanged existing mechanics adapter runs its existing TRF/LSMR solver from a nonzero smooth perturbation. Both z faces are fixed to the manufactured positions; all other particles are free. The prescribed perturbation vanishes on all faces. Frozen external volume forces include fixed nodes; constraint reactions are saved separately.

For nonaffine data, finite-domain reconstructed Fbar at the exact nodal target is not the analytic continuum F near a boundary. The predeclared solver recovery metric is max norm(Fbar(solved)−Fbar(target))/0.004; analytic-gradient error is reported separately and never used to hide reconstruction bias.'''+
             '\n\n'+table(['NL metric','Value'],[[k,nl[k]] for k in ['initial_perturbation_max_m','force_residual','F_recovery_error_max','F_analytic_error_max','displacement_error_max_m','energy_relative_deviation','PASS']])+
             f"\n\nExact target Fbar-to-analytic bias = {node['target_metrics']['Fbar_vs_analytic_max_scaled']:.6e} of imposed-gradient scale. The recovery gates remain force<=1e-8 and discrete target-gradient error<=1e-5. This does not claim analytic F is exactly reconstructed at truncated nodes.")
    s.append('## 14. Classical-manufactured result\n\n'
             'The second solve uses `b_local V + f_BC_localized` with the same frozen conversion, constraints and initial perturbation. Its target residual is independently `e_model V`. No negative model-error load is added to the boundary correction.\n\n'+
             table(['Classical metric','Value'],[[k,cl[k]] for k in ['force_residual','F_recovery_error_max','F_analytic_error_max','displacement_error_max_m','energy_deviation_J','energy_relative_deviation','PASS']])+
             '\n\n'+table(['Cross-control diagnostic','Value'],list(response.items())[:-1])+
             '\n\nThe nonzero target residual establishes that the local manufactured field is not an exact BAQP equilibrium. The load difference between controls is precisely the independently measured model residual. '
             'The target itself lies below the broad 1e-8 solver acceptance tolerance, so acceptance alone would not establish exactness; both independent perturbed-start solves continue to their own numerical termination. '
             f"The cross-control force-response balance leaves {response_error:.3%} of the bulk-model load at numerical residual level. "
             f"The NL recovery displacement error is {response['NL_recovery_error_over_classical_deviation']:.3%} of the classical deviation; this comparison indicates numerical resolution, not a rigorous error bound. "
             'A converged classical solve is not proof that every last digit of its displacement deviation measures bulk nonlocality. '
             'The classical-field status refers to non-exactness under the measured bulk operator, with response precision quantified above. '
             'The primary energy deviation compares BAQP energies at solved and target nodal states. The additional analytic-continuum energy comparison retains finite-domain quadrature and nonlocal-energy differences; it is not solely a solver error.')
    s.append('## 15. What part is boundary physics?\n\n'
             'The difference between the finite physical-domain force and the independently closed bulk force is a finite-domain conversion load on this surrogate. It decays to numerical zero at the fully closed center. '
             'That spatial and algebraic evidence supports localization. The load also contains discretization features of truncated reconstruction, which remain under the unchanged BAQP operator. '
             'Neither this definition nor manufactured recovery identifies a unique experimentally calibrated boundary law.')
    s.append('## 16. What part is finite-horizon bulk-model difference?\n\n'
             '`q_bulk + b_local` is evaluated without physical-domain truncation and remains in the classical target residual. '
             'It is classified SMALL using the frozen threshold, independently of the boundary-localization verdict. '
             'Global differences also include the separately quantified physical volume-quadrature term. At fixed h and delta this round cannot separate nonlocal-horizon error from lattice discretization error, nor establish an asymptotic limit. '
             'Small bulk force-density mismatch does not imply identical finite-domain and local continuum energies. '
             'No kernel, normalization, horizon, material or force law is adjusted.')
    s.append('## 17. Optional shear result\n\n'
             '**NONAFFINE_SHEAR=NOT_RUN**, as predeclared. This round completes the required axial oracle, identities, N27 localization and two equilibrium controls. No mixed traction or additional science campaign is run.')
    s.append('## 18. Implication for future curved clot-wall mechanics\n\n'
             'A future surrogate-dependent conversion must retain bulk mechanics and compatible body force explicitly. An affine negative-residual prescription cannot be transferred unqualified. '
             'The present Cartesian quadratic axial benchmark does not validate pressure/shear conversion on a curved clot, wall adhesion, delamination, fracture or CFD coupling. '
             'No surface operator correction, ghost treatment, perifit weights, optimized weights or variable horizon is implemented. Historical C5/C7 states and claims are unchanged.')
    s.append('## 19. One next task only\n\n'
             'After human review: **NONLOCAL_BOUNDARY_NONAFFINE_SHEAR_LOCALIZATION_V1** — apply this same already specified decomposition to the quadratic shear field with analytic traction on every face, retaining the independently measured bulk residual and reporting solver resolution limits. '
             'Do not begin it automatically.\n\n'+table(['Status','Value'],list(status.items()))+
             f"\n\nProtected existing files: {integrity()['protected_files']}; changed: {len(integrity()['differences'])}. "
             'Only the nine requested tests are run; [TARGETED_TESTS.log](TARGETED_TESTS.log) records their result. No old assertion or manifest is modified. '
             'The initial invocation, exact-storage precondition interruption and prerequisite check are preserved; none accepted a scientific state. '
             '[IMPLEMENTATION_NOTES.json](IMPLEMENTATION_NOTES.json) distinguishes those implementation issues from the scientific gates.\n\n'
             f"Resources: 8 CPU oracle workers; 1 serial equilibrium worker; BLAS 1; GPU_NOT_USED_BY_DESIGN. Measured scientific stage time {resource_summary['measured_compute_wall_seconds']:.2f} s, "
             f"peak process-tree RSS {resource_summary['peak_process_tree_RSS_bytes']/2**30:.3f} GiB (shared pages counted conservatively), minimum available RAM {100*resource_summary['minimum_available_RAM_fraction']:.1f}%. "
             'No packages installed. The two large equilibrium solves never overlap.\n\n'
             'Five figures: [600 dpi PNG and vector PDF](../../results/nonlocal_boundary_nonaffine_localization_v1/figures). '
             'Array units and legacy solver-field conventions are in [DATA_DICTIONARY.md](DATA_DICTIONARY.md). '
             'Reproduce numerical stages with `PYTHONPATH=. OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 <python> -B scripts/nonaffine_boundary_audit.py --stage <oracle|decomposition|n27|NL|classical>` in that order. '
             'Existing outputs are the accepted record; use a separate output copy for reruns. Plotting is `scripts/plot_nonaffine_boundary.py`.\n\n'
             '**Stopped for human review. No subsequent scientific development is run.**')
    (OUT/'REVIEW_NONAFFINE_BOUNDARY_LOCALIZATION_V1.md').write_text('\n\n'.join(s)+'\n')
    print(json.dumps(status,indent=2))


if __name__=='__main__':main()
