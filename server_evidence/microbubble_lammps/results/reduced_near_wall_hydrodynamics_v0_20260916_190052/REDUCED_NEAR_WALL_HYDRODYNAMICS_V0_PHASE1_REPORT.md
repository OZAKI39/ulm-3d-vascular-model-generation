# REDUCED_NEAR_WALL_HYDRODYNAMICS_V0 — Phase 1

**STATUS: BLOCKED — STOP_GCB_REFERENCE_INSUFFICIENT.**

The reference audit is complete. The standalone C++ kernel was not implemented or compiled. Algebra tests A–J were not run, and no epsilon domain is certified. This stop follows sections 5 and 25 of the current user request. The stage is an audit deliverable, not a qualified wall model.

## 1. Scope

The authorized target is a rigid sphere near an infinite rigid no-slip plane in Stokes flow, with wall excess resistance and analytic wall-tangential simple shear entering the RHS. No trajectories, LAMMPS integration, Palabos changes, WSS generation, arbitrary ambient strain, RBC, adhesion, BIE/QBX, fixed-multiblob rerun, or GPU rewrite were performed.

## 2. Source and provenance

Stage timestamp: `20260916_190052` (UTC).

- Local: `/home/lzy/projects/compre_output/reduced_near_wall_hydrodynamics_v0/20260916_190052`
- Vast results: `/workspace/microbubble_lammps/results/reduced_near_wall_hydrodynamics_v0_20260916_190052`
- Vast work: `/workspace/microbubble_lammps/work/reduced_near_wall_hydrodynamics_v0_20260916_190052`
- Build locator: `<Vast work>/build`; not created because implementation stopped.

`provenance/UPSTREAM_INPUTS.tsv` records absolute path, classification, SHA256 and size for all seven requested upstream inputs, their relevant local counterparts, and inherited contracts/state files. `BASELINE_BEFORE.json` and `BASELINE_AFTER.json` record the baseline checks and local/remote OS, compiler, CMake, MPI and timestamp. The former is the machine-readable snapshot before evidence finalization; initial read-only checks had already confirmed the same source identities before stage creation.

Initial checks: expected Git branch and commit matched; only the two permitted untracked environment reports were present. Vast was `f7c62a262077`, RTX 4090 24564 MiB; all 34 inherited key path locators existed. The selected production/wall source hashes, Frozen Flow hash and vascular STL hash matched the handoff. The local reference scripts also matched their Git archive copies.

The code added here consists only of reference-audit, plotting and provenance utilities. No historical runner was executed. Final source/payload synchronization evidence is in `provenance/REMOTE_LOCAL_SHA_CHECK.json`; source matching there refers to these audit utilities, not to a nonexistent C++ kernel.

## 3. Mathematical convention

The intended, unimplemented kernel equations remain:

`R_total = R_bulk + R_wall_excess`, where
`R_bulk = diag(6*pi*mu*a,6*pi*mu*a,6*pi*mu*a,8*pi*mu*a^3,8*pi*mu*a^3,8*pi*mu*a^3)`.

Velocity order is `[V_t1,V_t2,V_n,Omega_t1,Omega_t2,Omega_n]`; conjugate load order is `[F_t1,F_t2,F_n,T_t1,T_t2,T_n]`. Loads are externally applied forces and torques; torque is about the sphere center. Units are m, s, Pa*s, N and N*m.

## 4. Local frame

`d=a+h`, `epsilon=h/a`, and `n` points from the wall to the sphere center. `Q=[t1 t2 n]` is right-handed and orthonormal; `T=diag(Q,Q)`. Local matrices would transform as `R_global=T R_local T^T`. Rotations and tangent-basis sign invariance have not been tested on a new kernel.

## 5. Resistance table interpretation

RMBW input:

`/workspace/microbubble_lammps/results/wall_hydrodynamics_v0_20260916_130547/tables/RMBW_WALL_RESISTANCE_TABLE_V0.h5`

SHA256: `71e62ef6bc22edc66ae0c4d1dcc224378c2f6c5d3848578f3fbedd2cd9a43f98`; size: 4,897,730 bytes.

The 8,211-row table explicitly stores `R_wall_excess_scaled`. Across every row, `R_total_scaled - R_wall_excess_scaled - I` is exactly zero in the stored floating-point arrays. With `S=diag(sqrt(6*pi*mu*a) x3,sqrt(8*pi*mu*a^3) x3)`, SI excess is `S R_wall_excess_scaled S`. The independently checked stored total-SI reconstruction relative error is approximately `1.13e-16`. Thus bulk must be added once, and must not be subtracted again from the stored excess.

Reference metadata uses `a=9.683592065545495e-7 m` and `mu=0.001 Pa*s`; the table is dimensionless in epsilon. Consequently TT, TR and RR scale respectively as `mu*a`, `mu*a^2`, and `mu*a^3`. Dataset names, shapes, attributes and limitations are fully recorded in `reference/RMBW_TABLE_METADATA.json`.

The old lookup uses componentwise linear interpolation of the full matrix in log(epsilon), restores SI scaling, then rotates to global axes. It also clamps below 0.001 and sets excess to zero above 20. Those exterior behaviors were not adopted. A future reference kernel must reject nonpositive gaps and unsupported table queries explicitly.

Historical support is mode-dependent: Brenner supports normal translation; O'Neill provides discrete translation/induced-torque comparisons; parallel RR retains an unresolved 0.3817/0.3709 constant discrepancy; continuous finite-gap RR/TR accuracy remains unverified. The normal-axis rotational contact limit is supported, but no finite-gap accuracy bound is established. Table presence is not all-gap exactness.

## 6. GCB convention and blocking evidence

The requested source is [Goldman, Cox & Brenner 1967, Part II](https://doi.org/10.1016/0009-2509(67)80048-4). The publisher abstract confirms the Couette-flow problem, but does not provide the needed numerical closure. None of `classical_wall_reference.py`, `compare_oneill_table.py` or `rmbw_reference_adapter.py` implements free-sphere Couette translation and rotation. The `GCB_Ytt/Ytr/Yrr` keys are quiescent-fluid resistance asymptotics, not `FU/FOmega`.

The intended normalization is `FU=V/(d*S)` and `FOmega=Omega/(S/2)`, both tending to one far from the wall. In the primary paper's notation, center height `h` maps to this project's `d`, and surface gap `delta` maps to this project's `h`.

A [public transcription of the primary paper](https://www.scribd.com/document/684300240/Goldman1967-b), printed p657, presents expressions (4.11)/(4.12) with numerator constants 0.7431/0.8436 and denominator `0.6376-0.200*ln(epsilon)`. Its Table 3 transcription does not numerically agree:

| epsilon | quantity | printed-expression evaluation | transcribed Table 3 | relative difference to table |
|---|---|---:|---:|---:|
| 0.001 | FU | 0.3680259572 | 0.3966 | 7.204751% |
| 0.001 | FOmega | 0.4177993507 | 0.4268 | 2.108868% |

These are **UNVERIFIED_TRANSCRIPTION_DIAGNOSTIC** values. No authenticated original-page image or trusted numerical dataset was obtained. They establish a discrepancy requiring resolution, not proof that the original paper itself is erroneous. Table 3 is described as derived from the near-wall expressions and cannot serve as an independent numerical validation merely by repeating those expressions.

[Chaoui & Feuillebois 2003](https://doi.org/10.1093/qjmam/56.3.381), section 6, discusses accuracy problems in the earlier printed expression and a possible misprint; this was accessed through a [paper transcription](https://doczz.net/doc/7418969/creeping-flow-around-a-sphere-in-a-shear-flow-close-to-a-...). A [corrigendum](https://doi.org/10.1093/qjmam/hbs012) was identified but not checked. Their later approximants were not selected. [Tang et al. 2024, Table 1](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023WR035456) reproduces the familiar small-gap expressions; repetition does not resolve the discrepancy or independently certify their accuracy.

The user allows a sourced near-wall asymptotic result with PARTIAL physical validation. This stop is more specific: the candidate closure and reference evidence disagree, and the primary-page/numerical check needed to resolve the disagreement is missing. We therefore did not freeze even a near-wall closure, repair constants, choose another formula, or interpolate intermediate-gap values.

For a later resolved closure, shear must use `g_t=(I-n*n^T)g_wall`, `q_free=[d*g_t,0.5*(n cross g_t)]`, `q_GCB=[FU*d*g_t,0.5*FOmega*(n cross g_t)]`, and `b_wall_shear=R_total*q_GCB-R_bulk*q_free`. No solve or post-solve velocity replacement was performed here.

## 7. Gap domains

- RMBW data domain: `[0.001,20]`; inherited near-field/qualified/limited labels retain their mode-specific limitations.
- GCB candidate printed-expression scope: `epsilon << 1` only; no all-gap claim.
- Diagnostic plot range `[0.001,0.05]` is a display range, not an accepted kernel domain or active cutoff.
- Certified Phase 1 epsilon domain: **NONE**.

No HDF5 extrapolation, production cutoff, GCB interpolation or dynamics-layer gap correction was selected.

## 8. Tests and diagnostic artifacts

`validation/PHASE1_VALIDATION.json` marks all required kernel tests A–J **NOT_RUN_REFERENCE_STOP**. Negative/zero gap, NaN/Inf, SI scaling, dissipation, symmetry, frame invariance, shear identities, zero normal forcing and actual linear-solve recovery have not been qualified in a C++ implementation. `scripts/validate_phase1.py` was not created: there is no C++ raw output to validate.

Completed checks only inspect the frozen input table's identity, SI reconstruction, ordering and finite entries. The plotted resistance uses total input-table coefficients, with normalization explicit on the axes. Three static reference-only plots and corresponding CSVs were saved:

- `visualization/resistance_vs_epsilon.png` and `raw/resistance_reference.csv`
- `visualization/gcb_translation_vs_epsilon.png`
- `visualization/gcb_rotation_vs_epsilon.png`
- `raw/gcb_transcription_differences.csv` and `raw/gcb_printed_expression_diagnostic.csv`

GCB figures explicitly label the unqualified printed expression and unverified transcription points. `shear_rhs_components_vs_epsilon.png` was not generated: showing a zero RHS without an implemented, tested closure would be misleading.

## 9. Acceptance results

| Gate | Result |
|---|---|
| Git / selected baseline identities | Matched; before/after provenance retained |
| RMBW excess / units ambiguity | Resolved for the frozen table |
| Trusted continuous GCB closure | BLOCKED |
| ALGEBRA_PASS | NO — not run, not a failing numerical solve |
| PHYSICS_REFERENCE_STATUS | BLOCKED |
| Certified epsilon domain | NONE |

Historical wall-table interpolation action tolerance is `0.001` and reciprocity tolerance `1e-10`. The older reference audit explicitly has no universal 1% physical accuracy gate; its 2% agreement band is descriptive. No GCB physical tolerance was inherited or invented. Future tolerance should be set against a checked reference uncertainty budget before certification. Proposed algebra thresholds from the user's request are recorded as unexecuted requirements in `contracts/PHASE1_SCOPE_AND_STOP.json`.

## 10. Known limitations

No standalone kernel, resistance API, shear API, enums, C++ tests or kernel Python validator exists in this stage. Existing mathematical and reference capabilities remain upstream evidence only. The GCB transcription is not an authenticated data source; later-paper corrections still need verification. Existing finite-gap RR/TR limitations, pair twist pending, transient PBS/BSA field, unavailable validated WSS and other inherited project uncertainties are unchanged.

## 11. Physics reference status

**BLOCKED**, not PASS or a production promotion. `REFERENCE_CAPABILITY.json` distinguishes available resistance evidence from unavailable qualified free-shear references. The diagnostic expressions in audit scripts are not a selected hydrodynamic closure and must not be integrated into production.

## 12. Recommendation for Phase 2

Do not start Phase 2. First resolve Phase 1 by obtaining and visually checking the GCB equation/table pages or a trustworthy numerical reference with known normalization, checking relevant corrections, and explicitly selecting the allowed closure and validation bounds. If continuity requires interpolation or a different published approximation, present that choice for user review before implementation. Then resume the standalone A–J and independent-reference qualification. `next_authorized_stage=NONE`.

No build, simulation, GPU computation, historical source modification, Git add/commit/push, or Phase 2 work was performed.
