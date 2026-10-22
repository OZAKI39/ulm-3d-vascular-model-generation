# Physics provenance and decision record

## Current Stage 3 decisions

The fixed production chain and equations are in [Methods](METHODS_MODEL_DESCRIPTION.md). Water at 2.5 MHz uses a new matching four-element k-Wave basis, reused for all controls, and the unchanged Stage 2.5 adaptive APECSS table. No new production radial ODE or uniform table was created. Table range/medium rejection remains strict. Geometric separation warnings do not change the unregularized secondary force.

Classical secondary force is accepted by independent phase, distance, weak-pressure, reciprocity and zero-response checks. The correlation follows from the spherical monopole pressure and integration by parts over a periodic cycle; it is not accepted on the basis of Zhang reproduction.

Zhang 2011 was investigated but was not used as a validation basis because a key published trend could not be independently reproduced. Its evidence is [logically archived](../validation/deprecated/zhang2011/README.md); no downstream Zhang/high-fidelity pair gate remains. Lazarus 2017 is a trend benchmark, not parameter-identical reproduction.

History equations are implemented only in validation/history_force_sensitivity; full source/version distinctions and quantitative errors are in [the sensitivity report](HISTORY_FORCE_SENSITIVITY.md). No parameter or sign was changed to obtain desired trends. Initial surface tension, fixed shell manifold, absence of dissolution and the production translation approximation remain explicit assumptions.

The acceptance ledger separates numerical correctness, model scope and warnings. Neither CPU/GPU consistency, global consistency nor smooth plots prove local physical accuracy.

## Historical Stage 2 / 2.5 investigation


Stage2.5 development is `stage2.5-pair-physics`; the Stage2 branch below is retained
as a fixed baseline. New equation decisions, forensic ledger, cache identities,
source hashes and per-run evidence are indexed in
[STAGE2_5_VALIDATION_REPORT.md](STAGE2_5_VALIDATION_REPORT.md).
Independent residual agreement does not resolve the published Fig.3 discrepancy.
The adaptive database checks interpolation of the same APECSS model and fixed
parameters; no constitutive parameter was fitted to pass validation. Uncoated
pair, solid-shell translation, Lazarus paper and SonoVue scopes remain separate.
No N-body Doinikov sum, physical contact bond or patient model is introduced.

All runtime quantities use SI; complex fields use Re[P exp(+iωt)]. Scientific status is distinct from software-test status.
The original solver/results are preserved by `stage1-demo-baseline`; development is `stage2-validated-physics`.

| Code | Equations/source | Reuse level and scope | Status |
|---|---|---|---|
| acoustics/kwave_backend.py, basis.py | k-Wave first-order k-space pseudospectral propagation; Treeby2010; Wise2019 off-grid array | Existing k-wave-python0.6.2 / CuPy and kWaveArray, no solver rewrite | Software validated; two-grid sensitivity remains |
| bubbles/apecss.py | APECSS RPAR and interface.c; Marmottant2005 Eqs2–4, reversible σbreak=σwater specialization | Level1/2: unmodified C executable subprocess; v1.8 pinned commit | Official software comparison passed |
| bubbles/parameters.py | Tu2009 Tables1–3/Fig4, Eq6 | Direct transcription plus explicitly approximate digitization | Provenance established; nonlinear shell calibration incomplete |
| bubbles/response_table.py | Complex trilinear interpolation of APECSS harmonics | Offline CPU pool, NumPy/CuPy bounded batch lookup | Only valid nodes / midpoint-tested cells allowed |
| bubbles/response.py LinearResonantBubbleModel | Linear oscillator, explicitly assumed resonance/damping | Existing implementation retained only for regression/debug | Demo only, not SonoVue |
| bubbles/forces.py primary_bjerknes; translation.py | −⟨V∇p⟩; −1/2 Re[V1 conj(∇P1)] | Existing harmonic expression plus full-cycle quadrature | Monochromatic identity tested; no invented nonlinear discrepancy |
| validation/doinikov_dayton.py | DOI10.1121/1.2215228 Eqs6,31–33,36,40; small-signal Eq63/Eqs44–49 | Level3: no usable public author code found; minimal four-state thick-shell reference | Analytic/numerical check passed; not a SonoVue calibration |
| validation/zhang2011.py | DOI10.1088/1674-1056/20/11/114302 Eqs3–11 | Level3: publisher supplement absent, searches found no usable code | Partial trends; Fig3 sign inversion failed |
| validation/doinikov_bouakaz_pair.py | DOI10.1103/PhysRevE.92.043001 Eqs8–13,39–43,67–68 | Level3 fallback: author Mathematica supplement exists but inaccessible; direct finite Legendre system | Limited N=2 checks passed; Gaussian width assumed |
| bubbles/forces.py ClassicalFarFieldBjerknes | Classical periodic monopole identity: Fij=ρ⟨Vdot_i Vdot_j⟩eij/(4πd²) | Existing physical expression separated from old clipping; harmonics weighted h² | Guarded far-field kernel, CPU/GPU cross-check |
| bubbles/forces.py LegacyRegularizedBjerknes | Stage1 max-distance regularization | Compatibility alias NumpyPairBackend retained | DEPRECATED / DEMO ONLY |
| bubbles/translation.py, tracker.py | Quasi-steady slip F/(6πµR0) | Existing production approximation, compared with coupled reference | Approximation errors recorded; no history force |
| bubbles/contact.py, wall.py | Non-penetration projection / velocity feasibility | Existing geometric contact, no force calibration | Geometric only; no adhesion/coalescence |
| bubbles/wall.py DoinikovWallReference | DOI10.1088/0031-9155/60/20/7909 | Documented interface; explicit NotImplementedError | TODO, unavailable for physics |
| flow/fields.py | VTK cell interpolation, linear time interpolation | Existing PyVista/VTK dependency; generic wrapper retains velocity(x,t) | Synthetic VTU checks passed |
| validation/lazarus2017_clustering | DOI10.1121/1.5010170 Eq1, TablesI–IV, Figs3–7 | Full manuscript audit; simulation gated by failed Stage2G | Not reproduced; no manufactured curves |

## Search and reuse evidence

APECSS and PROTEUS official repositories were cloned before implementation; commits/licenses/tests are in
`third_party/THIRD_PARTY.md`, `APECSS_MODEL_AUDIT.md`, and `PROTEUS_REUSE_AUDIT.md`.
The original Marmottant paper, Tu2009 full text, DoinikovDayton2006 full text, Zhang2011 PDF,
DoinikovBouakaz2015 author full text and Lazarus2017 accepted manuscript were read.
APECSS/PROTEUS do not supply the particular coupled translation/pair equations used here.
For Doinikov2015, APS supplemental requests failed HTTP403 and HAL DOI search returned no record;
the manuscript explicitly says Mathematica coefficient-generation code exists. We do not claim it does not exist.
Only the necessary reference equations were implemented; APECSS remains the nonlinear radial production backend.
No upstream physics source was copied into the Python package.

## Operational validity is not a proof

Default policies: d/Ri and d/Rj≥10, d/(Ri+Rj)≥5, radial excursion≤0.2, kd≤0.3,
wall center distance/R≥5. These conservative configurable thresholds are operational choices, not sharp published universal cutoffs.
Table interpolation is allowed only in cells passing a 5% midpoint response check or at valid exact nodes.
Corner extrema form a conservative *node envelope*, not a mathematical bound on unsampled interior responses.
No extrapolation, hidden epsilon, field-medium mismatch or unresolved derivative-weighted harmonic tail is allowed in strict runs.
New table trajectories default to strict. Legacy configurations explicitly selecting the linear surrogate retain demo behavior and output labels.
Rejected states, last accepted HDF5 data and invalid reasons are saved before failure is propagated.
Evaluation-weighted and accepted-time-weighted pair fractions are distinguished.

`pair_reference` integrates a full uniform Gaussian N=2 problem; it refuses vascular flow, external PZT fields,
N≠2 and claims of coated physics. It cannot be inserted as a pairwise-summed many-body force.

## Contemporary cross-check

Cattaneo & Supponen is **Physical Review Fluids11,033606 (16March2026)**, DOI10.1103/p475-xggm.
The author preprint and CC-BY4.0 dataset DOI10.5281/zenodo.15692129 were inspected; the XLSX was downloaded,
while the 1.65GB video archive was not needed. No general public solver was found in that dataset.
Their history-force/thermal treatment and short-pulse experiments identify limitations of our quasi-steady model;
neither their fitted shell parameters nor translation scaling were substituted for SonoVue calibration.

## Reproducibility

Parent run SUMMARY/STATE files identify parameters, source hashes, hardware, seed or deterministic null seed, timing and status.
Individual APECSS jobs retain options, stdout/stderr, raw radius/speed trajectories and extracted harmonics.
The final source commit and Stage2 artifact manifest accompany the report; hashes captured during development
identify pre-formatting implementations and are not retroactively changed to the final commit.
