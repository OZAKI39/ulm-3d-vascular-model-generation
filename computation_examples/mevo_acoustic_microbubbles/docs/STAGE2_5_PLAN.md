# Stage 2.5 plan and predeclared gates

Baseline: Stage2 commit `6d659c395007a87741d997364867c9a5b7e9cdbe`;
work branch `stage2.5-pair-physics`. Stage1 tag and Stage2 branch remain fixed.
All new runs go to `results/stage2_5`; previous artifacts are read-only.

1. Read actual source/results and re-run server doctor; verify old manifests.
2. Zhang forensic: transcribe a parameter ledger; independently implement the
   four second-order equations as residuals, audit SI dimensions/sign/time origin,
   reproduce a predeclared matrix and retain full failed-case intermediate data.
   Investigate integration windows and documented parameter ambiguities without
   fitting parameters or changing a sign to match a curve. A documented unresolved
   discrepancy is an acceptable forensic outcome, not a physics PASS.
3. Reattempt author supplementary-code retrieval. Compare three pair descriptions
   on the same **uncoated** weak-field problem. Separate the formal incompressible
   far-separation limit from finite-sound-speed/retardation limits; never compare
   coated Zhang directly with uncoated Doinikov as if their radial physics matched.
   Record force definition, residual radial-model mismatch and observed scaling.
4. Build empirical pair and translation validity maps. Accuracy budgets are 5%
   (safe) and 20% (approximate); these are project budgets, while boundaries must
   come from measured errors. Unknown shell/parameter combinations stay outside
   the validated range; maps cannot silently extrapolate to SonoVue.
5. Adaptively subdivide the Stage2 focused domain R=1.6–1.96 um,
   f=2.3–2.7 MHz, P=5–40 kPa, using the unchanged cited shell mapping.
   Profiles: strict2%, production5%, exploratory10%; test complex fundamental,
   amplitude, phase, extrema and derivative-weighted volume energy. Use midpoint
   plus interior checks, shared deterministic crash-safe APECSS node cache,
   explicit depth/node limits and rejected leaves. Passing samples are empirical
   interpolation evidence, not proof of a cell-wide error bound or agent calibration.
   Benchmark 4/6/7 workers under the actual cgroup quota; never duplicate a
   completed deterministic node. Reserve disjoint deterministic node sets for
   worker benchmarks, report workload dependence.
6. Map quasi-steady versus Doinikov–Dayton displacement error at matched forcing
   and radius/pressure/frequency, keeping its shell distinct from SonoVue.
7. Add executable evidence-backed gate management and status watermarks.
   Lazarus requires G1 resolved, or a demonstrated disjoint validated regime;
   G2 measured pair map; G3 valid response domain or independent paper model.
   Missing evidence closes the gate. No override is used in this work.
   Cylinder requires pair/response/translation/Lazarus gates for its own config;
   frequency sweep additionally requires cylinder PASS. No patient/brain/skull run.
8. Independent contact-candidate benchmarks (N128..10000, multiple densities),
   explicit GPU fallback/selection tests, existing VTK regression plus SI scaling.
   Optional one production-size acoustic comparison if memory/time are suitable.
9. Preserve failures; complete tests/lint, baseline hash/trajectory regression,
   provenance, 21-section report and final source/result manifests.

Mandatory deliverables are forensic evidence, independent implementation,
cross-model comparison, empirical range guards, adaptive database/cache/resume,
translation map, provenance and passing software checks. Scientific outcomes
remain partial if any required physics gate is unresolved.

Execution refinement: the initial longest-axis tree failed an independent holdout
and remains rejected. Measured-axis curvature improved coverage to98.7427%; a
new child reused its tree and expanded the budget to24,000 nodes/depth20 without
changing tolerance or physics. It used12,902 nodes and passed256 fresh holdouts
with full declared-domain coverage. The final report records both failed and
successful attempts. G1 remains unresolved; all conditional downstream runs stay closed.
