# Stage 4 — predeclared execution and scope

The Stage 3 baseline is the annotated tag `stage3-minimal-validated`, commit
`dfe79a640fa54c164ce523dcdd8e129217148b5a`. Its scientific results and physics acceptance remain unchanged.
Development occurs on `stage4-bifurcation-control`. No patient geometry or Stage 5 work is authorized here.

Execute in order: A freeze/regression; B configurable smooth Y and crossing classifier;
C synthetic software-only flow; D existing SimVascular steady-flow workflow; E VTU import checks;
F new geometry-keyed k-Wave basis; G primary-only B0–B3; H paired multi-seed study;
I separate secondary-ON B2/B3; J time/grid/sampling/classification sensitivity;
K reference-only history analysis of deterministically selected critical particles.

Production constitutive laws, radial-table rejection and geometric contact remain Stage 3 physics.
Geometry and recording adapters may extend the application while retaining the existing
`ParticleTracker(config, acoustic_field, flow)` and `FlowField.velocity(x,t)` signatures.
No production history, streaming, wall hydrodynamics or feedback is introduced.

Acceptance must preserve negative evidence. Synthetic trajectories are `SOFTWARE_TEST_ONLY`.
CFD needs finite SI data, no-slip checks, inlet/outlet flux integration, near-symmetric branch split
and steady-state evidence. Flow import is checked independently at nodes, interior points,
near-wall points and the junction before coupling. Coarse/fine acoustic fields are compared at
identical coordinates for pressure, gradient, force and branch assignments.

Use paired initial IDs and positions across B0–B3 for at least 10 seeds and N=256.
The baseline sampling policy and geometry/flow/acoustic choices are recorded before formal trajectories.
Uniform-area sampling and N=512 are sensitivity conditions if cost permits.
First outward plane crossing determines the permanent branch; unresolved/wall-loss outcomes remain explicit.
The target effect is compared with paired bootstrap 95% intervals, between-seed variation,
dt/2, geometric classification tolerance and acoustic-grid differences; no arbitrary enrichment threshold.
Primary+secondary outcomes are analyzed separately, including close-pair validity and clustering.

History is allowed only under `validation/history_branch_robustness/`, on 10–20 particles chosen
by an automatic recorded decision-boundary rule. Its own numerical/model limits must be reported.
The requested agreement categories (≥95%, 80–95%, <80%) are reporting conventions, not physical laws.
An inconclusive or failed component cannot be silently marked PASS or replaced by synthetic evidence.
