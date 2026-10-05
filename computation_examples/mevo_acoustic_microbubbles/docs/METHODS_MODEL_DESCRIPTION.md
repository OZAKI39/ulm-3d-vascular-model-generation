# Methods: minimal PZT–microbubble simulator

The simulator addresses comparative changes in migration and spatial clustering under PZT control. Its numerical acceptance is scoped to the reported medium, response table and demonstrations; it is not an experimentally calibrated predictor of transient displacement.

A linear k-Wave solution is cached for each array element, H_i(x), with complex convention Re[P exp(iωt)]. The field is P=Σu_iH_i, u_i=A_i exp(iφ_i). A new propagation solve is needed for a changed medium, source geometry, grid or frequency, but not changed control amplitudes/phases. The water, 2.5 MHz four-element basis was computed once for Stage 3 because the old acoustic bases have a different medium.

The existing APECSS RPAR/Marmottant adaptive database supplies V_h(R0,f,|P|), rotated by h times the local pressure phase. Strict medium, domain, cell-validity and harmonic-resolution checks are retained. Production never calls a bubble ODE at each particle step. The database uses the Tu2009-derived fixed shell manifold and an explicitly assumed initial surface tension; these are not fitted to the Stage 3 trajectories. Effective radial viscosity and physical translation viscosity are separate declared quantities.

Primary force is −½ Re[V_1 conj(∇P)]. The reciprocal classical secondary force is ρ/(4πd²) times the harmonic derivative correlation times (x_j−x_i)/d; retained harmonics carry h² weights. Same-phase identical responses attract. Stage 3 has no distance regularization, pair-force cutoff or empirical force rescaling. Geometric contact is handled separately.

Particle velocity is u_flow+(F_primary+ΣF_secondary)/(6πμR0). The existing RK integration and displacement-limited geometric substeps advance this slow-time law. cKDTree finds candidates; sorted exact-distance filtering and the original mobility-weighted position/velocity constraints enforce nonpenetration. Brute-force search remains selectable. Contact does not imply adhesion, coalescence, fusion or accurate near-gap hydrodynamics.

FlowField.velocity(x,t) supplies zero, Poiseuille or static/time-indexed VTU velocity without changing bubble physics. Coupling is one-way. Synthetic space/time-linear VTU particle trajectories were compared with direct analytic trajectories. This verifies the interface, not the validity of any future CFD field.

Configurable pair separation categories and oscillation/retardation/wall checks report extrapolation. They do not alter forces. The Stage 2.5 evidence informs warning thresholds rather than providing a universal coated-pair error map. Results that enter near-contact carry `NEAR_CONTACT_PHYSICS_APPROXIMATED`.

The separate [history study](HISTORY_FORCE_SENSITIVITY.md) isolates instantaneous added-mass and history effects with matched inputs and reports the difference introduced by production averaging separately. Its discrepancies require explicit limits on quantitative transport times. No history-force reference module is a production dependency.

[Physics acceptance](../PHYSICS_ACCEPTANCE.json), [Stage 3 report](STAGE3_VALIDATION_REPORT.md) and [known limitations](KNOWN_LIMITATIONS.md) identify the evidence and scope. Zhang 2011 was investigated but was not used as a validation basis because a key published trend could not be independently reproduced.
