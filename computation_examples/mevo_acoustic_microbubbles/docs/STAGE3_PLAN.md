# Stage 3 implementation and acceptance plan

Baseline: ff6afaf, retained on stage2.5-pair-physics. Work branch: stage3-minimal-simulator.
No BraVa, old WSS, pressure-sensitivity, or stopped population jobs are resumed.

1. Preserve the Zhang investigation by a logical archive; replace downstream high-fidelity gates with scoped implementation/transport gates.
2. Keep the verified radial table/cache. Use its water medium and 2.5 MHz domain. The old acoustic bases use a different medium, so one new matching four-source basis is necessary. Amplitude and phase changes reuse it.
3. Add warning-only geometric pair regimes, deterministic cKDTree candidates and production profiling, retaining brute reference and exact contact resolution.
4. Establish secondary algebra + direct APECSS weak-field scaling, multi-bubble early-time pressure trend and timestep convergence.
5. Independently compare matched forcing/radius A/B/C translation in validation/history_force_sensitivity, using the primary variable-radius zero-Re no-slip history integral. Document averaging effects separately; do not attribute them to history.
6. Verify synthetic VTU trajectories through the existing FlowField API. Admit four formal cylinder controls only after prerequisite evidence exists.
7. Run regression, record scoped PHYSICS_ACCEPTANCE.json, generate scientific figures, update Methods/limitations and the final report. No patient geometry or frequency sweep.

Planned tolerances (fixed before runs): algebra/scaling synthetic slopes 1e-10; APECSS weak-field exponent within 0.05 of 2; early clustering timestep displacement relative difference < 2%, at least 80% of early pair checks inside configured classical policy; history successive-grid endpoint changes < 2%; analytic linear-VTU trajectory absolute error < 1e-10 m; four cylinder runs finite/nonpenetrating, C1 analytic drift, detectable C2 lateral migration and C3/C4 distribution change. These are numerical acceptance limits, not experimental uncertainty.
