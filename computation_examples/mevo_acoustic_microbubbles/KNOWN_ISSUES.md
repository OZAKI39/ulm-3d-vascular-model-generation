# Known issues and retained failures

Stage2.5: the independent printed-equation reference confirms the Fig.3 discrepancy
(UNRESOLVED). The adaptive focused domain1.6–1.96 µm /2.3–2.7 MHz /5–40 kPa now has
790 valid cells and256 passing independent holdouts. This numerical result does
not calibrate SonoVue or unlock clustering. The first adaptive attempt failed a
holdout and remains explicitly rejected. Pair evidence is sampled uncoated weak
field; the translation map uses a solid shell and20-cycle pulse.
See [Stage2.5 report](docs/STAGE2_5_VALIDATION_REPORT.md).
The numbered items below describe retained Stage2 results and limitations.

1. **Stage2G is not fully validated.** Zhang2011 Fig3 sign reversal at120 µm fails with the transcribed equations and χ=1 N/m. Attraction/repulsion and the2.9 µm resonance peak pass qualitative checks. Time refinement does not fix the sign. Do not tune parameters to hide it. H/J remain gated.
2. **No new Lazarus clustering or validated cylinder sweep.** Three pressure configurations are explicitly NOT_RUN. The source uses a near-wall2D agent/population and imposed cluster reorganization; required measured initial coordinates are unavailable numerically. FigureD is a gate-status panel, not clustering data.
3. **Response-table coverage is limited.** Broad coarse table:59/60 valid nodes,0/24 valid interpolation cells. Focused table:125/125 valid nodes,39/64 cells pass midpoint testing. Untested/failed cells and out-of-bounds queries are rejected. A midpoint pass does not certify a whole cell.
4. **SonoVue calibration is incomplete.** Fig4 points are approximate digitization; σinit=.036 N/m is a declared hypothesis. The40 kPa example is BR14. Continuous steady-state excitation differs from the diagnostic pulse. Table medium is water, not blood.
5. **Effective radial damping differs from translation viscosity.** Tu2009 uses .002 Pa s including thermal losses. `radial_effective_viscosity_Pa_s` can declare this separately from physical `medium.viscosity_Pa_s` used in Stokes drag. Reusing effective viscosity for drag needs its own justification.
6. **Near-field reference limitations.** Doinikov2015 is N=2, uncoated/spherical, finite Legendre projection, assumed Gaussian width. Author Mathematica supplement could not be downloaded. Numerical order convergence does not validate contact, lubrication or shape modes.
7. **Translation approximations.** Quasi-steady Stokes lacks added-mass transients/history force and finite-Re corrections. Reference examples include Re>1; no general accuracy claim is made.
8. **Wall hydrodynamics remains TODO.** Geometric projection is not the Doinikov wall model. Near-wall strict queries stop; no adhesion/coalescence law is provided.
9. **GPU path is hybrid.** Force kernels/table interpolation use CuPy; flow interpolation, tracker clock and NNLS/contact projection remain CPU. N² allocations are capped; a physically justified far-field cutoff is not implemented.
10. **C++ acoustic adapter is version-specific.** Missing serial HDF5 libraries were resolved by private package extraction.0.6.2 also lacks mass-source scaling/C↔F ordering and returns time×sensor full records. Our adapter corrects these for homogeneous isotropic3D pressure sources using the upstream scale function. Other package versions are refused pending re-audit. Small32³/64³ equivalence passes; a full production-grid speed/memory benchmark has not been performed. CuPy remains default.
11. **Retained setup/benchmark failures.** Early APECSS orchestration failures, failed table cells, Zhang gate failure, first C++ shape failure, unscaled C++ pressure error≈1.43e7 and a temporary scalar/array adapter error remain in timestamped runs/logs. Final successful checks do not erase them.
12. **Stage1 limitations remain.** Dense aggregation trajectories were not converged; geometric clusters are not fusion. The old100-bubble interrupted probe is unchanged. BraVa formal1565-contactfix trajectories retain141 observation cutoffs and2 numerical failures; no new BraVa injection-protocol batch was run.
