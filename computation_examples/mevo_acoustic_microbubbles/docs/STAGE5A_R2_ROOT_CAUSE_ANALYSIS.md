# Stage 5A-R2 root-cause analysis

Evidence root: `results/stage5/motion_solver_r2/diagnosis/`.
These findings are new read-only reconstructions, not reinterpretation of old
failures as passes. All three reconstructed A/G matrices match the immutable
failure fixtures exactly (zero entry differences). Row metadata includes stable
IDs, stages, normals, anchors, original/scaled lower bounds and wall exterior
feature indices. Physical velocity/resistance/force scales are recorded.

| Hypothesis | Finding | Evidence and limits |
|---|---|---|
| A: stacked start/predictor near-dependent tangents | CONFIRMED | `reconstruction_v1/REPLAY.json`, spectra and near-row lower differences; identical A/G reconstructions. Dependence is not infeasibility and no near-parallel rows were deleted. |
| B: guessed active set/dual polish fails | CONFIRMED | Re-executed old exceptions retain pre/post-polish residuals; B2/B3 guessed saddle polish has tiny linear residual but fails original KKT. Native ProxQP solves all rows to <=1e-12 without that polish. |
| C: native stopping differs from independent audit | CONFIRMED | Old B3 inaccurate candidate passes 1e-9 KKT but old status rule rejects; new native Clarabel can report Solved with KKT >1e-9. Neither native status alone nor a saddle linear residual is sufficient. |
| D: coefficient-lagged/time error limits old PC | CONFIRMED | All 20 original LOCAL_ERROR trials reproduced; responsible contacts, normals, activation/floor switches and coefficients saved. At B2 39.0625 ns, IDs 110/170 cross activation with gap 190.59->147.16 nm and error 1.7409. B3 IDs 206/221 have coefficient change 5.78e-9->2.70e-8 N s/m, error 2.1831. These are actual changes, not only algebraic failures. |
| E: QP exception bypasses interval recovery | CONFIRMED | Old interval catches only SubdivideInterval. Replays terminate with MotionSolverFailure. R2 uses bounded native alternatives and classifies subsequent interval recovery. |
| F: stage workspace thrashing | CONFIRMED | Old component key contains IDs only; start/corrector G patterns alternate in the same workspace. Reconstructed stages have different row sets. R2 stage+solver+component+pattern identities separate them and numeric updates are explicit. |
| G: known feasible predictor ignored | CONFIRMED | Reconstructed predictor minimum slack: -2.57e-16 / -4.61e-19 / -4.34e-18 on all original rows. Zero violates lower bounds by 1.69e-4 / 2.84e-7 / 2.02e-5 (normalized). Old OSQP ignores its feasible parameter. |
| H: continuous wall certificate absent | CONFIRMED | Old source checks wall endpoints only. New test has endpoint-clear path crossing a triangle, plus finite-radius edge graze. R2 uses full swept finite sphere–WALL certification. This confirms missing coverage, not proof an old accepted path actually crossed WALL. |

The initial checkpoints have strong transient motion and many near-contact
switches. The old 20 us PC trial changes hundreds of pair states and is far from
the accepted time-error budget. A faster frozen QP therefore cannot alone repair
continuous motion. R2 separates this from algebraic recovery.

The first new endpoint-contact S1 candidate introduced a distinct intermediate
path problem in initial negative-gap tolerance bands. Its rejected sweep is
retained; V2 fixes start tangents and iterates endpoint resistance. V3 uses bounded
Anderson iterations. No initial coordinates, force, radius, floor or geometric
tolerance were changed. Subsequent numerical/scientific/holdout results remain
separate from these root-cause findings.
