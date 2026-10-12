# F10F — bounded Small Steps result

Engineering PASS; all four SS2/SS4 scientific attempts fail before their
corresponding B0 failure windows. MORPHOLOGY_NOT_CONVERGED.
EXPLORATORY / PHYSICS_UNCALIBRATED / training_ready=false.

- analysis/FINAL_GATE.json: separate engineering, contact, morphology and timing gates.
- short_v1/: six original-initial-state attempts, including two exact B0 regressions.
- checks_v1/, failure_checks_v1/: CPU/GPU, graph, control and first-error checks.
- SUBSTEP_LEDGER.json.gz in each attempt: every actual substep and round scalar.
- CONTACT_TRACE.json.gz: last32 full attempted substep vector records.
- trajectory.h5: accepted saved frames; separate exact last_legal_state and rejected_interval.
- analysis/: actual contact, same-proposal PBD, morphology and all-cost yield.
- figures/: four English PNG/PDF figures, source hashes and visual review.
- progress_document/: unchanged-original backup evidence and minimal insertion diff.
- REMOTE_VERIFICATION.json / ARTIFACT_MANIFEST.json: hash-verified evidence.

No SS12, full-duration SS trial, formal speed benchmark or batch run was admitted.
No rejected state was used as an accepted continuation.
See docs/FAST_PZT_F10F_DEVELOPMENT.md and FAST_PZT_F10F_SMALL_STEPS_AUDIT.md.
