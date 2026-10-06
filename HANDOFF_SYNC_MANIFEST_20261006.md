# PD2 boundary handoff sync manifest — 2026-10-06

- Source branch: `clot-clearance/pfpd-progress-20261005`
- Source HEAD: `695e75166aa864809460ea9a1845fe72abce229c`
- New branch: `clot-clearance/pd2-boundary-handoff-20261006`
- Initial state recorded at: `2026-10-06T14:27:34.720324+00:00`
- Scientific status: **HIGHER_ORDER_BOUNDARY_MECHANICS=VERIFIED_FOR_INTACT_BENCHMARKS**.
- Boundary reintegration: **PASSED**; four manufactured recoveries plus exactly one classical-shear solve were already complete before handoff.
- Files added: **727**; updated: **2**; total added/updated: **729**; deleted: **0**.
- Total new staged content: **25099842 bytes (23.937075 MiB)**. This is the sum of complete added/updated blobs, including handoff metadata; Git pack compression and pre-existing ancestry are not counted.
- Largest synchronized file: `clot_clearance/results/higher_order_baqp_boundary_reintegration_v1/figures/c1_vs_pd2_boundary_field_quality.png`, **1145937 bytes**.
- Recommended next scientific task: **CURVED_CLOT_WALL_INTACT_MECHANICS_V1**, not started.

## Scope and inventories

The source working tree included untracked completed research. This branch was created directly from its local HEAD and selectively adds that state. There was no reset to an old remote, cleanup, scientific edit, historical deletion, history rewrite, PR, or new scientific campaign.

[Complete proposed staging paths](github_handoff_20261006/PROPOSED_STAGING_PATHS.txt) · [per-file staging size/type/reason](github_handoff_20261006/STAGING_INVENTORY.csv) · [scientific include/exclude inventory with SHA256](github_handoff_20261006/SCIENTIFIC_FILE_INVENTORY.csv). The staged Git index is audited independently before commit.

Source trees: `clot_clearance/pd_clot/pfpd_candidate/` including `higher_order_rk/` and `nonlocal_boundary/`; current related `scripts/`, `tests/`, and analysis source/data. Existing material, spectral split, geometry, regularization solver, configurations and the two-file `vendor/` helper remain inherited. No clone is copied.

Included verification stages:

- `clot_clearance/verification/baqp_near_boundary_kinematic_energy_audit_v1/` — final reports/configuration/compact evidence; raw-file exclusions below.
- `clot_clearance/verification/higher_order_baqp_boundary_reintegration_v1/` — final reports/configuration/compact evidence; raw-file exclusions below.
- `clot_clearance/verification/higher_order_baqp_energy_force_candidate_v1/` — final reports/configuration/compact evidence; raw-file exclusions below.
- `clot_clearance/verification/higher_order_baqp_rk_oracle_v1/` — final reports/configuration/compact evidence; raw-file exclusions below.
- `clot_clearance/verification/nonlocal_boundary_formulation_v1/` — final reports/configuration/compact evidence; raw-file exclusions below.
- `clot_clearance/verification/nonlocal_boundary_nonaffine_localization_v1/` — final reports/configuration/compact evidence; raw-file exclusions below.
- `clot_clearance/verification/nonlocal_boundary_nonaffine_shear_localization_v1/` — final reports/configuration/compact evidence; raw-file exclusions below.
- `clot_clearance/verification/pfpd_baqp_affine_coefficient_reference_audit_v1/` — final reports/configuration/compact evidence; raw-file exclusions below.
- `clot_clearance/verification/pfpd_baqp_finite_domain_boundary_consistency_v1/` — final reports/configuration/compact evidence; raw-file exclusions below.
- `clot_clearance/verification/pfpd_baqp_fully_closed_force_stencil_audit_v1/` — final reports/configuration/compact evidence; raw-file exclusions below.

## Representative figures

- `clot_clearance/results/higher_order_baqp_boundary_reintegration_v1/figures/c1_vs_pd2_boundary_field_quality.png`
- `clot_clearance/results/higher_order_baqp_boundary_reintegration_v1/figures/pd2_manufactured_recovery.png`
- `clot_clearance/results/higher_order_baqp_energy_force_candidate_v1/figures/c1_variational_force_regression.png`
- `clot_clearance/results/higher_order_baqp_energy_force_candidate_v1/figures/pd2_energy_force_directional_derivative.png`
- `clot_clearance/results/higher_order_baqp_rk_oracle_v1/figures/c1_vs_pd2_boundary_bias.png`
- `clot_clearance/results/higher_order_baqp_rk_oracle_v1/figures/pd2_polynomial_reproduction.png`

Only these six PNGs are newly synchronized: two per current stage, all under 3 MiB; no new PDFs or redundant raster variants. Their small JSON/CSV data and source plotting scripts are included. The polynomial heatmap has an additional exact 36-row table in `github_handoff_20261006/PD2_POLYNOMIAL_FIGURE_DATA.csv`; extracting this table did not rerun the oracle.

## Local-only artifacts and size policy

**305 relevant files, 398236550 bytes**, remain local. Each path, exact byte size, SHA256, exclusion reason and compact replacement is listed in [LOCAL_ARTIFACTS_NOT_SYNCED.md](LOCAL_ARTIFACTS_NOT_SYNCED.md), with machine-readable records in [EXCLUDED_RELEVANT_FILES.json](github_handoff_20261006/EXCLUDED_RELEVANT_FILES.json). This linked per-file manifest is part of this sync manifest.

Excluded categories: large N19 reference caches and geometry; large/repeated node/bond field arrays; duplicate solve inputs and historical state copies with retained metrics; detailed large per-node tables with compact summaries; non-selected historical/duplicate/high-resolution figures; paper PDFs and extracted full texts; superseded implementation-check output. Unrelated raw data, environments, cloned repositories, compiler products, old recovery indexes and scratch caches were outside selection and remain untouched.

No staged file exceeds 5 MiB; no staged raw array exceeds 2 MiB; total content is below 50 MiB; no large-file exception is needed. Twenty-one compact NPZ files are selected only as evidence/tests: sixteen energy-force fixtures (including a 613108-byte N7 reference-operator fixture needed by unchanged tests) and five final boundary accepted states. The small N7 file is intentional fixed test data, not a disposable runtime cache; no N19 cache is included. Exact quantities are checked from the staging inventory.

The repository-wide NPZ ignore rule was inspected and left unchanged. Only the explicitly listed NPZ paths are force-added individually through the allowlist; this does not widen ignore rules or stage a directory blindly.

Original scientific manifests describe full server outputs and are preserved byte-for-byte. Excluded files may still be named in those manifests and linked in original reports. Their absence on GitHub is declared here, and does not change the original science or assertions. Existing historical tracked binaries/papers in ancestry were not removed.

## Tests, integrity and limitations

- No active project writer was found before creating the branch; all five boundary solve outputs and final report were complete. Boundary FINAL_STATUS and FINAL_ARTIFACT_QA report success.
- Existing focused suites on server: **36 passed in 4.42 s** (RK 12, energy-force 12, boundary 12). [Log](github_handoff_20261006/SERVER_TARGETED_TESTS.log). No full pytest or simulation campaign.
- Isolated copy containing only selected/inherited required files: **20 passed** (11 energy-force, 9 boundary; 4 explicitly deselected). [Log](github_handoff_20261006/PORTABLE_TARGETED_TESTS.log). First copy omitted an already tracked vendor helper; adding that existing dependency to the copy resolved the import failure without editing tests or scientific code. No new branch dependency was missing.
- Full 36-test replay and all original plots require listed local-only artifacts. The portable subset checks live small mechanics and saved-result gates; it does not rerun N19 equilibria or certify every excluded array.
- All 1204 files frozen before handoff retain SHA256 identity, including omitted scientific artifacts; the expanded inventory additionally includes two inherited vendor files. No physics, production source, old scientific result, old test assertion, old manifest, or .gitignore is modified.
- New handoff documentation local links resolve within the selected/ancestral branch files. Historical reports retain original local-only links. Source import inspection and the isolated focused tests check required local module coverage.
- Credential-pattern scan found no matches; credentials and tokens were not printed. No environments, paper PDFs, clones or build artifacts are newly staged.

See [HANDOFF_CONTEXT_20261006.md](HANDOFF_CONTEXT_20261006.md) for scientific caveats and next-task scope; [HANDOFF_ENVIRONMENT.md](HANDOFF_ENVIRONMENT.md) for versions and commands. Current context was updated from accepted reviews; pre-handoff active context is archived for provenance and the historical index text remains intact.

## Commit and remote verification procedure

Create one handoff commit, push only the new branch to origin without force, then compare local HEAD with `git ls-remote origin refs/heads/clot-clearance/pd2-boundary-handoff-20261006`. Recheck branch nonexistence immediately before push; stop if it has appeared. Compare the three protected remote refs with the recorded initial values. Post-push command evidence is kept locally under `github_handoff_20261006/` to avoid a second self-referential commit; the final response provides the verified SHA. This pre-commit manifest does not pretend the push had already occurred.

Deliverable URL: [handoff branch](https://github.com/OZAKI39/ulm-3d-vascular-model-generation/tree/clot-clearance/pd2-boundary-handoff-20261006). Stop after remote verification; no PR or further science.
