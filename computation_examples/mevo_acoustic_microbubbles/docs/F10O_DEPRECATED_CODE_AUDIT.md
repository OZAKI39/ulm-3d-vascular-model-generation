# F10O bounded dead-code audit

## Actual deletion

Removed only the unused Warp device conversion function f64 (and its decorator/blank lines) from src/mevo_mb/fast_sim/vessel_kernels.py. It converted vec3 to vec3d, but no kernel, test, CLI, configuration, registry or dynamic loader selected or called it. This is a deliberately small deletion: no evidence supports wholesale removal of historical solvers.

Evidence: targeted AST definition/name inventory across tracked source/script/test files identified this sole unique definition in the active fast_sim/straight_tube scope. Reference checks cover tracked source, scripts, tests, configs, docs, pyproject and dynamic import/getattr/export patterns. All actual vessel kernel consumers use explicit attributes or explicit imports; none enumerates functions or depends on registering this helper. Removing an uncalled Warp function does not alter a launched kernel's arithmetic. The equivalent built-in vec3d constructor remains; no scientific algorithm or sole reproduction capability was removed. Historical original code is unchanged in the separate F10N and F4 worktrees and Git commit 15e8ce4de.

## Retained deliberately

- F10G/F10H CPU contact and lubrication references; F10I/J GPU diagnostics; F10K validation: unique historical comparison functions.
- Original Y-vessel and three-PZT modules: imported by inheritance and historical runners; inactive geometry is bypassed explicitly by the existing R3/F10N adapter, not deleted.
- Historical CLI scripts with no inbound Python imports: valid executable entrypoints, not dead functions.
- Marmottant/Doinikov refusal stubs: tested fail-closed interfaces, not dead code.
- All reports, fields, APECSS, raw trajectories, failures, execution snapshots, hash manifests and drawing inputs.
- No unused import findings from Ruff F401/F811/F841; no invented deletion of valid algorithms.

## Verification and diff

results/straight_tube_exp1/f10o/cleanup/ records baseline Git identity, exact source diff, reference search, original F10N verification and lint checks. Run full read-only F10N verifier before any new scientific work. The cleanup commit is a separate prerequisite to the B2 bridge. Scientific data in F10N is untouched.

This cleanup reduces a proven dead helper only; the concise active map supplies most of the reduction in development context. It is not a claim that hundreds of historical files are obsolete.
