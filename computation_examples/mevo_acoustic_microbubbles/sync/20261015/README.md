# MEVO scientific synchronization — 2026-10-15

Branch: `sync/mevo-acoustic-f10s-20261015`.
Complete local delivery: `31a44d3f342e5d0bca55f7f9b83f4d32b8189696`.
Published parent: `467b17411deeb13a77370565b7382f071b41a774` (F10Q sync).

The exact F10R–F10S scientific delta contains **451 retained files / 172,676,842 bytes**, plus small synchronization scripts, metadata and progress/handoff documents. **8 rebuildable PDF raster previews / 1,257,965 bytes** are excluded. The actual Git transfer is compressed and deduplicated; these are filesystem payload sizes.

All 36 F10R/S original trajectories and 245 returned server files are retained, including parameter controls, intermediate accepted states, all-step contact/force evidence, independent source tests, timing, historical failures and rejected/NOT_RUN branches. The largest new file is the original F10R half-step phase trajectory at 9,892,427 bytes. It and the F10S 3 ms half-step trajectory are essential convergence evidence. Final PNG/PDF research figures, plotting inputs, review records and immutable execution-source archives are retained.

Current scientific code/configuration/test/dependency parity covers **988 files**, including inherited files. Original four-disc geometry, APECSS response and both acoustic caches are inherited from the previous published snapshot; they are not regenerated. Current Chinese progress is copied byte-for-byte.

No new movies, environments, downloaded packages, compiled caches or credentials are included. Earlier intentionally omitted bulk remains omitted. No local scientific files are deleted; no development-history branch, default-branch merge, new solver run or server task is pushed/started by this task. Other repository projects are unchanged.

## Verification and recovery

Run from `computation_examples/mevo_acoustic_microbubbles`:

```bash
python scripts/verify_sync_20261015.py
```

The verifier checks copied SHA-256, exact code parity, all new registered inputs, 245 returned server files, plot provenance, stored success/failure/admission states, 32-ID zero-flow trajectory headers, all new HDF5/NPZ contents, Python syntax, source archives and credential patterns. It does not execute scientific calculations.

- `EXPORT_POLICY.json`: scope, identities and size policy.
- `COPIED_FILES.json` / `EXCLUDED_FILES.json`: source locations, sizes, Git blobs and SHA-256; reasons for each omission.
- `SCIENTIFIC_CODE_MANIFEST.json`: exact current scientific code/config/test/dependency hashes.
- `PROGRESS_PROVENANCE.json`: live Chinese progress copy and identity.
- `PACKAGE_MANIFEST.json`: actual new/updated Git payload, excluding itself and verification to avoid circular hashes.
- `VERIFICATION.json`: verification results and package-manifest SHA.

Historical full-delivery manifests and their original hashes are unchanged; some list omitted preview PNGs. This synchronization manifest defines the actual payload. Excluded originals remain at their recorded source paths and can be restored by SHA for a requested full audit. Their PDF source files are retained. Older sync validators describe frozen past packages and may not match intentionally updated landing documents.

Scientific execution SHAs need not be reachable as GitHub commits because bulk development history is not uploaded. Current source is byte-identical; original execution tarballs and input hashes preserve historical reproduction. Retain **EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED / training_ready=false**.

Start with [current project state](../../docs/CURRENT_PROJECT_STATE_20261015.md).
