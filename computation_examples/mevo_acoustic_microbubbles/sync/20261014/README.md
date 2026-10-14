# MEVO scientific synchronization — 2026-10-14

Branch: `sync/mevo-acoustic-f10q-20261014`.
Complete local source: `41ae8ed2b88f2cbbed834fffccaa0b5de6e91589`.
Published parent: `15c42e0e422ff5d1521e5cc8c12ede2688b3de59` (F10K sync).

This adds completed F10L–F10Q development and the subsequent four-PZT/per-bubble
force summary. Scientific source, results, thresholds and conclusions are copied
unchanged. No simulation, server job, history rewrite or default-branch merge is
performed. The uncurated local development commits are not made ancestors.

Copied scientific delta: **1,136 files / 199,048,839 bytes**, excluding small new
sync scripts, documentation, inventories and the progress copy. Omitted delta:
**50 files / 27,392,551 bytes**: replaced visual previews, three per-element F10L
intermediate fields (combined original basis retained), and the superseded
pre-port F10M material grid. Original locations, sizes and SHA-256 are recorded.
Previously published files are inherited; prior deliberately omitted bulk remains
omitted. Filesystem sizes differ from compressed/deduplicated Git transfer sizes.

The 24.6 MB fine water-field cache and 13.9 MB Y half-step trajectory are retained
because they are essential convergence evidence, including failed shape/ID
admission. All F10Q returned server files, all F10N–Q trajectories, controls,
failures, response/cache inputs, source execution archives, final PNG/PDFs and
numerical logs are retained. No new movies, environments or package downloads
are added.

## Evidence and verification

- `EXPORT_POLICY.json`: identities, scope, size rule and essential larger files.
- `COPIED_FILES.json`: exact copied scientific files, Git blobs and SHA-256.
- `EXCLUDED_FILES.json`: omissions, reasons and original recovery paths.
- `SCIENTIFIC_CODE_MANIFEST.json`: 950 current code/config/test/validation/
  dependency files, including inherited files.
- `PROGRESS_PROVENANCE.json`: unchanged live Chinese progress-document copy.
- `PACKAGE_MANIFEST.json`: actual new/updated payload, excluding this manifest
  and verification report to avoid circular hashes.
- `VERIFICATION.json`: byte identity, input/figure provenance, dataset reads,
  archive inspection, Python parsing and credential-pattern scan.

From `computation_examples/mevo_acoustic_microbubbles`, with Python, NumPy and h5py:

```bash
python scripts/verify_sync_20261014.py
```

This verifier does not invoke a solver. Historical full-delivery manifests remain
unchanged and can reference omitted previews/intermediates; they do not claim
that the curated branch contains every original artifact. Older sync validators
describe frozen past packages and may not match intentionally updated files.

One original F10O record, `results/straight_tube_exp1/f10o/cleanup/F10N_VERIFY.json`,
contains a Python dictionary representation despite its extension. It is preserved
byte-for-byte and safely parsed with `ast.literal_eval`, explicitly recorded by
the verifier. Original Markdown whitespace and CSV line endings are preserved too.

## Recovery

Excluded originals remain at paths listed in `EXCLUDED_FILES.json`. Verify their
SHA-256 before restoring them for a separately requested full audit. Intermediate
F10L fields can also be regenerated with retained scripts, without assuming
deterministic byte identity. The old full migration backup deleted at the user's
request is not a recovery source.

Scientific source/execution SHAs are provenance identities and need not be
GitHub-reachable commits. Current source is copied exactly; original immutable
source-only tarballs are under `results/straight_tube_exp1/*/execution/`, with
their input/control and execution records. No bulky source-history push is needed.

See [current project state](../../docs/CURRENT_PROJECT_STATE_20261014.md).
Retain **EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY /
PHYSICS_UNCALIBRATED / training_ready=false**. Simulated target success and
engineering checks do not establish experimental calibration.
