# Current-project Git synchronization · 2026-10-10

Branch: `sync/mevo-acoustic-f5-migration-20261010`.
Parent delivery: `66d4cebec18555a26bc9a8035cc3082756d0e76f`.
See [current state](../../docs/CURRENT_PROJECT_STATE_20261010.md) before using older
Stage 4 handoffs, absolute provenance paths or legacy remote launchers.

The parent commit already contains the delivered F1–F5 numerical code, configs,
tests, core inputs, results, and compact figures/videos. This update adds only the
remaining local visualization work and new-server migration/environment evidence.
It preserves original model code, parameters, trajectories and failed-attempt logs.

## Verification and contents

- `COPY_PROVENANCE.json`: source path, byte count and SHA-256 for each copied
  artifact; explicit metadata projection where unrelated filesystem inventory is omitted.
- `PACKAGE_MANIFEST.json`: hashes for the files actually introduced or updated
  by this sync, excluding the manifest itself and the final verification report.
- `VERIFICATION.json`: copy/hash, source preservation, structured-data, size and
  credential-pattern checks performed for this synchronization.
- `EXCLUDED_CURRENT_FILES.json`: per-file reasons for omitted current exports.
- `EXPORT_POLICY.json`: additional classes excluded from the sync.
- `SPARSE_CHECKOUT.txt`: optional patterns for the active project, selected
  results and new handoff, without checking out unrelated repository projects.

The migrated GPU checks are dated 2026-10-09 and copied with their original
status and source identity. This sync's verification is a new read-only packaging
check, not a rerun or extension of the scientific simulations.

## Deliberate limits

No virtual environments, managed Python binaries, downloaded wheels, private keys,
credential files, system fonts, CUDA caches, raw response-node cache, large event
streams or unrelated projects are added. Duplicate SVG/TIFF exports, previews,
decoded video frames and derived duplicate trajectory arrays are omitted from the
new visualization folders; original HDF5 sources, compact delivered videos,
PNG/PDF figures and essential metadata/logs remain. Existing tracked historical
files are inherited unchanged, not purged or re-uploaded as new artifacts.

Some original export inventories reference files intentionally excluded here.
Their bytes are retained as historical evidence and are not rewritten to make
the curated tree appear to be a full original archive. The exact included scope
is described by this directory's manifests. External Arial and other-project
style-reference videos remain external; no font/video license is assumed.

The complete local migration backup was deleted at the user's request on
2026-10-10. Statements about it in older reports describe their original dates.
