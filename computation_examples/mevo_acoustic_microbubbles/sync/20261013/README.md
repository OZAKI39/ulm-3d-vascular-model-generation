# MEVO scientific synchronization — 2026-10-13

Branch: `sync/mevo-acoustic-f10k-20261013`.
Complete scientific source: `850fd2df3225351f4fe9564f39cbe72a38662b30`.
Published parent: `c8ecf9deece1bb0e86b0d0f1c574b984a30ebbd0`.

This is a read-only scientific packaging task. No solver, control, quality
threshold, trajectory or prior scientific conclusion was changed or rerun.
The new commit is based directly on the published synchronization so that the
intervening bulky local development results are not uploaded through Git history.

Copied from F10K: 2,951 new/updated files, 149,984,961 bytes. Omitted delta:
475 files, 1,862,197,578 bytes. These figures exclude the small new sync documentation,
manifests, research-progress copy, and five source-only execution archives
(about4.3 MB total). The already-public parent contains about844.5 MB of project
files, inherited without purging; this is not an estimate of additional upload size.
Git transport compression/deduplication further differs from filesystem size.

## Evidence inventory

- `EXPORT_POLICY.json`: exact scope, size rules, selected original HDF5 and source commit.
- `COPIED_FILES.json`: every copied delta artifact, byte count, source path, Git blob and SHA-256.
- `EXCLUDED_FILES.json`: every omitted delta artifact, reason, local recovery path and original hash.
- `EXECUTION_SNAPSHOTS.json`: exact source tarballs for the five F10K immutable execution commits;
  file hashes are cross-checked against the original remote source manifest.
- `PROGRESS_PROVENANCE.json`: unchanged live Chinese progress-document copy.
- `PACKAGE_MANIFEST.json`: actual new/updated files in this Git package, excluding
  this manifest and its verification report to avoid circular hashes.
- `VERIFICATION.json`: actual packaging checks, including byte identity, dataset
  readability, scientific status preservation, source archives and credential-pattern scan.

The 12 explicitly selected original HDF5 files contain F10A OFF/S/O10/5 µs,
F10D prepare1101/1102, historical F10H B2 dt10/5 µs, and both F10K B2 windows.
Additional tiny historical HDF5 files are kept when below256 kB. Compact full
CPU/GPU benchmark records, latest arrays, failure stderr and recovery records,
all current source/configuration/reports and publication PNG/PDFs are retained.
No new movie, environment, downloaded package, raw radial waveform or full
per-step contact-matrix stream is added. Original large files remain locally.

## Verification and recovery

From `computation_examples/mevo_acoustic_microbubbles`:

```bash
python scripts/verify_sync_20261013.py
```

The packaging verifier requires Python with h5py/NumPy. It does not invoke GPU,
CFD, acoustics or trajectory integration. It verifies this committed subset;
historical full-delivery validators still require their complete original files.
The original full manifests and scientific gate bytes are not rewritten to make
the curated subset look complete. Missing full raw audits are explicitly listed.

Recovery sources are the recorded complete local archive paths and retained
scientific server stage directories. The original local commits are provenance
identities and need not be GitHub-reachable on this curated branch. For F10K's
five numerical execution commits, exact source archives are supplied here.
Check each omitted file's SHA-256 before restoring it for a separately requested
full audit; this sync does not start that audit or resume old tasks.
