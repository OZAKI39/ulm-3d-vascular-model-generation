# GitHub synchronization manifest

This directory contains the formal `anchor_003274` Ultraliser smoke and final
run selected for repository synchronization. The help-only diagnostic `r1`
and all unrelated datasets/runs are intentionally excluded.

## Console logs

Ultraliser emits a progress line for a very large number of operations. The
uncompressed stdout files therefore exceed GitHub's per-file size limit. They
remain unchanged in the local run directory and are ignored by Git; the
losslessly compressed streams below are tracked instead.

| Run | Tracked file | Original bytes | Original SHA-256 | Compressed SHA-256 |
| --- | --- | ---: | --- | --- |
| smoke | `smoke/stdout.log.gz` | 248568413 | `ff738fdcd46fdbe83fbec8794ebd651bc3c5f49389d086005bb195840b2d1906` | `22577a56b137a6b0593a0f295933d7d0b527ac0f268dcf971b43b33e53032df6` |
| final | `final/raw_ultraliser/stdout.log.gz` | 791183232 | `9760a89a55c803a2a4e7a4d467ecce5bd7c72765921ef9cac0139108fa96a0d7` | `5414491bc3264820b30b533ed05ac0e91c8bc1cafc72241bb5fcfbbbc851caa1` |

Both archives were decompressed and checked byte-for-byte by SHA-256 before
commit. Use `gzip -dk stdout.log.gz` to restore each original stream.

## Included material

- canonical ROI SWC, CUT_PORT mapping, metadata, and official H5 adapter;
- Ultraliser version/help and the exact smoke/final commands;
- raw official meshes and mesh statistics from both permitted runs;
- final micrometre VTP/STL and metre STL;
- topology, resolution, radius-fidelity, and provenance QC;
- exactly four validation figures and the final decision report.

The upstream Ultraliser source is referenced as a submodule pinned to commit
`3e4b0eee685adbf513e40720a68fd92e66a34b44`; local build products and the
paper PDF are not vendored.
