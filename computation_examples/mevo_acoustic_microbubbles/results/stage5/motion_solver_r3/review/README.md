# Stage 5A-R3 lightweight scientific review

This directory contains real full-system B2 states, candidate QP/geometry arrays,
the actual open WALL mesh, and complete losslessly compressed original accepted
and rejected ledgers for both R2 failed runs. No particles are removed from the
fixtures. Local repaired trajectories start at the actual R2 terminal state;
they must not be mistaken for the separately recorded original-3-s formal runs.

Read the project docs/STAGE5A_R3_ROOT_CAUSE_ANALYSIS.md and final gate first.
The finite-feature geometric repair recovers B2 continuation. Independent
all-active global trajectory accuracy is not established; no promotion follows.
The early +30..50 us fixtures use distinct actual nominal/refined states with
explicit cold caches. A local difference does not identify either as truth.

REVIEW_SOURCES.json supplies hashes and exact copy/compression transformations.
SOURCE_DIFF_FROM_R2.patch records all project source changes through its named
parent commit. The containing Git archive also includes later packaging scripts.
The branch contains the required audit JSON, PNG/PDF figures and actual-state
NPZ data. Large trajectories, videos and every formal failed-candidate array are
indexed with actual paths/hashes by final/ARTIFACT_MANIFEST.json. No old failure
or gate is overwritten. V1 and V2 execution identities remain separate.

To read a retained original ledger: gzip -dc original_B2/nominal/REJECTED.jsonl.gz
For executable reproduction commands, see scripts/stage5a_r3/README.md.
