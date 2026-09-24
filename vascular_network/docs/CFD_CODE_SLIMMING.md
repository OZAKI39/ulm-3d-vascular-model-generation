# CFD code slimming

The research branch accumulated many one-shot solver runners, forensic
harnesses, duplicate constants, and tests tied to experiments that are now
frozen evidence. This phase reduced the active context needed for production
integration without changing `cfd_flow.py`, `configs/cfd_flow.yaml`, production
Lua rendering, production QC semantics, or any scientific artifact.

## What changed

- Added one small, I/O-free validated Tau=1 contract and one compact evidence
  loader.
- Consolidated reusable restart decoding, Full V2 replay, V3/MLS physical flux,
  steady auditing, field reconstruction, boundary parsing, and grid-analysis
  math into readable solver-free modules.
- Reduced each historical bug to one canonical regression plus immutable
  evidence.
- Removed 127 superseded paths: 33 one-shot research modules, 34 research CLI
  wrappers, and 60 duplicated/runner-bound tests.
- Kept three reusable research utilities outside the default production
  context: dimensionless geometry, solver-free three-grid math, and analytic
  wall-q geometry.
- Kept all tracked historical evidence and all critical binaries; no output
  directory or scientific payload was deleted.

## Frozen conclusions

Base Tau=1 is accepted, Coarse-to-Base sensitivity is partial two-grid evidence,
Fine steady is incomplete because of the resource budget, Full V2 is the only
final conservation identity, and physical Q uses V3 planes with continuous
aperture Gauss/MLS V2. The adaptive inlet denominator uses 375 active
populations, not 376. Formal three-grid GCI is not complete.

## Audit and recovery

The QC run at
`outputs/cfd_flow/code_slimming_for_production_20260902_002647/qc` contains the
before/after import graphs, context metrics, full classification inventory,
retired-component map, evidence hashes, static production snapshots, and the
machine-readable integration surface. Production snapshots and all critical
evidence hashes must report `PASS`.

Deleted source is recoverable from Git at starting commit
`9e66e902b70e15cff17441c915ef4990ffe327d1`; use normal `git show`/history for
forensic retrieval. Do not restore retired code as a second active scientific
implementation.

For the next phase, read `CFD_PRODUCTION_INTEGRATION_CONTEXT.md` first and then
follow `context_priority_order` in `production_integration_surface.json`.
