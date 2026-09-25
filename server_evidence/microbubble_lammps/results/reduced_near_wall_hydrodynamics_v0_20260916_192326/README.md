# Phase 1R review archive

Read `REDUCED_NEAR_WALL_HYDRODYNAMICS_V0_PHASE1R_REPORT.md`, `CURRENT_STATE.json`, and `validation/PHASE1_VALIDATION.json` first. The CF2003 reference and standalone C++ algebra pass; the combined resistance/shear physics status remains PARTIAL. Phase 2 is not authorized.

- `reference/`: frozen CF2003 artifact, checked page imagery, corrigendum and coefficient audit.
- `src/`: standalone kernel and generated frozen-reference constants.
- `tests/`: C++ A–K tests with a true scaled linear solve.
- `raw/`: all input/output CSVs, including the entire 10,009-point recovery grid.
- `validation/`: source, reference, C++ and independent NumPy results; inherited table-transition diagnostic.
- `visualization/`: six static scientific figures with CSV sources.
- `provenance/`: before/after upstream checks, build/run logs and local/remote SHA256 receipt.

This is a completed evidence archive. For a replay, copy it into a new isolated stage and update STAGE_PATHS.json before running any stage-writing scripts. Do not execute stage-writing audit runners in immutable historical directories.

The completed build used CMake Release and GCC 13.3 on Vast. `scripts/run_remote_phase1.py` records the standalone command and sources; it never invokes LAMMPS or Palabos. The independent validator and plots use existing `/usr/bin/python3` packages. Required PDF pages were rendered with the existing `wallref_pystokes` environment's PyMuPDF 1.24.10, without installation.

The earlier files with PENDING/UNVERIFIED in their names are acquisition history only. Current coefficient authority is the fully verified primary-PDF extraction CSV and `CF2003_FREE_SHEAR_REFERENCE_V0.json`.
