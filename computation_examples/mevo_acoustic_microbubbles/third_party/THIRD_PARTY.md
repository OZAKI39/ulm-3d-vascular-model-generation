# Audited dependencies

## APECSS

- Official source: https://github.com/polycfd/apecss ; paper DOI https://doi.org/10.21105/joss.05435.
- Latest release verified against GitHub releases/tags API on 2026-10-01: **v1.8**. Pinned commit `a0ee071ca247b5489e997d90e83b5f53c40e4ec7`, release source dated 2025-07-21. Development HEAD at discovery was `41392784d3c7196e6ee877e9b61055da0acbcb0c`; not used.
- Full clone retained under `third_party/src/apecss`, detached at v1.8, ignored by project Git. API discovery response retained in `APECSS_DISCOVERY.json`. No upstream patches.
- Actual remote build: gcc 13.3.0, CMake Release `-O3 -Wall -Werror -march=native -DNDEBUG`, double (`sizeof(APECSS_FLOAT)=8`); no fast-math. Out-of-tree static library and official executables under `third_party/build/apecss`. Native binaries must be rebuilt on a different CPU architecture.
- Run `scripts/stage2/apecss_build_validate.py` from project root. Actual passing official run: `results/stage2/apecss_validation/20261001T195618Z`. Seven ultrasound cases and Rayleigh emissions case executed; eight upstream plot scripts completed. Four supplied reference arrays pass upstream tolerances (rtol=1e-6, atol=3e-16), relative L2 2.74e-14 to 4.05e-12. The lipid example has no digitized experimental reference array; its execution is not an experiment reproduction.
- First setup attempts failed because a source sync was still in progress, then an assumed example filename was wrong. Logs retained; corrected before the passing run. These are orchestration failures, not hidden physical validations.
- Integration: subprocess around the unmodified official ultrasound executable; no independent Python radial solver. Pressure phase correction and exact shell-law scope documented in `docs/APECSS_MODEL_AUDIT.md`.
- License: MPL-2.0, full notice in `licenses/APECSS-MPL-2.0.txt`. Copyright The APECSS Developers. Do not omit source availability/license obligations when redistributing binaries.

## PROTEUS

- Official source https://github.com/PROTEUS-SIM/PROTEUS ; DOI https://doi.org/10.1109/TUFFC.2024.3427850.
- Audited commit `fc24fcbcf9623006e8f79989c0126c45a76bc23b`, 2024-07-21. Source clone retained under `third_party/src/proteus`; source inspection only, not a MATLAB execution claim.
- Root MIT copyright (c) 2023 NathanBlanken. Bundled mesh voxelisation has BSD-3-Clause license; bundled vtkToolbox has MPL-2.0. All three notices retained in `licenses/`.
- No PROTEUS source copied into runtime. Component decisions and source-level limitations: `docs/PROTEUS_REUSE_AUDIT.md`.

## Retrieved literature and analysis tools

Open original manuscripts/PMC pages are retained in ignored `third_party/downloads/` for scientific reference, not represented as project-owned content. Extracted parameter values and digitization provenance are stored separately in configurations. pypdf 6.19.0 and fonttools 4.66.1 were installed only into ignored `third_party/build/pdf_tools` for local PDF extraction; not production dependencies. Upstream notices remain in their distributions.

The actual k-Wave installation and its backend dependencies are recorded separately by the server doctor and backend benchmark. A package version constraint is not evidence of an installed version.

## k-Wave and runtime repair

Installed/latest stable k-wave-python0.6.2; private Ubuntu HDF5 packages and the tested version-specific C++ input/output adapter are documented in `docs/KWAVE_BACKEND_BENCHMARK.md`. LGPL notices and actual installed license files are retained. No upstream files were patched. PyMuPDF was additionally used only for local paper equation rendering; its distribution license remains with the ignored reading-tool installation.
