# Third-party licenses and use

| Dependency | Actual version/identity | License evidence | Integration |
|---|---|---|---|
| APECSS | v1.8, a0ee071ca247b5489e997d90e83b5f53c40e4ec7 | MPL-2.0, full notice `third_party/licenses/APECSS-MPL-2.0.txt` | Unmodified C executable subprocess; clone/source retained |
| PROTEUS | fc24fcbcf9623006e8f79989c0126c45a76bc23b | MIT, copyright2023 NathanBlanken; bundled voxelisation BSD3; vtkToolbox MPL2 | Audited only, no source copied/executed |
| k-wave-python | Installed0.6.2; latest stable verified0.6.2 | Installed distribution/LICENSE and repository LGPL3; README citation's GPL label differs | Python dependency, no upstream patch |
| k-Wave MATLAB | ucl-bug/k-wave repository | LGPL3, retained repository notice and source SHA | Reference only; no MATLAB runtime executed |
| k-Wave CUDA binary | Download metadata v1.4.1; banner v1.3/buildv3.6, git9ca8c89bc5b833adfafbef73d2e22122f7243bd4 | Repository notice LGPL3-or-later; copyright2012–2020 SC@FIT/Brno | Existing packaged binary; actual SHA in benchmark |
| HDF5 serial runtime | Ubuntu1.10.10+repack-3.1ubuntu4 | Both package copyright files retained | .deb extracted under ignored project build directory; no system library replacement |
| NumPy/SciPy/PyVista/VTK/CuPy/Matplotlib/h5py | Actual lock/server doctor | Upstream package notices remain in installed distributions | Runtime dependencies, no copied solver source |
| pypdf/fonttools/PyMuPDF | Local extraction tools only | Distribution notices retained in ignored pdf_tools; PyMuPDF is AGPL/commercial dual-licensed | Used to read/render papers, not included in runtime or redistributed as project code |
| Cattaneo2026 XLSX data | Zenodo15692129 | Dataset CC-BY4.0 metadata retained | Read-only contemporary cross-check, not fitted to project outputs |

The project wrappers call external implementations; they do not replace license notices or claim ownership of them.
The C++ compatibility adapter calls the installed upstream `scale_pressure_source_uniform_grid` function;
no copied implementation or patched site-package is distributed. Conversion and index mapping are project code.
Paper PDFs/HTML are scientific reading material in ignored downloads, not code under this project's license.
When redistributing binaries/dependencies, include the respective source/notices and satisfy their upstream terms.
Exact source URLs/SHAs for retrieved license files are beside the notices in `third_party/licenses/`.
