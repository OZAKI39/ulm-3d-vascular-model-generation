# k-Wave backend benchmark

Installed and latest stable verified: **k-wave-python0.6.2** (GitHub release and PyPI checked2026-10-01).
Python/CuPy and unified `backend="cpp", device="gpu"` were both actually run on RTX4090.
Packaged binary metadata is v1.4.1; executable banner says v1.3/buildv3.6, May17,2026,
git9ca8c89bc5b833adfafbef73d2e22122f7243bd4. The actual binary SHA is in KWAVE_STATE.json.

## Dependency and interface repair

`ldd` found missing `libhdf5_serial_hl.so.100` and `libhdf5_serial.so.103`.
Ubuntu24.04 packages `libhdf5-hl-100t64` and `libhdf5-103-1t64`, version1.10.10+repack-3.1ubuntu4,
were downloaded and extracted under `third_party/build/kwave_runtime`. No incompatible SONAME symlink,
global package replacement or site-package patch was used. The benchmark sets process-local LD_LIBRARY_PATH.

Actual installed source inspection and failed comparisons exposed additional0.6.2 interface defects:
raw additive pressure is written without the density-source scaling required by the binary;
source masks use F order while project source rows use C order; output uses F-order sensor rows and time×sensor layout;
record_start_index is ignored. The project adapter calls the existing upstream
`scale_pressure_source_uniform_grid`, maps mask row order, transposes the output and crops the requested absolute-time window.
It is guarded to0.6.2 and isotropic homogeneous3D sources. Python/CuPy calculation is unchanged.

## Same-case results

Run: `results/stage2/performance/kwave_20261001T205318768781Z`.
Identical disc/source/medium/CFL/PML/sensor/CW controls;1MHz;32³ and64³ physical grids.
Only new backend fixtures were computed; existing1/2MHz production basis caches were reused untouched.

| Grid | CuPy seconds | C++ CUDA seconds | Pressure relativeL2 | Gradient relativeL2 | Sampled VRAM CuPy/C++MiB |
|---|---:|---:|---:|---:|---:|
|32³|0.992|0.578|3.265e-7|4.114e-7|434/434|
|64³|1.142|0.635|4.374e-7|5.341e-7|578/1029|

These are single measured end-to-end runs, including initialization and C++ file I/O. They do not establish a universal1.7–1.8×speedup.
VRAM is0.5s whole-device sampling, not isolated exact peak allocation. C++ full-time recording may substantially increase RAM/I/O
for the production ROI; no production-grid C++ benchmark is claimed. **CuPy remains the default.**

Reproduce: `OPENBLAS_NUM_THREADS=1 .venv/bin/python scripts/stage2/kwave_benchmark.py` from server project root.
Failed timestamped runs are retained, including the≈1.43e7 unscaled-pressure discrepancy.
See [official release](https://github.com/waltsims/k-wave-python/releases/tag/v0.6.2),
[official C++ documentation](https://www.k-wave.org/documentation/example_cpp_running_simulations.php)
and `third_party/KWAVE_DISCOVERY.json` for source/version evidence.
