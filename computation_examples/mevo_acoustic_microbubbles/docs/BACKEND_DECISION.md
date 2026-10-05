# Backend verification, 2026-10-01

The server already ran an unrelated BraVa local P2 study. It was not stopped or modified.
The supplied endpoint matches the existing `vast4090` SSH alias; that alias supplies the user's existing key.
No server driver or global Python environment was changed. The project uses its own Python 3.12.3 venv.

The initial main venv had neither k-wave-python nor CuPy. Installed project versions are recorded in
SERVER_DOCTOR.json and requirements.lock.txt. CUDA driver/toolkit capability was 13.2; RTX 4090 has 24564 MiB.
Host RAM is about 61 GiB, but the cgroup limit is only 45156925440 bytes (~42.06 GiB); CPU quota is 7.68 cores.

The preferred C++/CUDA backend was actually attempted in the GPU smoke run. Its binary failed to load
`libhdf5_serial_hl.so.100`. The full subprocess error is retained in the run's CPP_FAILURE.json.
Installed k-wave-python 0.6.2 also warns that the unified C++ branch ignores sensor.record_start_index.
We therefore use the explicitly allowed Python/CuPy GPU fallback, without replacing system libraries or drivers.
The fallback's sensor output shape is asserted to equal ROI sensor count by exactly four CW periods.
The CPU smoke pressure peak was 405.963806 Pa; the GPU fallback was 405.963684 Pa.
These peaks are a smoke consistency check, not a full solver equivalence certificate.

Source geometry uses the installed, supported kWaveArray 3D disc API. Its BLI weights are combined in
explicit C mask order, matching the unified API. PML is outside the physical Cartesian grid, retaining
the requested physical domain. Expanded solver dimensions are recorded separately in run metadata.

References checked against installed source:

- https://k-wave-python.readthedocs.io/en/latest/get_started/new_api.html
- https://k-wave-python.readthedocs.io/en/latest/kwave.kspaceFirstOrder.html
- https://k-wave.org/documentation/kWaveArray.php
- https://k-wave.org/documentation/extractAmpPhase.php

Coherent integer-period Fourier demodulation implements the same single-frequency extraction idea,
with an explicit absolute-time phase correction and no window leakage. No complete 4D domain time history
is produced. Only the rectangular vessel ROI plus a gradient halo is sampled at the end.
