# Stage SV1.3H native probes

`cuda_kernel_smoke.cu` was compiled and run three times on the RTX 4090 with CUDA 12.6.3. The native executable is mirrored under `outputs/sv1_3h/native/`; runtime and library evidence is in `reports/sv1_3h/cuda12_runtime.json`. Process time includes CUDA initialization and is not a vascular-flow benchmark.

`petsc_cuda_smoke.c` is an unchanged copy of the prepared SV1.3G sparse-system test. It was not compiled or executed in SV1.3H because the original PETSc 3.19.6 CUDA12 build failed. CPU1 / CPU4 / GPU1 benchmark measurements do not exist; `status.json` explicitly records NOT_RUN.

The CUDA environment wrapper and native binary use the remote server's installation paths. No CUDA toolkit was installed into WSL.
