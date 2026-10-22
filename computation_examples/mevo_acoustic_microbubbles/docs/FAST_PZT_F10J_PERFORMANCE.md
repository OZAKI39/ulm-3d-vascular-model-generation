# F10J fixed-state performance

**Stop GPU B2 production integration for this B=1, N=32 workload.** Numerical admission passed first; three matched measurements show optimized GPU 22.048 ms versus CPU 11.867 ms. GPU cost is 1.858 times CPU. This is a measured lack of advantage, not an extrapolated long-trajectory speed ratio.

Hardware: RTX 4090, Ryzen 9 7950X3D; Python 3.12.3, Warp 1.18.0+cu12, NumPy 2.2.6, SciPy 1.15.3. Existing isolated environment and driver unchanged. OPENBLAS/OMP/MKL/NUMEXPR threads=1 via the inherited wrapper. GPU ownership checked before each job; process monitoring archived. No unrelated GPU job was interrupted.

Same saved x/free, 10 us, R, viscosity, kc, discrete A, geometric forces, original quality budgets. CPU uses unchanged `geometric_contact` source with atol=1e-15, which additionally satisfies the original GPU 1e-14 gradient-bound target (verified), not merely its looser original CPU residual stop. Both methods meet the same independent position/force/geometry gates. CPU matrix assembly is **inside** timing. GPU coefficient/normal preparation, validation and complete final vector/operator/residual readback are inside. Same independent quality checks and same JSON byte payload fields are added to both. Disk evidence writing is outside both per-solve timers. Reference state loading is outside both. No matched output is skipped.

One warmup per method, then exactly three rotated-order pairs. First compile and reusable graph construction excluded from resident cost and recorded separately. One additional instrumented profile per GPU version is diagnostic only. Original F10I still fails and is not a valid production baseline. Historical 35.457 ms is not reused as a new result.

|Method|Full median ms|Range ms|Prep+solve median ms|Independent quality ms|Matched output ms|
|---|---:|---:|---:|---:|---:|
|original|22.789198|22.716498–24.268950|22.371730|0.174589|0.236709|
|corrected|42.255502|39.056687–42.413241|41.840444|0.171979|0.241709|
|optimized|22.047892|21.798132–23.477144|21.634534|0.177209|0.242389|
|cpu|11.867071|11.847171–12.165530|11.476153|0.152009|0.241039|

Independent component medians do not necessarily sum to the median whole measurement. All raw paired timing objects include CPU assembly/contact and GPU buffer preparation/nonlinear-loop/scalar waits/result readback. Construction records include buffer allocation, inherited graph build and replacement graph build; those inherited graphs remain part of this minimal subclass prototype's startup cost.

|Diagnostic profile|PCG graph device spans ms|Armijo evaluation wall ms|Scalar wait/read wall ms|
|---|---:|---:|---:|
|original|18.102272|2.321169|20.695629|
|corrected|36.072448|0.712787|37.384156|
|optimized|19.232768|0.711197|20.125436|

**These categories overlap.** In particular scalar read wall time includes waiting for earlier PCG work. It is not all PCIe copying or CPU overhead that can disappear by deleting a read. The solve-wide CUDA event span contains host decision gaps and is not pure kernel time. PCG-only events enclose graph execution without host decision gaps, but still include device graph/kernel scheduling overhead. The measured PCG work dominates this very small 96-DOF GPU solve; reducing 35 to 12 intermediate reads is insufficient.

New processes (one each, existing compiler cache, new CUDA context/buffers/graphs): original 0.650504 s, corrected 0.746669 s, optimized 0.661756 s, cpu 0.445735 s. Includes common harness imports, input loading, one solve/check/output, artifact write and shutdown. CPU imports the common Warp-containing harness but initializes no CUDA context. These are diagnostic startup costs, not a cold compiler benchmark or repeated throughput statistic.

All six registered remote wrappers total 11.082999 s; this is process wall time, not all development time or pure GPU occupancy. Independent child costs are already included in process_v1 and are not added again. No B64, 10 ms trajectory, AI dataset, pressure scan, new CFD or new acoustic solve was run.

Keep opt-in code for reproducibility and the CPU B2 as a finite-state numerical research reference. Do not continue production optimization merely because the old GPU failed or because scalar-read counts improved. A future architecture change would require a separate task, new admission and demonstrated workload need.
