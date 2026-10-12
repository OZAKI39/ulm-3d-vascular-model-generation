# F10G reproduction

Use the existing isolated project environment; NumPy/SciPy/h5py/Matplotlib suffice
for offline analysis. Warp/CUDA is required only for the GPU integration adapter
and the original-kernel microbenchmarks. Set OPENBLAS_NUM_THREADS=1 and PYTHONPATH=src:.

Completed fixed runs (do not rerun them to regenerate figures):

1. baseline: verify original F10F artifacts, original sources and 13 input hashes.
2. lubrication --label <unique>: saved same-state 2/4/8/Cholesky+sparse comparison.
3. prototype --label <unique>: two/four independent references, conditional 32-state proposals.
4. integration --label <unique>: original full-environment static gates, deliberate rejection,
   then two original-initial-state short windows only. Hard cap: 5 and 2.35 ms.
5. checks --label <unique>: two-bubble adapter oracle and retained fail-stop proposal.
6. performance --label <unique>: current F4 masked / F3 ordered kernels, fixed inputs,
   no time advance. performance_v1 is the retained unmasked early microbenchmark;
   performance_v2 is used in the final current-kernel cost comparison.

The remote module accepts a module and unique --label and deploys a committed
Git archive to /workspace/mevo_acoustic_microbubbles/f10g/.execution/<sha>.
It checks the GPU and concurrent tasks, does not stop another project, and creates
one non-restarting Supervisor job with a fixed timeout. archive verifies every
remote source against local Git archive and input/result SHA-256 before reporting
successful local backup.

Offline: python -m scripts.fast_pzt_f10g.analyze, then figures, then delivery,
then verify. delivery is idempotent only if the external progress document has not
changed since its byte-preserving insertion; review a changed document manually.
Figures preserve true radius and coordinates. A failure record is never promoted
to an accepted frame. VISUAL_REVIEW.json records manual inspection of PNG and PDF.

Production remains original M0. The F10G CPU dense contact adapter is neither
GPU production code nor a calibrated SonoVue contact model. See both F10G reports.
