# Runtime and resources

Server: AMD Ryzen9 7950X3D,16physical/32logical cores; one RTX4090 24GiB. NVIDIA550.120, PyTorch2.14.0+cu126 / TorchVision0.29.0+cu126, BF16 operations and gradients verified. Actual environment, free RAM/disk, driver query and bounded64MiB fsync/warm-read probe are in `evidence/environment.json`. Other Supervisor jobs were not stopped. Isolated env passed dependency checks; a missing setuptools dependency was copied into the new env from an existing installed wheel. Slow package download and one interrupted parallel transfer were resolved; transfer hashes were verified.

Formal suite: 12 trained runs, 258 epochs and 11094 optimizer steps; suite wall 705.00s (11.75min). The individual run total is 703.25s. These are substantive repeated full-dataset passes, not smoke batches. Each run early-stopped by the declared validation rule, with best/last/periodic checkpoints saved.

| Configuration / seed | Epochs / best | Optimizer steps | Wall s | Stop reason |
|---|---:|---:|---:|---|
| direct_gray / 1729 | 16 / 2 | 688 | 33.19 | VALIDATION_EARLY_STOP |
| direct_fused / 1729 | 20 / 12 | 860 | 42.52 | VALIDATION_EARLY_STOP |
| direct_heat / 1729 | 28 / 20 | 1204 | 60.80 | VALIDATION_EARLY_STOP |
| direct_fused_robust / 1729 | 19 / 11 | 817 | 71.28 | VALIDATION_EARLY_STOP |
| geo_fused / 1729 | 16 / 1 | 688 | 56.34 | VALIDATION_EARLY_STOP |
| direct_gray / 2718 | 22 / 14 | 946 | 46.20 | VALIDATION_EARLY_STOP |
| direct_fused / 2718 | 19 / 11 | 817 | 39.88 | VALIDATION_EARLY_STOP |
| direct_heat / 2718 | 26 / 18 | 1118 | 56.38 | VALIDATION_EARLY_STOP |
| direct_fused_robust / 2718 | 25 / 17 | 1075 | 91.01 | VALIDATION_EARLY_STOP |
| geo_fused / 2718 | 16 / 1 | 688 | 57.39 | VALIDATION_EARLY_STOP |
| direct_fused_2frame / 1729 | 21 / 13 | 903 | 44.80 | VALIDATION_EARLY_STOP |
| geo_no_local_supervision / 1729 | 30 / 22 | 1290 | 103.47 | VALIDATION_EARLY_STOP |

Training GPU utilization (1s device samples) averaged 59.7%; peak allocated CUDA memory 5087.4MiB. Parent CPU core-equivalent range 0.66–1.20; DataLoader child CPU is not included in that measure. Device-wide utilization includes any concurrent device work; GPU process audit showed no other compute job at launch. No unsupported kernel-time utilization total is claimed.

| Method | Pure model p50/p95 ms | Full pipeline p50/p95/max ms | Consumed /120 | Fresh six-parameter Hz |
|---|---:|---:|---:|---:|
| PIV single pass | N/A | 34.35, 44.29, 47.32 | 114 | 25.24 |
| Direct6 gray | 2.45, 3.25 | 19.72, 25.16, 32.23 | 120 | 0.00 |
| Direct6 + heat | 3.24, 4.39 | 16.08, 21.13, 22.98 | 120 | 0.00 |
| GeoMotion | 2.94, 3.42 | 18.92, 23.56, 26.06 | 120 | 0.00 |

Cold new-process probes: Direct6 gray: 2.910s (imports 2.659s, checkpoint/CUDA 0.101s, first inference 0.149s); Direct6 + heat: 3.389s (imports 3.059s, checkpoint/CUDA 0.128s, first inference 0.201s); GeoMotion: 4.019s (imports 3.557s, checkpoint/CUDA 0.097s, first inference 0.363s). These include Python module imports and first CUDA inference, excluding interpreter spawn/SSH/video decode. Offline batch128 cached-input throughput: Direct6 gray: 4356.8 triples/s; Direct6 + heat: 4374.3 triples/s; GeoMotion: 1131.6 triples/s; this excludes detection and is not streaming fresh-state Hz. CUDA-synchronized training-step wall time sums to 491.13s (includes host/optimizer overhead; GPU kernel time was not separately measured).

Model timing uses60 synchronized warmed repetitions; cold checkpoint load is separate. Full pipeline includes detection, heat creation if used, resize/transfer/model, geometric conversion, distribution, schema and consumer send. Processing is paced to actual recording PTS with a single latest-frame slot. PIV consumed114/120, dropped6, and formed100 valid fresh windows; AI consumed120/120, formed118 windows and had zero valid outputs due to frozen rejection. Empty counts are not presented as low-risk successes. History spans66.6ms, separate from software age; interval reference age includes half that span. No-new-frame TTL checks passed.

Explicitly timed stage wall seconds: `{"training_suite": 704.998767557001, "validation_and_calibration_predictions": 11.310612702989602, "final_predictions_including_H3": 16.765406527025334, "PIV_CPD_accuracy_runs": 10.931113819002348, "causal_replays_monitored": 25.2097787169987, "server_data_generation": 15.356769906997215, "supplemental_demo": 27.06455608199758, "new_process_cold_probes": 10.31809217800037, "video_render_full_decode_local": 32.56709431199124}`. Their sum is 854.52s (14.24min). This ledger excludes installation/download, some setup/smoke and manual review, and must not be described as whole engineering elapsed time or pure GPU kernel time. Local rendering is listed separately from server computation.

Local independent reproduction uses PyTorch2.11.0+cu128 / TorchVision0.26.0+cu128 and RTX4060. Gray/heat images and detector centers reproduce byte-for-byte; label differences are <=2.71e-12. Checkpoint inference differs by at most0.087 frozen error-tolerance units in the tested128 validation triples (BF16/hardware/version variation), not bit-identical. Optimizer state restores and changes weights on a temporary in-memory train step; formal files remain unchanged.
