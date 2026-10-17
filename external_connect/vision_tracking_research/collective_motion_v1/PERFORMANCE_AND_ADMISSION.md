# Measured runtime and admission

| Replay pipeline | Processed / input | Software p50 / p95 / max ms | Fresh motion Hz, active span | Fresh motion Hz, startup included |
| --- | --- | --- | --- | --- |
| PIV 2-pass | 64 / 120 | 6.406/148.3/157.7 | 0 | 0 |
| PIV candidate 1-pass | 104 / 120 | 39.84/48.02/57.74 | 17.51 | 13.16 |

Active acquisition spans 3.9961/3.9969 s. Startup-inclusive raw rates include decoder preroll; both denominators are reported explicitly. M4 message rates during active replay are 16.016/26.020 Hz, while fresh six-parameter rates are 0.000/17.513 Hz. Candidate distributions can remain available during motion warmup. Skipped frames: 56/16. Consumer rejection includes both expired processing and frame-discontinuity resets, not only stale processing. Latest slots/queues are bounded at one; there is no accumulated FIFO lag.

The 2-pass p50 includes cheap warmup/reset iterations, **not** a 6 ms full multipass solve. Detailed PIV and fit-only timing is in `causal_replay.json`. Single-pass is a secondary latency probe selected from already tested development options, not a new independently tuned accuracy winner. It did not undergo a new independent held-out accuracy trial after this speed probe. No claim of 30/60 Hz complete real-time six-parameter processing.

Transport p95 0.976/1.018 ms; CPU core equivalents 0.970/1.000, process RSS samples retained. PIV/CPD run on CPU. SEA-RAFT tensor peak 102.56 MiB, measured device-wide utilization and process RSS retained; includes other desktop GPU users. AI times are forward+reverse flow plus mask/transfer; not a single-forward-only benchmark.

Measured independent synthetic test processing 38.384 s, real 120-frame all-window offline analysis 56.901 s. Benchmark timers exclude synthetic generation/setup; AI resource interval 2.510 s excludes weight loading/warmup. No invented total GPU-hour/whole-development estimate. Install/plot/video time is not online performance.

Admission: software/math/interface PASS with explicit limits; clear-case synthetic recovery LIMITED_PASS; difficult-case false acceptance and calibrated uncertainty FAIL; real absolute motion accuracy NOT_VERIFIED; live camera unavailable; hardware/PZT/AFG NOT_RUN; CAMERA_TO_CONTROLLER_END_TO_END_LATENCY NOT_MEASURED. `training_ready=false`. Freshness does not establish scientific accuracy.
