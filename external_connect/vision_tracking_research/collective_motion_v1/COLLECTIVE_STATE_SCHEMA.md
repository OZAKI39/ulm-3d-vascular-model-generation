# M4 collective observation interface

Machine-readable schema: `configs/collective_state.schema.json`. Actual message: `evidence/state_example.json`. Implementation: `src/state.py`.

The payload includes source/session/frame ID, optional capture timestamp, host receive/processing timestamps on one monotonic clock, source PTS, ROI, calibration and image quality. Real recording PTS are not camera exposure timestamps. Live receiver time, if available, is explicitly host decode time, never silently an exposure time. Camera-to-controller end-to-end latency remains NOT_MEASURED.

Distribution: visible count, center, covariance, eigenvalues, principal axis and orientation (null if near isotropic), Rg, axis ratio, density and boundary fraction. Area validity currently means geometric ROI coverage, not a validated in-focus/PDMS mask. Visible point statistics are not segmented physical bubble-group boundaries. Density bandwidth 12 px and candidate/px² units are explicit. Density/count changes can reflect ingress/egress or missed detections.

Motion: six values plus six booleans and six reasons. Invalid = null, never old value or zero. Support size, hull coverage, normalized condition number, residual, uncertainty method and interval are retained when available. CPD intervals are multi-start empirical ranges, not confidence intervals. Bootstrap absent in timed live/replay mode is explicitly null, not a fabricated precision estimate. `valid=true` means algorithmic support/quality/freshness only; `scientific_accuracy_verified=false` and `training_ready=false` remain.

Future actions: channel, actual state, amplitude and unit, phase degrees, frequency Hz, command/execution time, synchronization quality, pressure-calibration flag. Every run used `action_event=null`; there is no actuator connection. Acoustic Pa-equivalent is never assigned to an AFG voltage field.

Local ZeroMQ queues have capacity one, no unbounded backlog. Duplicate/old/PTS-reversed/stale messages are rejected. When no new image arrives, `ConsumerGate.current()` expires the previous state after 100 ms. M4 freshness invalidates the distribution as well as motion. Consumer poll watchdog records the expiration. Resets require causal warmup and do not reuse pre-discontinuity frames. `tensor_input()` returns float32 values (NaN for missing) and a separate boolean mask; it never silently fills missing values with zero. PyTorch callers can use `torch.as_tensor(values)` and `torch.as_tensor(mask)`; no model is trained.

Unit and actual local-message tests cover partial field validity, null preservation, duplicate/time reversal, stale observations, dropped frame warmup and no-new-frame expiry. `results/mock_fault_tests.json` records the transport test.
