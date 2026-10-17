# Causal 3/5/10-frame windows

Endpoints span N-1 actual PTS intervals. Recorded cadence is ~30.030 Hz, not assumed 60 Hz. Synthetic time uses explicitly registered 0.0333 s intervals. All frames available to the method are at or before the current endpoint. No future frames or bidirectional temporal smoothing.

| Frames | PTS span ms | Output fraction | MAE vx / vy (px/s) | MAE divergence / vorticity (1/s) | 90% step response ms |
| --- | --- | --- | --- | --- | --- |
| 3 | 66.6 | 0.9167 | 2.401/1.131 | 0.01887/0.07251 | 66.6 |
| 5 | 133.2 | 0.625 | 1.182/0.8785 | 0.01664/0.02351 | 133.2 |
| 10 | 299.7 | 0.5556 | 1.08/0.823 | 0.02561/0.09421 | RIGHT_CENSORED at 266.4 ms |

Errors condition on accepted outputs; lower error when many cases are rejected is not by itself an improvement. Detailed bias, SD, signs, CI coverage, false accepts, availability and processing times for each method/window are in `results/metrics.json`. Figure 4 shows the accuracy/availability tradeoff and a causal sudden translation step.

The 10-frame window has not reached 90% by the short test horizon, 266.4 ms after onset; this is right-censoring, not a measured infinite lag. The interval-equivalent velocity mixes old and new motion as it should. Window midpoint ages are 33.3/66.6/149.85 ms, in addition to software processing delay. Longer windows can also lose correspondence or average different local motions. No adaptive-window policy was tuned on the final test.
