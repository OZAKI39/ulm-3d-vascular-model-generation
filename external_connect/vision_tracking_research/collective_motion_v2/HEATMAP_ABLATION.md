# Heatmap ablation

H0 gray, H1 gray+actual detector heat, H2 heat-only, H4 separately trained robust gray+heat each have two seeds and the same data/cap/early-stop policy. H3 applies fixed dropout/false-center/jitter/severe corruption to the validation-selected H1, H2 and H4 checkpoints without retraining or recalibration. Underlying independent test gray images and truth are identical. Full numeric tables are in `results/analysis.json`; failures remain in `results/all_failures.jsonl`.

| Method | MAE denominator | Six-parameter MAE (vx, vy, div, vort, s1, s2) | Accepted / 1536 | Failed accepted |
|---|---:|---|---:|---:|
| PIV single pass | 1114 | 1.894, 1.698, 0.020, 0.039, 0.009, 0.010 | 1111 | 85 |
| Direct6 gray | 1320 | 11.420, 7.994, 0.057, 0.327, 0.045, 0.030 | 0 | 0 |
| Direct6 + heat | 1320 | 11.394, 8.011, 0.061, 0.332, 0.040, 0.031 | 0 | 0 |
| Direct6 heat only | 1320 | 13.416, 9.526, 0.059, 0.300, 0.044, 0.031 | 0 | 0 |
| Direct6 + heat robust | 1320 | 11.554, 7.937, 0.057, 0.288, 0.043, 0.026 | 0 | 0 |
| GeoMotion | 1320 | 10.833, 7.597, 0.057, 0.282, 0.040, 0.023 | 0 | 0 |

No reliable six-parameter advantage is established by adding heatmaps. H1 slightly changes shear error but its overall recovery remains far below PIV. H2 is especially sensitive: vx MAE rises from 13.42 to 16.65 px/s under severe heat corruption. H4 is less sensitive, but low sensitivity is not evidence of useful robustness when outputs can collapse toward a prior.

The predicted vx standard deviation / truth standard deviation is 0.0385 for H1 and 0.0046 for GeoMotion. This post-hoc failure diagnostic was not used to retrain or select a new model. Three-frame versus two-frame and local-loss ablations also fail the joint risk gate; no temporal or geometric reliability improvement is established.
