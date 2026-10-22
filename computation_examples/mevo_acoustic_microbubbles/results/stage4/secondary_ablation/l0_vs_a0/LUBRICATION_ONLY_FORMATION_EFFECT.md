# Lubrication-only formation effect

This comparison isolates the effect of the Stage 4L translation change before any secondary ablation.

Exact shared HDF5 saved frame: 0.08 s. No simulation was run. Initial positions, IDs and radii match bitwise.

| Case | Model | Largest cluster | Cluster count | Near-contact fraction | Median NN (µm) |
|---|---|---:|---:|---:|---:|
| B2 | L0 | 19 | 74 | 0.8125 | 3.560 |
| B2 | A0 | 19 | 73 | 0.8438 | 3.560 |
| B3 | L0 | 19 | 72 | 0.8281 | 3.560 |
| B3 | A0 | 20 | 72 | 0.8789 | 3.560 |

Cluster threshold: 10 µm center distance; near-contact particle fraction: minimum surface gap ≤0.5 µm, unchanged from the stored Stage 4 configuration. Dense same-ID diagnostics separately use the inherited 2.1R graph definition. Full distributions, centroids, gyration radii and membership counts are in the JSON.

L0→A0 is the translation effect. A0→A2 will be the incremental ablation effect. L0→A2 combines both changes.
