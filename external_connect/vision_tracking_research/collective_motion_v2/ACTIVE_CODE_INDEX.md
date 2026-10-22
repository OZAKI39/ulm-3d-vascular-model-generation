# Active code index

| Area | Entry | Purpose |
|---|---|---|
| Source reuse | src/legacy.py | Read-only F10U adapters |
| Data | src/dataset.py; scripts/prepare_data.py | Real patches, continuous triples, frozen split |
| Heatmap | src/heatmaps.py; src/learning.py | Actual detector input and H3/H4 perturbations |
| Geometry | src/geometry.py; tests/test_math.py | exp/log, fixed center, full resize conversion |
| Networks | src/networks.py | Official MobileNet blocks, Direct6, GeoMotion |
| Training | scripts/train.py; scripts/smoke_profile.py | Bounded full training and hardware profile |
| Selection | scripts/evaluate.py | Validation gate, independent calibration, locked test |
| Comparators | scripts/baselines.py | Frozen single/multipass OpenPIV and CPD |
| Metrics | src/metrics.py; scripts/analyze.py | Explicit denominators, risk curves, grouped analysis |
| Failure analysis | scripts/posthoc_diagnostic.py | Post-hoc prior-collapse description only |
| State | src/state_v2.py; scripts/infer.py | New schema and old-consumer adapter |
| Runtime | scripts/replay_runtime.py | Actual-PTS causal latest-slot replay |
| Visualization | scripts/figures.py; scripts/synthetic_demo.py; scripts/videos.py | Real figures and chronological MP4 |
| Reproduction | scripts/independent_check.py; tests/ | Independent checkpoint, gradients, math, MockConsumer |
| Audit | scripts/dependency_audit.py; scripts/report.py | Dependencies and result-derived reports |

Historical collective_motion_v1 and smart_online_tracking_v2 remain frozen. No Chronos receiver, PZT or AFG implementation is modified.
