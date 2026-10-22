# Deformation and recovery

First action: X pair 700 Pa-equivalent, 180 degrees, 2 ms; numerical regression precedes recovery. Thirty-two 4 ms candidate protocols were attempted; 31 completed and one failed before the endpoint: 16 X-pair second actions, then 16 bounded X/Y combinations. No reverse playback or invented restoring force. Selected active candidate: `recovery_05`.

| Metric | At switch error | Active endpoint error | OFF endpoint error | Recovery fraction |
| --- | --- | --- | --- | --- |
| Rg | 0.00693763 | 0.00745502 | 0.00693763 | -0.07457622888322768 |
| axes | 0.028231 | 0.0299055 | 0.028231 | -0.059311082459390985 |
| covariance | 0.0125199 | 0.0131998 | 0.0125199 | -0.05430455499052611 |
| Kabsch_m | 2.39408e-07 | 2.49338e-07 | 2.39408e-07 | -0.04147921746753891 |
| pair_distance_m | 2.27263e-07 | 2.37357e-07 | 2.27263e-07 | -0.04441587672446623 |
| neighbor_error | 0.146667 | 0.146667 | 0.146667 | 0.0 |

Recovery levels: {"RG_RECOVERY_PASS": false, "PRINCIPAL_AXIS_RECOVERY_PASS": false, "COVARIANCE_RECOVERY_PASS": false, "FIXED_ID_CONFIGURATION_RECOVERY_PASS": false, "TOPOLOGY_RECOVERY_PASS": false}

Rg recovery is never substituted for structure or topology recovery. At a zero initial topology error, unchanged edges are maintenance, not evidence of a restoration process.

EXPLORATORY · HOMOGENEOUS_WATER_CONTROL_PROXY · PHYSICS_UNCALIBRATED · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
