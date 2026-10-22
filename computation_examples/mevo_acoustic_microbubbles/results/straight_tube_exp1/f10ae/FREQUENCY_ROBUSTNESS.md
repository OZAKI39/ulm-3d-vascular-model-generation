# Frequency robustness

| Case | Status | Aspect change (%) | Rg change (%) | COM (nm) |
| --- | --- | --- | --- | --- |
| freq_T1_2400000 | PASS | 5.4531 | -2.7598 | 3.0669 |
| freq_T1_2500000 | PASS | 2.9974 | -0.7359 | 2.1325 |
| freq_T1_2600000 | PASS | 7.0153 | -4.2029 | 0.32061 |
| shell_T1_0.8 | PASS | 3.6672 | -0.82634 | 1.979 |
| shell_T1_1.2 | PASS | 2.156 | -0.50317 | 1.9897 |
| position_T1_-5e-06 | PASS | 3.0016 | -0.73766 | 2.1492 |
| position_T1_5e-06 | PASS | 2.9965 | -0.73548 | 2.126 |
| freq_T2_2400000 | PASS | 0.34095 | -0.15498 | 1.2637 |
| freq_T2_2500000 | PASS | 0.40956 | -0.18955 | 1.2628 |
| freq_T2_2600000 | PASS | 0.45084 | -0.21204 | 1.2076 |
| shell_T2_0.8 | PASS | 0.48433 | -0.23553 | 1.2683 |
| shell_T2_1.2 | PASS | 0.29808 | -0.13159 | 1.1372 |
| position_T2_-5e-06 | PASS | 0.40954 | -0.18955 | 1.2627 |
| position_T2_5e-06 | PASS | 0.40954 | -0.18955 | 1.2628 |
| Ynode_-5e-06 | PASS | 0.019705 | -0.012395 | 0.12268 |
| Ynode_0 | PASS | 0.010057 | -0.0056781 | 0.057222 |
| Ynode_5e-06 | PASS | 0.021052 | -0.012419 | 0.031877 |

Each 2.4/2.5/2.6 MHz test uses its own admitted field, gradient, response and F1/F2. Same frozen actions, no new optimization at each frequency. Half-step and fine-cache comparisons are in ROBUSTNESS_CONVERGENCE.json; accepted integration alone does not establish numerical admission. New small frequency offsets: NOT_RUN. Three discrete frequencies cannot establish a universal frequency ranking.

EXPLORATORY · HOMOGENEOUS_WATER_CONTROL_PROXY · PHYSICS_UNCALIBRATED · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
