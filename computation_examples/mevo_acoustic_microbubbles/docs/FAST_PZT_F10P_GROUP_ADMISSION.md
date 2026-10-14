# F10P layered admission

Registered before numerical results in commit 36343716a. Each layer has independent fields in analysis/FINAL_GATE.json; no particle-ID boolean enters the COM or shape gate.

| Layer | Registered decision |
|---|---|
| BASIC_NUMERICS_PASS | Complete declared window with strict acoustic/response/flow, sphere, resistance equation, complementarity, correction and real-wall checks. Failed intervals are never admitted; legal prefixes retain their actual coverage. |
| GROUP_POSITION_NUMERICALLY_ADMITTED | Complete 2 ms dt1/dt0.5/cache comparison, maximum COM difference ≤0.02 µm, control displacement greater than measured numerical difference. OFF is an explicitly stationary null reference. The 0.02 µm gate is 26.7% of the unchanged 0.075 µm target tolerance. |
| GROUP_SHAPE_NUMERICALLY_ADMITTED | Rg difference ≤0.02 µm; X/Y absolute ratio difference ≤0.01; Frobenius Q difference / reference Q norm ≤1%; global aspect relative difference ≤1%; identity-free centered Hungarian matching RMS/R ≤0.02; direction difference ≤2 degrees only where both top eigenvalue relative gaps exceed0.001. |
| GROUP_TOPOLOGY_NUMERICALLY_ADMITTED | Initial 32-member single network survives complete 2 ms in all three comparisons. Component membership and event ledgers retained. This bounded gate does not automatically admit arbitrary later split/merge events. |
| PARTICLE_ID_NUMERICALLY_ADMITTED | Independently maximum same-ID displacement difference ≤0.05 µm for complete step and cache comparisons. |
| TARGET_REACHED | Group-position admitted AND actual endpoint distance to preregistered X/Y ≤0.075 µm. |
| SHAPE_TARGET_REACHED | Group-shape admitted AND original X/Y ratio target1.25 reached; no new shape-optimization campaign is required here. |
| EXPERIMENTALLY_CALIBRATED / TRAINING_READY | NOT_VERIFIED / false, regardless of numerical layers. |

All positions use fixed original32 IDs. Centered shape removes translation; Hungarian matching removes identity only, not orientation, so it cannot hide real directional deformation. Principal directions at degeneracy are NOT_IDENTIFIABLE, not assigned a spurious rotation. A split would make the global COM a center of the original population, not one intact cluster; membership is retained explicitly.

Instantaneous edge changes use the unchanged 2.1R geometric threshold and a separate10 µs event-time comparison. Persistent events require2 ms observation; later events are right-censored. Initial single-network survival is a genuine2 ms statement but not evidence of a new aggregation event, bubble fusion or long-term persistence.

The additional simulation_group_position_eligible, simulation_group_shape_eligible, simulation_topology_eligible and simulation_particle_id_eligible fields describe numerical consistency inside the current model only. If the ID gate fails it does not change the other gates. Conversely, COM consistency cannot admit shape, topology, contact physics or experimental labels.

B2/hard-sphere/M0 contrasts are saved separately from within-model dt/cache differences. Model-dependent shape responses cannot establish which fixed-sphere contact model matches real coated, deformable bubbles. EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / training_ready=false apply throughout.
