# Restricted slow-time model

Original F1/F2 are evaluated by Warp at each RK stage. In the existing resistance implementation A is exactly identity when all surface gaps exceed 0.1 R; there is no extra far-field hydrodynamic term silently dropped here. SciPy RK45 uses rtol=1e-7, atol=1e-12 m, max step=0.5 ms. Its dense-polynomial coefficient bound certifies the entire accepted interval for all pair gaps, central-bore clearance and cache interior. Near contact, the method returns to the unchanged 1 us F10AA solver. A rejected RK stage never becomes an accepted trajectory point; a mandatory HF step prevents guard/retry cycling without physical progress. Failed pressure/domain checks remain fail-stop. Only a fixed action is admitted.

## Evidence

- `scripts/straight_tube_exp1/f10af/slow.py`
- `configs/straight_tube_exp1/f10af/protocol.json`
- `lf/*/STEP_RECORDS.json`

EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED / training_ready=false
