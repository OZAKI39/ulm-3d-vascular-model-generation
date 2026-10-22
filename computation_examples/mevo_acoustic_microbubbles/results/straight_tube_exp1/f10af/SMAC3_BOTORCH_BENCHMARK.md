# SMAC3 and BoTorch benchmark scope

This delivery distinguishes an actual physics experiment from an optimizer API exercise. MultiFidelityFacade/Hyperband replays the finite action catalog with fidelity budgets 1/2 at the SAME 2 ms duration. SingleTaskMultiFidelityGP + inverse-cost qMFKG fits measured records and evaluates acquisitions. BoTorch inner hypothetical optimization is relaxed; no proposed control is executed as a scientific trajectory. Initial SMAC catalog exhaustion is preserved; the retry is bounded to eight trials. Replay overhead cannot establish saved physics time. A formal superiority comparison with Sobol or CMA-MAE is INCONCLUSIVE / NOT_RUN.

## Evidence

- `OPTIMIZER_API_REPLAY.json`
- `OPTIMIZER_API_REPLAY_V2.json`
- `SOBOL_RESULTS.json`

EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED / training_ready=false
