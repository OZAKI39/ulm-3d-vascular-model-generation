# Low-pressure APECSS audit

Official APECSS commit `a0ee071ca247b5489e997d90e83b5f53c40e4ec7`, unchanged RPAR / reversible Marmottant model. 27 actual radial attempts: pressure nodes and independent midpoints, 48/96/192 cycles, 512/1024/2048 samples per cycle and tightened tolerances. The 0.0001 Pa finest attempt timed out at 120 s; its raw output is preserved and this domain is excluded.

Validated diagnostic suffix: 0.001–1 Pa. Complex harmonic and volume-rate energy errors are used; small noisy phase angles are not independently interpolated. Genuine F1/F2 values and CPU/GPU relative comparisons are in `LOW_PRESSURE_FINAL.json`. The final explicit low-pressure trajectories use source 0.003 Pa-equivalent, with actual incident pressure below 0.1 Pa. The older nominal 0.1 Pa-equivalent tests produced about 0.37 Pa locally and are not used alone as low-incident-pressure proof.

The earlier three failures are independently reconstructed at the last accepted state, not mislabeled as the failed kernel force.

EXPLORATORY · HOMOGENEOUS_WATER_CONTROL_PROXY · PHYSICS_UNCALIBRATED · ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED · training_ready=false
