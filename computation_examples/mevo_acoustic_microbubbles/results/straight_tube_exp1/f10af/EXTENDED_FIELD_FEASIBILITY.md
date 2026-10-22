# Extended field feasibility

3 frequencies; bounds X +/-620 um, Y +/-520 um, Z +/-100 um; 10 / 5 um caches; strict cubic interpolation interior.

48-point FNM quadrature was checked against order 96, independent Rayleigh integration, and pressure finite differences. Every pressure and gradient derives from the same complex basis. Sampled aggregate-L2 coarse gradient relative errors (128 independent points) are 0.074–0.086%; fine 0.019–0.022%. Full trajectory cache refinements are separately recorded. The prescribed near-top state is outside this central field and not simulated.

## Evidence

- `FIELD_BUILD.json`
- `fields/*/FIELD_GATE.json`
- `SCIENTIFIC_RESULTS.json`

EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED / training_ready=false
