# Known limitations

- PZT controls specify ideal acoustic-source pressure and phase. No voltage-to-pressure calibration or internal electromechanical model is implied.
- Stage 3 accepts the minimal implementation and tested engineering trends. It does not establish experimental calibration, patient-specific predictions or complete local physical accuracy.
- The independent history study finds substantial quantitative differences in transient displacement, even for mild radial excursion. Averaging, added mass and memory must not be conflated. The zero-Re reference itself omits finite-Re corrections and is not a universal truth model; H3 includes Reynolds excursions beyond the cited experimental mean-Re range.
- The water demonstration cannot be relabelled as blood transport; changing the medium requires matching radial-response and acoustic evidence.
- The radial table covers only its declared R/f/P box and fixed Tu2009-derived shell manifold. Effective radial viscosity incorporates the existing damping approximation; it is not the viscosity used for translational drag. Initial surface tension and reversible shell dynamics remain assumptions. No dissolution, shedding or nonspherical instability is modelled.
- The two-grid acoustic diagnostic has nonzero pressure/gradient differences. Superposition accuracy is an implementation check, not proof of spatial convergence or local physical accuracy. Relative phase is undefined at zero acoustic pressure.
- Classical secondary forces neglect mutual radial coupling, propagation delay and detailed near-field interactions. Separation categories only warn. Many widely separated pairs in the cylinder exceed the configured kd limit. Dense clusters are geometrical/topological approximations.
- No near-wall drag, lubrication, adhesion, fusion or coalescence is predicted. Wall/contact projection only enforces geometry within the configured tolerance.
- No acoustic streaming, two-way flow/acoustic coupling or patient geometry. Synthetic VTK equivalence does not validate future imported CFD.
- GPU and CPU agreement verifies implementation consistency. Runtime crossover measurements are specific to this server, workload and software stack.
- Lazarus 2017 is a qualitative pressure/clustering trend benchmark, not parameter-identical reproduction. Zhang 2011 is logically archived and excluded from development admission gates.

See [the quantitative report](STAGE3_VALIDATION_REPORT.md) and [history sensitivity](HISTORY_FORCE_SENSITIVITY.md). The next scope is simple bifurcation/SimVascular integration, with these limitations retained.
