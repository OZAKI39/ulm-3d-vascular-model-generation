# Literature and physics scope

Only primary papers and official software documentation support the method choices. They do not calibrate this chip.

- [Kelly and McGough, circular-piston FNM, Eq. 9](https://www.egr.msu.edu/~fultras-web/publications/tsd-circ.pdf): finite circular-source pressure integral underlying the retained WaterDiscs implementation; current derivatives and quadrature are independently checked.
- [Polymer-coated microbubble acoustofluidics study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6187510/): measured primary/secondary forces and sub-wavelength organization motivate keeping shell response and pair interactions. Its different bubble population and frequency cannot calibrate 2.5 MHz, 1.78 um bubbles here.
- [Acoustic vortex trapping study](https://pmc.ncbi.nlm.nih.gov/articles/PMC7354944/): microbubble trapping in a designed vortex motivates examining phase and forces; trapping is not proof that this four-source cloud rotates.
- [Barnkob et al., radiation and streaming measurements](https://arxiv.org/abs/1208.6534): particle transport can involve both radiation and streaming. It is not a quantitative microbubble calibration.
- [SciPy integration API](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html): mature adaptive integration; this implementation uses its RK45 class to interpose conservative geometric acceptance and original contact fallback.
- [SMAC3 official multi-fidelity interface](https://automl.github.io/SMAC3/latest/api/smac/facade/multi_fidelity_facade/) and [SMAC3 paper](https://www.jmlr.org/beta/papers/v23/21-0888.html): existing mixed-space optimization and budget allocation, not a newly implemented optimizer.
- [BoTorch official discrete-fidelity tutorial](https://botorch.org/docs/next/tutorials/discrete_multi_fidelity_bo): existing multi-fidelity GP and knowledge-gradient acquisition with cost accounting. Offline API testing is separate from scientific admission.
- [SciPy Sobol](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.qmc.Sobol.html): bounded, reproducible baseline; [pyribs](https://github.com/icaros-usc/pyribs) CMA-MAE is retained as prior morphology exploration, not rerun as an unrelated contest.

Acoustic streaming, buoyancy, wall hydrodynamics, spin, PDMS elastic scattering, size/shell dispersion and long-time bubble stability remain uncalibrated or absent. A spatial phase singularity, nonzero curl, local F1 equilibrium or a few-micrometre simulated displacement is not by itself an experimental pattern, sustained vortex or camera-observable effect.

## Evidence

- `PREREGISTRATION.json`
- `ENVIRONMENT.json`

EXPLORATORY / HOMOGENEOUS_WATER_CONTROL_PROXY / PHYSICS_UNCALIBRATED / ELECTRICAL_TO_ACOUSTIC_UNCALIBRATED / training_ready=false
