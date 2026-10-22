# F10N water-filled, zero-flow control prototype

The baseline is F10M-R3 commit `e008a420b8bc386ac20454c556955b389be75fec`. F10N is a separate opt-in adapter. Historical F1–F10M code, interface admission and results remain unchanged.

The 32.5 × 9 × 9 mm block, central 1.6 mm bore, two mirror 1.5 mm open top holes and sealed X ends retain the R3 geometry. The 5.4 mm end-plane-to-hole-centre reference remains provisional. All cavity volumes relevant to motion are **assumed water filled** for this control experiment; this is an explicit user-authorized simplification, not a new measurement. The R3 numerical interior scope ceiling remains a fail-stop, not a free-surface wall. All present initial states and goals are in the central bore, far from that ceiling and from the holes. No flow, injection, droplets, free-surface dynamics, residual flow or streaming are simulated.

## Four source fields

Four 7 mm circular source faces at Y± and X± retain the measured-layout dimensions and inward normals. Each is represented by an independent uniform baffled piston in lossless homogeneous water, c = 1500 m/s, at 2.5 MHz. The fields superpose linearly. The approximation includes each disc's diffraction; it omits PDMS reflections, mutual source scattering, elastic waves and electrical/transducer transfer. The ideal baffle belongs to the isolated source approximation; it does not claim that a real apparatus seal is acoustically rigid.

The model label is **HOMOGENEOUS_WATER_CONTROL_PROXY**. Control amplitude in Pa is the plane-wave-equivalent quantity rho*c*v_normal, not measured face pressure and not applied voltage. Phase follows Re[P exp(+i omega t)]. Four amplitudes, four phases and four enable flags are accepted; each event has a duration. Optimization fixes the unobservable common phase to zero, leaving seven independent real search variables.

The implementation uses the circular-piston single-integral expression in [Kelly and McGough (2006), equation 9](https://pmc.ncbi.nlm.nih.gov/articles/PMC2117898/), with analytic spatial differentiation under the same integral. The CPU preflight compares it to a separately implemented two-dimensional Rayleigh area integral, the on-axis closed form, doubled quadrature order and finite differences. This is numerical verification of the proxy, not experimental calibration.

Four complex basis pressures are cached on a local 3 µm grid, with a separately built 1.5 µm refinement. During control updates only linear combinations are computed. A C1 Catmull–Rom pressure interpolant and its exact polynomial derivative are evaluated on GPU; this avoids a one-sided trilinear gradient at the Z symmetry plane. The new kernel retains the established APECSS interval/knot lookup and pressure phase rotation. No clipped pressures or table extrapolation are admitted. Cache borders have a strict interior guard.

## Frozen initial states and goals

The response table supports R = 1.78 µm, f = 2.5 MHz and nonzero amplitudes 0.1–40000 Pa. Exact zero drive is the quiescent branch. A smaller registered actual-particle pressure ceiling of 12000 Pa bounds this campaign. The source amplitude bound is 1500 Pa per disc.

Two fixed 32-ID states are stored before any trajectory: a seeded separated cloud and a connected 4×4×2 demonstration lattice with a positive surface gap of 0.0001R. The latter is manually prescribed, physically nonoverlapping geometry with a near-contact graph; it is not the output of an acoustic aggregation experiment and does not imply adhesion or fusion.

The position goal is fixed at (3, 2, 1) µm relative to the initial origin, with a 0.75 µm endpoint tolerance after 2 ms. The shape goal is sqrt(Qxx/Qyy) ≥ 1.25, with an identifiable major axis within 15° of X. The shape objective additionally penalizes centroid displacement. Q is the centred fixed-ID second moment; global aspect, Rg, eigenvectors and connected components are saved independently. An unidentifiable axis is not reported as a meaningful orientation.

## Search and numerical checks

The campaign begins with OFF and fixed BASELINE. Two small bounded static primary-force optimizations propose controls; their predictor ignores deformation and is never used as trajectory evidence. Actual full-M0 rollouts select candidates by registered endpoint error plus failure, connectivity and contact/wall penalties. Four bounded coordinate perturbations provide further actual-trajectory comparisons. There are at most eight navigation and four shape candidates, 22 total trajectories, 45 s per trajectory and 600 s for the campaign. Each static optimization is capped at 160 function evaluations. A 700 s external Supervisor timeout bounds process failure as well.

Selected position and shape controls are replayed from fresh simulators, compared at dt = 1 and 0.5 µs and compared with the refined acoustic cache. The original F2 near-field force candidate, local lubrication, 12-pass PBD and budgets remain unchanged: accepted overlap ≤0.01R, per-pass correction ≤0.1R, wall correction ≤0.1R. No new bulk hydrodynamics or calibrated contact law is introduced. Failed prefixes and proposals are retained. Source power cost is reported only as integral(sum amplitude² dt), in Pa² s, not physical electrical energy.

With only X/Y sources the proxy is even in Z. A symmetric initial cloud cannot establish independent signed Z translation through these controls. Local force-response singular values and actual Z motion are saved; no claim of arbitrary 3D controllability is assumed.

All outcomes remain **EXPLORATORY**, **PHYSICS_UNCALIBRATED**, `training_ready=false`. Short-run step/cache sensitivity is measured, not substituted for asymptotic convergence or experimental validation. Sampled graph connectivity and instantaneous contact events retain their existing time-sensitivity limits.
