# F10F Small Steps audit

## Registered basis and limits

The frozen protocol is configs/fast_pzt/f10f/protocol.json, committed before GPU
checks or scientific trials. The baseline is F10E 741ad7d8683bfbb2fc9a270d9e92a218dead5e54.
Historical results are read-only and individually hashed in BASELINE_AUDIT.json.

[Small Steps, Macklin et al. 2019](https://mmacklin.com/smallsteps.pdf),
DOI [10.1145/3309486.3340247](https://doi.org/10.1145/3309486.3340247),
studies inertial XPBD. Equations 3–7 and Algorithm 1 motivate fresh predictions
with fewer iterations per substep. Its displacement from external forces scales
with dt squared. Here overdamped Stokes displacement scales with dt, and the
two-pass local lubrication and central Jacobi projection are retained. Neither
the paper's convergence behavior nor speed factors are inherited. Its amortized
collision detection is not adopted: distances and secondary force are fresh at
every physical step.

The [PBD author manuscript](https://matthias-research.github.io/pages/publications/posBasedDyn.pdf)
and [2007 DOI](https://doi.org/10.1016/j.jvcir.2007.01.005) discuss sequential
nonlinear projection. The present code instead reads all pair positions before
writing simultaneous relaxed central corrections. The [XPBD paper](https://matthias-research.github.io/pages/publications/XPBD.pdf),
DOI [10.1145/2994258.2994272](https://doi.org/10.1145/2994258.2994272), introduces
compliance and accumulated multipliers; neither is added to the microbubble model.

[Warp 1.18 runtime documentation](https://nvidia.github.io/warp/v1.18/user_guide/runtime.html)
describes graph capture and replay. This implementation captures one substep,
checks its quality, and only then launches the next. Python decisions made
during capture are not runtime conditional branches. PZT enable changes use
the existing graph invalidation; controls on this fixed grid execute before
the substep starting at the requested physical time. Nonaligned controls are
explicitly rejected, never delayed or crossed.

## Source audit

engine.step performs fields, pairs, free velocity, two local mobility passes,
component mean balance, prediction, Jacobi PBD, acceptance and clock update.
F2 owns pair feedback, F3 ordered reductions, F4 wall passes after prediction
and every third contact round, F5 actual CFD queries at step input and final
proposal plus exact open-boundary processing, and F6 original strict budgets.
F9 ledger provides the precedent for passive observation. F10C remains the
default full-particle batch path; F10D postprocessing is reused unchanged.

B0 is 10us/12, SS2 is 5us/6, SS4 is 2.5us/3. All keep relaxation 0.5,
two mobility iterations, field and flow update stride 1. Original wall pass
schedule means respectively 5, 3, 2 wall passes per physical step. SS12 is not
enabled. Both stored-input and final-proposal CFD queries remain present.

Two recording/acceptance details need explicit opt-in handling:
- Odd rounds end in tmp2; the final acceptance argument, not a hardcoded tmp1,
  is saved as the rejected proposal.
- The original host pass budget is also checked before GPU acceptance. Its
  threshold stays 178nm and contact error bit stays 2048. Clock advancement
  occurs only without error. Positions and active flags must equal the prior
  legal state on a rejected substep. No retry of a rejected object is allowed.

Full traces retain the last 32 attempted steps, while the scalar ledger retains
every substep, round residual, correction, force extrema, velocity and wall
correction. Rejected proposals are separate from accepted saved frames.
No engine copy, frozen force, soft tolerance or new physical force is introduced.

## Experimental outcome

All four SS2/SS4 attempts failed earlier than the corresponding B0 failure.
Original contact budgets, correct odd-round proposal capture, exact B0 prefixes,
CPU/GPU oracles and graph first-error checks passed engineering validation.
The contact and morphology gates did not admit a candidate for production.
This two-case negative result is not a refutation of inertial XPBD or an
experimentally calibrated claim about microbubble shells. Actual tables and
cost limits are in FAST_PZT_F10F_DEVELOPMENT.md and analysis JSON files.
