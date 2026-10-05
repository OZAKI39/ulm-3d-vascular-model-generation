# Validation-only moving-trajectory history sensitivity

Production must never import this directory. The production flow, radial table,
force law and geometric contact remain frozen Stage 3 models.

The intended reference separates two timescales. A periodic variable-radius
translation problem uses actual cached APECSS radial cycles and the conservative
added-mass / no-slip zero-Re history equations from [Garbin et al. (2009)](https://pages.jh.edu/aprospe1/publications/PapersPublished/GasBubbles/GarbinHistory.pdf)
and [Cattaneo & Supponen, author preprint v2 (2026)](https://arxiv.org/html/2506.16413v2),
Eqs. (11), (14). In transformed time theta = integral dt/R², the infinitely
periodic Abel integral has Fourier multiplier sqrt(pi) sqrt(i omega_theta).
This allows the fast radial oscillations to be resolved without billions of
moving-trajectory microsteps. The local pressure and complex gradient are updated
along the moving reference path, using the same k-Wave field and SimVascular flow.

An additional equilibrium-radius added-mass/Abel kernel models the slow slip
envelope with zero initial slip. This combination is an **adiabatic two-scale
reference**, not a resolved moving variable-radius history solution. Startup,
fast/slow coupling, unmodelled fluid-acceleration terms and wall hydrodynamics
remain limitations. No empirical displacement correction is applied.

The reference with history disabled retains the same periodic radial response
and added mass; its comparison with the history-enabled reference isolates the
history terms within this approximation. Comparing production with the enabled
reference additionally includes microcycle averaging and inertia differences.
A reduced-population production control must identify any change caused by
omitting the other particles' geometric contacts before attributing a branch
change to history. Particle selection is automatic and recorded.

The production adaptive radial table and the reference's independently sampled
periodic closure have different interpolation errors. Production-vs-reference
differences therefore also include that numerical approximation; B-vs-C retains
the same radial nodes and closure construction. Neither comparison validates
the missing wall hydrodynamics.

Required checks before scientific use: constant-radius harmonic response,
Fourier resolution convergence, interpolation holdouts, slow-step refinement,
and production-import exclusion. Passing these numerical checks does not make
this reference an experimental validation or remove the Stage 3 displacement
limitations. Branch agreement is reported only for the selected sample.
