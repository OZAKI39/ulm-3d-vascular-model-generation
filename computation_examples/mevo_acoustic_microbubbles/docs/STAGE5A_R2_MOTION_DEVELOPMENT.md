# R2 development ledger

Base 7c39828d9; plan 5877e733a. Source before formal validation is development,
not promoted. All frozen force and lubrication classes remain unchanged.

1. QP strategy 1: native ProxQP 0.7.2 at eps_abs=1e-12. All three original
   failed QPs pass unchanged, KKT 1.5e-14 / 3.5e-14 / 9.1e-13 (dense).
   Sparse also passes. Native comparison includes OSQP and Clarabel failures.
   34 new manufactured/interface/real-QP tests passed (tests_dev01.xml).

2. Integrator V1, source 18c7fa38f: converged endpoint resistance AND endpoint
   contact tangent. B2 prefix exposed SWEEP_NOT_CERTIFIED: initial negative gaps
   close to -1e-11 m can pass endpoint conditions but fail intermediate straight
   segment clearance. E.g. dt=1e-8 s minimum swept gap -1.2411e-11 m.
   Steps shrank below 1 ns, not acceptable progress. This diagnostic was stopped
   by the developer at a checkpoint (not resumed old campaign, not a user stop).
   Preserve run `20261006T215835450341Z_B2_S1_1e-06s_prefix_v1`.

3. Integrator V2 hypothesis: a SINGLE start-tangent normal set ensures a convex
   pair distance remains above its initial tangent lower bound along a straight
   segment; retain endpoint self-consistency only in A(x1). Equation:
   A(x1)U = F(x0,t0)+D*u_flow(x0,t0)+J(x0)^T lambda, x1=x0+dt U.
   This is first-order semi-implicit contact/forcing, nonlinear implicit normal
   lubrication resistance. End wall geometry still needs full WALL sweep.
   No new force term, projection, gap floor or tolerance. New source/tests/run
   identity required; V1 evidence is retained in its source snapshot.

4. V2 completed B2 1 us with 60 accepted half steps, 15 nonlinear rejections,
   33.43 s, KKT 9.04e-10, linear residual 3.59e-14, minimum swept pair gap
   -9.99757e-12 m and no projection. Its 10 us continuation was stopped for
   diagnosis; repeated Picard iterations were the bottleneck (1888 QP calls for
   the 1 us prefix). The solver-map velocity difference hovered at 1e-9 while
   independent nonlinear equilibrium can already be substantially better.

5. V3: bounded 6-history Anderson acceleration of the nonlinear solve, with a
   trial-iterate step safeguard; acceptance remains unregularized KKT <=1e-9.
   Velocity-map delta <=1e-7 triggers evaluation of the final original-equation
   residual, it does NOT itself accept a step. Fix the candidate contact set for
   each substep. 34 regressions including whole-step squeeze/wall refinement pass.
