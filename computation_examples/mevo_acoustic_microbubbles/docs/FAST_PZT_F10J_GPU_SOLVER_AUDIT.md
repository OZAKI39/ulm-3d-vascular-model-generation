# F10J GPU solver audit

The original saved 32-bubble failure is reproduced exactly. Both the base PCG and **the first contact-direction PCG** exceed their own targets. The last direction is accurate and descending; total-energy roundoff rejects all 20 Armijo candidates. The original mathematical stopping rule was retained and met after the correction.

## Actual equation and units

Let R=1.78e-6 m, dt=1e-5 s, zeta=6 pi mu R, mu=0.001 Pa s, and kc=0.04641319518378658 N/m (F10G 2kc, exploratory). The initial-step resistance is dimensionless A=I+sum c_ij B_old^T B_old. The code uses c=0.5[phi(max(h/R,0.001))-phi(0.1)] for h/R<0.1; phi(e)=0.5/e-0.45 log(e)-(3/56)e log(e). It is fixed during the contact solve. Free velocity already includes the archived external acoustic/background inputs; no new force or flow is generated.

With x_s=(x-origin)/R, b=dt*free/R, A u=b, pred=x_s+u, q dimensionless, alpha=dt kc/zeta, y=pred+q:

`E(q)=0.5 q^T A q + 0.5 alpha sum_pairs max(2-|y_j-y_i|,0)^2`.

`g=Aq-alpha f(y)`; physical contact force is `kc R f(y)`. The physical equation is `x_new-x_old = dt A^-1 [free + F_contact/zeta]`. The global-coordinate residual in metres is evaluated independently from actual returned positions and real pair geometry.

Every outer iteration updates pair distances, normals and the unilateral active set. A stays frozen. The direction solves `H d=-g`, `H=A+alpha B_active^T B_active`, using a block Jacobi inverse of each diagonal 3x3 block. H is a positive Gauss-Newton metric, **not** the full nonlinear Hessian and not a proof of global convexity or unique physical equilibrium. Measured symmetry errors are at most 2.842e-14; min eigenvalues are approximately 1; all block inverses are positive definite; no used direction has nonpositive PCG curvature. A's maximum eigenvalue is 1513.655982. CPU and GPU A differ by 2.359e-10 in absolute matrix entries from coordinate evaluation roundoff; this is not a physical model change.

## Three different residuals

1. PCG recursive r is updated by the recurrence; raw norms are dimensionless and compared with `max(1e-14,1e-12 ||b||_2)` for that solve's RHS.
2. True linear r is recomputed on the GPU as b-Ax or -g-Hd using the actual operator and GPU solution. The preconditioned diagnostic is sqrt(r^T M r), **not** the unpreconditioned stopping norm. Relative residual divides by RHS norm. All arrays and independent dense-reference linear errors are archived.
3. Full nonlinear position residual is independently recomputed in metres from the actual final contact forces and A inverse. It is not a PCG residual.

For exact base prediction, `R ||A^-1 g||_2 <= R ||g||_2` since lambda_min(A)>=1. With approximate base solve, include the base residual: `||position residual||_2 <= R (||g||_2+||r_A||_2)` up to coordinate rounding. The old quantity R||g|| is a conservative **equation-residual** bound, not the true position error relative to a nonlinear root. F10J does not reinterpret the failed 1.850e-14 m bound as actual position error, and does not change the stopping rule. Independent root comparison supplies separate observed error evidence on this fixed state; no global a priori error theorem is claimed.

|Original linear solve / outer|Recursive norm|True norm|Original target|Decision|
|---|---:|---:|---:|---|
|A / None|1.254988e-12|1.253710e-12|5.189764e-13|FAIL|
|H / 0|5.930708e-13|5.929474e-13|2.063341e-13|FAIL|
|H / 1|6.171999e-15|6.172066e-15|1.000000e-14|PASS|
|H / 2|6.175399e-15|6.175486e-15|1.000000e-14|PASS|
|H / 3|6.802414e-15|6.802416e-15|1.000000e-14|PASS|
|H / 4|8.369456e-15|8.369456e-15|1.000000e-14|PASS|

Base recursive/true residuals are 1.254987782e-12 / 1.253709625e-12. They closely agree: this is an insufficient bounded iteration budget, not a catastrophic recursive-residual drift. Its independent position-equivalent linear error is only 7.217151e-20 m, which also shows why a residual norm cannot be equated to solution error.

## Armijo diagnosis and minimal correction

The inherited rule is `E(q+t d) <= E(q)+1e-4 t g.d`, at most 20 halvings, at most 24 accepted updates. Final original g.d=-2.547748395e-17<0. All 20 failed candidates retain the same 38 active contacts. For the full candidate, direct float64 total-energy subtraction gives **+1.023486851e-16**; independent extended-precision stable difference gives **-1.328128981e-17**; Armijo requires at most -2.547748395e-21. The rejected full step has actual equation residual 1.107076461e-16 m and original gradient bound 7.885565382e-16 m. Thus the stopping rule itself need not be loosened: the original code never accepted the available sufficient step.

`f10j_gpu.py` subclasses the entire original `FixedStateB2.solve`; it does not duplicate or replace the F10I solver. It reuses resistance, contact forces, normals, block preconditioning, PCG and original Armijo control flow. The bounded PCG budget becomes 256; true and recursive residual checks both gate admission. Conditional execution actually needs 132 base iterations and 93 for the final direction. No residual tolerance is loosened. A deliberately one-iteration base solve is rejected with TRUE_OR_RECURSIVE_LINEAR_RESIDUAL.

The merit is translated so current energy is zero and trial energy is the algebraically identical increment:

`Delta E = q^T A s + 0.5 s^T A s + alpha/2 sum (delta_new-delta_old)(delta_new+delta_old)`.

Distances use the rationalized increment `Delta d=(2 dx.ds+ds.ds)/(|dx+ds|+|dx|)`. For contacts staying active, `Delta delta=-Delta d`. This avoids subtracting two nearly equal total energies and two nearly equal distances. The inherited acceptance inequality then compares Delta E directly with 1e-4 t g.d. All 24 original candidates agree with the independent long-double acceptance calculation in the actual GPU kernel check. Maximum absolute energy-increment discrepancy across the different-sized increments is 1.708e-18; no acceptance decision differs.

## Numerical evidence and scope

Original 18 small fixtures: 12 two-bubble and 6 four-bubble tests pass, including two expected overlap-budget rejections. The real 32-bubble case passes, as do nonfinite input, coincident-center, severe-overlap and incomplete-PCG rejection checks. Final conservative bound 7.885584867e-16 m; base true residual 1.191054933e-13; overlap 17.497568596 nm; maximum one-update correction 1.931690282 nm. Original budgets 17.8/178 nm are unchanged. Contact/edge sets match CPU; reconstructed force mismatch 1.615687e-20 N; net reaction 8.507188e-23 N. Nonnegative pair forces and zero force for separated pairs follow the unchanged geometric unilateral kernel and are audited against returned geometry.

GPU vs original CPU (atol=1e-14): position 2.605195605e-15 m; force 1.319722056e-18 N. GPU vs independent force-loop/Powell reference: position 1.155703174e-16 m; force 5.851729339e-20 N. Actual position equation residual 1.107182785e-16 m. Both references solve this same exploratory discrete model; neither is experimentally calibrated truth.

## Synchronization semantics

Old intermediate reads: initial outer checks (5), after-direction checks (5), candidates (24), final check (1): 35. Corrected: 6 outer checks + 5 after-direction + 5 trials + final =17. Optimized reuses the **accepted** candidate residual for five immediately following outer checks:12. Rejected trials are never reused. Positions are bitwise identical across all 19 fixtures. Final host reads remain 12 old /14 new, including full vectors/operator evidence and both PCG checks; explicit prep/final synchronizations also remain.

Warp 1.18's actual `_run_capturable_loop` supports `capture_while`. Conditional bounded PCG graphs replace fixed 256-iteration captured loops without changing their criterion. Zero-RHS fixtures do not move; converged states are not subsequently changed. The original fixed-loop masked CG kernels also froze x after convergence, but still executed unused work. No new graph framework, GPU-resident nonlinear loop, batch integration or production default was added. Failed PCG is rejected at final admission; this fixed-state module is not a streaming accepted-step integrator.

Read-only instrumentation and ordinary original execution give exactly identical positions, also identical to archived F10I v2. Raw traces: `audit_v1/AUDIT.json`; independent audit: `analysis/ORIGINAL_AUDIT.json`; numerical runs: `corrected_v1`, `optimized_v2`; negative and merit checks: `guards_v1`.

Primary implementation references: [Warp 1.18 CG source](https://github.com/NVIDIA/warp/blob/v1.18.0/warp/_src/optim/linear.py), [SciPy CG stopping rule](https://docs.scipy.org/doc/scipy/reference/generated/scipy.sparse.linalg.cg.html). Actual installed Warp source and SHA are archived in `audit_v1/WARP_CAPABILITIES.json`; project equations and independent checks remain the primary evidence.
