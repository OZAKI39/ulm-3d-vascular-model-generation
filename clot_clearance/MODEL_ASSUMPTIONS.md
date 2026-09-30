# Model, units and numerical limits

This is a local **research/software prototype**, not an experimentally calibrated model of ultrasound thrombolysis. The supplied paper studies aspiration thrombectomy. See [paper mapping](references/PAPER_READING.md) for the exact differences.

## A. Established framework used here

Positions are m, time s, velocity m/s, density kg/m³, stress/energy density Pa, particle force N, volume m³. The fluid stress is `sigma = -p I + mu (grad(u) + grad(u).T)`. With the **outward solid normal**, positive pressure produces inward traction. At each exposed voxel face, `f_i += traction * h²`; divide by `rho*h³` for acceleration. Multiple faces at corners are integrated separately. There is no arbitrary force amplification.

The particles have uniform volume `V=h³`. Neighbor pairs are unique unordered `i<j`, with `xi=Xj-Xi`, `eta=xj-xi`, horizon `delta=3.01 h`, and influence `w=(1-|xi|/delta)²`. The active influence is `q=w*integrity`.

`K_i = sum_j q_ij (xi tensor xi) V_j`; `F_i = [sum_j q_ij (eta tensor xi) V_j] inv(K_i)`. Damage does **not** remove the identity part of F in an undeformed solid: both numerator and K use the same active weights. K is recomputed after integrity changes; positions do not change the reference neighborhoods.

The initial constitutive model is compressible Neo-Hookean:

```
lambda = bulk_modulus - 2*shear_modulus/3
W(F) = G/2 (tr(F.T F)-3) - G log(J) + lambda/2 log(J)^2
P(F) = G (F-F^-T) + lambda log(J) F^-T
```

Main values are **synthetic** `G=1,000 Pa`, `K=50,000 Pa`, corresponding to small-strain `nu≈0.49007`. They are not copied from the paper's generalized Ogden fit. Stress is degraded by `g_i=1-D_i`, where `D_i` is the volume-weighted mean damage over the initial neighbor family. This is the explicit choice `alpha(D)=D` in the spirit of the paper's Eq.14.

The continuum contribution from an unordered pair to force on i is

```
f_ij = Vi Vj q_ij [g_i P_i inv(K_i) + g_j P_j inv(K_j)] xi_ij
```

The force on j is `-f_ij`. These signs follow the variation of the discrete energy `sum_i Vi*g_i*W(F_i)` with integrity held fixed. Tests check numerical energy derivatives, force balance, torque balance and frame objectivity.

## B. Deliberate modeling approximations

- Straight empty pipe: radius 1.25 mm, length 6 mm, `Q=18 mL/min`, BraVa's fluid `rho=1056`, `mu=0.00345312`. The background is **analytic Poiseuille**, exported as a tetrahedral P1 snapshot. No new Navier–Stokes solve was performed. The clot does not block or change this background.
- A `1.5 × 1.0 × 0.5 mm` block lies on an idealized flat, supported wall patch inside the cylinder. Its bottom cell layer is fixed exactly. This flat support approximates wall attachment; it is not a conforming cylindrical adhesion/contact model.
- All exposed initial faces except the bottom have reference quadrature locations, reference area and fixed normals. Traction is a **dead load**, evaluated on this reference surface. Large-deformation follower traction (Nanson mapping), newly exposed crack-face loading and moving-domain resampling are not implemented.
- The local synthetic traction is a Gaussian around a prescribed bubble center. The axial direction is projected onto each face tangent; normal traction is separately configurable. `traction_scale_Pa` is an independent test parameter. When it is null, the optional dimensional scale is `mu*streaming_velocity_scale/streaming_length_scale`.
- Bubble radius and ultrasound frequency do not solve or predict the traction magnitude; radius is metadata for a future resolved provider. There is no acoustic wave, coated bubble dynamics, cavitation, real microbubble cloud, RBC, clot chemistry or thrombolytic drug model.
- Fixed lower particles constrain displacement and velocity, and can accumulate incident bond damage; their support condition never delaminates.
- `FileStreaming` supports direct traction CSV/point VTP, or volumetric CSV and linear tetrahedral VTU/legacy VTK velocity+pressure. Tetra import uses the **original** containing cells, original P1 interpolation and copied BraVa gradient algebra. Outside-domain queries fail. Scattered CSV/point VTP uses a bounded convex sample domain; it does not infer vessel topology or curved-surface tangential derivatives. No silent nearest-neighbor extrapolation is used. Surface velocity alone is insufficient for the 3-D fluid gradient and is rejected.
- Imported fields are **static SI snapshots**. Time sequences and arbitrary SV element types are future work. Field names are configurable (`Velocity`, `Pressure`, `traction_Pa` by default). Pressure must be a consistent applied gauge relative to the supported clot; offsets are not automatically removed.

## C. Unvalidated components and numerical safeguards

### Provisional stabilization

Define the residual `r_ij=eta_ij-F_i*xi_ij`. Add the objective non-affine energy density

```
H_i = (c_i/2) sum_j q_ij |r_ij|² Vj
c_i = stabilization_alpha * G * g_i / trace(K_i)
```

Its force-state contribution is `q*c_i*r_ij`; assemble `Vi Vj q (c_i*r_ij + c_j*r_ji_signed)` with the latter residual expressed using the same i→j orientation. The residual's weighted first moment vanishes, so the derivative terms through F cancel. The penalty is zero for affine deformation and positive for tested non-affine perturbations. This full-vector energy penalty differs from the supplied paper's Eq.10. `stabilization_alpha=3` is provisional, **not calibrated or proven converged**. Setting it to zero disables it; this is not recommended for production. Stabilization energy is exported separately and is appreciable in this coarse example, so it must not be concealed as material energy. Tests of several motions are not a complete stability proof.

### Damage and two clocks

`instantaneous_stretch` uses stretch ratio `s=|eta|/|xi|`, thresholds `s1<s2`, a cubic smoothstep integrity transition and `new_integrity=min(old,candidate)`. Updates occur at accepted explicit steps. It intentionally differs from the paper's linear transition.

`cyclic_accumulation` measures `Q=(max(s)-min(s))/2` over resolved representative cycles, then applies

```
D_new = min(1, D_old + DeltaN*C_damage*max(Q/Q_ref - 1,0)^m_damage)
```

Only the stretch-amplitude driver is implemented; unsupported driver names fail. Damage is held fixed while measuring a representative cycle, then updated at the macro boundary. For the main run the first mechanical cycle is warmup and excluded from Q. Subsequent cycles are transient responses of the current damaged state; periodic steady response is not assumed.

The main run resolves a **500 Hz mechanical proxy** with `dt=2 microseconds`, while counting 1,000 **1 MHz carrier cycles** per macro step. Eight macro steps give N=8,000 and formal `N/f_carrier=0.008 s`, while the resolved proxy dynamics total 0.018 s including warmup. These clocks are intentionally distinct. The substitution does not justify transferring 500 Hz strain amplitudes to MHz dynamics, and neither clock establishes a clinically meaningful clearance time. `C_damage=2e-5`, `m=2`, `Q_ref=0.014`, and `DeltaN=1000` are artificial software parameters.

More time cannot heal any individual bond. The tested stronger-loading control increases overall damage, but no universal ordering of separate nonlinear trajectories is asserted. Rate zero and disabled damage preserve all integrity variables. Energy released by a damage jump is not balanced against a calibrated fracture energy; jump-size convergence and fatigue calibration remain open.

### Numerical protection

- Damped velocity Verlet uses half exponential damping steps with rate `3500 s^-1`. This is artificial damping, not measured clot viscoelasticity.
- A conservative **estimate**, not an eigenvalue proof, is `dt_limit = 0.2*h/sqrt((K+(4/3+2*alpha)*G)/rho)`. Main limit is approximately 3.39 microseconds; actual dt is 2 microseconds. A small 96-particle two-step comparison with 1 microsecond is reported separately; it is not full mesh/horizon convergence.
- Nonzero families with fewer than four neighbors or shape condition number above `1e8` fail loudly. There is no pseudoinverse, silent shape regularization or automatic deletion of damaged particles/bonds.
- A truly isolated particle has no continuum F or material stress. The output has `deformation_gradient_valid=0` and a zero sentinel for J. It remains a particle with mass and any prescribed force. Connected but rank-deficient fragments **stop the run**. This can limit highly damaged simulations.
- Nonfinite states, `J<=0.2`, acceleration above `1e7 m/s²`, or displacement above 0.4 mm stop with `FAILED.json` and `failed_state.npz`. Existing output directories are refused.
- Fluid/fragment transport, contact, collision, coagulation, reconnection and physical particle removal are absent. “Erosion” can only mean loss of bond connectivity; it is not a measured removed mass.

### Connectivity, not appearance

Fragments are graph components with edges `integrity > connectivity_threshold`; default threshold is exactly zero. All particles, including isolated and fixed points, count. The original 384-particle demo has **one** component and no broken bonds. A separate deliberately imposed cleavage test produces two components; it is not a hydrodynamic fracture result. Largest-fragment volume and detached-from-base volume are computed from particle reference volumes, with no deletion of mass.

## D. Next work

1. Import a small **real microbubble-resolved** tetrahedral SI snapshot containing velocity and pressure, or a clot-surface traction field; verify normal orientation, pressure reference, spatial coverage, field gradients and net force before running PD. Existing BraVa frozen flow contains no resolved acoustic microstreaming and cannot simply be renamed as such.
2. Add timestamped field sequences and conservative moving-surface transfer with deformation-aware normals/areas.
3. Implement and verify a first-order/generalized Ogden material, then experimental clot/fatigue calibration, horizon/particle/jump convergence and robust fragment treatment.
4. Only subsequently add clot-wall adhesion, fragment-fluid transport, coated-bubble/acoustic physics and two-way CFD–PD coupling.
