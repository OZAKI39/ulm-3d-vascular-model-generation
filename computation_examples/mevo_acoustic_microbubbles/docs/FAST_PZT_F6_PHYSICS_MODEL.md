# F6 physical model and evidence boundaries

## Stored acoustics and pressure gradient

`acoustics/kwave_backend.py` solves linear water acoustics with k-Wave and disc BLI
source weights, extracts the coherent positive-time phasor, and stores unit-source
basis fields. P(x)=sum a_k exp(i phi_k) P_k(x); p(t)=Re[P exp(i omega t)].
Actual F5/F6 input is the Stage-3 four **3 mm discs**, 2.5 MHz, original grid
256 x 128 x 128, 100 micrometre spacing, stored ROI 205 x 37 x 37. The unrelated
Stage-4 6 mm source configuration is not its identity. Input hashes are preserved.
Source amplitudes are ideal acoustic pressures, not hardware voltages.

Legacy `kernels.field_response` interpolates pressure and stored nodal gradients
separately. `consistent_field.query_numpy/consistent_primary` instead differentiate
the same trilinear pressure shape functions: dP/dx=sum P_abc (2a-1) w_b w_c / dx.
Cell faces use the positive-side cell, the domain maximum the final cell; the
coordinate quotient's floating-point floor defines the exact computational tie.
One-sided derivatives generally jump at faces. Invalid/nonfinite/outside positions
are rejected, never extrapolated. This resolves an interpolation inconsistency,
not acoustic mesh convergence, measured pressure accuracy or spatial resolution.

## Radial oscillator and response table

The pinned APECSS v1.8 source (a0ee071ca247b5489e997d90e83b5f53c40e4ec7) solves RPAR:

rho(R Rddot+3 Rdot^2/2) = pg(1-3 gamma Rdot/c)-p0-p_ac
                         -2 sigma(R)/R-4 mu_rad Rdot/R-4 kappa_s Rdot/R^2,
pg=(p0+2 sigma0/R0)(R0/R)^(3 gamma).

The reversible Marmottant law uses Rb=R0/sqrt(1+sigma0/chi),
Rr=Rb sqrt(1+sigma_water/chi), and sigma=0, chi[(R/Rb)^2-1], sigma_water
below/between/above those radii. It omits irreversible rupture history and
nonspherical deformation. APECSS's -P sin(omega t) is phasor iP; stored volume
harmonics are phase-rotated to real cosine input. Harmonic n rotates by n arg(P).
`ResponseSection` restricts actual valid cells to R0=1.78 um and f=2.5 MHz,
retaining eight harmonics. Nonzero out-of-table pressures are invalid; P=0 gives
exactly zero oscillatory response. No new radial table is generated for F6.

Tu2009's fitted SonoVue parameters describe particular measured bubbles and
conditions; the figure-derived manifold is not a universal blood/in-vivo calibration.
Sigma0=0.036 N/m is an explicit existing modeling assumption. Radial effective
viscosity 0.002 Pa s includes damping approximations; translation uses 0.001 Pa s.
See SONOVUE_PARAMETERS.md and APECSS_MODEL_AUDIT.md for exact provenance.

## Primary and secondary acoustic forces

For a sinusoidal external field, the cycle average is
Fp=-Re[V1 conj(grad P)]/2. The sign follows the actual complex APECSS response,
not a generic below-resonance statement. An in-phase sum at a selected point
maximizes coherent amplitude there for fixed individual magnitudes; it does not
establish a spatial extremum or stable translational equilibrium.

Far-field Fs_ij = rho omega^2/(8 pi) sum(n^2 Re[V_in conj(V_jn)]) dx_ij/d^3.
Positive signed volume-velocity correlation attracts; negative correlation repels.
F2 (`nearfield.py`) uses fundamental pair monopole feedback:
Vi + K chi_i Vj = Vi0, Vj + K chi_j Vi = Vj0,
K=rho omega^2/(4 pi d), chi_i=Vi0/Pi. Its signed correlation blends smoothly
back to the uncoupled value by d=4(Ri+Rj). Other harmonics are unchanged.
The separate soft core replaces d^3 by (d^2+a^2)^(3/2),
a=(Ri+Rj)(1-q)^2, q=min[d/(2(Ri+Rj)),1]. These are exploratory near-field
closures, not a collective scattering solution or fitted attraction multiplier.
Different pairs do not define a unique radial state for a bubble. Singular
closures and invalid effective pressures fail. Primary remains uncoupled.

## Flow, translation, lubrication and geometry

`flow_gpu.py` interpolates the unchanged 146230-tetrahedral prescribed water CFD
using strict P1 barycentric coordinates. Gradients are exact within each stored
tetrahedron, discontinuous across faces. No outside extrapolation. This does not
validate the CFD near-wall resolution. The slow affine flow experiment is explicitly
synthetic, not this CFD and not an independent blood-flow prediction.

F5 free velocity is uf+Fac/(6 pi mu R). Two local relaxed Jacobi normal-resistance
passes plus common-component mean balance approximate (D+Lpp)U=D uf+Fac.
Lubrication activates at h/R<0.1 and floors h/R at0.001; this is not a globally
converged many-body solution. F6 retains it as `f5_local`, offers independent
`diagonal_stokes`, and `rpy_free_space` (next section).
Euler prediction is followed by 12 central PBD nonpenetration passes and geometric
wall checks. PBD is a numerical diagonal-mobility constraint, not acoustic attraction,
adhesion, gas fusion or physical stress. WALL is the original open Y geometry;
INLET/OUTLETs are legal exits; acoustic support exit is a distinct data-domain event.
There is no wall hydrodynamic drag, streaming, two-way flow feedback or inertia
in the production model. `vessel_kernels.py` retains finite-triangle/swept checks.

## Independent RPY option

`rpy.py` implements equal-radius free-space no-slip translational RPY:
Mii=I/(6 pi mu R); for d>=2R,
Mij=[(1+2R^2/(3d^2))I+(1-2R^2/d^2)nn]/(8 pi mu d).
For 0<d<2R it uses [(1-9d/(32R))I+3d nn/(32R)]/(6 pi mu R).
The self term is included once; coincident centers fail. GPU full-pair multiplication
computes uf+M Fac. Old local lubrication is bypassed to avoid double counting;
no combined lubrication/RPY solver is claimed. No wall images, shear stresslets,
rotational mobility or finite-Re effects are included. A coated bubble's no-slip
sphere closure and diagonal PBD with RPY remain physical/contact approximations.
No group speed multiplier is applied and no acceleration is presupposed.

## Translation reference and evidence

Cycle-averaged production, instantaneous F/(6 pi mu R(t)), conservative added mass,
and no-slip history are distinct models. New local references use the existing
independent history-force solver and actual low navigation pressures. The old
22–163% result at 5–150 kPa cannot correct F5 or define its universal error.
Unsupported conditions reject correction requests: F6 correction_mode='none',
translation_calibration='UNCALIBRATED' until independent evidence supports a bounded
correction. Flow and contact displacements must never receive acoustic correction.

Formula derivation, software agreement, numerical convergence and physical validation
are separate. All production datasets remain EXPLORATORY; software tests do not
promote to MODEL_CHECKED, and no PHANTOM_VALIDATED claim is made.

## Control and time integration

`navigation.BoundedFeedback` chooses only four nonnegative source amplitudes and
three relative phases. The amplitude cap is the inherited 291.237687721856 Pa per
source, not a voltage limit. Every 25 ms it minimizes the mismatch between the
centroid's table-based primary Stokes drift plus the prescribed background flow
and 2 s^-1 times the target offset. A soft local pressure preference above100 Pa
helps avoid removing attraction; it is neither an added force nor a cohesion
guarantee. All32 separate bubble coordinates are still integrated with pair forces,
mobility and contact. Cohesion is assessed from those coordinates after integration.
The cheap controller does not predict the full constrained many-body trajectory.
The same controller model is retained in the RPY and legacy-gradient ablations;
they test robustness of this controller, not the best possible controller in each
alternative model. Static preflight additionally probes the target and six2µm
offsets for restoring drift. Neither fixed phase alignment nor a local optimizer
proves global reachability or Lyapunov stability.

`adaptive.StepDoubling` optionally compares one full Euler/contact step with two
half steps, accepts the two half steps only within a maximum same-ID position
error budget (2 nm +1e-4 times Rg by default) and unchanged contact budgets, and
rolls back positions, diagnostics, clocks, boundary-event count and flow caches
on rejection. Trials and accepted physical steps are counted separately. The
inherited control-event segmentation prevents crossing an unhandled switch.
The minimum accepted half step is1.25µs; exhaustion fails explicitly. This host
controller is not a claimed speed improvement and is independently admitted;
its combination with feedback is deliberately rejected until verified.
Formal navigation uses fixed10µs, with three controls preregistered at2.5µs after
10/5µs contact-budget pilot failures. Exact10/5µs full-horizon feedback comparison
is independent of this admission. PBD iterations never count as physical steps.

Local translation references use imposed pressure305 and2038.7896 Pa at2.5MHz,
with ON/OFF/ON at48/72 acoustic cycles,96cycles total (38.4µs). They resolve local
switching transients; they do not reproduce full25ms feedback or0.3–1s navigation
memory. The single-bubble reference uses a fixed standing-wave gradient kP;
the pair reference uses identical prescribed radial responses at fixed10R spacing,
where the production near-field blend is zero. Mutual radial scattering, walls
and moving pressure sampling are absent. `reference_audit.py` separately compares
the resolved-cycle mean with the actual unchanged stored response table; this
prevents calling a newly resolved radial baseline identical to table interpolation.
No correction is inferred outside these narrow conditions.

## Primary sources

- [kWaveArray](https://k-wave.org/documentation/kWaveArray.php), source representation.
- [APECSS v1.8](https://github.com/polycfd/apecss/tree/v1.8), pinned radial source.
- [Marmottant 2005](https://doi.org/10.1121/1.2109427), shell model.
- [Tu 2009](https://pmc.ncbi.nlm.nih.gov/articles/PMC2803720/), measured fit scope.
- [Mettin 1997](https://doi.org/10.1103/PhysRevE.56.2924), far-field pair force.
- [Wajnryb 2013 author PDF](https://www.fuw.edu.pl/~piotrek/publications/Rotne_Prager_Yamakawa.pdf), Eq.3.12 equal-sphere translation and overlap regularization.
- [Garbin 2009](https://pages.jh.edu/aprospe1/publications/PapersPublished/GasBubbles/GarbinHistory.pdf), variable-radius history reference.
- [Cattaneo and Supponen author preprint v2](https://arxiv.org/html/2506.16413v2), translation reference scope; its experiment is not reproduced here.
