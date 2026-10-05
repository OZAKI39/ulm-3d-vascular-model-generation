# APECSS lipid-model audit

Source: official `polycfd/apecss`, latest GitHub release/tag verified as **v1.8**, pinned commit
`a0ee071ca247b5489e997d90e83b5f53c40e4ec7` (2025-07-21). The discovered development HEAD was
`41392784d3c7196e6ee877e9b61055da0acbcb0c`; it is not the production dependency.
Repository/API discovery is retained in `third_party/APECSS_DISCOVERY.json`.

## Actual implemented choices

| Option | Source | Mathematical identity | Use |
|---|---|---|---|
| `Marmottant` | `src/interface.c`, `apecss_interface_surfacetension_marmottant` | Piecewise zero / area-elastic / clean-water tension | Selected reproducible nonlinear baseline |
| `Gompertz-Marmottant` | `src/interface.c`, `apecss_interface_surfacetension_gompertzmarmottant` | Smooth Gompertz approximation of the piecewise law | Distinct optional model, never silently substituted |
| `RPAR` | `src/rayleighplesset.c`, `apecss_rp_rayleighplessetacousticrationvelocity_ode` | Rayleigh–Plesset plus gas-pressure acoustic-radiation damping | Matches the radial equation used by the simple coated example |
| `KM`, `Gilmore` | same file, separately named functions | Different compressibility approximations | Must retain separate metadata; not automatically equivalent to RPAR |

For `Marmottant`, `src/bubble.c` computes
`Rb = R0/sqrt(1 + sigma0/chi)` and `Rr = Rb*sqrt(1 + sigma_water/chi)`.
The tension is `0` below Rb, `chi*(R²/Rb²-1)` between Rb and Rr, and `sigma_water` above Rr.
The shell viscous pressure is `4*kappa_s*Rdot/R²`.
Units: R in m, sigma and chi in N/m, kappa_s in kg/s, pressures in Pa.

The ideal-polytropic gas is initialized with `pg0 = p0 + 2*sigma0/R0`.
With `RPAR`, the term `R*pg_dot/c` on the radial-equation right-hand side produces
the factor `(1 - 3*kappa*Rdot/c)` multiplying the polytropic pressure.
This maps to **Marmottant 2005 Eq. (3)**, with the monolayer normal stress in **Eq. (2)**.

## Is this the complete 2005 model?

No blanket equivalence is claimed. It implements the widely used reversible three-branch
specialization of **Eq. (4)** with break-up tension equal to the clean-water value.
The paper also allows a distinct, higher break-up tension and a history-dependent transition
after rupture. APECSS's stateless R-only law has no independent sigma_breakup parameter
or irreversible rupture memory. Its simple Fig. 5b example does not require implementing
that broader resistant-bubble scenario. A solver of that broader scenario is not added here.

Gompertz–Marmottant is explicitly attributed by upstream to Guemmer et al. (2021).
It is continuously differentiable and does not have an exactly flat zero-tension branch.
It requires a usable strictly positive sigma0; the log(sigma0/sigma_water) mapping cannot
be blindly used for the buckled sigma0=0 official simple case.

## SonoVue-like suitability and parameter mapping

This supplies a nonlinear lipid-monolayer radial backend, not a universal calibrated SonoVue model.
`Elasticity` ↔ chi; `DilatationalViscosity` ↔ kappa_s; `SigmaInit` ↔ sigma(R0);
`SurfaceTensionCoeff` ↔ clean-water tension; `PolytropicExponent` ↔ gas exponent.
Radius, fluid properties, ambient pressure and gas law remain explicit inputs.
Literature-fitted values, their size/pressure/temperature conditions and uncertainty must be
mapped separately in `SONOVUE_PARAMETERS.md`; choosing the model does not validate its parameters.

## Numerical and phase details

The original C library uses adaptive Dormand–Prince RK5(4). This project wraps the unmodified
official ultrasound executable and standardized options/results files; it does not rewrite the ODE.
The default APECSS_FLOAT is double (verified at build time), Release uses -O3 with -march=native,
without -ffast-math. Every subprocess has its own working directory.

`src/bubble.c:apecss_bubble_pressure_infinity_sinexcitation` uses
`p_inf - p0 = -pa*sin(omega*t)`. Under our `Re[P exp(+i omega t)]` convention, its input
phasor is **i*pa**. A response measured from this waveform must be divided by i before being
stored for a positive-real cosine pressure; all higher harmonics require the corresponding
harmonic phase shift. Ignoring this gives a 90-degree force error.

Two documentation issues are treated as documentation, not copied into a new solver:
the RPAR discussion writes a small Mach number followed by approximately 1;
the Gompertz derivative printed in the TeX uses R in a denominator where the implemented
chain rule uses Rb. The source formulas and direct checks are authoritative for the wrapper.
For the piecewise law, tension itself is continuous at the two breakpoints, while its derivative jumps.

## Verification evidence

`examples/ultrasound/lipidcoated_simple.apecss`: R0=0.975 µm, 2.9 MHz, 130 kPa,
chi=1 N/m, kappa_s=15e-9 kg/s, sigma0=0, SF6 exponent 1.095, water viscosity 1 mPa s.
Upstream identifies this with Marmottant Fig. 5b, but supplies no digitized radius reference.
Execution and tolerance convergence are software checks; figure agreement must be reported separately.

The Rayleigh-collapse emissions example does supply four reference arrays and an official check
at rtol=1e-6, atol=3e-16. Exact commands, stdout, run times and comparison outcomes reside in
`results/stage2/apecss_validation/STATE.json` and its timestamped run directory.

Sources: [APECSS v1.8](https://github.com/polycfd/apecss/tree/v1.8),
[APECSS paper](https://doi.org/10.21105/joss.05435),
[Marmottant 2005 author manuscript](https://liphy-annuaire.univ-grenoble-alpes.fr/pages_personnelles/philippe_marmottant/Publications_files/Marmottant2005_JASACoatedBubble.pdf).
