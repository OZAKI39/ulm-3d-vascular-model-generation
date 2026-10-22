# Lazarus 2017 — extraction completed; simulation gated

Full 38-page accepted manuscript read, including equation images:
[Imperial accepted manuscript](https://spiral.imperial.ac.uk/server/api/core/bitstreams/02e45b2e-b3bb-4dc7-a1ee-02f0d3c85649/content),
[DOI 10.1121/1.5010170](https://doi.org/10.1121/1.5010170).
Repository contains manuscript and thumbnail; supplementary videos are referenced, not a public numerical solver.
No usable author code was found in title/DOI/repository searches. No contact with authors was made.

The experimental agent is **in-house perfluorobutane/DPPC-DPPA-PEG lipid**, not SonoVue.
The experiment is a 2D population resting at the top glass wall of a 0.6 mm ID capillary.
1 MHz, 50/75/100 kPa, 20 ms; optical field 230×100 µm, 40 kfps. Most diameters 3–6 µm;
concentration ≈10^7/ml. Measured positions and radii from individual video frames provide initial conditions,
which are not supplied as machine-readable N>40 populations. Fig5 has N=47, not an arbitrary random 3D cloud.

Equation1 uses Fij=(2πρ/9)(ωPa)² β² Ri³Rj³/dij² eij, drag −6πµRi vi,
added mass (2π/3)ρRi³ and gas mass (4π/3)ρPFB Ri³.
Table1: β=5e-7 Pa^-1, µ=.001 Pa s, ρ=1000, ρPFB=9.7 kg/m³, c=1500 m/s,
contact gap100 nm, empirical rearrangement time0.5 ms, Matlab ode45 time increment0.1 µs.
It does **not integrate radial dynamics**: fixed mean radii, all oscillations in phase, below resonance.
The paper explicitly extends the far-field formula beyond its validity, then imposes inelastic collisions,
rearranges smaller bubbles around the largest and immobilizes cluster members in the cluster center frame.
These empirical operations are not validated near-field hydrodynamics and are not our geometric contact solver.

Cluster metric: ≥2 bubbles with surface gap<100 nm; cluster position is volume-weighted center.
Number of clustered bubbles is fit with A(1−exp(−Bt)); normalization distance is sqrt(S/N).
This differs from nearest-neighbor spacing and from reporting only the largest connected component.
The paper's reported normalized B units include a bubble count although B in that exponential has units1/time;
retain this ambiguity in any future digitization instead of silently changing the metric.

50/75/100 kPa configs and NOT_RUN status are recorded in `results/stage2/lazarus2017/`.
Reason: Stage2G Zhang Fig3 validation fails, so the user's explicit gate forbids large aggregation.
Additional reproduction limitations: missing measured population coordinates/radii, wall/history drag,
shell/compressibility difference, planar dimensionality and imposed cluster rearrangement.
No pressure curves or fitted clustering rates have been fabricated. Figure D contains gate status only.
