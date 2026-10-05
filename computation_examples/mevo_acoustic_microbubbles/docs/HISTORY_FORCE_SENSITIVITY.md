# Hydrodynamic history force sensitivity

## 1. Why this comparison was performed

Quantify the approximation in slow-time Stokes transport without adding memory to ParticleTracker. This is a controlled model sensitivity study, not reproduction or calibration of experimental displacement.

## 2. Models

Aprod uses the steady APECSS cycle-averaged primary force and equilibrium radius, matching the production translation law (table interpolation error is excluded). Aresolved uses instantaneous F(t)/(6πμR(t)). B adds conservative added mass to Aresolved. C adds the no-slip zero-Re history force to B. Aresolved, B and C share the same radius, forcing, time grid, medium and zero initial velocity. Aprod is algebraic and has no independently specifiable initial velocity; its startup/averaging difference is explicitly separate. Only B versus C isolates history.

## 3. Equations

For prescribed R(t), V=4πR³/3 and M=ρV/2, B solves d(Mv)/dt+6πμRv=F. C adds 6πμ/√(πν) ∫[dq/dτ]/√(θ(t)−θ(τ)) dτ to the left side, where q=Rv and θ(t)=∫R⁻²dt. Gas inertia is neglected in both B and C. F=−V ∂p/∂x. The initial history is zero. No finite-Re drag correction is used.

## 4. Literature sources

[Garbin et al. (2009), DOI 10.1063/1.3227903](https://pages.jh.edu/aprospe1/publications/PapersPublished/GasBubbles/GarbinHistory.pdf), full seven-page paper, Eq. (1), provides the variable-radius no-slip balance and transformed-time singular-endpoint treatment. Its experiment used BR-14 pairs and an eight-cycle, 2.25 MHz tapered pulse; its radius-based Reynolds definition differs by a factor of two from ours. Its model excludes close-gap boundary-layer interactions.

[Cattaneo & Supponen (2026), DOI 10.1103/p475-xggm](https://journals.aps.org/prfluids/abstract/10.1103/p475-xggm) was identified as PRFluids 11, 033606, published 16 March 2026. The fully read equation source is [author preprint v2, 7 February 2026](https://arxiv.org/html/2506.16413v2), Eqs. (11), (14); identical final typesetting is not asserted. The paper studies isolated SonoVue bubbles (1.8–4.9 μm), 20-cycle 1.5 MHz pulses, and initial-to-final straight-line displacement. It reports comparable zero-/finite-Re history predictions up to mean diameter-based Re about 2. This supports selecting the simpler reference, not extending its experimental validation to our cases. Public dataset metadata were accessible; restricted files were not downloaded.

## 5. SonoVue conditions

H1 and H2 use the existing Tu2009 Fig. 4 parameter manifold at R0=1.78 μm and 2.5 MHz. H3 uses the published Table 2 SonoVue values χ=0.30 N/m and κs=3.2e−9 kg/s at 150 kPa. All use σinitial=0.036 N/m as the existing model assumption, not a new fit. The APECSS radial effective viscosity is 0.002 Pa s; translational physical viscosity is 0.001 Pa s. Water density/sound speed are 1000 kg/m³ and 1500 m/s. These are not the exact Garbin or Cattaneo experimental conditions.

| Case | Pressure (kPa) | max R/R0 | min R/R0 |
|---|---:|---:|---:|
| H1 | 5 | 1.010940 | 0.988912 |
| H2 | 20 | 1.043206 | 0.955253 |
| H3 | 150 | 1.512558 | 0.551530 |

## 6. Exposure duration

20 cycles (8 μs) and 0.1 ms (250 cycles). APECSS supplies the complete 48-cycle startup. Only its final verified periodic cycle is repeated afterwards. The local standing-wave pressure is −P sin(ωt), with gradient −kP sin(ωt), representing a fixed-position standing-wave point with gradient/pressure ratio k. There is no pulse-off interval, dissolution or moving-position feedback. Longer exposure is a numerical sensitivity scenario, not a prediction of long-pulse bubble survival.

## 7. Numerical convergence

Uniform 64, 128 and 256 samples per cycle. BDF2 advances added momentum, using backward Euler at the first step. Exact interval integration of piecewise-linear q in θ resolves the Abel singularity without a cutoff. Independent manufactured-variable-radius kernel and fixed-radius mass/force tests verify the discretization. Endpoint change threshold: 2%; it was fixed before running. Both short and long endpoints of every model are checked.

| Case/window | Largest final-refinement displacement change (%) |
|---|---:|
| H1/short | 0.1021 |
| H1/long | 0.1700 |
| H2/short | 0.1427 |
| H2/long | 0.1776 |
| H3/short | 0.1790 |
| H3/long | 0.3300 |

## 8. Displacement comparison

| Case/window | Aprod (μm) | Aresolved (μm) | B (μm) | C (μm) | Aprod–C error (%) |
|---|---:|---:|---:|---:|---:|
| H1/short | 0.00362741 | 0.00235124 | 0.00540388 | 0.00250951 | 44.55 |
| H1/long | 0.0453426 | 0.0301609 | 0.0447421 | 0.0213095 | 112.78 |
| H2/short | 0.0578846 | 0.0375041 | 0.059187 | 0.0220101 | 162.99 |
| H2/long | 0.723557 | 0.481174 | 0.68801 | 0.316549 | 128.58 |
| H3/short | -4.46197 | -2.50965 | -4.91131 | -1.83671 | 142.93 |
| H3/long | -55.7747 | -32.9872 | -68.3109 | -45.8147 | 21.74 |

All mean/peak velocities, mean/max Reynolds numbers, half-reference-displacement first-crossing times, and shared input hashes are in the linked machine-readable summary. First crossings can be affected by intra-cycle oscillation; they are not sustained-arrival times.

## 9. Force comparison

| Case | C primary RMS (nN) | C drag RMS (nN) | C added-mass RMS (nN) | C history RMS (nN) | C mean/max Re |
|---|---:|---:|---:|---:|---:|
| H1 | 0.875075 | 0.0726796 | 0.400104 | 0.509912 | 0.00692/0.01523 |
| H2 | 3.52642 | 0.309547 | 1.63777 | 2.0372 | 0.02847/0.0615 |
| H3 | 51.0596 | 16.1029 | 47.5648 | 49.2392 | 1.538/3.254 |

## 10. Estimated approximation error

| Case/window | Averaging/startup error (%) | Added-mass comparison (%) | History-only B–C (%) | Classification of total Aprod–C |
|---|---:|---:|---:|---|
| H1/short | 54.28 | 56.49 | 115.34 | quasi-steady approximation requires explicit limitation |
| H1/long | 50.34 | 32.59 | 109.96 | quasi-steady approximation requires explicit limitation |
| H2/short | 54.34 | 36.63 | 168.91 | quasi-steady approximation requires explicit limitation |
| H2/long | 50.37 | 30.06 | 117.35 | quasi-steady approximation requires explicit limitation |
| H3/short | 77.79 | 48.90 | 167.40 | quasi-steady approximation requires explicit limitation |
| H3/long | 69.08 | 51.71 | 49.10 | quasi-steady approximation requires explicit limitation |

Errors use the reference named by each comparison as denominator, floored at 1e−15 m only for reporting; they are not additive. The classification thresholds 5%, 10%, 20% are reporting conventions. A small radial excursion alone does not guarantee small translational approximation error.

## 11. Implications for production simulator

The tested errors require an explicit limitation on quantitative transient displacement and arrival-time predictions. Production remains a slow-time engineering model for control-induced migration and clustering trends. No history term is imported into the tracker, and no history-enabled production config is supported. This sensitivity does not certify finite-Re motion, near-wall drag, dense clusters or blood transport. It does not supply an empirical correction factor.

## 12. Methods / Limitations wording

> Bubble translation was approximated by one-way flow advection plus cycle-averaged acoustic force divided by equilibrium-radius Stokes drag. An independent sensitivity analysis used identical prescribed APECSS radial histories and instantaneous acoustic forces to separate added-mass and no-slip zero-Re hydrodynamic-history contributions. Production averaging/startup effects were evaluated separately. The tested discrepancies limit quantitative interpretation of transient displacement and arrival time; no experimental calibration or universal error bound is claimed. History forces were not included in the population simulator.

Evidence: [results/stage3/history/20261002T120701030569Z/SUMMARY.json](../results/stage3/history/20261002T120701030569Z/SUMMARY.json). Figures: [HF1](../results/stage3/history/20261002T120701030569Z/figures/HF1.png), [HF2](../results/stage3/history/20261002T120701030569Z/figures/HF2.png), [HF3](../results/stage3/history/20261002T120701030569Z/figures/HF3.png), [HF4](../results/stage3/history/20261002T120701030569Z/figures/HF4.png).
