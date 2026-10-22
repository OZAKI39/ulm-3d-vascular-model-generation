# SonoVue parameters and limits

Primary source: Tu, Guan, Qiu & Matula, JASA 126, 2954–2962 (2009), DOI [10.1121/1.3242346](https://doi.org/10.1121/1.3242346), [full text PMC2803720](https://pmc.ncbi.nlm.nih.gov/articles/PMC2803720/). Full HTML and figure 4 were downloaded and read; numerical values below are transcribed from the article, not inferred from its abstract.

| Case | Agent | Radius (µm) | chi (N/m) | kappa_s (kg/s) | Pressure / frequency | Provenance |
|---|---|---:|---:|---:|---|---|
| `br14_table1` | **BR14** | 1.70 | 0.25 | 4.0e-9 | 40 kPa / 2.5 MHz | Table 1, Fig. 1, Gaussian tapered 8 cycles, data originally Van der Meer 2007 |
| `sonovue_table2` | SonoVue | 1.78 | 0.30 | 3.2e-9 | 150 kPa / 2.5 MHz | Table 2, Fig. 3, measured Ultramark 4 Plus M-mode pulse |
| `sonovue_table3` | SonoVue | 1.72 | 0.22 | 2.5e-9 | 150 kPa / 2.5 MHz | Table 3, radius included as fitting parameter; alternative fit to the same observation |

The fitted equation is their Eq. (6), a **linearized shell constitutive law without buckling or rupture**, retaining nonlinear radial inertia/gas terms. The paper reports increasing elasticity/viscosity with radius, not a universal SonoVue constant. The Table 2 and 3 differences quantify a calibration/model choice sensitivity; they are not two independent population measurements or confidence limits. The 40 kPa example must not be renamed SonoVue.

`configs/stage2/sonovue_parameters.json` contains SI values with DOI, table/figure, conditions, and extraction method for every parameter. `SonoVueParameterSet` validates SI units and positive values. `SonoVueParameterDistribution` interpolates selected SonoVue markers of Fig. 4 (elasticity linearly; viscosity linearly in log10). The 400×672 source image, SHA256, selected pixel centres and axis calibration are retained. Hollow-square BR14 viscosity markers are excluded. This is an approximate digitization, **not the authors' raw data**. Radius bounds are the intersection of the two sets of digitized markers; extrapolation is rejected. It does not justify a 0.75–4 µm parameter range.

Digitization uncertainty is explicitly an implementation resolution envelope: ±0.02 N/m and ±0.06 log10(kg/s). It is not the paper's standard deviation bars and does not represent biological variability. Seeded draws within this envelope support sensitivity tests; no normal/lognormal population law has been inferred. Figure 4's reported standard deviations and radius calibration uncertainties would require better raw data for a population-grade statistical model. Radius sampling remains an independent input.

The source's simulation constants include p0=101300 Pa, rho=1000 kg/m³, c=1500 m/s, sigma_water=0.072 N/m, gamma=1.07, and effective liquid viscosity=0.002 Pa s (includes thermal damping). Keeping that effective viscosity and separately adding thermal damping would double count it. Vapor pressure 2330 Pa is listed for the general model but not subtracted in their Eq. (6). APECSS ideal gas equilibrium is mapped explicitly to the selected equation.

## Mapping to the nonlinear backend

The paper does **not** identify a unique equilibrium surface tension or buckling/rupture radii for the full piecewise shell. Eq. (6)'s water-tension term cannot simply be inserted as sigma_init into the clipped law while preserving its elastic derivative. Our initial nonlinear study declares sigma_init=0.036 N/m as a **modeling hypothesis**, with 0 and 0.072 N/m sensitivity cases; it is not a measured SonoVue value. APECSS's reversible rupture-at-water-tension law is separately audited in `APECSS_MODEL_AUDIT.md`. Exported responses remain partially validated/software-benchmarked, rather than experimentally calibrated SonoVue predictions.

Continuous sinusoidal steady-state tables differ from the paper's finite measured pulses. Neither a fit to Figure 3 nor experimental waveform reproduction is claimed when the original excitation waveform samples are unavailable. Frequencies/pressures away from 2.5 MHz/150 kPa are a constitutive-model extrapolation, even when strictly inside the numerical lookup table. Both distinctions must travel with table metadata.

The trajectory adapter accepts `medium.radial_effective_viscosity_Pa_s` to match the radial table while retaining `medium.viscosity_Pa_s` for physical translation drag. The Stage2 template declares .002 and .001 Pa s respectively. Radial-only damping changes do not invalidate the acoustic cache.
