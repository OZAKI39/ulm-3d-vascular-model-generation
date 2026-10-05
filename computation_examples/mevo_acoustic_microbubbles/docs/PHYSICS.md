# Scientific scope and sources

The implemented force equations are those explicitly requested in the task. The harmonic sign convention
is positive-time complex exponent. Linearization of an inertial spherical radial oscillator gives
the negative pressure-to-volume response below resonance: increasing pressure compresses a bubble.
Multiplication and cycle averaging give the factor 1/2 in both forces. Thus a subresonant bubble
in a standing wave moves toward pressure antinodes, and identical in-phase volume oscillations attract.
Tests include global phase invariance; changing every pressure and volume phase together cannot change force.

The pair force assumes a dilute, quasistatic far-field monopole interaction with prescribed oscillations.
It omits retardation, multiple scattering and mutual modification of radial response. Its use at short
separations is an explicit numerically regularized surrogate. Spatial contact/cluster results should not
be interpreted as experimentally validated aggregation rates. Stokes mobility also omits added mass,
unsteady drag, inter-bubble hydrodynamic mobility and lubrication. These limitations become stronger in a dense cloud.

Primary sources checked:

- k-Wave Python unified API: https://k-wave-python.readthedocs.io/en/latest/get_started/new_api.html
- k-Wave array representation: https://k-wave.org/documentation/kWaveArray.php
- Harmonic extraction: https://k-wave.org/documentation/extractAmpPhase.php
- Bjerknes forces between small cavitation bubbles in a strong acoustic field (1997):
  https://doi.org/10.1103/PhysRevE.56.2924
- Marmottant et al. (2005), coated-bubble buckling and rupture model:
  https://doi.org/10.1121/1.2109427
  Author manuscript: https://liphy-annuaire.univ-grenoble-alpes.fr/pages_personnelles/philippe_marmottant/Publications_files/Marmottant2005_JASACoatedBubble.pdf
- Acoustofluidic measurements of coated-bubble forces (2018): https://pmc.ncbi.nlm.nih.gov/articles/PMC6187510/

These sources motivate the scope and deferred shell model. They do not establish the YAML resonance/damping
as SonoVue constants. No nonlinear shell equation or material constant was guessed.

`k-wave-python 0.6.2` installed source was inspected directly for constructor signatures, C-order source mapping,
1-based sensor.record_start_index conversion and backend dispatch. The Python/GPU implementation records only
the final selected samples. The incompatible C++ binary failure and its recording API limitations are preserved.

The initial parameter drafting used 25 kPa per source for migration. Before any particle production,
the solved unit fields showed that this could violate the small radial-oscillation premise. Migration
source pressure was therefore set to 2 kPa per element, with no change to field solutions or force formulae.
Aggregation retains the task's 75 kPa local-pressure scenario and reports its larger surrogate radial amplitude explicitly.

The initial 100-bubble dense-contact development probe was interrupted after 1.941884 ms of its
20 ms horizon because of measured computational cost. All 7935 complete saved states are retained,
with an explicit interrupted status. The accepted integration demonstration uses 32 bubbles for
20 ms, with unchanged cloud radius, pressure target, force equations and geometric thresholds.
It is a smaller initial population, not a continuation or relabeling of the 100-bubble probe.
The matched secondary-off control quantifies the primary-force-only behavior with the same initial cloud.
