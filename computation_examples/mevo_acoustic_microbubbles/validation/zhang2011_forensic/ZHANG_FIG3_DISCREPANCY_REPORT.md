# Zhang Figure 3 discrepancy — UNRESOLVED, localized evidence

**The nonlinear Fig3 reproduction still fails. No force sign or source parameter
was changed to make it pass.** The printed equations have now been checked by two
separate numerical implementations. This closes the software-sign investigation
substantially, but does not resolve the missing information behind the published curve.

## 1. Errors excluded within the tested scope

The independent SI residual implementation shares no production bubble classes,
force backend, APECSS code or acceleration-matrix assembly. It evaluates all four
second-order residual equations before extracting an acceleration system.
Seven matrix cases agree with the Stage2 implementation: maximum relative radius
difference below6e-12 and center differences below1.5e-18m under the same high
accuracy settings. Random non-equilibrium residual checks supplement trajectory
comparisons. SI dimensional tests cover every family of radial, translation,
shell and coefficient terms. The total derivative of Eq9 includes both radius
and radial-velocity derivatives. The1/D² factor in Eqs3/4 is outside the derivative.

The paper explicitly says positive fB means attraction. Its drive is minus sine;
the exp(+iwt) convention and harmonic phase conversion are consistent. Periodic
time-origin and reflected-coordinate tests preserve physical force. Eq11 gives
N m²; the Fig3 axis is fB/(1e-18 N m²), not force in N or a dimensionless parameter.
There is no overall sign mistake, radius/diameter switch, surface-gap substitution,
or an omitted common sin-to-cos phase factor supported by these checks.

## 2. Where the disagreement occurs

At R10=3um,R20=6um,D0=120um,Pa100kPa,f2.25MHz,chi1N/m,
the cycle-average Vdot product is negative after transients. The coefficient is
−6.143989e-18 N m² over cycles56–64 and−6.131573e-18 over248–256.
Approximate image extraction gives about+6.6e-18 (image uncertainty about0.22e-18).
The discrepancy is already in the radial correlation, before vector construction
or division by distance². The full HDF5 records the pressure, shell, scattered,
coupling, radii, velocities, volume derivatives, product and force averages.

The first cycle yields+5.44574e-18, but cycles0–8 give−2.49594e-18,
8–16 give−5.28161e-18, and24–32 give−6.02012e-18. Selecting the initial
positive cycle would be an unjustified change of averaging protocol: the50um
case's first-cycle coefficient is also positive, contrary to Fig3's negative value.
No selected time window is presented as a recovered publication protocol.

## 3. Unexcluded possibilities and paper ambiguity

Section3 specifies chi≤1N/m without a figure-specific value. Figure3 does not
state its averaging window, transient removal, integration horizon or solver.
The model uses the same boundary radius in the gas law and shell terms; a separate
coated/gas radius conversion is not defined. These ambiguities remain unresolved.
A transcription omission in the publication, use of an unstated approximation,
or a different numerical procedure remain possible; none is established as fact.

Eq10 is an approximate resonance expression that omits a surface-tension stiffness
term present in the derivative of Eq9. The rounded2.9um resonance cannot uniquely
recover chi. The separate Fig4 prose has apparent time-unit typos, but this does
not establish a typo in the Fig3 force or justify correcting its sign.

## 4. Small-signal diagnostic, not a fitted correction

Linearizing the independent printed-equation residual gives+5.248832e-18 at120um.
An actual100Pa nonlinear calculation scaled by(100000/100)² gives+5.248786e-18.
Their agreement supports the linearization calculation. The small-signal curve
crosses zero near70um, whereas the image crosses near60um; its30um value is
about−14.88e-18 instead of the image's roughly−9.3e-18. Thus it has a similar
sign-reversal trend but does **not** reproduce the published quantitative curve.
The100Pa run is clearly identified as a diagnostic and is never renamed100kPa.
This evidence suggests investigating an unstated linear/weak-amplitude treatment;
it does not prove that the authors used one.

## 5. Other literature and independent reproduction

Doinikov–Bouakaz2015 supports separate coupling/translation and far-limit tests
for uncoated bubbles. Agreement there is not an independent reproduction of this
coated100kPa Fig3 case. Searches by exact title/DOI/Figure3/Mathematica, APS,
author/institutional sources, HAL and public archive lookup did not retrieve an
independent exact Fig3 solver or original numerical values. A2018 finite-element
study appeared in search, but the linked full page returned404, so its search
excerpt is not treated as a validated reproduction. We make no claim that no
independent reproduction exists.

## 6. Author code retrieval

The Doinikov2015 paper explicitly refers to a Mathematica supplement. This turn's
APS supplemental endpoint returned401; the DOI supplemental link returned the
abstract page whose supplemental section is marked data-auth-required. HAL DOI
lookup returned zero and the public archive query returned an empty list. An
author-institution candidate timed out. Raw responses/statuses are saved under
results/stage2_5/audit/literature_retrieval. This is inaccessible/unretrieved code,
not nonexistent code. No access restriction was bypassed and no author contacted.

## 7. Current interpretation and gate decision

The most defensible interpretation is that both implementations solve the printed
nonlinear equations consistently for the explicitly declared chi1N/m assumption,
while the paper's Fig3 cannot yet be reproduced from its stated information.
Zhang is **partially validated; Fig3 UNRESOLVED**. Low-amplitude uncoated asymptotic
comparisons do not prove that coated50–100kPa clustering avoids this uncertainty.
G1 therefore remains closed for the proposed production aggregation domains.
Lazarus, SonoVue clustering, cylinder and frequency sweep remain NOT_RUN_GATE.

Sources: [Zhang publisher article](https://cpb.iphy.ac.cn/en/article/doi/10.1088/1674-1056/20/11/114302),
[Doinikov2015](https://doi.org/10.1103/PhysRevE.92.043001).
