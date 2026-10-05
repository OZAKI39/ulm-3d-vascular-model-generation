# Matched pair asymptotics

The comparison sets chi, shell viscosity and shell thickness to zero in Zhang;
Doinikov is uncoated throughout. Both use R1=1.5um,R2=2um,rho998,mu.001,
sigma.072,p0101300,gamma1.4 and the same continuous -Pa sin(wt) excitation.
Only the Gaussian envelope of the Doinikov pressure input is replaced; mass,
radial damping and translation equations are unchanged. This is an explicit
benchmark boundary condition, not reproduction of its published Gaussian pulse.

Two distinct forces are retained:

* rho mean(Vdot1 Vdot2)/(4pi d0²), evaluated on coupled solutions;
* the antisymmetric interaction reconstructed from each bubble's isolated
  added-mass momentum rate minus its explicit drag, using the full ODE acceleration.

The latter is not merely Eq11 reused on both outputs. Terms involving the other
bubble's velocity in Zhang drag remain in its specified drag; Doinikov Eq68 has
a different drag closure. Finite-c radial damping formulations also differ. Those
model differences are retained rather than fitted away.

Classical reference uses analytic isolated small-signal uncoated volume phasors
with the Doinikov pressure damping. Formal c=1e9 m/s removes finite-compressibility
model mismatch and isolates the hydrodynamic far-distance limit. It is not water.
At physical c1482, increasing d eventually violates the quasistatic kd policy;
such samples are plotted as formal comparisons and excluded from the regime map.

Errors at sign-zero crossings need care: the relative denominator is the actual
classical magnitude with a1e-40 safety floor; a changed sign always rejects safe
classification. Fits report observed slopes, not a imposed d^-n claim.
The independent momentum residual can suffer cancellation for extremely weak,
very distant pairs; radial coefficient and momentum residual are reported separately.

PairPhysicsRegimeMap accepts only exact sampled radii, pressure, frequency,
separation, medium sound speed and uncoated shell. It checks BOTH models, BOTH
force error measures, sign, kd, radial excursion and center drift.5%/20% are
project accuracy budgets; the locations meeting them come from data. Untested
intervals are rejected. A valid uncoated point never certifies SonoVue or Lazarus.
