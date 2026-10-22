# Phase and sign conventions

All project acoustic phasors mean Re[P exp(+i omega t)]. Zhang2011 after Eq8
and APECSS use -Pa sin(omega t), so their input phasor is +i Pa.
A volume Fourier coefficient measured against exp(-i h omega t) must be
multiplied by i^(-h) to tabulate the response to positive cosine pressure.
For a local pressure phase theta the hth volume coefficient gains exp(i h theta).

Time-origin shifts multiply both matching phasors by the same harmonic phase.
Consequently Re[V1h conj(V2h)] and Re[Vh conj(gradPh)] are unchanged.
Changing sin to cos while keeping a suddenly started, resting bubble's initial
state is a DIFFERENT transient. The invariance test shifts complete periodic
signals (or shifts their initial state as well), not unrelated startup windows.

Define e12=(x2-x1)/d. Eq11's fB=rho mean(Vdot1 Vdot2)/(4pi) is positive for
attraction according to the text directly following Eq11. Thus F_on1=+fB e12/d²,
F_on2=-F_on1 at classical order. Defining e21=-e12 changes the written component
coefficient to -fB; it does not change the physical vector or attraction/repulsion.
A global coordinate reflection reflects both vectors and positions.

For real periodic volumes, mean(Vdot1 Vdot2)=omega²/2 sum_h h² Re[V1h conj(V2h)].
Its unit is m6/s2; multiplying rho/(4pi) gives N m2; dividing d² gives N.
Zhang Fig3's 1e-18 N m2 is an axis unit, not a dimensionless force definition.
Primary force remains -0.5 Re[Vhat conj(gradP)] and is not rewritten.

The forensic artifact keeps the product, pressure/coupling terms, instantaneous
force and multiple average windows separately. No negative prefactor or arbitrary
phase adjustment is introduced to repair a literature discrepancy.
