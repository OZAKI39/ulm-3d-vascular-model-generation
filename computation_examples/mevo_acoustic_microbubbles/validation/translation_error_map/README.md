# Translation approximation error map

Source Doinikov–Dayton2006, DOI10.1121/1.2215228, Eqs6,31–33,36,40.
Use existing four-state reference and identical reference-cycle mean forces in
quasi-steady6pi mu R0 balance.3radii×3pressures×3frequencies,20cycles from rest;
9exact Stage2 cases are reused,18new cases are solved. No postpulse drift is added.

The error is |x_fast(T)-x_reference(T)|/max(|x_reference(T)|,1e-14m).
The raw signed displacements and absolute errors are retained. In weak short
pulses, the endpoint includes oscillatory/startup motion while radiation drift
scales as pressure². A small reference displacement can therefore produce a very
large relative error, even at low Re. This is a result, not a justification to
change the force or hide the negative displacement. This finite-pulse map does
not certify20ms CW or SonoVue translation.

FAST_MODEL_SAFE≤5%, FAST_MODEL_APPROXIMATE≤20%, otherwise REFERENCE_REQUIRED.
These are project error budgets. TranslationValidityMap rejects unsupported
radius/pressure/frequency/shell/pulse length; it does not infer a continuous
accuracy surface from sparse nodes. Existing analytic Eq63 test remains intact.
