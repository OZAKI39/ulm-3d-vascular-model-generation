# Doinikov–Bouakaz 2015 N=2 reference

[PRE 92,043001](https://doi.org/10.1103/PhysRevE.92.043001); author-uploaded full text was read on ResearchGate.
The authors explicitly provide Mathematica code for inverse-distance expansion coefficients in the supplement.
APS PDF/supplement requests returned HTTP403 and the exposed supplemental link required subscription;
HAL DOI lookup returned zero results. This is **inaccessible author code**, not an assertion that no code exists.

`validation/doinikov_bouakaz_pair.py` implements the finite Legendre boundary system Eqs8–13,
evaluates kinetic energy Eq41, then applies Lagrange Eqs39–43, pressure Eq67 and drag Eq68.
It does not copy the appendix or claim numerical identity with the authors' d^-9 expansion.
Complex-step derivatives of the symmetric mass matrix avoid cancellation. Exchange/reflection symmetry,
isolated-bubble added mass, positive mass eigenvalues, and multipole order convergence are tested.
M=8 versus M=12 maximum distance difference before contact is 2.37e-10 m in the strong-pulse example.
The two weak-field motion directions agree with Fig3, using a shortened 60-cycle pulse.

Limitations: uncoated gas, spherical boundaries, potential-flow interaction, spatially uniform pressure,
approximate radial acoustic damping, thin-boundary-layer drag; no lubrication/contact/adhesion/coalescence.
Gaussian full duration is assumed to be 6 sigma because the paper does not specify its width convention.
Fig2 uses its stated 30 cycles; quantitative agreement with the original waveform is not established.
Contact terminates the ODE. Multipole convergence does not validate overlapping viscous layers or nonspherical deformation.
Only N=2 is supported by `near_field_policy=pair_reference`; vascular flow and PZT fields are explicitly refused.
