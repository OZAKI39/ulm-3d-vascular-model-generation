# CURVED_MURAL_MULTILOBE_CLOT_V3
Single deterministic engineering morphology. Not segmented patient anatomy.
Coordinates: a=x-x0, b=Rv atan2(y,z), u=a/0.68mm, v=b/0.62mm.
U=u-.10u^2+.05; V=[v-.10(1-u^2)-.055u]/[1+.10sin(pi u)]; q^2=U^2+V^2.
Height raw=max(1-q^2,0)^1.15 * [0.44mm*(1-.09 tanh(2u)) + 0.43mm*sum(w_j exp(-((a-a_j)^2+(b-b_j)^2)/(2sigma_j^2)))].
Three actual lobe centers/weights/sigmas are in CONFIG; sigma130-140um, FWHM306-330um (6+ dx). A broad base preserves connectivity. Height and derivative vanish at the footprint edge. Thicker upstream platform and unequal lobe weights form an asymmetric shoulder and downstream taper.
One predetermined scalar height normalization sets max height0.65mm. No random draws. Native WALL triangles are clipped/subdivided; height goes inward from actual wall points. Fresh global50um cubic mask, 27 support points and generic degree<2 whisker removal; new typed bonds.
The analytic definition is frozen before assembly. At most one explicit pre-dynamics morphology correction is allowed if the normal-scale preview has no visible lobes; never change after observing fracture.
