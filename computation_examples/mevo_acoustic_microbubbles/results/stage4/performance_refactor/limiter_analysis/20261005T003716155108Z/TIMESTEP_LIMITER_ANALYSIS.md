# Timestep limiter analysis

Global particle IDs are retained. Macro-boundary remainders are counted separately from displacement limits.
Categories overlap; a contact/cluster limiter can also be near-wall and outside the far-field approximation.

| Segment | Steps | Particle-limited | IDs covering 90% | Median stiff/active | Median contact component |
|---|---:|---:|---:|---:|---:|
| B2_S2 | 3162 | 3122 | 65 | 171.0/180.0 | 55.0 |

Local-subcycling eligibility is evidence-based: few limiting IDs alone are insufficient if the contact/near-contact buffer covers much of the population. Long-range active–inactive forces must still be evaluated; freezing the remaining bubbles is not allowed.
