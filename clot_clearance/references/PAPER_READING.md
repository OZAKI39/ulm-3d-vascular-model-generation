# Supplied paper and implementation boundary

Karmakar, Burgreen & Antaki (2026), *Peridynamics coupled computational fluid dynamics as a paradigm for simulating problems involving fracture of thrombi under fluid flow: Demonstration in aspiration thrombectomy*, Applications in Engineering Science 27, 100347. DOI: https://doi.org/10.1016/j.apples.2026.100347. The supplied 12-page PDF is copied alongside this note, with extracted text.

The paper concerns **aspiration thrombectomy**, not ultrasound or microbubble fatigue. It uses a two-way PD/OpenFOAM formulation and presents 2-D bifurcation demonstrations. Its equations do not validate the empirical cycle-jump law requested in the separate prompt.

| Paper location | This prototype |
|---|---|
| p.2, Eqs.4–8: shape tensor, nonlocal F, correspondence force state | Implemented in 3-D with pairwise force/energy verification. |
| p.3, Eqs.9–11: Silling stabilization | Use a separately documented objective, full-vector non-affine energy penalty, **not an exact implementation of Eq.10**. Configurable and explicitly provisional. |
| p.3, Eqs.12–15: irreversible bond integrity and particle stress degradation | Implement irreversible integrity; average incident integrity defines particle damage; choose alpha(D)=D. Smooth cubic instantaneous transition rather than the paper's linear transition. |
| p.3, Eqs.16–19: damped explicit dynamics | Use damped velocity Verlet with exact exponential damping split, not the paper's particular gamma update. |
| p.4, Eq.26 and Eq.28: fluid stress and surface force | Implement sigma=-pI+mu(grad u+grad u.T), outward **solid** normal, and reference surface quadrature t*A. |
| pp.3–4: two-way conservative pMesh coupling | Deferred. FileStreaming is one-way only; direct traction, linear tetra interpolation and bounded scattered interpolation are supported. |
| p.5, Eq.33: generalized Ogden; Eq.34 material constants | Deferred. Start with compressible Neo-Hookean and explicit synthetic material values; do not relabel these as fitted paper parameters. |
| p.5: s1=2.0, s2=2.5; p.9 limitations | The paper itself states its fracture thresholds were not calibrated to clot composition. Our software tests use separately disclosed artificial thresholds. |
| Cycle jump / MHz microbubble streaming | Not in this paper. Added solely as the requested unvalidated prototype interface. |

The original PDF remains untouched. The preserved copy is for this user's local research project; no third-party source code is copied from the paper. Do not describe the straight-pipe demo as reproducing the paper's aspiration outcomes.
