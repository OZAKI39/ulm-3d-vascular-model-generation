# CPD independent comparison

Official pycpd 2.0.0 AffineRegistration; 60 iterations, tolerance 1e-4, w=0.1; center/scale-normalized points; initial rotations 0, -0.18 and +0.18 rad. Previous point set is source Y; current is target X. Transform converted consistently to the same six parameters. Independent implementation source convention is checked against a known affine map.

Require >=8 detected points, count change <=25%, convergence, residual <=3 px, axis ratio >=1.15 and nondegenerate design. Multi-start parameter spread above preregistered tolerances invalidates the corresponding field. Isotropy is a conservative warning, not a proof that every isotropic-looking irregular point set is ambiguous. Empirical ranges are not calibrated confidence intervals.

| Method | Endpoint tests | All-six output fraction | Accepted failures | Failure / accepted truth outputs | Non-affine correct rejection |
| --- | --- | --- | --- | --- | --- |
| PIV_multipass | 216 | 0.6991 | 21 | 0.1429 | 0.6667 |
| CPD | 216 | 0.5417 | 12 | 0.1081 | 0.5 |

CPD sometimes has lower errors on the subset it accepts, and rejects more cases. In the matched AI subset its 10 accepted cases have no tolerance failures; that does not prove it dominates PIV on all 18 inputs. On the full test it has 12 erroneous accepted truth estimates and admits 6/12 non-affine opposing-motion cases. Static candidates can produce a wrong correspondence even when point counts stay constant; high dropout often breaks convergence/residual tests. Small nearly symmetric sets invalidate the material rotation/strain interpretation. This is an independent diagnostic, not a mandatory real-time ID tracker.

Full convergence, residual, count, start sensitivity, time and null reasons: `results/benchmark_test.json`. Additional local defocus and extra candidates: `results/supplement_stress.json`. No persistent IDs are used.
