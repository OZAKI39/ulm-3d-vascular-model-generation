# Stage 4I.5 backend decision

recommended_backend_for_JK: **P1_ordered**

Correctness takes priority. Candidates need all four completed S1/S2 regressions and the strict trajectory/contact gate. Aggregate speedup uses summed matched solver times. A component-timed gain greater than 5% is a screening criterion; adoption additionally requires the same four windows with per-component wrappers removed, strict timer-removal/state regression and a confirmed gain greater than 5%. This controls the API-dependent cost of millions of per-point VTK timer calls. The 2× engineering target is not a scientific acceptance threshold.

| Backend | Completed segments | Strict gate | Aggregate speedup |
|---|---:|---|---:|
| P1_exact_cache | 4 | True | 0.9689189782766292 |
| P1_neighbor | 4 | True | 0.9380730252192133 |
| P1_combined | 4 | True | 0.9145765313608052 |
| P1_static | 4 | True | 0.9830988829429476 |
| P2_warp | 2 | False | Not comparable |
| P1_ordered | 4 | True | 1.5170516370178377 |

## Minimal-timer confirmation

Evidence: results/stage4/performance_refactor/timing_confirmation_report/20261005T025613151336Z/SUMMARY.json. Common accepted-step diagnostic recording, dt and classifier timers remain; per-component and per-point wall wrappers are disabled. These short replay timings are not full-production runtime measurements.

| Backend | Strict confirmation | Aggregate speedup with component wrappers removed |
|---|---|---:|
| P1_ordered | True | 1.5833487944633318 |

## Independent wall fixture

| Backend | Geometry gate | Query speedup | Max distance error [m] | Max normal angle [deg] | Inside mismatches | Projection decisions | Velocity decisions |
|---|---|---:|---:|---:|---:|---:|---:|
| vtk | True | 1.4704 | 0 | 0 | 0 | 0 | 0 |
| warp | False | 106.43 | 1.3509e-09 | 92.036 | 2 | 1019 | 1017 |

vtk full projection: `{"max_position_error_m": 0.0, "rms_position_error_m": 0.0, "hit_mismatches": 0, "maximum_penetration_m": 9.99999999392329e-13, "projection_failures": 0}`. Evidence: results/stage4/performance_refactor/p1_vtk_contact/20261005T010705538792Z/SUMMARY.json.


warp full projection: `{"projection_failures": 1, "error": "RuntimeError('Y lumen projection did not converge')"}`. Evidence: results/stage4/performance_refactor/p2_warp/20261005T023057899660Z/SUMMARY.json.


## Original physics and terminal boundaries

All 42 original production source files match sealed gold byte-for-byte. New optional modules do not change force laws, constraints or physical tolerances.
Existing assigned-particle removal is an absorbing model boundary, not proven terminal-outlet physics. No additional deletion or endpoint shortcut is enabled. J retains the 6 s observation with objective-policy validation; a plateau with unassigned particles is insufficient to stop.

## Local subcycling

No experimental local solver or full optimized P3 runs launched. Global all-pair forces and widespread stiffness make a small active subset unsupported by these samples.
[
  {
    "id": "B2_S1",
    "ids_covering_90_percent": 67,
    "median_stiff": 127.0,
    "active": 256.0,
    "median_contact_component": 1.0
  },
  {
    "id": "B2_S2",
    "ids_covering_90_percent": 65,
    "median_stiff": 171.0,
    "active": 180.0,
    "median_contact_component": 55.0
  },
  {
    "id": "B3_S1",
    "ids_covering_90_percent": 51,
    "median_stiff": 109.0,
    "active": 256.0,
    "median_contact_component": 1.0
  },
  {
    "id": "B3_S2",
    "ids_covering_90_percent": 78,
    "median_stiff": 155.0,
    "active": 168.0,
    "median_contact_component": 75.0
  }
]

## Scope of gold regression

The original B2/B3 are the unique full gold runs and were never rerun. P0 checkpoints are independently reconstructed and compared with gold endpoints. P1 tests compare every accepted state, dt, contact hash, wall flags and assignments with matched P0 short replays. No measured full-6-s optimized branch agreement is claimed when no such run was performed. Most dense windows contain no new outlet crossing; their strongest tests are state/contact equivalence.

## Warp and SDF

Warp remains optional/benchmark-only unless its explicit mesh, projection, trajectory and runtime gates all pass. Float32 geometry is never made acceptable by loosening the original tolerances. SDF is not substituted: the new profile attributes major dense-case cost to contact projection/NNLS, so another approximate wall representation does not address that dominant cost.
