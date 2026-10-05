# Stage 4I.5 baseline profile

P0 includes the per-particle exact interior-wall certificate already present in gold. These costs must not be compared to the older unaccelerated Stage 4 profile as if it were current gold.
Initialization, scientific stepping and HDF output are reported separately. No VTK frame output occurs in the original scientific loop. Python overhead includes unwrapped work and instrumentation overhead, not an independently exact physical component.
VTK already builds its original locator once and reuses generic cells. Value and gradient calls each perform a closest-point search.

## B2_S1

Accepted steps: 1755; physical interval 0.0–0.04 s. Instrumented wall time: 85.606 s.
Gold endpoint maximum position difference: 0 m. This reconstruction check is separate from candidate-vs-P0 agreement.

| Exclusive component | Seconds | Calls |
|---|---:|---:|
| wall_closest_distance | 16.686493 | 898560 |
| NNLS | 12.806526 | 3509 |
| wall_normal_including_repeated_closest_point | 12.594313 | 898560 |
| physics_validity | 10.754263 | 3510 |
| contact_position_projection | 6.021219 | 3510 |
| contact_candidate_search | 5.285207 | 63022 |
| wall_cache_and_batch | 5.030737 | 10530 |
| branch_classification | 4.071334 | 1755 |
| flow_interpolation | 3.085124 | 3510 |
| pair_force | 2.236917 | 3510 |
| wall_projection | 2.098355 | 59512 |
| contact_velocity_assembly | 1.392031 | 3510 |
| response_lookup | 1.384171 | 3510 |
| acoustic_interpolation | 1.052893 | 7020 |
| wall_velocity_constraint | 0.274464 | 7020 |
| diagnostic_recording | 0.187484 | 3510 |
| timestep_computation | 0.037326 | 1755 |
| Unclassified Python/timer overhead | 0.606910 | — |
| HDF output (outside loop) | 0.686482 | 1 |

## B2_S2

Accepted steps: 3162; physical interval 5.94–5.96 s. Instrumented wall time: 566.323 s.
Gold endpoint maximum position difference: 0 m. This reconstruction check is separate from candidate-vs-P0 agreement.

| Exclusive component | Seconds | Calls |
|---|---:|---:|
| contact_position_projection | 238.099265 | 6324 |
| NNLS | 104.537611 | 6324 |
| wall_cache_and_batch | 59.420954 | 678489 |
| contact_candidate_search | 40.808360 | 446826 |
| wall_projection | 37.665762 | 440502 |
| wall_closest_distance | 24.503384 | 8283681 |
| wall_normal_including_repeated_closest_point | 21.351906 | 8283681 |
| physics_validity | 11.104570 | 6324 |
| pair_force | 9.588349 | 6324 |
| contact_velocity_assembly | 5.613646 | 6324 |
| branch_classification | 5.190649 | 3162 |
| flow_interpolation | 2.237879 | 6324 |
| response_lookup | 2.011120 | 6324 |
| acoustic_interpolation | 1.812376 | 12648 |
| diagnostic_recording | 0.630098 | 6324 |
| wall_velocity_constraint | 0.498874 | 12648 |
| timestep_computation | 0.068821 | 3162 |
| Unclassified Python/timer overhead | 1.179283 | — |
| HDF output (outside loop) | 0.895243 | 1 |

## B3_S1

Accepted steps: 1828; physical interval 0.0–0.04 s. Instrumented wall time: 92.589 s.
Gold endpoint maximum position difference: 0 m. This reconstruction check is separate from candidate-vs-P0 agreement.

| Exclusive component | Seconds | Calls |
|---|---:|---:|
| wall_closest_distance | 17.377418 | 935936 |
| NNLS | 15.102391 | 3649 |
| wall_normal_including_repeated_closest_point | 13.090080 | 935936 |
| physics_validity | 11.449491 | 3656 |
| contact_position_projection | 6.596756 | 3656 |
| contact_candidate_search | 5.733306 | 65805 |
| wall_cache_and_batch | 5.423200 | 10968 |
| branch_classification | 4.304184 | 1828 |
| flow_interpolation | 3.456608 | 3656 |
| pair_force | 2.393314 | 3656 |
| wall_projection | 2.270681 | 62149 |
| contact_velocity_assembly | 1.571992 | 3656 |
| response_lookup | 1.487596 | 3656 |
| acoustic_interpolation | 1.131261 | 7312 |
| wall_velocity_constraint | 0.294659 | 7312 |
| diagnostic_recording | 0.203210 | 3656 |
| timestep_computation | 0.041026 | 1828 |
| Unclassified Python/timer overhead | 0.662275 | — |
| HDF output (outside loop) | 0.712289 | 1 |

## B3_S2

Accepted steps: 3487; physical interval 4.88–4.9 s. Instrumented wall time: 680.611 s.
Gold endpoint maximum position difference: 0 m. This reconstruction check is separate from candidate-vs-P0 agreement.

| Exclusive component | Seconds | Calls |
|---|---:|---:|
| contact_position_projection | 317.208722 | 6974 |
| NNLS | 108.961732 | 6974 |
| wall_cache_and_batch | 67.820177 | 822825 |
| contact_candidate_search | 49.671821 | 558238 |
| wall_projection | 46.466829 | 551264 |
| wall_closest_distance | 26.297389 | 9216041 |
| wall_normal_including_repeated_closest_point | 22.779846 | 9216041 |
| physics_validity | 11.019083 | 6974 |
| pair_force | 9.282433 | 6974 |
| contact_velocity_assembly | 6.283629 | 6974 |
| branch_classification | 5.325759 | 3487 |
| flow_interpolation | 2.692458 | 6974 |
| response_lookup | 2.147677 | 6974 |
| acoustic_interpolation | 1.908242 | 13948 |
| diagnostic_recording | 0.752837 | 6974 |
| wall_velocity_constraint | 0.553724 | 13948 |
| timestep_computation | 0.077078 | 3487 |
| Unclassified Python/timer overhead | 1.361442 | — |
| HDF output (outside loop) | 0.957384 | 1 |

