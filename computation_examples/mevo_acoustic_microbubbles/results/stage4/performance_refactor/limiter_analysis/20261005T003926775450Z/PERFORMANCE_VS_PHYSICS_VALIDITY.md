# Performance versus physics validity

These are four short checkpoint replays, not an instrumented wall-time census of the full gold runs.
The original evaluation-pair fractions include RK stages. Accepted-time and measured-step wall-time fractions below use the **start** of each accepted interval.
Pair fractions and the fraction of steps containing any near/contact pair answer different questions; neither establishes local physical accuracy.

## B2_S2
- accepted_physical_time: {"FAR_FIELD": 0.7633823363645228, "TRANSITION": 0.15345717985717025, "NEAR_CONTACT": 0.07464567040333425, "CONTACT": 0.008514813374972752, "any_NEAR_CONTACT_or_CONTACT": 1.0}
- measured_step_wall_time: {"FAR_FIELD": 0.7633766359843683, "TRANSITION": 0.1534425923788317, "NEAR_CONTACT": 0.07466889933754534, "CONTACT": 0.008511872299254642, "any_NEAR_CONTACT_or_CONTACT": 1.0}
- accepted_state_pair_samples: {"FAR_FIELD": 0.7633780998833526, "TRANSITION": 0.153443612482337, "NEAR_CONTACT": 0.0746650655616765, "CONTACT": 0.008513222072633943, "any_NEAR_CONTACT_or_CONTACT": 1.0}

Large cost in NEAR_CONTACT/CONTACT reflects the numerical stiffness of advancing an approximate far-field model in close configurations. Smaller timesteps do not by themselves make that force model physically valid. Finite-Re and near-wall limitations remain.
