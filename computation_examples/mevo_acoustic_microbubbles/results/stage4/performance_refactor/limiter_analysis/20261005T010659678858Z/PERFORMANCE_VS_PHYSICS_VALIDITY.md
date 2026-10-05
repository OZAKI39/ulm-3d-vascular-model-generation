# Performance versus physics validity

These are four short checkpoint replays, not an instrumented wall-time census of the full gold runs.
The original evaluation-pair fractions include RK stages. Accepted-time and measured-step wall-time fractions below use the **start** of each accepted interval.
FAR/TRANSITION/NEAR/CONTACT values are pair fractions averaged with the named weights, not a partition of CPU time between force kernels. The any_NEAR_CONTACT_or_CONTACT value measures time/step coverage by at least one such pair.
Pair fractions and the fraction of steps containing any near/contact pair answer different questions; neither establishes local physical accuracy.

## B2_S1
- accepted_physical_time: {"FAR_FIELD": 0.9949524076056655, "TRANSITION": 0.0018891833027209607, "NEAR_CONTACT": 0.0021268805070899623, "CONTACT": 0.001031528584523572, "any_NEAR_CONTACT_or_CONTACT": 0.9888033577984812}
- measured_step_wall_time: {"FAR_FIELD": 0.9929136621636325, "TRANSITION": 0.0025952700110201624, "NEAR_CONTACT": 0.003177128051074193, "CONTACT": 0.0013139397742731973, "any_NEAR_CONTACT_or_CONTACT": 0.9995880771958948}
- accepted_state_pair_samples: {"FAR_FIELD": 0.9932979826546003, "TRANSITION": 0.0024678439751969162, "NEAR_CONTACT": 0.002967362158538629, "CONTACT": 0.0012668112116641529, "any_NEAR_CONTACT_or_CONTACT": 0.9994301994301994}

## B2_S2
- accepted_physical_time: {"FAR_FIELD": 0.7633823363645228, "TRANSITION": 0.15345717985717025, "NEAR_CONTACT": 0.07464567040333425, "CONTACT": 0.008514813374972752, "any_NEAR_CONTACT_or_CONTACT": 1.0}
- measured_step_wall_time: {"FAR_FIELD": 0.7633766359843683, "TRANSITION": 0.1534425923788317, "NEAR_CONTACT": 0.07466889933754534, "CONTACT": 0.008511872299254642, "any_NEAR_CONTACT_or_CONTACT": 1.0}
- accepted_state_pair_samples: {"FAR_FIELD": 0.7633780998833526, "TRANSITION": 0.153443612482337, "NEAR_CONTACT": 0.0746650655616765, "CONTACT": 0.008513222072633943, "any_NEAR_CONTACT_or_CONTACT": 1.0}

## B3_S1
- accepted_physical_time: {"FAR_FIELD": 0.9949751746358211, "TRANSITION": 0.0018712232049851278, "NEAR_CONTACT": 0.0020553688554870095, "CONTACT": 0.001098233303706757, "any_NEAR_CONTACT_or_CONTACT": 0.9893039713039545}
- measured_step_wall_time: {"FAR_FIELD": 0.9925406296194981, "TRANSITION": 0.0027238451283347334, "NEAR_CONTACT": 0.003263673056237009, "CONTACT": 0.0014718521959302226, "any_NEAR_CONTACT_or_CONTACT": 0.9996374387211616}
- accepted_state_pair_samples: {"FAR_FIELD": 0.9929447664596474, "TRANSITION": 0.002588462559531471, "NEAR_CONTACT": 0.0030498985015231473, "CONTACT": 0.001416872479298065, "any_NEAR_CONTACT_or_CONTACT": 0.99945295404814}

## B3_S2
- accepted_physical_time: {"FAR_FIELD": 0.6653938618483559, "TRANSITION": 0.23907102236801378, "NEAR_CONTACT": 0.085443100325494, "CONTACT": 0.010092015458136204, "any_NEAR_CONTACT_or_CONTACT": 1.0}
- measured_step_wall_time: {"FAR_FIELD": 0.6653826582670046, "TRANSITION": 0.23906655159102266, "NEAR_CONTACT": 0.0854610808181952, "CONTACT": 0.010089709323777472, "any_NEAR_CONTACT_or_CONTACT": 1.0}
- accepted_state_pair_samples: {"FAR_FIELD": 0.6653837844406234, "TRANSITION": 0.23906760611269576, "NEAR_CONTACT": 0.08545821626442719, "CONTACT": 0.010090393182253625, "any_NEAR_CONTACT_or_CONTACT": 1.0}

Large cost in NEAR_CONTACT/CONTACT reflects the numerical stiffness of advancing an approximate far-field model in close configurations. Smaller timesteps do not by themselves make that force model physically valid. Finite-Re and near-wall limitations remain.
