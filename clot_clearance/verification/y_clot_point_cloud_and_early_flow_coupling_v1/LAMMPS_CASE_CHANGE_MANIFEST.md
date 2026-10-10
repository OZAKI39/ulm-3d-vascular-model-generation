# One permitted actual mechanics case

Old C1_GAUSSIAN_A200 is read-only. New C1_GAUSSIAN_A200 resides only under this task. No baseline rerun, calibration, amplitude scan, grid study, CFD or material upgrade.

Unchanged: body.data,1927 clot nodes,1561 anchors,157637 initial bonds,3 anchor layers,V3/native Y,K33333.333333333336Pa,G940Pa,rho1000kg/m3,bulk s00=.15,interface=.140624859375,alpha0,dx50um,delta3.01dx,dt.5us,4660 steps,Gaussian reference200Pa,pulse1.1605ms,1.0dx post-break contact,native-wall callback and original CFD. Input file hashes are recorded.

Necessary changes: remove unconditional fix viscous; add earlyflow fix external after nativewall with event-triggered COM fluid force and attached-only original numerical damping. Add separate d_drag* / d_num* diagnostics. Increase output float precision15->17digits for faithful terminal handoff. Two MPI ranks and existing shared library retained. Interface run0 advances zero steps and is not a new mechanics case. MECHANICAL_CASE_LAUNCHED.json enforces at most one actual run. No changes to the existing install.
