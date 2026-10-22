# Zhang Figure 3 parameter ledger

Source: [publisher PDF](https://cpb.iphy.ac.cn/en/article/pdf/preview/10.1088/1674-1056/20/11/114302.pdf), DOI10.1088/1674-1056/20/11/114302.
The unaltered downloaded PDF and rendered source are preserved locally under third_party/downloads.

| Quantity | Value / definition | Evidence type and location |
|---|---|---|
| Figure / curve | Figure3, one curve | Figure caption, p4 |
| R10 / R20 | 3e-6 / 6e-6 m | Figure3 caption |
| Radius meaning | Rj is spherical bubble boundary in liquid; same R enters gas law and shell term | Eqs7–9; surrounding theory. Distinct inner gas vs outer shell radius is not specified; no subtraction of epsilon from R0 in gas law is authorized |
| D0 | 30–120 um; failed endpoint120 um | p4 surrounding paragraph |
| Distance | D=x2−x1>0, center-to-center | Fig1 and initial conditions p3; not surface gap |
| Frequency | 2.25e6 Hz | Fig3 caption |
| Pressure amplitude | 1e5 Pa (peak, not RMS) | Fig3 caption; drive definition after Eq8 |
| Ambient pressure | 1.01e5 Pa | Section3 first paragraph |
| Density | 998 kg/m3 | Section3 first paragraph |
| Liquid viscosity | .001 Pa s | Section3 first paragraph |
| Surface tension | .051 N/m | Section3 first paragraph |
| Shell elasticity | chi ≤1 N/m; exact Fig3 value absent | Section3 first paragraph |
| Implemented chi | 1 N/m | DECLARED ASSUMPTION: endpoint of stated range; consistent with rounded2.9um resonance, not directly specified |
| Shell viscosity | 1.27 Pa s | Section3 first paragraph; not surface viscosity kg/s |
| Shell thickness | 2e-9 m | Section3 first paragraph |
| Gas exponent | 1.07 | Section3 first paragraph |
| Sound speed | 1500 m/s | Section3 first paragraph |
| Initial radii / velocities | equilibrium radii, all four derivatives zero | Section3 first paragraph |
| Initial centers | x1=0, x2=D0 | Section3 first paragraph |
| Drive phase | −Pa sin(2pi f t), zero at t=0 | after Eq8 |
| Average | angle brackets in Eq11; window unspecified | Eq11. 64/256cycles and window comparisons are OUR NUMERICAL AUDIT, not paper values |
| Force coefficient | fB=rho/(4pi) mean(Vdot1 Vdot2), units N m2 | Eq11 |
| Plot scaling | fB/(1e-18 N m2) | Fig3 y-axis. This is dimensional unit scaling, not a physical nondimensionalization |
| Sign | positive attraction; negative repulsion | explicit text after Eq11 and Fig3 discussion |
| Vector on1 | e12=(x2−x1)/D; F1=(fB/D2)e12 at classical order | Inference from Eq11 plus attraction convention; not a sign change to match the graph |
| Volume | 4pi Rj3/3; derivative4pi Rj2 Rdotj | spherical geometry and Eq11 |
| Shell model | elastic law, no buckling/rupture branch | Eq9 and elastic-state paragraph |
| Solver/tolerance | not specified | OUR DOP853 rtol2e-10, vector SI atol and max stepT/96 |

No table is supplied for Fig3. A separately recorded approximate image digitization is not author data.
The Eq10 resonance approximation omits the surface-tension stiffness contribution that is present
when Eq9 is differentiated. That analytical approximation must not replace Eq9 in a reproduction.
The paper's Fig4 prose has apparent time-unit typos (um); the plotted horizontal units are ms.
These observations do not prove Fig3 is a typo.
