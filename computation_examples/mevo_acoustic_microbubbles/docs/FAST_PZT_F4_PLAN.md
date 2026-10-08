# F4 frozen geometry-only integration plan

Base F3 archive: 3eeb94eb885f5df4f0e214ba9b7904dbdecec744.
Branch: dev/fast-pzt-f4-y-geometry-20261008. No F1/F2/F3 mutation.

First confirm the server's original labelled Stage-4 VTP and its requested hash,
audit closed outward orientation and tags, and measure the actual lumen/field
intersection. Initialization uses 32 nonoverlapping spheres per actual parent and
junction wall patch, real 1.78 um radii and original SI mesh coordinates. Rejection
sampling only; no invalid initial positions corrected into the vessel.

Freeze F2 feedback/central F1 PBD, two local mobility and twelve pair-contact passes,
10 us actual step and one substep. Default stable reduction. Wall additions are
geometric normal corrections only, after prediction and after pair passes 3/6/9/12.
No flow, wall hydrodynamics, scattering, friction, adhesion, or group locking.
Wall and pair overlap budgets remain 0.01R; per-pass correction remains 0.1R.

Static Warp mesh BVHs remain on GPU. Float32 broad-phase geometry is independently
checked against the original triangle coordinates; local finite-triangle distances
are refined in float64 without loosening the budgets. Inside/outside comes from the
closed original surface, physical distance from the WALL subset. Caps are open.
Candidate physical displacement and final wall-corrected paths get a cheap
distance-bound fast path or local swept finite-triangle checks. Unresolved paths,
excess correction, pair/wall quality failures must not be accepted silently.
The field boundary emits FIELD_DOMAIN_EXIT and ends that particle's tracking;
it is not a reflecting wall. True finite-cap crossings emit TERMINAL_EXIT with ID,
label and physical time. Inactive IDs remain recorded but leave force/contact loops.

Bounded geometry tests distinguish imposed GEOMETRY_UNIT_TEST velocities from
actual acoustics. Cover signs, triangle/edge precision, slide, near-wall pairs,
unsafe path rejection, three finite open caps, rim safety, field exit and exact
F3 equivalence when vessel is disabled. Then run the two initialized patches for
120 ms with real baseline PZT and paired WALL OFF/ON at identical inputs. Record
whether actual wall contact and aggregation occur, without tuning forces to get it.
Choose the parent patch as the predeclared continuation region if it has a connected
cluster, otherwise the junction if it does; if neither forms a cluster, report that
and use the largest connected state without declaring maturity. Continue matched
unchanged and F3 amplitude/enable controls for 120 ms; retain all original IDs.
No new control search. Inherited F3 initial/model sensitivity limits still apply.

If 10 us fails in a representative near-wall trajectory, only that representative
may be compared at 5 us. No broad scan. One 32-bubble matched wall-cost comparison
uses the scientific pair; after functional gates pass, one N200/5ms regression
per wall mode quantifies overhead and GPU kernel costs. No old R3 task resumes.

Three English figures and one actual-frame movie show the original global Y and
a synchronized local view at true sphere radii. The original mesh is parametric,
not patient-specific. All claims remain EXPLORATORY and GEOMETRIC_WALL_ONLY.
Execution SHAs, source/mesh/input hashes, failures, logs and closed evidence are
preserved separately from the final archive commit and pushed on the new branch.
