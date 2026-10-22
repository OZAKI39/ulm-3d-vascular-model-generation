# F10M-R3 geometry and boundary contract

The 32.5 × 9 × 9 mm PDMS block contains one centered X-directed bore of radius
0.8 mm and two Z-directed holes of radius 0.75 mm. The geometric void is their
union. Hole axes are provisionally X = ±10.85 mm, Y = 0, using a 5.4 mm distance
from the adjacent end-PZT face plane to the hole center. Center-versus-edge
measurement remains unresolved. Holes connect downward to the bore without
perforating the bottom; no unreported blind extension is added.

Four 7 mm source faces are at (0, ±4.5, 0) and (±16.25, 0, 0) mm with inward
normals. There are no Z source faces. Planar source faces are distinct from
unconfirmed solid PZT packaging and seal stacks. They support four independent
amplitudes, phases and enable flags through the existing PZTArray / kWaveArray.

Both X ends are always mechanically sealed. The surfaces at ±16.25 mm exclude
finite-radius bubbles. Material, physical thickness, bonding, insulation and
acoustic transmission/reflection of these seals are independently unknown.
No PDMS plug or rigid acoustic wall is inferred from mechanical closure.

Both top mouths are geometrically open. Fill height, phase, meniscus, wetting and
external coupling medium are unknown. Region labels describe void topology,
not proof that a hole is liquid-filled. Unresolved medium arrays carry NaN and
are deliberately rejected by the acoustic adapter. Boundary IDs distinguish
sealed X interfaces (1/2) and open top mouths (3/4) from volume labels.

The implicit union membership is not a true signed distance. PDMSCavity queries
exposed cylinder patches, closed X disks and cylinder-intersection curves.
PDMSCavitySimulator applies the same geometry after prediction and PBD, using
the original M0 force/contact path and 0.1R wall-correction budget. Legacy open-X
code remains only for earlier experiments; R3 does not inherit that adapter.

Only the declared interior z+R < 1.5 mm scope is admitted in this revision.
This is NOT an assumed physical liquid level. Crossing it is a scope failure,
without projection or reflection. OFF checks at interior port locations assume
those test regions are wet and are explicitly numerical DEMO checks. The
geometry query at junctions is tested; free-surface bubble physics is not.

All full-apparatus fields and ON/navigation cases remain NOT_RUN until the
unknown physical interfaces, driving frequency, response table and resource
admission are resolved. Zero active flow does not model injection residuals,
ultrasound streaming, buoyant/free-surface effects or evaporation.
