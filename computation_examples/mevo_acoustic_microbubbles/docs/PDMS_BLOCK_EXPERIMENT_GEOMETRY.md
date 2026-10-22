# F10M confirmed geometry and unresolved mounting

The user confirms a centered 32.5 × 9 × 9 mm PDMS block, with a needle-cast through-channel along X, radius 0.8 mm. The block is centered at the origin. PDMS is exactly the rectangular solid minus the lumen; there is no independent tube wall, outer diameter or 0.1 mm shell. The side PDMS path is 4.5 − 0.8 = 3.7 mm. Earlier F10L dimensions remain historical demonstration data.

Exactly three SMD07T05R412WL ceramic faces have diameter 7 mm. Centers in mm and inward normals:

| ID | Center | Normal |
|---|---|---|
| side_y_minus | (0, −4.5, 0) | (0, 1, 0) |
| side_y_plus | (0, 4.5, 0) | (0, −1, 0) |
| end_x_plus | (16.25, 0, 0) | (−1, 0, 0) |

Side faces contact the block. No Z source exists. The centered end face covers the entire main-channel opening, area π(0.8 mm)² = 2.010619 mm²; 5.22449% of the disc face. The remaining PDMS contact area is 36.473891 mm².

During F10M the user clarified an open top injection port, diameter 1.5 mm, connected downward to the center channel, located 5.4 mm along the length from the end PZT; X+ is always closed. These supersede the initial unresolved-port snapshot, preserved under `history/pre_top_port/`. The new model has an open vertical cylindrical branch unioned with the main lumen and a closed X+ end. The provisional face-to-hole-center interpretation gives x=16.25−5.4=10.85 mm; y=0 is an explicit lateral-centering assumption. The drilled port is modeled from z=0 to z=4.5 mm, without an unconfirmed blind extension below the main-channel center. Actual lateral position, reference for 5.4 mm and blind depth should be confirmed; diameter, connection and open/closed state were confirmed.

X− remains `open_fail_stop` as an unconfirmed original computational endpoint. The top port is open and X+ closed based on the user's reply, never inferred from rendering. Acoustic sealing, bonding and the solid PZT mechanical transfer remain uncalibrated even though the particle boundary state is now known. Full fields remain NOT_RUN because of resources and model limitations.

`CylindricalLumen` reuses the F10L coordinate, gap and containment operations. The GPU cylinder projection and complete M0 solver are unchanged. A new admission guard stops particles before the top-port axial neighborhood; junction transport is NOT_ADMITTED rather than approximated with a false cylindrical wall. This bounded OFF experiment only validates the central main bore. Initial cloud coordinates are the post-injection demo, not measured positions. Every particle retains its identity.

Whole-apparatus translation and proper rotation are supported. The actual three faces remain parameterized; a position outside the confirmed setup must be re-audited before field generation. Rendering is in millimeters with true aspect ratios. Figure opacity is a display setting only.
