# Postprocessing notes

- Added q/b/e Euclidean norm columns to BULK_BAQP_SHEAR_REFERENCE.csv directly from its saved component values, preserving those component strings and row order. The orchestration CSV writer now emits the same columns on a future independent reproduction. No oracle or force evaluation was repeated.
- The first kinematic-bias plot used a lower gradient-error axis limit of 1e-13. Visual inspection showed some interior RMS points below that display limit. Lowered the plot limit to 1e-14 and re-exported only that figure; data and thresholds did not change.
- Normal/tangential splitting uses diagonal/off-diagonal actual discrete Piola stress at fixed geometry. No unique node-to-face ownership is inferred.
- Reference-coordinate continuum moment retains its analytically nonzero source. Current-position angular momentum is evaluated separately.
