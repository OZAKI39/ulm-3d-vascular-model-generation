# Mechanical load requirement bracketing

Entry: `campaign.py` runs the initial sequential shear decades with evidence-based early stopping. `campaign.execute(mode, peak)` submits one additional frozen-amplitude case to the native subprocess. All cases use the same 2 ms half-cosine rise and 0.5 ms hold, original geometry/material/anchors, and the accepted native library.

`prescribed_traction.py` defines nominal reference-surface dead loads. Local inward/outward patches remain material-attached, and the normal/tangent stays frozen. Thus the high-deformation response is **not** a validated current-surface follower-pressure experiment. Load area, force, moment and work interpolation are explicitly audited.

`native_case.py` reuses the accepted Native/Contact classes and material model without editing them. Only the external traction, evidence retention and topology bookkeeping are new. The full specimen reference is stored once, with relative symlinks per case. Key states retain full precision, while histories are scalar CSV and compressed damage updates. Failures retain the unaccepted trial and last accepted state separately.

`analyze.py`, `report.py`, and `validate_outputs.py` produce tables, six figures, the 13-question review, eight final statuses and evidence checks. Run only on the frozen task run; new physical protocols require a new run directory. No historical production state or native source is overwritten.
