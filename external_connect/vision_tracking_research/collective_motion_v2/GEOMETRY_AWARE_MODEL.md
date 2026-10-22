# Geometry-aware model

See `src/networks.py`, `src/geometry.py`, `src/inference.py` and `tests/test_math.py`. Coordinate conventions and full affine log conversion are shared with F10U. The Torch matrix-log series is restricted to the normalized near-identity training domain; the official SciPy log is the inference authority. Invalid/non-real transforms are rejected, not silently replaced by zeros.

Three-frame local matching produces endpoint and composed adjacent displacements. Learned confidence and endpoint/chain agreement weight a robust affine fit. Border locations whose correlation windows require padding are excluded. Heatmap support is a weak weighting term, not geometric ground truth. A soft support floor allows visual evidence where the detector misses points.

Local supervision does not ensure material identity on real microscopy. Stationary source dots, repeated microbubbles, dropout and finite matching range can dominate matching. Real outputs remain unverified. `geo_no_local_supervision` is a trained, same-data ablation, not a synthetic invented comparator. The final analysis records both its failures and the full model's failures.
