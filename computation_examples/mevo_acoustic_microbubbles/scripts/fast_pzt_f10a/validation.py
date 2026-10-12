"""Independent CPU verification of the GPU ledger and actual accepted samples."""

import numpy as np
from .common import ROOT, js


def validate_ledger(ledger, positions, sample_times):
    q = js(ROOT / "configs/fast_pzt/f9/protocol.json")["ledger"]
    d, meta, events = [ledger[k] for k in ["displacement_m", "step_meta", "wall_contact_events"]]
    assert not ledger["rejected_intervals"], "Rejected audit intervals are not accepted training samples"
    assert np.isfinite(d).all() and np.isfinite(meta).all() and np.isfinite(events).all()
    assert np.all(meta[:, 2] == 1)
    assert meta[:, 3].max(initial=0) < q["per_particle_step_residual_m"]
    assert np.linalg.norm(d[:, 8].sum(0)) < q["cumulative_group_residual_m"]
    # Higher interior advection changes floating-point cancellation scale, not physics budgets.
    # CPU summation vs GPU reductions: conservative 32-epsilon bound on summed operands.
    arithmetic_bound = 32 * np.finfo(float).eps * (np.abs(d[:, 7]) + np.abs(d[:, :7]).sum(1)) + 1e-23
    assert np.all(np.abs(d[:, 7] - d[:, :7].sum(1) - d[:, 8]) <= arithmetic_bound)

    cumulative = np.vstack([np.zeros((1, 9, 3)), np.cumsum(d, axis=0)])
    ends = np.r_[0, meta[:, 0] + meta[:, 1]]
    idx = np.searchsorted(ends, sample_times - 1e-12)
    assert np.max(np.abs(ends[idx] - sample_times)) < 1e-11
    c = positions[:, 0].mean(1)
    np.testing.assert_allclose(c - c[0], cumulative[idx, 7], atol=5e-17, rtol=0)
    wall = np.zeros_like(d[:, 5])
    wall_lengths = np.zeros(len(d))
    if len(events):
        step = events[:, 0].astype(int)
        np.add.at(wall, step, events[:, 11:14] / positions.shape[2])
        np.add.at(wall_lengths, step, events[:, 14] / positions.shape[2])
        np.testing.assert_allclose(events[:, 8:11] - events[:, 5:8], events[:, 11:14], atol=0, rtol=0)
        assert set(events[:, 4].astype(int)).issubset({0, 3, 6, 9, 12})
        assert events[:, 19].min() >= 0
        assert events[:, 23].max() < q["wall_query_reconstruction_error_m"]
        expected = np.minimum(events[:, 15], 0)[:, None] * events[:, 16:19]
        np.testing.assert_allclose(
            events[:, 11:14], expected, rtol=0, atol=q["wall_query_reconstruction_error_m"]
        )
    np.testing.assert_allclose(wall, d[:, 5], atol=1e-22, rtol=1e-10)
    np.testing.assert_allclose(wall_lengths, meta[:, 5], atol=1e-22, rtol=1e-10)
    return dict(
        status="PASS",
        steps=len(d),
        wall_event_count=len(events),
        max_particle_step_residual_m=float(meta[:, 3].max(initial=0)),
        cumulative_group_residual_m=float(np.linalg.norm(d[:, 8].sum(0))),
        max_wall_metadata_reconstruction_m=float(meta[:, 6].max(initial=0)),
        mean_wall_vector_m=d[:, 5].sum(0).tolist(),
        mean_wall_length_m=float(meta[:, 5].sum()),
        semantics="Within-trajectory exact ledger, not independent physical causal forces or a no-wall counterfactual",
    )
