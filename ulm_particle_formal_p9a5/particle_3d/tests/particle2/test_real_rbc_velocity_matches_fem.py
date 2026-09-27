import numpy as np


def test_every_saved_velocity_and_gradient_matches_unchanged_p0(real_cases,center_samples):
    for rows,summary in real_cases:
        samples=center_samples[summary["dt_index"]]
        for prefix in ["V","u"]:
            np.testing.assert_array_equal([[r[f"{prefix}_{a}_m_s"] for a in "xyz"] for r in rows],samples.velocity_m_s)
        np.testing.assert_array_equal([r["tetra_id"] for r in rows],samples.tetra_id)
        np.testing.assert_array_equal([[r[f"G_{i}{j}_s_inv"] for i in range(3) for j in range(3)] for r in rows],samples.velocity_gradient_s_inv.reshape(-1,9))
        dt=summary["validation_dt_s"]
        for a,b in zip(rows[:-2],rows[1:-1]):
            expected=np.array([a[f"{s}_m"] for s in "xyz"])+dt*np.array([a[f"V_{s}_m_s"] for s in "xyz"])
            np.testing.assert_array_equal(expected,[b[f"{s}_m"] for s in "xyz"])
