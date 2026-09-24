import numpy as np


def test_gradient_is_constant_on_each_side_and_jump_is_recorded(face_case, real_field):
    rows, checks = face_case
    for row in rows:
        actual = np.array([[row[f"gradient_{i}{j}_s_inv"] for j in range(3)] for i in range(3)])
        np.testing.assert_array_equal(actual, real_field.gradients_s_inv[row["tetra_id"]])
    assert all(np.isfinite(c["gradient_jump_s_inv"]) for c in checks)
    # A discontinuity is a measured property, never a continuity failure gate.
    print("measured gradient jumps (1/s)", [c["gradient_jump_s_inv"] for c in checks])
