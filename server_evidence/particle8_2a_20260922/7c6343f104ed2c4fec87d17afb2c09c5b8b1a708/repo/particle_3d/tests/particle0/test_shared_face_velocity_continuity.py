import numpy as np


def test_two_traces_equal_independent_face_mean(face_case, real_field, faces):
    rows, checks = face_case
    for check in checks:
        assert check["velocity_trace_error"] <= check["velocity_bound"], check
    for face in faces:
        sample = real_field.sample(face["center"])
        assert sample.inside_lumen and sample.tetra_id == min(face["ids"])
        for fraction, cell in [(-1., face["ids"][0]), (1., face["ids"][1])]:
            position = face["center"] + fraction * face["distance"] * face["normal"]
            sample = real_field.sample(position)
            assert sample.tetra_id == cell
            expected = real_field.velocity_nodes_m_s[face["nodes"]].mean(0) + real_field.gradients_s_inv[cell] @ (position - face["center"])
            bound = 8 * real_field.geometry.weight_tolerance[cell] * np.max(np.abs(real_field.velocity_nodes_m_s[real_field.tetra[cell]]))
            np.testing.assert_allclose(sample.velocity_m_s, expected, atol=bound, rtol=0)
