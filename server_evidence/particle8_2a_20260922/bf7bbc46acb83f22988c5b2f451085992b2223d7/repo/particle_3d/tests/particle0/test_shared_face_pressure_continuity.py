import numpy as np


def test_pressure_two_traces_are_continuous(face_case, real_field, faces):
    for check in face_case[1]:
        assert check["pressure_trace_error"] <= check["pressure_bound"], check
    for face in faces:
        value = real_field.sample(face["center"]).pressure_pa
        expected = real_field.pressure_nodes_pa[face["nodes"]].mean()
        bound = max(c["pressure_bound"] for c in face_case[1])
        np.testing.assert_allclose(value, expected, atol=bound, rtol=0)
