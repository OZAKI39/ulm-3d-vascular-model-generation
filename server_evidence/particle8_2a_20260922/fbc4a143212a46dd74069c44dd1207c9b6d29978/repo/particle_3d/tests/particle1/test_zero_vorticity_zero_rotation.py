
def test_zero_rotation_at_every_uniform_step(uniform_rows):
    assert all(r['angular_velocity_norm_s_inv']==0. for r in uniform_rows)
