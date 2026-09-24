from particle_3d.particle9a_diagnostics import triangle_change

def test_triangle_identity_and_normal_angle_are_independent():
    a=dict(nearest_wall_triangle_id=1,wall_normal_xyz=[0.,0.,1.])
    b=dict(nearest_wall_triangle_id=2,wall_normal_xyz=[0.,0.,1.])
    c=dict(nearest_wall_triangle_id=2,wall_normal_xyz=[0.,1.,0.])
    assert triangle_change(a,a)==dict(wall_triangle_changed=False,normal_change_angle_from_previous_deg=0.)
    assert triangle_change(a,b)['wall_triangle_changed']
    assert triangle_change(a,b)['normal_change_angle_from_previous_deg']==0.
    assert not triangle_change(b,c)['wall_triangle_changed']
    assert triangle_change(b,c)['normal_change_angle_from_previous_deg']==90.
