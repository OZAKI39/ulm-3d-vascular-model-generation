from particle_3d.particle9a1_audit import instrument

def test_no_default_observer_installation():
 obj=object();assert instrument(obj) is obj
