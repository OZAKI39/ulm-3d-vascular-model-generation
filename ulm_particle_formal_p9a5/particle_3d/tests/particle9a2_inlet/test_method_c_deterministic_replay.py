from method_c_test_helpers import source
from particle_3d.injection_method_c import canonical_bytes

def test_replay_is_byte_identical_and_order_independent(distribution):
    a=source(distribution);b=source(distribution)
    forward={i:canonical_bytes(a.event(i)) for i in range(1,25)}
    backward={i:canonical_bytes(b.event(i)) for i in range(24,0,-1)}
    assert forward==backward
