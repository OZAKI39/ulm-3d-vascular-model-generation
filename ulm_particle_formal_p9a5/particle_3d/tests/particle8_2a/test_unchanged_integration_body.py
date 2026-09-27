import inspect
from particle_3d.particle81_simulation import integrate_one
from particle_3d.particle82a_integration import integrate_admitted


def test_dynamics_and_guard_body_is_literal_original():
    old=inspect.getsource(integrate_one);new=inspect.getsource(integrate_admitted)
    # Everything from flow provider through termination, accepted samples and
    # cap crossing audit is byte-identical to P8.1. Only birth initialization differs.
    start='        def provider('
    assert old[old.index(start):]==new[new.index(start):]
