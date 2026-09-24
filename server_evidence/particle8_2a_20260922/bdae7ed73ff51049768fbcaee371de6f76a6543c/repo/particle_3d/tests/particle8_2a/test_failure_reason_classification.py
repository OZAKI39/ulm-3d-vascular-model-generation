import pytest
from particle_3d.particle82a_geometry import classify_failure


@pytest.mark.parametrize('extra,expected', [
    ({'perimeter_distance': .5e-6}, 'INLET_PERIMETER_CLEARANCE_FAIL'),
    ({'wall_distance': .5e-6}, 'WALL_INTERSECTION'),
    ({'wall_distance': 1.001e-6}, 'NEARFIELD_HANDOFF_VIOLATION'),
    ({'full_domain': True, 'cap_distance': 0.}, 'OPEN_CAP_ONLY_OUTSIDE_DOMAIN'),
    ({'owner': False}, 'NO_OWNER_TETRA'),
    ({'finite': False}, 'NUMERICAL_GEOMETRY_UNRESOLVED'),
    ({'pair_conflict': True}, 'PAIR_CONFLICT'),
    ({}, 'ACCEPTED')])
def test_distinct_causes(extra, expected):
    args = dict(radius=1e-6, wall_distance=5e-6); args.update(extra)
    assert classify_failure(**args) == expected
