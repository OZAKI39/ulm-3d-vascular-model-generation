from dataclasses import fields
from particle_3d.microbubble import MicrobubbleState


def test_inertial_quantities_and_orientation_are_not_state_variables():
    names={f.name for f in fields(MicrobubbleState)}
    assert not names.intersection({'mass','mass_kg','acceleration','linear_acceleration','angular_acceleration','force','quaternion'})
