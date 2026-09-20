from sv_validation.validation import fields_finite

def test_actual_velocity_and_pressure_finite(measured):
    assert fields_finite(measured[1],measured[2])
