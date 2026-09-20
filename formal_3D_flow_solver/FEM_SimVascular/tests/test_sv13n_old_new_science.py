from sv13n_support import *
def test_old_new_native_science():science_gate(accepted('old_new_cpu_science'))
@pytest.mark.parametrize('key',list(SCIENCE_LIMITS))
def test_each_scientific_threshold_enforced(key):
 d=dict(pressure_shift_applied=False,errors={k:0 for k in SCIENCE_LIMITS});d['errors'][key]=SCIENCE_LIMITS[key]*1.01
 with pytest.raises(GateError):science_gate(d)
def test_pressure_shift_forbidden():
 d=dict(pressure_shift_applied=True,errors={k:0 for k in SCIENCE_LIMITS})
 with pytest.raises(GateError):science_gate(d)
