from sv13h_support import *
def test_actual_cpu_gpu_scientific_equivalence():
    science_gate(actual('science_equivalence'),policy())
@pytest.mark.parametrize('key',['velocity_relative_L2','pressure_relative_L2','Qin_relative',
    'Qout_relative','mass_error_difference','max_velocity_relative'])
def test_non_equivalent_metric_rejected(key):
    d=equivalent_fixture()
    if key=='Qout_relative':d[key]['OUTLET_02']=1.
    else:d[key]=1.
    with pytest.raises(GateError,match='GPU_NUMERICAL_DIFFERENCE'):science_gate(d,policy())
def test_frozen_equivalence_thresholds():
    p=policy()
    for k in ('velocity_relative_L2','pressure_relative_L2','Qout_relative','max_velocity_relative'):assert p[k]==1e-5
    assert p['Qin_relative']==p['mass_error_difference']==1e-6
    science_gate(equivalent_fixture(),p)

