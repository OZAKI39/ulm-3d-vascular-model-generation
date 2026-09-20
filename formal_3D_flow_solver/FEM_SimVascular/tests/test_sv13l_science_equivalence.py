from sv13l_support import *
def test_actual_gpu_cpu_equivalence():
    science_gate(actual('science_equivalence'),policy())
@pytest.mark.parametrize('key',['velocity_relative_L2','pressure_relative_L2','Qin_relative','mass_error_difference','max_velocity_relative'])
def test_changed_scientific_answer_rejected(key):
    d=dict(velocity_relative_L2=0.,pressure_relative_L2=0.,Qin_relative=0.,Qout_relative={'OUTLET_01':0.,'OUTLET_02':0.,'OUTLET_03':0.},mass_error_difference=0.,max_velocity_relative=0.)
    d[key]=1.
    with pytest.raises(GateError,match='NUMERICAL_DIFFERENCE'):science_gate(d,policy())
def test_proof_input_changes_limited_to_start_and_stop():
    d=load('flow_input_manifest')
    assert d['vascular_mesh_physics_BC_dt_nonlinear_method_unchanged']
    assert {x['tag'] for x in d['vascular_XML_differences']}=={'Continue_previous_simulation','Number_of_time_steps'}

