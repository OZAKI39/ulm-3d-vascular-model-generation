from sv13n_support import *
def test_official_cpu_fields_and_convergence():flow_gate(accepted('svmp_cpu_smoke'),gpu=False)
@pytest.mark.parametrize('key,value',[('overlap',1),('fill_level',3),('rtol',1e-8),('restart',30),('ordering','rcm')])
def test_runtime_solver_override_rejected(key,value):
 d=good_flow();d['runtime_semantics'][0][key]=value
 with pytest.raises(GateError):flow_gate(d,gpu=False)
