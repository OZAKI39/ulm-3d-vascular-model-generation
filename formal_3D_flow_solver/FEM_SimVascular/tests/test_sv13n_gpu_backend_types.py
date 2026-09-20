from sv13n_support import *
def test_actual_official_cuda_types():flow_gate(accepted('svmp_gpu_smoke'),gpu=True)
@pytest.mark.parametrize('key,value',[('mat_type','seqaij'),('vec_type','seq')])
def test_cpu_masquerading_as_gpu_rejected(key,value):
 d=good_flow();d[key]=value
 with pytest.raises(GateError):flow_gate(d,gpu=True)
