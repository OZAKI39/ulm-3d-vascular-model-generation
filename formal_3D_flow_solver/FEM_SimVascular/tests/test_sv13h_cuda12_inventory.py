from sv13h_support import *
def test_actual_cuda12_paths_and_preserved_default13():
    d=actual('cuda12_runtime')
    cuda_environment_gate(d,d['prefix'])
    before=load('pre_install_environment');after=load('final_environment')
    assert before['cuda_default_link']=='/usr/local/cuda-13.2'
    unchanged_driver_and_cuda13(before,after)
@pytest.mark.parametrize('change',[
    {'wrapper_enabled':False},
    {'nvcc_path':'/usr/local/cuda-13.2/bin/nvcc','nvcc_version':'release 13.2'},
    {'libcudart_path':'/usr/local/cuda-13.2/lib64/libcudart.so'},
])
def test_wrong_or_inactive_cuda12_environment_rejected(change):
    d=load('cuda12_runtime');d.update(change)
    with pytest.raises(GateError):cuda_environment_gate(d,d['prefix'])

