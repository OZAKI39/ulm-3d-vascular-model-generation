from sv13h_support import *
def test_three_real_cuda_kernels_and_log_hashes():
    d=actual('cuda12_runtime')
    assert len(d['runs'])==3 and d['compile']['exit_code']==0
    for r in d['runs']:
        assert r['exit_code']==0 and not r['timeout'] and r['correct']
        assert 'runtime=12060' in r['stdout'] and 'RTX 4090 arch=8.9' in r['stdout']
        assert 'elements=1024 correct=1 last_error=0' in r['stdout']
        assert sha256(ROOT/'logs/sv1_3h/remote'/Path(r['log']).name)==r['sha256']
def test_runtime_link_and_native_binary_mirror():
    d=load('cuda12_runtime')
    assert d['prefix']+'/lib64/libcudart.so.12' in d['ldd'] and 'not found' not in d['ldd']
    assert sha256(ROOT/'outputs/sv1_3h/native/cuda_kernel_smoke')==d['binary_sha256']
    assert sha256(ROOT/'scripts/use_cuda12_gpu_env.sh')==d['wrapper_sha256']

