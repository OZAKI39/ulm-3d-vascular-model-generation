from sv13j_support import *
def test_actual_cuda123_three_kernels():
    d=actual('cuda123_runtime');assert len(d['runs'])==3
    assert '-ccbin' in d['compile']['command'] and '/usr/bin/g++-12' in d['compile']['command']
    for r in d['runs']:
        assert r['exit_code']==0 and not r['timeout']
        assert 'runtime=12030' in r['stdout'] and 'RTX 4090 arch=8.9' in r['stdout']
        assert 'correct=1 last_error=0' in r['stdout']
        assert sha256(ROOT/'logs/sv1_3j/remote'/Path(r['log']).name)==r['sha256']
def test_actual_thrust_generation_and_toolchain():
    d=actual('cuda123_versions');assert d['THRUST_VERSION']==200200
    toolchain_gate('12.3.2',12,d['nvcc_version'],d['nvcc_path'],d['prefix'],load('petsc_cuda123_build')['configure_command'])
    assert sha256(ROOT/'scripts/use_cuda123_env.sh')==d['wrapper_sha256']
def test_cuda_environment_mpi_regression():
    d=actual('cuda123_MPI');assert len(d['rank1'])==5 and d['MPI_wrapper_unchanged'] and not d['MPI_rebuilt']
    for r in d['rank1']:mpi_process_gate(r,1)
    mpi_process_gate(d['rank2'],2)
def test_wrapper_pointing_to_wrong_cuda_is_rejected():
    with pytest.raises(GateError):toolchain_gate('12.3.2',12,'release 12.6','/cuda126/bin/nvcc','/cuda123',[])

