from sv13g_support import *
def test_actual_environment_fingerprint():
    d=load('mpi_environment')
    assert all(k in d for k in ('uid','paths','realpaths','environment','probes'))
    assert all(k in d['environment'] for k in ('PATH','LD_LIBRARY_PATH','CUDA_VISIBLE_DEVICES'))
    assert any('OpenRTE' in p['stdout'] for p in d['probes'])
