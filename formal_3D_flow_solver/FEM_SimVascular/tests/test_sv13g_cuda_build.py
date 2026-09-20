from sv13g_support import *
def test_actual_cuda_build_and_selftest():
    d=actual('cuda_build');assert d['self_test']=='PASS' and len(d['library_sha256'])==64
