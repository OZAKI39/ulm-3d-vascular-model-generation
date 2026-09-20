from sv13h_support import *
def test_actual_configure_cuda12_sm89_real_double_cxx17():
    d=load('petsc_cuda12_build')
    assert d['configure_exit']==0 and d['cuda_enabled']
    c=d['configure_command']
    for x in ('--with-cuda=1','--with-cuda-arch=89','--with-cxx-dialect=C++17',
              '--with-cuda-dialect=C++17','--with-debugging=0','--with-precision=double',
              '--with-scalar-type=real','--with-64-bit-indices=0'):
        assert x in c
    flags=d['compiler']
    assert d['CUDA_prefix'] in flags['CUDAC']
    assert d['MPI_prefix'] in flags['CC'] and d['MPI_prefix'] in flags['CXX']
    for k in ('CXX_FLAGS','CUDAC_FLAGS'):
        assert '-std=c++17' in flags[k] and '-std=c++20' not in flags[k]
    assert '--download-hypre' not in ' '.join(c)
def test_actual_generated_header_and_logged_configure():
    d=load('petsc_cuda12_build')
    assert '#define PETSC_HAVE_CUDA 1' in (R/'remote/build_failure_evidence/petscconf.h').read_text()
    r=d['steps'][0]
    assert sha256(ROOT/'logs/sv1_3h/remote'/Path(r['log']).name)==r['sha256']

