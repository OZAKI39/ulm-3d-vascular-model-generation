from sv13l_support import *
def test_actual_mpi_wrappers_share_prefix():
    d=actual('mpi_fortran_build')
    for w in d['wrappers'].values():assert Path(w['path']).is_relative_to(d['prefix'])
    assert all('sv1_3g/external/gpu_mpi' not in s for s in d['linked_libraries'].values())
def test_mpi_source_archive_is_same():
    d=actual('openmpi_source');ref=load('reference_manifest')
    assert d['source']['sha256']==ref['MPI_source']['sha256']
