from sv13l_support import *
def test_actual_all_bindings_build():
    mpi_config_gate(actual('mpi_fortran_build'))
    assert 'or usempif08 (or all, build mpifh' in (ROOT/'logs/sv1_3l/remote/openmpi_configure_help.log').read_text()
def test_disabled_bindings_rejected():
    d=load('mpi_fortran_build');d['configure_command'].append('--disable-mpi-fortran')
    with pytest.raises(GateError,match='DISABLED'):mpi_config_gate(d)
def test_missing_wrapper_rejected():
    d=load('mpi_fortran_build');d['wrappers'].pop('mpifort')
    with pytest.raises(GateError,match='WRAPPER_MISSING'):mpi_config_gate(d)
def test_mixed_compiler_major_rejected():
    d=load('mpi_fortran_build');d['compiler_environment']['CC']='/usr/bin/gcc-13'
    with pytest.raises(GateError,match='MAJOR_MISMATCH'):mpi_config_gate(d)
