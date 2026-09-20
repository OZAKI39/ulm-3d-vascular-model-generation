from sv13g_support import *
from sv_validation.sv13g import *
from sv_validation.provenance import sha256
def test_single_working_wrapper():
    d=actual('mpi_resolution');assert d['rank1_passed']==5 and d['rank2']=='PASS'
    assert sha256(ROOT/'scripts/run_gpu_mpi.sh')==d['wrapper_sha256']
    assert d['system_MPI_modified'] is False
    assert d['fallback_count']<=1
def test_root_is_not_silently_overridden():
    with pytest.raises(GateError,match='ROOT_OVERRIDE'):explicit_root_gate('Open MPI',0,['mpiexec','-n','1'],{})
    explicit_root_gate('Open MPI',0,['mpiexec','--allow-run-as-root'],{})
