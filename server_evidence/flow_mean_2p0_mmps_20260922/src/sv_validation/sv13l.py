"""MPI application compatibility and GPU flow acceptance, with explicit failure gates."""
from pathlib import Path
from .sv13j import (GateError, require, cuda_types_gate, proof_gate, science_gate,
                   median_benchmark, performance_classification, restart_gate,
                   mpi_process_gate, source_integrity_gate)

FORTRAN_TYPES={'MPI_INTEGER','MPI_DOUBLE_PRECISION','MPI_CHARACTER','MPI_LOGICAL'}

def mpi_config_gate(record):
    flags=record['configure_command']
    require('--disable-mpi-fortran' not in flags,'MPI_FORTRAN_DISABLED')
    require('--enable-mpi-fortran=all' in flags,'FULL_FORTRAN_BINDINGS_NOT_ENABLED')
    require(record['version']=='4.1.6','MPI_VERSION_CHANGED')
    require('mpifort' in record['wrappers'],'FORTRAN_WRAPPER_MISSING')
    require('gfortran-12' in record['wrappers']['mpifort']['showme'],'FORTRAN_COMPILER_MISMATCH')
    for key,name in [('CC','gcc-12'),('CXX','g++-12'),('FC','gfortran-12')]:
        require(Path(record['compiler_environment'][key]).name==name,'MPI_COMPILER_MAJOR_MISMATCH')

def datatype_rows_gate(rows,ranks,sizes):
    for rank in range(ranks):
        observed=[r for r in rows if r['rank']==rank and r['name'] in FORTRAN_TYPES]
        require({r['name'] for r in observed}==FORTRAN_TYPES and len(observed)==4,'FORTRAN_DATATYPE_COVERAGE_MISSING')
        for r in observed:
            require(r['ranks']==ranks,'WRONG_PROBE_RANK_COUNT')
            require(r['size']>0 and r['size']==sizes[r['name']],'MPI_FORTRAN_DATATYPE_SIZE_INVALID')
            require(r['type_size_rc']==0 and r['bcast_rc']==0 and r['error_class']==0,'MPI_FORTRAN_BCAST_FAILED')
            require(r['data_ok']==1,'MPI_FORTRAN_PAYLOAD_MISMATCH')

def mpi_application_gate(record):
    require(record['C_rank1_passed']==5 and record['C_rank2']=='PASS','MPI_C_BASIC_FAILED')
    require(record['Fortran_rank1']==record['Fortran_rank2']=='PASS','MPI_FORTRAN_BASIC_FAILED')
    require(record['datatype_repeats_passed']==3,'MPI_FORTRAN_PROBE_REPETITIONS_FAILED')

def stack_link_gate(resolved,expected,label):
    require(Path(resolved)==Path(expected),'WRONG_'+label+'_LINKAGE')

def official_flow_gate(record):
    require(record['startup_pass'] and record['exit_code']==0,'OFFICIAL_FLOW_STARTUP_OR_EXIT_FAIL')
    require(not record.get('PETSc_error_detected',False),'OFFICIAL_FLOW_PETSC_ERROR')
    require(record['linear_solves']>0 and record['linear_failures']==record['nonlinear_failures']==0,'OFFICIAL_FLOW_CONVERGENCE_FAIL')
    require(record['VTU_count']>0 and record['velocity_finite'] and record['pressure_finite'] and record['reload_pass'],'OFFICIAL_FLOW_FIELD_FAIL')
    cuda_types_gate(record['mat_type'],record['vec_type'])
    require(record['KSP']=='gmres' and record['PC']=='asm','OFFICIAL_FLOW_SOLVER_SEMANTICS_CHANGED')
