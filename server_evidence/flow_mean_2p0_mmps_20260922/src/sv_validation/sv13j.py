"""Fail-closed compatibility and GPU bring-up gates for SV1.3J."""
from pathlib import Path
from .sv13h import (GateError,require,cuda_types_gate,mpi_process_gate,linkage_gate,
                   restart_gate,proof_gate,science_gate,median_benchmark,performance_classification)

CANDIDATES=('12.3.2','12.2.2','12.1.1')

def toolchain_gate(version,host_major,nvcc_version,nvcc_path,prefix,flags):
    require(version in CANDIDATES,'CANDIDATE_NOT_AUTHORIZED')
    require(host_major==12,'HOST_COMPILER_NOT_GCC12')
    require('release '+'.'.join(version.split('.')[:2]) in nvcc_version,'WRONG_CUDA_VERSION')
    require(Path(nvcc_path).is_relative_to(Path(prefix)),'WRONG_CUDA_PREFIX')
    normalized=' '.join(flags).lower().replace('_','-')
    require('allow-unsupported-compiler' not in normalized and '-override' not in normalized,'UNSUPPORTED_COMPILER_OVERRIDE')

def source_integrity_gate(record):
    require(record['fresh_source'] and not record['previous_objects_reused'],'CANDIDATE_OBJECT_REUSE')
    require(record['modifications']==[],'PETSC_SOURCE_MODIFIED')

def isolation_gate(records):
    for key in ('source','PETSC_ARCH','prefix'):
        values=[r[key] for r in records]
        require(len(set(values))==len(values),'CANDIDATE_ISOLATION_VIOLATION: '+key)

def compatibility_gate(record):
    require(record['kernel']=='PASS' and record['configure']=='PASS' and record['make']=='PASS','UNPATCHED_CANDIDATE_NOT_COMPATIBLE')

def matrix_order_gate(records):
    ran=[r for r in records if r['status']!='NOT_REQUIRED']
    require([r['CUDA_version'] for r in ran]==list(CANDIDATES[:len(ran)]),'CANDIDATE_ORDER_VIOLATION')
    winner=None
    for i,r in enumerate(ran):
        require(winner is None,'SEARCH_CONTINUED_AFTER_SUCCESS')
        if r['make']=='PASS':compatibility_gate(r);winner=r['CUDA_version']
        elif i<len(ran)-1:require(r['status']=='FAIL','PREVIOUS_CANDIDATE_NOT_FAILED')
    return winner

def application_mpi_gate(rows):
    needed={'MPI_INTEGER','MPI_DOUBLE_PRECISION','MPI_CHARACTER'}
    require(needed <= {r['name'] for r in rows},'APPLICATION_MPI_TYPES_UNTESTED')
    require(all(r['size']>0 and r['bcast_rc']==0 for r in rows if r['name'] in needed),'APPLICATION_MPI_DATATYPE_UNAVAILABLE')
