"""Reproduce datatype failure without changing MPI, PETSc, or the solver."""
import json,re,shutil
from pathlib import Path
from runner_remote import *
mpi=json.loads((BASE/'configs/mpi_resolution.json').read_text());w=load('compatibility_winner')
prefix=Path(mpi['prefix']);binary=B/'mpi_datatype_probe'
extra={'OMPI_CC':'/usr/bin/gcc-12','OMPI_CXX':'/usr/bin/g++-12'}
compile=run([prefix/'bin/mpicc',B/'mpi_datatype_probe.c','-o',binary],'mpi_datatype_compile',cuda=w['candidate_wrapper'],extra_env=extra)
assert okay(compile)
runs=[]
for n in (1,2):
    r=run([mpi['working_launcher'],'-n',str(n),w['candidate_wrapper'],binary],'mpi_datatype_rank'+str(n),timeout=10,cuda=w['candidate_wrapper'])
    r['stdout']=text(r);r['rows']=[]
    for line in text(r).splitlines():
        if line.startswith('DATATYPE '):
            row=dict(field.split('=') for field in line.split()[1:])
            r['rows'].append({k:v if k=='name' else int(v) for k,v in row.items()})
    runs.append(r)
info=run([prefix/'bin/ompi_info','--parsable','--all'],'mpi_datatype_ompi_info',cuda=w['candidate_wrapper'])
evidence=R/'mpi_datatype_evidence';evidence.mkdir()
sources=[BASE.parent/'sv1_3g/external/openmpi-4.1.6/ompi/datatype/ompi_datatype_module.c',
         BASE.parent/'sv1_3g/external/openmpi-4.1.6/opal_config.h',
         BASE/'external/solver-inputs/svMultiPhysics/Code/Source/solver/CmMod.h']
mirrors=[]
for p in sources:
    if p.exists():shutil.copyfile(p,evidence/p.name);mirrors.append({'source':str(p),'sha256':digest(p),'mirror':str((evidence/p.name).relative_to(BASE))})
c_names={'MPI_INT','MPI_DOUBLE','MPI_CHAR','MPI_CXX_BOOL'}
fortran={'MPI_INTEGER','MPI_DOUBLE_PRECISION','MPI_CHARACTER','MPI_LOGICAL'}
rows=[row for r in runs for row in r['rows']]
reproduced=all(not r['timeout'] and len(r['rows'])==8*r['rows'][0]['ranks'] for r in runs)
reproduced=reproduced and all(row['bcast_rc']==0 and row['type_size_rc']==0 for row in rows if row['name'] in c_names)
reproduced=reproduced and all(row['error_class']==3 and row['size']==0 for row in rows if row['name'] in fortran)
write('mpi_datatype_diagnosis',{'status':'CONFIRMED' if reproduced else 'UNRESOLVED','compile':compile,'runs':runs,
    'source':str(B/'mpi_datatype_probe.c'),'source_sha256':digest(B/'mpi_datatype_probe.c'),
    'MPI_wrapper_sha256':digest(mpi['working_launcher']),'MPI_changed':False,'solver_changed':False,
    'working_C_datatypes':sorted(c_names),'failed_Fortran_datatypes':sorted(fortran) if reproduced else None,
    'source_evidence':mirrors,'ompi_info':info,'reason':'SVMULTIPHYSICS_REQUIRES_MPI_FORTRAN_PREDEFINED_TYPES'})
print('MPI datatype diagnosis: '+('CONFIRMED' if reproduced else 'UNRESOLVED'),flush=True)
