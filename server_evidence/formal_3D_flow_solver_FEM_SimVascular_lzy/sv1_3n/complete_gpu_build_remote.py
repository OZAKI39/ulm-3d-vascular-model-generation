"""Finish provenance after a case-sensitive log parser error; no rebuild or rerun."""
import re,shutil
from runner_remote import *
d=load('petsc_gpu13_build');assert d['status']=='BUILDING'
assert all(okay(r) and digest(BASE/r['log'])==r['sha256'] for r in d['steps'])
check=d['steps'][-1];assert check['name']=='gpu13_configfix_check'
raw=text(check);observed={k:s in raw.lower() for k,s in [('CPU','run successfully with 1 mpi process'),('MPI2','run successfully with 2 mpi processes'),('CUDA','run successfully with cuda')]}
assert all(observed.values()) and not re.search(r'Possible (?:error|problem)',raw)
source=Path(d['source']);prefix=Path(d['prefix']);integrity=load('petsc_gpu13_configfix_source_integrity')
assert integrity['source']==str(source)
changes=[n for n,h in integrity['files'].items() if digest(source/n)!=h];assert not changes
evidence=R/'gpu13_configfix_build_evidence';evidence.mkdir();conf=source/d['PETSC_ARCH']/'lib/petsc/conf'
for p in (conf/'configure.log',conf/'make.log',conf/'petscvariables',source/d['PETSC_ARCH']/'include/petscconf.h'):shutil.copyfile(p,evidence/p.name)
ldd=run(['ldd',prefix/'lib/libpetsc.so'],'gpu13_adoption_ldd',cuda=d['candidate_wrapper'],extra_env={'LD_LIBRARY_PATH':str(prefix/'lib')})
assert okay(ldd) and 'not found' not in text(ldd)
resolved=re.search(r'libmpi[^\s]*\s+=>\s+(\S+)',text(ldd));assert resolved and Path(resolved[1]).resolve()==(Path(d['MPI_prefix'])/'lib/libmpi.so').resolve()
d['steps'].append(ldd);d.update(status='PASS',selftest=check,selftest_observed=observed,source_unmodified=True,source_files_verified=len(integrity['files']),source_modifications=[],library_sha256=digest(prefix/'lib/libpetsc.so'),MPI_link=resolved[1],parser_correction='Self-test prints CUDA uppercase; case-insensitive observation. Original self-test exit0 and all3 lines retained; not rerun.')
write('petsc_gpu13_build',d);integrity.update(status='PASS',phase='after_build',modifications=[]);write('petsc_gpu13_configfix_source_integrity',integrity)
print('GPU13 original configure/make/install/CPU-MPI2-CUDA checks all PASS; provenance completed without rerunning.',flush=True)
