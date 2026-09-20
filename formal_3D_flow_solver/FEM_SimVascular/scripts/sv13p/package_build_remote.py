"""Isolated, unmodified PETSc and official GPU package build; no CFD or CPU runs."""
import sys,tarfile,hashlib,shutil,re,urllib.request
from runner_remote import *
kind=sys.argv[1];assert kind in ('hypre','amgx')
ref=json.loads((BASE/'configs/reference_freeze.json').read_text());old=ref['PETSc_build']
E=BASE/('external/petsc325_cuda13_'+kind);E.mkdir(parents=True,exist_ok=False)
archive=BASE.parent/'sv1_3n/external/petsc325/petsc-3.25.5.tar.gz';assert digest(archive)==old['source_archive_sha256']
source=E/'source';prefix=E/'install';arch='arch-sv13p-'+kind
original={}
with tarfile.open(archive) as t:
 for m in t:
  if m.isfile():original[str(Path(m.name).relative_to('petsc-3.25.5'))]=hashlib.sha256(t.extractfile(m).read()).hexdigest()
 temp=E/'extract';temp.mkdir();t.extractall(temp,filter='data')
(temp/'petsc-3.25.5').rename(source)
assert all(digest(source/n)==h for n,h in original.items())
package=source/('config/BuildSystem/config/packages/'+kind+'.py')
package_url='https://github.com/hypre-space/hypre/archive/85b779557005b2eb94c231c1b516e988b87f4e53.tar.gz' if kind=='hypre' else 'https://web.cels.anl.gov/projects/petsc/download/externalpackages/amgx-2.4.0.tar.gz'
assert ('85b779557005b2eb94c231c1b516e988b87f4e53' if kind=='hypre' else package_url) in package.read_text()
package_archive=E/(kind+'.tar.gz')
# Fetch the precise archive supported by the frozen PETSc package mechanism.
urllib.request.urlretrieve(package_url,package_archive)
cmd=[a for a in old['configure_command'] if not a.startswith(('--prefix=','PETSC_ARCH='))]+['--prefix='+str(prefix),'PETSC_ARCH='+arch,'--download-'+kind+'='+str(package_archive),'--with-make-np=8']
if kind=='hypre':cmd+=['--download-hypre-configure-arguments=--disable-gpu-aware-mpi --disable-unified-memory --enable-cusparse','--download-hypre-openmp=0']
d=dict(old,status='BUILDING',source=str(source),prefix=str(prefix),PETSC_ARCH=arch,configure_command=cmd,steps=[],package=kind,package_source_URL=package_url,package_archive_sha256=digest(package_archive),package_definition_sha256=digest(package),package_version='3.1.0' if kind=='hypre' else '2.4.0',package_commit='85b779557005b2eb94c231c1b516e988b87f4e53' if kind=='hypre' else None)
for k in ('selftest','selftest_observed','library_sha256','MPI_link'):d.pop(k,None)
extra={'LD_LIBRARY_PATH':str(prefix/'lib'),'PETSC_OPTIONS':'-skip_petscrc -use_gpu_aware_mpi 0'}
def execute(command,name,timeout=5400):
 r=run(command,kind+'_'+name,cwd=source,timeout=timeout,cuda=old['candidate_wrapper'],extra_env=extra);d['steps'].append(r);write('petsc_'+kind+'_build',d)
 if not okay(r):
  d.update(status='FAIL',failed_step=name);write('petsc_'+kind+'_build',d)
  evidence=R/(kind+'_build_evidence');evidence.mkdir(exist_ok=True)
  for p in (source/arch/'lib/petsc/conf/configure.log',source/'configure.log',source/arch/'lib/petsc/conf/make.log'):
   if p.is_file():shutil.copyfile(p,evidence/p.name)
  raise SystemExit(1)
 return r
execute(['/usr/bin/python3','-B','configure','--help'],'configure_help',60)
execute(cmd,'configure');execute(['make','-j8','V=1','PETSC_DIR='+str(source),'PETSC_ARCH='+arch,'all'],'make');execute(['make','PETSC_DIR='+str(source),'PETSC_ARCH='+arch,'install'],'install')
changed=[n for n,h in original.items() if not (source/n).is_file() or digest(source/n)!=h];assert not changed
conf=source/arch/'lib/petsc/conf';evidence=R/(kind+'_build_evidence');evidence.mkdir(exist_ok=True)
for p in (conf/'configure.log',conf/'make.log',conf/'petscvariables',source/arch/'include/petscconf.h'):
 if p.exists():shutil.copyfile(p,evidence/p.name)
if kind=='hypre':
 h=prefix/'include/HYPRE_config.h';shutil.copyfile(h,evidence/h.name);contents=h.read_text()
 d['hypre_features']={f:bool(re.search(r'^#define '+f+r'\b',contents,re.M)) for f in ['HYPRE_USING_CUDA','HYPRE_USING_GPU','HYPRE_USING_CUSPARSE','HYPRE_USING_UNIFIED_MEMORY','HYPRE_BIGINT','HYPRE_MIXEDINT','HYPRE_USING_GPU_AWARE_MPI']}
 assert all(d['hypre_features'][f] for f in ['HYPRE_USING_CUDA','HYPRE_USING_GPU','HYPRE_USING_CUSPARSE'])
 assert not any(d['hypre_features'][f] for f in ['HYPRE_USING_UNIFIED_MEMORY','HYPRE_BIGINT','HYPRE_MIXEDINT'])
 for file in ['solvers-ilu.rst','solvers-boomeramg.rst']:
  matches=list((source/arch/'externalpackages').glob('**/'+file));assert matches;shutil.copyfile(matches[0],evidence/file)
ldd=execute(['ldd',prefix/'lib/libpetsc.so'],'ldd',30);assert 'not found' not in text(ldd)
resolved=re.search(r'libmpi[^\s]*\s+=>\s+(\S+)',text(ldd));assert resolved and Path(resolved[1]).resolve()==(Path(old['MPI_prefix'])/'lib/libmpi.so').resolve()
d.update(status='PASS',source_unmodified=True,source_files_verified=len(original),source_modifications=changed,library_sha256=digest(prefix/'lib/libpetsc.so'),MPI_link=resolved[1],selftest='PENDING_SINGLE_GPU_STANDALONE')
write('petsc_'+kind+'_build',d);write(kind+'_source_integrity',dict(status='PASS',files=original,modifications=changed,source=str(source)))
print(kind+' isolated PETSc build PASS',flush=True)
