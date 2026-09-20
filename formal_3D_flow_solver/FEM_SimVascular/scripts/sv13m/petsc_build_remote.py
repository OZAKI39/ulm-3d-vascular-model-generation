"""Fresh archive, narrowly allowlisted backport, identical Stage L configuration."""
import hashlib,shutil,tarfile,sys,re
from runner_remote import *
iteration=sys.argv[1];assert iteration in ('repair_01','repair_02','repair_03')
old=load('baseline_L_compatibility_winner');E=BASE/'external/compat_cuda';E.mkdir(exist_ok=True)
suffix='' if iteration=='repair_01' else '-'+iteration
source=E/('petsc-3.19.6-cuda123-ghostfix'+suffix);prefix=E/('petsc-cuda123-ghostfix'+suffix);arch='arch-sv13m-cuda123-'+iteration
archive=BASE.parent/'sv1_3/petsc-3.19.6.tar.gz';assert digest(archive)==old['source_archive_sha256']
assert not source.exists() and not prefix.exists()
original={}
with tarfile.open(archive) as t:
 for m in t:
  if m.isfile():original[str(Path(m.name).relative_to('petsc-3.19.6'))]=hashlib.sha256(t.extractfile(m).read()).hexdigest()
 temp=E/('extract_'+iteration);temp.mkdir();t.extractall(temp,filter='data')
(temp/'petsc-3.19.6').rename(source)
assert all(digest(source/n)==h for n,h in original.items())
patch=BASE/'patches'/iteration/'petsc319_cuda_ghost_backport.patch'
expected=json.loads((BASE/'configs'/('patch_scope_'+iteration+'.json')).read_text())['modified_files']
actual_files=re.findall(r'^\+\+\+ b/(.+)$',patch.read_text(),re.M);assert sorted(actual_files)==sorted(expected)
apply=run(['patch','--batch','--fuzz=0','-p1','-i',patch],iteration+'_patch_apply',cwd=source);assert okay(apply)
after={n:digest(source/n) for n in original};changed=[n for n in original if original[n]!=after[n]];assert sorted(changed)==sorted(expected)
hash_tree=lambda d:hashlib.sha256(json.dumps(d,sort_keys=True,separators=(',',':')).encode()).hexdigest()
write(iteration+'_source_integrity',dict(status='PASS',archive=str(archive),archive_sha256=digest(archive),source=str(source),fresh_source=True,previous_objects_reused=False,original_files=original,patched_files=after,before_tree_sha256=hash_tree(original),after_tree_sha256=hash_tree(after),modifications=changed,expected_changed_files=expected,patch_sha256=digest(patch),source_files=len(original)))
command=[('--prefix='+str(prefix)) if a.startswith('--prefix=') else ('PETSC_ARCH='+arch) if a.startswith('PETSC_ARCH=') else a for a in old['configure_command']]
diff=[dict(before=a,after=b) for a,b in zip(old['configure_command'],command) if a!=b];assert len(diff)==2
wrapper=old['candidate_wrapper'];extra={'OMPI_CC':'/usr/bin/gcc-12','OMPI_CXX':'/usr/bin/g++-12','PETSC_OPTIONS':'-skip_petscrc -use_gpu_aware_mpi 0','LD_LIBRARY_PATH':str(prefix/'lib')}
def execute(args,name,cwd=source,timeout=3600):
 r=run(args,iteration+'_'+name,cwd=cwd,timeout=timeout,cuda=wrapper,extra_env=extra);assert okay(r),name;return r
cfg=execute(command,'configure');make=execute(['make','-j8','V=1','PETSC_DIR='+str(source),'PETSC_ARCH='+arch,'all'],'make')
install=execute(['make','PETSC_DIR='+str(source),'PETSC_ARCH='+arch,'install'],'install')
check=execute(['make','PETSC_DIR='+str(prefix),'PETSC_ARCH=','check'],'selftest',timeout=600)
observed={k:s in text(check) for k,s in [('CPU','run successfully with 1 MPI process'),('MPI','run successfully with 2 MPI processes'),('CUDA','run successfully with cuda')]}
assert all(observed.values()) and 'Possible error' not in text(check) and 'Possible problem' not in text(check)
assert all(digest(source/n)==h for n,h in after.items()),'UNEXPECTED_SOURCE_CHANGE_DURING_BUILD'
ldd=execute(['ldd',prefix/'lib/libpetsc.so'],'ldd');assert 'not found' not in text(ldd)
resolved=re.search(r'libmpi[^\s]*\s+=>\s+(\S+)',text(ldd));assert resolved and Path(resolved[1]).resolve()==(Path(old['MPI_prefix'])/'lib/libmpi.so').resolve()
evidence=R/iteration;evidence.mkdir(exist_ok=True);conf=source/arch/'lib/petsc/conf'
for p in (conf/'configure.log',conf/'make.log',conf/'petscvariables',source/arch/'include/petscconf.h'):shutil.copyfile(p,evidence/p.name)
d=dict(old,status='PASS',version='3.19.6 + upstream CUDA ghost compatibility backport',stack_label='PETSc 3.19.6 + upstream CUDA ghost compatibility backport',source=str(source),prefix=str(prefix),PETSC_ARCH=arch,configure=cfg,make=make,install=install,selftest=check,selftest_observed=observed,configure_command=command,configuration_diff_from_L=diff,source_unmodified=False,source_expected_changes_only=True,source_modifications=changed,repair_iteration=iteration,library_sha256=digest(prefix/'lib/libpetsc.so'),MPI_link=resolved[1])
d.pop('configuration_diff_from_J',None)
write(iteration+'_petsc_build',d);write('compatibility_winner',d)
print(iteration+' CLEAN BUILD / CPU / MPI / CUDA SELFTEST PASS',flush=True)
