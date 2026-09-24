"""Isolated MPI4.1.6 with all Fortran bindings and native application checks."""
import re,shutil
from runner_remote import *
assert load('compiler_environment')['status']=='PASS'
src=Path(load('openmpi_source')['source_directory']);prefix=BASE/'external/gpu_mpi_fortran'
helptext=(L/'openmpi_configure_help.log').read_text()
assert 'or usempif08 (or all, build mpifh' in helptext
hwloc=src/'opal/mca/hwloc/hwloc201/configure.m4'
for n in ('gl','opencl','cuda','nvml'):assert 'enable_'+n+'=no' in hwloc.read_text()
shutil.copyfile(hwloc,R/'embedded_hwloc_configure.m4')
assert not prefix.exists()
compiler_env={'CC':'/usr/bin/gcc-12','CXX':'/usr/bin/g++-12','FC':'/usr/bin/gfortran-12'}
command=['./configure','--prefix='+str(prefix),'--enable-mpi-fortran=all','--disable-oshmem',
         '--with-hwloc=internal','--with-libevent=internal','--with-pmix=internal',
         '--without-cuda','--without-ucx','--without-ofi','--without-verbs']
d={'status':'BUILDING','version':'4.1.6','prefix':str(prefix),'source':load('openmpi_source'),
   'compiler_environment':compiler_env,'configure_command':command,'steps':[],'bindings_requested':'all'}
for args,name in [(command,'mpi_configure'),(['make','-j8'],'mpi_make'),(['make','install'],'mpi_install')]:
    result=run(args,name,cwd=src,timeout=1800,extra_env=compiler_env)
    d['steps'].append(result);write('mpi_fortran_build',d)
    if not okay(result):
        d.update(status='FAIL',reason='MPI_FORTRAN_DATATYPE_REPAIR_FAIL');write('mpi_fortran_build',d);raise SystemExit(1)
wrapper=BASE/'scripts/run_gpu_mpi_fortran.sh';wrapper.parent.mkdir(exist_ok=True)
wrapper.write_text('#!/usr/bin/env bash\nset -euo pipefail\n# Isolated MPI4.1.6 with Fortran bindings; no scientific options.\nexport PATH="'+str(prefix/'bin')+':${PATH:-/usr/bin:/bin}"\nexport LD_LIBRARY_PATH="'+str(prefix/'lib')+'${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"\nroot_args=()\nif [[ "$(id -u)" == 0 ]]; then root_args=(--allow-run-as-root); fi\nexec "'+str(prefix/'bin/mpiexec')+'" "${root_args[@]}" "$@"\n')
wrapper.chmod(0o755)
cuda=BASE.parent/'sv1_3j/external/compat_cuda/cuda-12.3.2'
cw=BASE/'scripts/use_cuda123_mpif_env.sh'
cw.write_text('#!/usr/bin/env bash\nset -euo pipefail\nexport CUDA_HOME="'+str(cuda)+'"\nexport PATH="'+str(cuda/'bin')+':'+str(prefix/'bin')+':/usr/bin:/bin"\nexport LD_LIBRARY_PATH="'+str(cuda/'lib64')+':'+str(prefix/'lib')+'${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"\nexport CC=/usr/bin/gcc-12\nexport CXX=/usr/bin/g++-12\nexport FC=/usr/bin/gfortran-12\nif [[ "$#" -gt 0 ]]; then exec "$@"; fi\n')
cw.chmod(0o755)
def execute(args,name,timeout=60):return run(args,name,timeout=timeout,cuda=cw)
def launch(binary,ranks,name,args=()):
    r=execute([wrapper,'-n',str(ranks),cw,binary,*args],name,timeout=20)
    r['stdout']=text(r);return r
info=execute([prefix/'bin/ompi_info','--parsable','--all'],'mpi_ompi_info')
wrappers={}
for n in ('mpicc','mpicxx','mpifort','mpif90'):
    p=prefix/'bin'/n
    if p.exists():
        show=execute([p,'--showme'],n+'_showme');assert okay(show)
        wrappers[n]={'path':str(p),'realpath':str(p.resolve()),'sha256':digest(p),'showme':text(show)}
assert 'mpifort' in wrappers and 'gfortran-12' in wrappers['mpifort']['showme']
linkage={}
for n in ('mpiexec','mpifort'):
    x=execute(['ldd',prefix/'bin'/n],n+'_ldd');linkage[n]=text(x)
    assert okay(x) and 'not found' not in text(x) and 'sv1_3g' not in text(x)
d.update(status='PASS',bindings_requested='all',wrappers=wrappers,ompi_info=info,
    mpi_library_sha256=digest(prefix/'lib/libmpi.so'),binary_sha256=digest(prefix/'bin/mpiexec'),
    linked_libraries=linkage,wrapper=str(wrapper),wrapper_sha256=digest(wrapper),cuda_wrapper=str(cw))
write('mpi_fortran_build',d)
for compiler,source,flags in [('mpicc','mpi_hello.c',[]),('mpifort','mpi_fortran_smoke.f90',[]),('mpicxx','mpi_datatype_probe.cpp',['-std=c++17'])]:
    x=execute([prefix/'bin'/compiler,B/source,*flags,'-o',B/Path(source).stem],Path(source).stem+'_compile')
    assert okay(x),'NATIVE_MPI_SMOKE_COMPILE_FAIL'
rank1=[launch(B/'mpi_hello',1,'mpi_c_rank1_'+str(i)) for i in range(1,6)]
rank2=launch(B/'mpi_hello',2,'mpi_c_rank2')
cpass=all(okay(x) and x['stdout'].strip()=='rank=0 size=1' for x in rank1) and okay(rank2) and set(rank2['stdout'].splitlines())=={'rank=0 size=2','rank=1 size=2'}
write('mpi_c_basic',{'status':'PASS' if cpass else 'FAIL','rank1':rank1,'rank2':rank2})
f_runs=[launch(B/'mpi_fortran_smoke',r,'mpi_fortran_rank'+str(r)) for r in (1,2)]
for x in f_runs:
    x['rows']=[{'rank':int(a),'ranks':int(b),'sizes':list(map(int,(c,d,e,f))),'failures':int(g)} for a,b,c,d,e,f,g in re.findall(r'FORTRAN rank=(\d+) ranks=(\d+) sizes=(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+failures=(\d+)',x['stdout'])]
fpass=all(okay(x) and len(x['rows'])==i+1 and all(r['failures']==0 and all(s>0 for s in r['sizes']) for r in x['rows']) for i,x in enumerate(f_runs))
write('mpi_fortran_basic',{'status':'PASS' if fpass else 'FAIL','runs':f_runs})
if not cpass or not fpass:
    write('mpi_application_gate',{'status':'FAIL','reason':'MPI_FORTRAN_DATATYPE_REPAIR_FAIL'});raise SystemExit(1)
sizes=f_runs[0]['rows'][0]['sizes'];assert all(row['sizes']==sizes for x in f_runs for row in x['rows'])
probes=[]
for repeat in range(1,4):
    runs=[]
    for ranks in (1,2):
        x=launch(B/'mpi_datatype_probe',ranks,f'mpi_datatypes_repeat{repeat}_rank{ranks}',list(map(str,sizes)))
        x['rows']=[dict(zip(('rank','ranks','name','size','expected_size','type_size_rc','bcast_rc','error_class','data_ok'),
            [int(v) if i!=2 else v for i,v in enumerate(m)])) for m in re.findall(r'DATATYPE rank=(\d+) ranks=(\d+) name=(\S+) size=(-?\d+) expected_size=(\d+) type_size_rc=(\d+) bcast_rc=(\d+) error_class=(\d+) data_ok=(\d+)',x['stdout'])]
        x['accepted']=okay(x) and len(x['rows'])==8*ranks and all(r['size']==r['expected_size'] and r['size']>0 and r['type_size_rc']==r['bcast_rc']==r['error_class']==0 and r['data_ok']==1 for r in x['rows'])
        runs.append(x)
    probes.append({'repeat':repeat,'runs':runs,'accepted':all(x['accepted'] for x in runs)})
passed=all(x['accepted'] for x in probes)
write('mpi_fortran_datatypes',{'status':'PASS' if passed else 'FAIL','native_Fortran_sizes':dict(zip(('MPI_INTEGER','MPI_DOUBLE_PRECISION','MPI_CHARACTER','MPI_LOGICAL'),sizes)),'repetitions':probes,'size_semantics':'C++ datatype sizes checked against independently measured Fortran storage_size; both ranks also check actual broadcast payloads.'})
write('mpi_application_gate',{'status':'PASS' if passed else 'FAIL','reason':None if passed else 'MPI_FORTRAN_DATATYPE_REPAIR_FAIL','C_rank1_passed':5,'C_rank2':'PASS','Fortran_rank1':'PASS','Fortran_rank2':'PASS','datatype_repeats_passed':sum(x['accepted'] for x in probes),'prefix':str(prefix),'working_launcher':str(wrapper),'cuda_wrapper':str(cw),'wrapper_sha256':digest(wrapper)})
if not passed:raise SystemExit(1)
print('MPI application compatibility PASS: C basic, native Fortran, C++ datatypes 3 complete rank1/rank2 repetitions',flush=True)
