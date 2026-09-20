"""Install only the authorized official Fortran12 compiler and necessary dependencies."""
import re, shutil
from runner_remote import *
before=load('pre_install_environment');commands=[]
probe={name:shutil.which(name) for name in ('gcc-12','g++-12','gfortran-12','gfortran','gfortran-13')}
write('compiler_probe',probe)
assert probe['gcc-12'] and probe['g++-12']
if not probe['gfortran-12']:
    policy=run(['apt-cache','policy','gfortran-12'],'gfortran12_package_policy')
    assert okay(policy) and 'archive.ubuntu.com/ubuntu' in text(policy)
    args=['apt-get','--no-install-recommends','--no-upgrade','install','gfortran-12']
    simulation=run(args[:1]+['--simulate']+args[1:],'gfortran12_simulation')
    assert okay(simulation) and re.search(r'0 upgraded, \d+ newly installed, 0 to remove',text(simulation))
    changes=[line.split()[1] for line in text(simulation).splitlines() if line.startswith('Inst ')]
    allowed={'gfortran-12','libgfortran-12-dev','libgfortran5','libquadmath0'}
    assert set(changes)<=allowed,changes
    install=run(args[:1]+['--yes']+args[1:],'gfortran12_install',timeout=600,extra_env={'DEBIAN_FRONTEND':'noninteractive'})
    assert okay(install)
    commands=[simulation,install]
compilers={};smokes=[]
for name in ('gcc-12','g++-12','gfortran-12'):
    path=Path(shutil.which(name));tag=name.replace('+','p')
    v=run([path,'--version'],tag+'_version');n=run([path,'-dumpfullversion'],tag+'_fullversion')
    assert text(n).strip().split('.')[0]=='12'
    compilers[name]={'path':str(path),'realpath':str(path.resolve()),'sha256':digest(path),'version':text(v),'number':text(n).strip()}
sources={
 'gcc-12':('compiler_c.c','#include <stdio.h>\n#if __GNUC__ != 12\n#error expected GCC12\n#endif\nint main(void){int a=21;printf("C result=%d\\n",2*a);return 2*a!=42;}\n',['-std=c17']),
 'g++-12':('compiler_cpp.cpp','#include <iostream>\n#include <tuple>\nstatic_assert(__cplusplus==201703L);\nstatic_assert(__GNUC__==12);\nint main(){auto [a,b]=std::tuple<int,int>{19,23};std::cout<<"CXX17 result="<<a+b<<"\\n";return a+b!=42;}\n',['-std=c++17']),
 'gfortran-12':('compiler_fortran.f90','program smoke\nimplicit none\ninteger :: a\na=21\nif (2*a /= 42) stop 1\nprint *, "FORTRAN result=42"\nend program\n',[])}
for name,(filename,source,flags) in sources.items():
    path=B/filename;path.write_text(source);binary=path.with_suffix('')
    compile=run([compilers[name]['path'],path,*flags,'-o',binary],binary.name+'_compile')
    assert okay(compile)
    result=run([binary],binary.name+'_run');result['stdout']=text(result)
    assert okay(result) and 'result=42' in text(result)
    smokes.append({'compiler':name,'compile':compile,'run':result,'binary_sha256':digest(binary)})
libs={}
for compiler,name in [('g++-12','libstdc++.so'),('gfortran-12','libgfortran.so'),('gfortran-12','libquadmath.so')]:
    p=run([compiler,'-print-file-name='+name],name.replace('.','_')+'_location')
    path=Path(text(p).strip());assert path.is_file()
    libs[name]={'path':str(path),'realpath':str(path.resolve()),'sha256':digest(path)}
for name,old in before['default_compilers'].items():
    path=Path(old['path']);assert str(path.resolve())==old['realpath'] and digest(path)==old['sha256']
write('compiler_environment',{'status':'PASS','compilers':compilers,'libraries':libs,'install_commands':commands,'smokes':smokes,'default_compilers_unchanged':True})
print('GCC12 / G++12 / GFortran12 native compiler smokes PASS',flush=True)
