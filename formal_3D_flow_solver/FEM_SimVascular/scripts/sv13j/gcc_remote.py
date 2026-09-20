"""Official Ubuntu GCC12 only; preserve default compiler paths and versions."""
import json, re, shutil
from pathlib import Path
from runner_remote import *
before=load('pre_install_environment')
commands=[]
if not shutil.which('gcc-12') or not shutil.which('g++-12'):
    policy=run(['apt-cache','policy','gcc-12','g++-12'],'gcc12_package_policy')
    assert 'archive.ubuntu.com/ubuntu' in text(policy)
    args=['apt-get','--no-install-recommends','--no-upgrade','install','gcc-12','g++-12']
    simulation=run(args[:1]+['--simulate']+args[1:],'gcc12_install_simulation')
    assert okay(simulation)
    assert re.search(r'0 upgraded, \d+ newly installed, 0 to remove',text(simulation))
    changes=[line.split()[1] for line in text(simulation).splitlines() if line.startswith('Inst ')]
    allowed=re.compile(r'^(gcc-12(-base|-x86-64-linux-gnu)?|g\+\+-12(-x86-64-linux-gnu)?|cpp-12(-x86-64-linux-gnu)?|libgcc-12-dev|libstdc\+\+-12-dev|libasan8|libtsan2|libcc1-0|libitm1|libatomic1|liblsan0|libubsan1|libquadmath0|libgomp1)$')
    assert all(allowed.fullmatch(n) for n in changes),changes
    command=run(args[:1]+['--yes']+args[1:],'gcc12_install',timeout=600,extra_env={'DEBIAN_FRONTEND':'noninteractive'})
    commands=[simulation,command]
    assert okay(command),'GCC12_OFFICIAL_INSTALL_FAILED'
compilers={}
for name in ('gcc','g++','gcc-12','g++-12'):
    p=Path(shutil.which(name));v=run([p,'--version'],'version_'+name.replace('+','p'))
    major=run([p,'-dumpfullversion'],'version_full_'+name.replace('+','p'))
    compilers[name]={'path':str(p),'realpath':str(p.resolve()),'sha256':digest(p),'version':text(v),'number':text(major).strip()}
    if name.endswith('-12'):assert compilers[name]['number'].split('.')[0]=='12'
    else:assert digest(p)==before['default_compilers'][name]['sha256'] and str(p.resolve())==before['default_compilers'][name]['realpath']
c=B/'gcc12_smoke.c';cpp=B/'gcc12_smoke.cpp'
c.write_text('#include <stdio.h>\n#if __GNUC__ != 12\n#error expected GCC12\n#endif\nint main(void){int a=21;printf("C_OK gcc=%d.%d result=%d\\n",__GNUC__,__GNUC_MINOR__,2*a);return 2*a!=42;}\n')
cpp.write_text('#include <iostream>\n#include <tuple>\nstatic_assert(__cplusplus==201703L);\nstatic_assert(__GNUC__==12);\nint main(){auto [a,b]=std::tuple<int,int>{19,23};std::cout<<"CXX17_OK gcc="<<__GNUC__<<"."<<__GNUC_MINOR__<<" result="<<a+b<<"\\n";return a+b!=42;}\n')
smokes=[]
for compiler,source in (('gcc-12',c),('g++-12',cpp)):
    binary=source.with_suffix('')
    compile=run([compilers[compiler]['path'],str(source),'-std=c17' if compiler=='gcc-12' else '-std=c++17','-o',str(binary)],'smoke_compile_'+compiler.replace('+','p'))
    assert okay(compile)
    result=run([binary],'smoke_run_'+compiler.replace('+','p'));result['stdout']=text(result)
    assert okay(result) and 'result=42' in text(result)
    smokes.append({'compile':compile,'run':result,'binary_sha256':digest(binary)})
lib=run(['g++-12','-print-file-name=libstdc++.so'],'gcc12_libstdcpp')
library=Path(text(lib).strip())
write('gcc12_environment',{'status':'PASS','compilers':compilers,'install_commands':commands,'smokes':smokes,
    'libstdcpp_path':str(library),'libstdcpp_realpath':str(library.resolve()),'libstdcpp_sha256':digest(library),
    'default_compilers_unchanged':True,'unsupported_compiler_override':False,
    'support_policy':'NVIDIA CUDA 12.3 host compiler policy checks major version and explicitly supports newer minor versions of listed GCC majors; exact Ubuntu24.04 distro is not listed in CUDA12.3 qualified OS table.'})
print('Official GCC12 C / C++17 smoke PASS; default GCC13 unchanged',flush=True)
