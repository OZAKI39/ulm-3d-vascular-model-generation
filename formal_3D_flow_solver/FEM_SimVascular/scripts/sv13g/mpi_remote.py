"""Bounded five-branch system MPI diagnosis, executed natively on the server."""
import hashlib,json,os,shutil,signal,subprocess,time
from pathlib import Path
BASE=Path(__file__).resolve().parent
R=BASE/'reports';L=BASE/'logs';B=BASE/'benchmarks'
for p in (R,L,B):p.mkdir(exist_ok=True)
ENV={k:os.environ[k] for k in ('HOME','USER','LOGNAME','LANG') if k in os.environ}
ENV.update(PATH='/usr/bin:/bin:/usr/local/cuda/bin',LC_ALL='C',OMP_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='0')
def write(name,data):
    (R/(name+'.json')).write_text(json.dumps(data,indent=2)+'\n')
def run(name,args,env=None,limit=10):
    start=time.monotonic();timeout=False
    # A singleton MPI daemon may create its own session and inherit stdout.
    # Regular files keep timeout bounded even when a daemon retains an fd.
    with (L/(name+'.stdout')).open('w') as fout,(L/(name+'.stderr')).open('w') as ferr:
        p=subprocess.Popen(args,cwd=B,env=ENV if env is None else env,stdin=subprocess.DEVNULL,stdout=fout,stderr=ferr,text=True,start_new_session=True)
        try:p.wait(timeout=limit)
        except subprocess.TimeoutExpired:
            timeout=True;os.killpg(p.pid,signal.SIGKILL);p.wait(timeout=2)
    out=(L/(name+'.stdout')).read_text();err=(L/(name+'.stderr')).read_text()
    rec={'name':name,'command':args,'exit_code':p.returncode,'timeout':timeout,'timeout_s':limit,'wall_time_s':time.monotonic()-start,'stdout':out,'stderr':err}
    write(name,rec);print(json.dumps({k:rec[k] for k in ('name','exit_code','timeout','wall_time_s')}),flush=True)
    return rec
def okay(r):return r['exit_code']==0 and not r['timeout']
def fingerprint():
    cmds=[['hostname'],['uname','-a'],['cat','/etc/os-release'],['id'],['id','-u'],['whoami']]
    cmds += [[p,'--version'] for p in ('mpiexec','mpirun')]
    cmds += [[p,'--show'] for p in ('mpicc','mpicxx')]
    cmds += [[p,'--version'] for p in ('ompi_info','prte_info','pmix_info') if shutil.which(p,path=ENV['PATH'])]
    probes=[run('environment_'+str(i),a) for i,a in enumerate(cmds)]
    paths={p:shutil.which(p,path=os.environ.get('PATH')) for p in ('mpiexec','mpirun','mpicc','mpicxx','ompi_info','prte_info','pmix_info','strace')}
    selected={k:os.environ.get(k) for k in ('PATH','LD_LIBRARY_PATH','CUDA_VISIBLE_DEVICES','TMPDIR')}
    environment={'uid':os.getuid(),'root':os.getuid()==0,'paths':paths,'realpaths':{k:str(Path(v).resolve()) if v else None for k,v in paths.items()},'environment':selected,'clean_test_environment':ENV,'probes':probes}
    write('mpi_environment',environment)
    openmpi=any('OpenRTE' in p['stdout'] or 'Open MPI' in p['stdout'] for p in probes)
    assert openmpi,'Unexpected MPI implementation: inspect fingerprint before proceeding'
    assert all(str(Path(paths[p]).resolve()).startswith('/usr/bin/') for p in ('mpiexec','mpicc','mpicxx')),'MPI_STACK_MISMATCH'
    launcher=run('launcher_ldd',['ldd',paths['mpiexec']])
    compile=run('mpi_hello_compile',[paths['mpicc'],str(B/'mpi_hello.c'),'-o',str(B/'mpi_hello')],limit=30)
    assert okay(compile)
    binary=run('mpi_hello_ldd',['ldd',str(B/'mpi_hello')])
    write('mpi_stack_consistency',{'status':'PASS','family':'Open MPI','prefix':'/usr','paths':paths,'realpaths':environment['realpaths'],'wrapper_show':[p for p in probes if '--show' in p['command']],'launcher_ldd':launcher,'binary_ldd':binary,'binary_sha256':hashlib.sha256((B/'mpi_hello').read_bytes()).hexdigest(),'compile':compile})
    run('gpu_snapshot',['nvidia-smi'])
    run('cuda_version',['/usr/local/cuda/bin/nvcc','--version'])
    return environment
def main():
    e=fingerprint();tests=[]
    for name,args in [('direct_true',['/bin/true']),('direct_hostname',['hostname']),('direct_singleton',[str(B/'mpi_hello')])]:tests.append(run(name,args))
    launch=e['paths']['mpiexec']
    for name,args in [('default_true',['/bin/true']),('default_hostname',['hostname']),('default_hello',[str(B/'mpi_hello')])]:tests.append(run(name,[launch,'-n','1']+args))
    working=None;recipe=[launch]
    if all(okay(t) for t in tests[-3:]):working='SYSTEM_MPI_DEFAULT_PASS'
    elif e['root']:
        recipe=[launch,'--allow-run-as-root']
        for name,arg in [('root_true','/bin/true'),('root_hello',str(B/'mpi_hello'))]:tests.append(run(name,recipe+['-n','1',arg]))
        if all(okay(t) for t in tests[-2:]):working='MPI_ROOT_LAUNCH_RESTRICTION_RESOLVED'
    help=run('mpi_help',[launch]+(['--allow-run-as-root'] if e['root'] else [])+['--help','all'])
    if working is None:
        assert '--bind-to' in help['stdout']
        recipe+=[ '--bind-to','none'];t=run('no_binding',recipe+['-n','1',str(B/'mpi_hello')]);tests.append(t)
        if okay(t):working='MPI_BINDING_CONFIGURATION_RESOLVED'
    if working is None:
        components=run('mpi_components',['ompi_info','--parsable','--all'])
        assert all(s in components['stdout'] for s in ('mca:pml:ob1:','mca:btl:self:','mca:btl:tcp:'))
        recipe+=['--mca','pml','ob1','--mca','btl','self,tcp'];t=run('minimal_transport',recipe+['-n','1',str(B/'mpi_hello')]);tests.append(t)
        if okay(t):working='MPI_LOCAL_TRANSPORT_RESOLVED'
    if working is None:
        tmp=BASE/'tmp_mpi';tmp.mkdir(exist_ok=True);(tmp/'writable_probe').write_text('writable\n')
        diagnostics=[['ls','-ld','/tmp',str(tmp),'/dev/shm'],['df','-h','/tmp',str(tmp),'/dev/shm'],['cat','/proc/1/cgroup'],['mount'],['bash','-c','ulimit -a'],['ls','-la','/sys/fs/cgroup'],['ls','-la','/tmp'],['ls','-l','/proc/self/ns']]
        for i,args in enumerate(diagnostics):run('restricted_environment_'+str(i),args)
        tmpenv=dict(ENV,TMPDIR=str(tmp));t=run('tmpdir',recipe+['-n','1',str(B/'mpi_hello')],env=tmpenv);tests.append(t)
        if okay(t):working='MPI_TMPDIR_RESOLVED'
        if working is None and e['paths']['strace']:
            run('launcher_strace',[e['paths']['strace'],'-f','-tt','-o',str(L/'mpi_launcher_strace.log'),launch]+(['--allow-run-as-root'] if e['root'] else [])+['-n','1','/bin/true'])
    write('system_mpi_tests',{'status':'PASS' if working else 'FAIL','classification':working or 'SYSTEM_MPI_RUNTIME_UNUSABLE','tests':tests,'working_command':recipe if working else None,'branches_allowed':5,'fallback_needed':working is None})
    if working:
        repeats=[run('rank1_'+str(i),recipe+['-n','1',str(B/'mpi_hello')]) for i in range(1,6)]
        rank2=run('rank2',recipe+['-n','2',str(B/'mpi_hello')])
        write('mpi_hard_gate',{'status':'PASS' if all(okay(t) and t['stdout'].strip()=='rank=0 size=1' for t in repeats) and okay(rank2) and set(rank2['stdout'].splitlines())=={'rank=0 size=2','rank=1 size=2'} else 'FAIL','rank1':repeats,'rank2':rank2,'command':recipe})
if __name__=='__main__':main()
