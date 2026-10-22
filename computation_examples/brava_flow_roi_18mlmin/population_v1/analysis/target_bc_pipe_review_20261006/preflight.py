from common import *
import datetime
import subprocess

def main():
    old=read(IDENTITY/'PRODUCTION_SOURCE_IDENTITY.json')
    src=read(IDENTITY/'EXACT_DEPLOYED_SOURCE_MANIFEST.json')
    assert all(sha(RETAINED/n)==h for n,h in src.items())
    assert sha(old['actual_binary'])==old['actual_binary_sha256']
    assert sha(old['PETSc']['path'])==old['PETSc']['sha256']
    delivered=read(PREVIOUS/'OUTPUT_SHA256.json')
    assert all(sha(PREVIOUS/n)==h for n,h in delivered.items())
    inputs={str((PREVIOUS/n).relative_to(BASE)):h for n,h in delivered.items()}
    for n,h in read(PREVIOUS/'INPUT_SHA256.json').items():
        assert sha(BASE/n)==h
        inputs[n]=h
    put('INPUT_SHA256.json',inputs)
    resources={k:(Path('/sys/fs/cgroup')/k).read_text().strip() for k in
        ['cpu.max','cpuset.cpus.effective','memory.max','memory.current','memory.stat','memory.events']}
    commands={'GPU':['nvidia-smi','--query-gpu=name,driver_version,memory.total,memory.used,memory.free,utilization.gpu','--format=csv'],
        'GPU_processes':['nvidia-smi','--query-compute-apps=pid,process_name,used_memory','--format=csv'],
        'CPU_RAM':['free','-b'],'disk':['df','-h',str(BASE)],
        'other_processes':['ps','-eo','pid,ppid,pcpu,rss,comm','--sort=-pcpu']}
    for k,cmd in commands.items():resources[k]=subprocess.check_output(cmd,text=True)
    git=lambda *a:subprocess.check_output(['git','-c',f'safe.directory={BASE}','-C',str(BASE),*a],text=True).strip()
    resources.update(UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),project_root=git('rev-parse','--show-toplevel'),
        branch=git('branch','--show-current'),HEAD=git('rev-parse','HEAD'),preexisting_status=git('status','--short'))
    put('RESOURCES.json',resources)
    ident=read(PREVIOUS/'SOURCE_IDENTITY.json')
    ident.update(verified_now_UTC=resources['UTC'],previous_artifacts_verified=len(delivered),production_modified=False)
    put('SOURCE_IDENTITY.json',ident)
    cuda=read(PREVIOUS/'CUDA_LIBRARY_IDENTITY.json')
    # Preserve the previous identity structure and independently verify every library path/hash.
    checks=[]
    def walk(v):
        if isinstance(v,dict):
            if 'path' in v and 'sha256' in v:
                assert sha(v['path'])==v['sha256'];checks.append(v['path'])
            for x in v.values():walk(x)
        elif isinstance(v,list):
            for x in v:walk(x)
    walk(cuda)
    for path,digest in cuda.get('current_libraries_sha256',{}).items():
        assert sha(path)==digest
        checks.append(path)
    put('CUDA_LIBRARY_IDENTITY.json',dict(previous=cuda,verified_paths=checks,UTC=resources['UTC']))
    print('Production identity PASS:',len(src),'source files;',len(delivered),'prior artifacts;',len(inputs),'inputs')

if __name__=='__main__':main()
