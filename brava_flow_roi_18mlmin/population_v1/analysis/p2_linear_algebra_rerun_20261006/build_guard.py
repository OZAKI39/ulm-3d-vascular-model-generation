from common import *
import subprocess
rt=HERE/'runtime';out=rt/'cuda_admission.so'
cmd=['gcc','-shared','-fPIC','-O2','-Wall','-Wextra',str(HERE/'cuda_admission.c'),'-ldl','-o',str(out)]
subprocess.run(cmd,check=True)
smoke=rt/'guard_smoke.c';smoke.write_text('#include <stddef.h>\n#include <stdlib.h>\nextern int cudaMalloc(void**,size_t);extern int cudaFree(void*);int main(int c,char**v){void*p=0;size_t n=c>1?((size_t)1<<42):4096;int e=cudaMalloc(&p,n);if(e)return e;return cudaFree(p);}\n')
subprocess.run(['gcc',str(smoke),'-L/usr/local/cuda/lib64','-Wl,-rpath,/usr/local/cuda/lib64','-lcudart','-o',str(rt/'guard_smoke')],check=True)
e=native_env();e.update(LD_PRELOAD=str(out),BRAVA_CUDA_ADMISSION_LOCK=str(rt/'guard.lock'),BRAVA_CUDA_ADMISSION_LOG=str(rt/'guard_smoke.jsonl'))
r=[]
for args,expected in [([],0),(['reject_without_allocation'],86)]:
 p=subprocess.run([str(rt/'guard_smoke'),*args],env=e,text=True,capture_output=True);r.append(dict(args=args,expected=expected,returncode=p.returncode,stderr=p.stderr,PASS=p.returncode==expected))
put('CUDA_ADMISSION_GUARD.json',dict(status='PASS' if all(x['PASS'] for x in r) else 'FAIL',source_sha256=sha(HERE/'cuda_admission.c'),binary_sha256=sha(out),compile_command=cmd,tests=r,coverage=['cudaMalloc','cudaMallocAsync','cudaMallocManaged'],limitations='Runtime entrypoints only; not a universal bound on driver/library allocations. CPU-only PC removes the previously failing CUDA SpGEMM product. Fresh launch also requires explicit resident budget.',reserve_MiB=2048))
assert all(x['PASS'] for x in r)
