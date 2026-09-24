"""Run the repository's N004 manufactured case with the baseline GPU backend."""
from pathlib import Path
import json,subprocess,shutil,xml.etree.ElementTree as ET,os,time,hashlib,sys,urllib.request
root=Path(sys.argv[1]);b=json.loads((root/'audit/remote_protection_manifest.json').read_text())['build'];s=Path(b['source'])
test=s/'tests/cases/stokes/manufactured_solution/P2P1'
out=root/'official_smoke_lfs';out.mkdir();(out/'N004').mkdir();(out/'mesh').mkdir();(out/'bforce').mkdir()
shutil.copytree(test/'mesh/N004',out/'mesh/N004');shutil.copytree(test/'bforce/N004',out/'bforce/N004')
commit=json.loads((root/'audit/remote_protection_manifest.json').read_text())['source_git']['commit']['stdout'].strip()
downloads=[]
for local in list((out/'mesh').rglob('*'))+list((out/'bforce').rglob('*')):
    if not local.is_file():continue
    pointer=local.read_bytes()
    if not pointer.startswith(b'version https://git-lfs.github.com/spec/v1'):continue
    info=pointer.decode().splitlines();digest=info[1].split(':',1)[1];size=int(info[2].split()[1])
    original=test/local.relative_to(out)
    url='https://media.githubusercontent.com/media/SimVascular/svMultiPhysics/'+commit+'/'+str(original.relative_to(s))
    with urllib.request.urlopen(url,timeout=60) as response:payload=response.read()
    assert len(payload)==size and hashlib.sha256(payload).hexdigest()==digest
    local.write_bytes(payload);downloads.append(dict(path=str(local.relative_to(out)),url=url,sha256=digest,bytes=size))
(out/'lfs_downloads.json').write_text(json.dumps(downloads,indent=2)+'\n')
original=test/'N004/solver.xml';shutil.copy2(original,out/'solver_original.xml')
tree=ET.parse(original);ls=tree.find('.//LS');ls.set('type','GMRES');la=ls.find('Linear_algebra');la.set('type','petsc');la.find('Preconditioner').text='petsc-jacobi'
tree.write(out/'N004/solver.xml',encoding='utf-8',xml_declaration=True)
opts=(root/'baseline_PETSC_OPTIONS.txt').read_text().strip()
env={k:os.environ[k] for k in ('HOME','USER','LOGNAME','LANG') if k in os.environ}
env.update(PATH='/usr/bin:/bin:/usr/local/cuda/bin',LC_ALL='C',OMP_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='0',LD_LIBRARY_PATH=b['runtime_library_path'],PETSC_OPTIONS=opts)
w=b['PETSc_build'];cmd=[w['candidate_wrapper'],w['launcher'],'-n','1',w['candidate_wrapper'],b['executable'],'solver.xml'];start=time.time()
with (out/'solver.log').open('w') as f:
    p=subprocess.run(cmd,cwd=out/'N004',env=env,stdout=f,stderr=subprocess.STDOUT,timeout=180)
record=dict(command=cmd,exit_code=p.returncode,elapsed_seconds=time.time()-start,original_case=str(test/'N004'),xml_changes=['LS GMRES','Linear_algebra petsc','petsc-jacobi'],PETSC_OPTIONS=opts,outputs=[str(x.relative_to(out)) for x in out.rglob('result*')],binary_sha256=hashlib.sha256(Path(b['executable']).read_bytes()).hexdigest())
(out/'execution.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record));print((out/'solver.log').read_text()[-3500:])
