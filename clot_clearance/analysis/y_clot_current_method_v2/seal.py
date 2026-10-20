"""Seal this application's files after processes and notebook updates finish."""
from common import *
import sys,platform,datetime,subprocess

def main():
    assert (O/'FINAL_STATUS.json').exists() and (O/'WSL_DOCUMENT_UPDATE.json').exists()
    frozen=json.loads((O/'SOURCE_FREEZE.json').read_text())
    for p,want in frozen['implementation_sha256'].items():assert sha(R/p)==want,p
    for p,want in frozen['input_sha256'].items():assert sha(O/p)==want,p
    for p,want in json.loads((O/'POSTPROCESS_SOURCE_FREEZE.json').read_text()).items():assert sha(R/p)==want,p
    import scipy,pyvista,vtk
    dump(O/'ENVIRONMENT.json',dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__,pyvista=pyvista.__version__,
        vtk=vtk.vtkVersion.GetVTKVersion(),platform=platform.platform(),native_library_sha256=sha(LIB),
        native_execution='Single MPI rank and one thread',environment='/workspace/y_clot/env.sh'))
    excluded={'DELIVERY_MANIFEST.json','LOCAL_DELIVERY_READBACK.json','SERVER_FINAL_SEAL_READBACK.json'};files={}
    for folder in [A,O,F]:
        for p in sorted(folder.rglob('*')):
            if p.name in excluded or '__pycache__' in p.parts or not p.is_file():continue
            q=dict(bytes=p.stat().st_size,sha256=sha(p))
            if p.is_symlink():q['symlink']=os.readlink(p)
            files[str(p.relative_to(R))]=q
    manifest=dict(task=TASK,run=RUN,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=files,
        file_count=len(files),total_bytes=sum(q['bytes'] for q in files.values()),
        excluded_attempt='verification/y_clot_current_method_v1/20261019T125103Z remains separate; exclusion record is included',
        native_library_external_dependency=str(LIB.relative_to(R)),native_library_sha256=sha(LIB),GitHub_pushed=False)
    dump(O/'DELIVERY_MANIFEST.json',manifest)
    bad=[p for p,z in files.items() if sha(R/p)!=z['sha256']];assert not bad
    dump(O/'SERVER_FINAL_SEAL_READBACK.json',dict(status='PASS',files=len(files),bytes=manifest['total_bytes'],mismatches=bad))
    print('Y_APPLICATION_SEALED',len(files),manifest['total_bytes'],flush=True)

if __name__=='__main__':main()
