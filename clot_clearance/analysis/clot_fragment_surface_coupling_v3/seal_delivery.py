"""Seal completed evidence after numerical, visual, document and commit records exist."""
from common import *
import datetime,argparse

def link(path,target):
    path=Path(path);target=Path(target);path.parent.mkdir(parents=True,exist_ok=True)
    if path.is_symlink():assert path.resolve()==target.resolve();return
    assert not path.exists(),path
    path.symlink_to(os.path.relpath(target,path.parent))

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--prepare',action='store_true');args=parser.parse_args()
    for name in ['ARTIFACT_VALIDATION.json','ARTIFACT_VISUAL_REVIEW.json','WSL_DOCUMENT_UPDATE.json']:
        assert json.loads((O/name).read_text())['status']=='PASS',name
    if not args.prepare:
        commit=json.loads((O/'COMMIT_RECORD.json').read_text());assert len(commit['commit'])==40
    status=json.loads((O/'FINAL_STATUS.json').read_text())
    if not status['final_delivery_complete']:
        assert args.prepare
        status.update(final_delivery_complete=True);dump(O/'FINAL_STATUS.json',status)
    core=['SOURCE_FREEZE.json','TARGET_COMPONENT_SURFACE_DIAGNOSIS.json','P0_TARGET_SURFACE_DIAGNOSIS.json','SURFACE_QUADRATURE_AUDIT.csv',
          'FRAGMENT_EXPOSURE_HISTORY.csv','FRAGMENT_FORCE_AND_POWER.csv','PHYSICAL_DOMAIN_SCALE_ASSESSMENT.md','DOMAIN_SCALE_PROVENANCE.json',
          'NATIVE_V3_RESULT.json','FRAGMENT_SIZE_HISTORY.csv','FRAGMENT_LINEAGE.json','FRAGMENT_MORPHOLOGY.csv','ENERGY_LEDGER.csv',
          'COMPUTE_COST.json','FINAL_STATUS.json','REVIEW_CLOT_FRAGMENT_SURFACE_COUPLING_V3.md','REPRODUCE.md','COMMIT_RECORD.json','ARTIFACT_VALIDATION.json']
    for name in core:link(A/name,O/name)
    for name,target in [('OPEN_RESULTS.html',F/'OPEN_RESULTS.html'),('figures',F/'figures'),('videos',F/'videos'),('states',O/'cases/biphasic/keys')]:link(A/name,target)
    for p in sorted(A.iterdir()):
        if p.is_file() and not p.is_symlink() and p.suffix in ('.py','.cpp','.sh'):link(A/'scripts'/p.name,p)
    sources={str(p.relative_to(R)):sha(p) for p in sorted(A.iterdir()) if p.is_file() and not p.is_symlink()}
    if args.prepare:
        dump(O/'POSTPROCESS_SOURCE_FREEZE.json',dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),sha256=sources))
        print('DELIVERY_PREPARED_FOR_GIT_COMMIT',flush=True);return
    for p,h in json.loads((O/'POSTPROCESS_SOURCE_FREEZE.json').read_text())['sha256'].items():assert sha(R/p)==h,p
    exclusions={'DELIVERY_MANIFEST.json','LOCAL_DELIVERY_READBACK.json'};files={};links={}
    for directory in [A,O,F]:
        for p in sorted(directory.rglob('*')):
            if '__pycache__' in p.parts or p.name in exclusions:continue
            key=str(p.relative_to(R))
            if p.is_symlink():links[key]=os.readlink(p)
            elif p.is_file():files[key]=dict(bytes=p.stat().st_size,sha256=sha(p))
    dump(O/'DELIVERY_MANIFEST.json',dict(task=TASK,run=RUN,utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        files=files,symlinks=links,file_count=len(files),total_bytes=sum(v['bytes'] for v in files.values()),
        intentionally_excluded=sorted(exclusions),native_dependencies='SOURCE_FREEZE.json and V2_DEPENDENCY_FREEZE.json; historical V2 archive remains separately sealed'))
    print('DELIVERY_SEALED',len(files),sum(v['bytes'] for v in files.values()),flush=True)

if __name__=='__main__':main()
