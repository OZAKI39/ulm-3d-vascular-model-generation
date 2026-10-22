"""One-time, read-only-source migration audit. Never deletes or overwrites source data."""
from pathlib import Path
import ast, collections, hashlib, json, os, shutil, subprocess, zipfile

PROJECT = Path(__file__).resolve().parents[1]
BASE = PROJECT.parent
SNAP = BASE/'temp_storage/mouse_flow_github_sync_20261022'
AUDIT = PROJECT/'docs/audit'
CE = BASE/'computation_examples'
CASE_NAMES = {
    'mean-2p0-mmps-A-ROI-only-balanced-pressure-v1':'roi_only',
    'mean-2p0-mmps-A-best-feasible-balance-v1':'best_feasible',
    'mean-2p0-mmps-A-H0-pressure-v1':'network_h0',
    'mean-2p0-mmps':'mean_2mmps',
    'fem_reference_stage_q':'stage_q',
}

def sha(p):
    with Path(p).open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()

def save(name,obj):
    (AUDIT/name).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')

def inventory():
    included=json.loads((SNAP/'sync_20261022/INCLUDED.json').read_text())
    rows=[]
    for r in included:
        p=Path(r['source']); actual=sha(p) if p.is_file() else None
        rows.append(dict(**r,local_sha256=actual,local_matches_snapshot=actual==r['sha256']))
    save('local_vs_published.json',rows)
    for n in ['INCLUDED','EXCLUDED','LINKS','FROZEN_FIELDS_CHECK','CASE_INDEX','UPSTREAM_TEST_RESOURCE_REFERENCES']:
        shutil.copyfile(SNAP/f'sync_20261022/{n}.json',AUDIT/f'published_{n}.json')
    git=[]
    for rel in ['computation_examples/_shared/particle_workflow',
                'computation_examples/_shared/geometry_network',
                'computation_examples/_shared/active_network_development',
                'computation_examples/_shared/formal_3D_flow_solver/FEM_SimVascular',
                'vascular_printing','temp_storage/mouse_flow_github_sync_20261022']:
        p=BASE/rel
        def cmd(*a):
            r=subprocess.run(['git','-C',str(p),*a],text=True,capture_output=True)
            return r.stdout.strip() if r.returncode==0 else None
        top=cmd('rev-parse','--show-toplevel');head=cmd('rev-parse','HEAD')
        status=cmd('status','--porcelain=v1','--untracked-files=all') or ''
        filename='git_'+rel.replace('/','_')+'.txt';(AUDIT/filename).write_text(status+'\n')
        git.append(dict(path=str(p),toplevel=top,commit=head,status_file=filename,status_lines=len(status.splitlines())))
    save('git_baseline.json',git)
    print('LOCAL_BASELINE',len(rows),'differences',sum(not r['local_matches_snapshot'] for r in rows),flush=True)
    return rows

def destination(path):
    """Only scientific data/evidence here; current Python sources are handled separately."""
    p=Path(path); parts=p.parts
    if path.startswith('geometry_inputs/'):
        if parts[1] in ('outputs','vessel_model'):
            return Path('inputs/geometry')/Path(*parts[1:])
        return None
    if parts[0]!='computation_examples':return None
    if parts[1] in CASE_NAMES:
        case=CASE_NAMES[parts[1]]; tail=Path(*parts[2:])
        if parts[2]=='CASE_README.md':return Path('cases')/case/'ORIGINAL_CASE_README.md'
        if parts[2]=='particle_3d' and len(parts)>4 and parts[3]=='reports':
            return Path('results')/case/'particles'/Path(*parts[4:])
        return Path('results')/case/tail
    if parts[1]=='validation_cases':return Path('validation/recorded_cases')/Path(*parts[2:])
    if parts[1]!='_shared':return None
    area=parts[2];tail=Path(*parts[3:])
    if area=='particle_workflow':
        if str(tail).startswith('formal_3D_flow_solver/FEM_SimVascular/'):
            return Path('inputs/reference/fem')/Path(*tail.parts[2:])
        if str(tail).startswith('particle_3d/contracts/'):
            return Path('configs/particles/contracts')/Path(*tail.parts[2:])
        if str(tail).startswith('particle_3d/reports/'):
            return Path('inputs/reference/particles')/Path(*tail.parts[2:])
        if len(tail.parts)==1 and tail.suffix in ('.md','.json','.txt'):
            return Path('archive/workflow_records')/tail
        return None
    if area=='sonovue_size_distribution_v0':
        if tail.parts[0]!='src' and tail.parts[0]!='tests':return Path('inputs/reference/sonovue')/tail
        return None
    if area=='active_network_development' and tail.parts[0]=='ulm_3D_vascular':
        category=tail.parts[1]
        if category in ('outputs','reports'):
            return Path('inputs/network')/Path(*tail.parts[1:])
        if category=='configs':return Path('configs/geometry/original')/Path(*tail.parts[2:])
        return None
    if area=='geometry_network':
        # Identical/older inputs are recorded in the manifest, not a second code tree.
        return None
    if area=='formal_3D_flow_solver' and tail.parts[0]=='FEM_SimVascular':
        category=tail.parts[1]
        if category=='external':
            sub=Path(*tail.parts[3:])
            if any(s.lower() in ('tests','test','testcases','testdata','.github','docs') for s in sub.parts):return None
            return Path('external/svMultiPhysics/source')/sub
        if category in ('inputs','outputs'):
            return Path('inputs/fem')/Path(*tail.parts[1:])
        if category=='configs':return Path('configs/flow/original')/Path(*tail.parts[2:])
        if category in ('reports','logs'):return Path('archive/fem')/Path(*tail.parts[1:])
        if category=='solver_support' and tail.suffix in ('.h','.json','.sh','.md'):
            return Path('external/solver_support')/Path(*tail.parts[2:])
        return None
    if area in ('wss_audit','wss_validation_v2'):
        if tail.suffix=='.zip':return None
        return Path('validation')/area/tail
    if area=='fem_reference_contracts':return Path('archive/fem_reference_contracts')/tail
    return None

def migrate(rows):
    files=[]; skipped=[]; evidence=[]; used={}
    for r in rows:
        dest=destination(r['path']);source=Path(r['source'])
        if dest is None:
            skipped.append(dict(source=str(source),reason='out_of_selected_data_scope_or_source_handled_separately'));continue
        if not source.is_file():raise FileNotFoundError(source)
        # Original executable scripts are evidence only; no runtime imports from results.
        if source.suffix in ('.py','.sh','.ps1') and not str(dest).startswith('external/'):
            evidence.append((source,r['path']));continue
        if source.name.endswith('.zip'):
            skipped.append(dict(source=str(source),reason='review_bundle_unpacked_evidence_retained'));continue
        target=PROJECT/dest;actual=r['local_sha256']
        if dest in used and used[dest]!=actual:raise ValueError(('collision',dest))
        used[dest]=actual;target.parent.mkdir(parents=True,exist_ok=True)
        if not target.exists():shutil.copyfile(source,target)
        check=sha(target)
        if check!=actual:raise ValueError(('copy hash mismatch',source,target))
        files.append(dict(source=str(source),published_path=r['path'],path=str(dest),size=target.stat().st_size,sha256=check,status='PASS'))
    with zipfile.ZipFile(PROJECT/'archive/original_execution_scripts.zip','w',zipfile.ZIP_DEFLATED) as z:
        for src,name in evidence:z.write(src,name)
    save('data_manifest.json',files);save('selection_decisions.json',skipped)
    print('COPIED',len(files),sum(r['size'] for r in files),'bytes',flush=True)

if __name__=='__main__':
    AUDIT.mkdir(parents=True,exist_ok=True)
    (PROJECT/'archive').mkdir(exist_ok=True)
    rows=inventory();migrate(rows)
