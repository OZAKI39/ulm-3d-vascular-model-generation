"""Capture byte hashes, classify existing trajectories, and prepare remote inputs."""
import json,hashlib,subprocess,tarfile
from pathlib import Path
from collections import Counter

ROOT=Path(__file__).resolve().parents[1];REPO=ROOT.parents[1];PROJECTS=REPO.parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def dump(p,x):p.write_text(json.dumps(x,indent=2,ensure_ascii=False)+'\n')

def main():
    src=REPO/'particle_3d/src';old=REPO/'formal_3D_flow_solver/FEM_SimVascular'
    new=PROJECTS/'ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps'
    sets={'P81':REPO/'particle_3d/outputs/particle8_1/trajectories',
          'PPT_B':REPO/'particle_3d/outputs/particle8_2a_ppt/data/B/trajectories',
          'P82A':ROOT/'data/imported_P82A'}
    records=[];seen={};allhash={};counts={}
    for name,folder in sets.items():
        c=Counter()
        for f in sorted(folder.glob('*.npz')):
            meta=f.with_suffix('.json');m=json.loads(meta.read_text());h=sha(f)
            assert h==m['samples_sha256'],f
            key=h+':'+str(m['radius_m']);c[m['end_reason']]+=1
            r=dict(dataset=name,id=m['particle_id'],radius_m=m['radius_m'],sample_count=m.get('sample_count',0),
                   completed=m['completed'],outlet=m['exit_outlet'],end_reason=m['end_reason'],sha256=h,
                   metadata_sha256=sha(meta),local_path=str(f),remote_path=f'trajectories/{name}/{f.name}',
                   duplicate_of=seen.get(key))
            if key not in seen:seen[key]=f'{name}/{f.name}'
            records.append(r);allhash[str(f)]=h;allhash[str(meta)]=sha(meta)
        counts[name]=dict(c)
    for p in src.rglob('*.py'):allhash[str(p)]=sha(p)
    for p in (old/'frozen_reference').rglob('*'):
        if p.is_file():allhash[str(p)]=sha(p)
    dependencies=[]
    for line in (old/'frozen_reference/SHA256SUMS.txt').read_text().splitlines():
        expected,relative=line.split('  ',1)
        if not relative.startswith('frozen_reference/'):
            p=old/relative;assert sha(p)==expected
            dependencies.append((p,relative));allhash[str(p)]=expected
    for p in [new/'frozen_flow/steady_flow_mean_2p0_mmps.vtu',new/'run/solver.xml',new/'reports/physics_validation.json']:
        allhash[str(p)]=sha(p)
    son=PROJECTS/'sonovue_size_distribution_v0'
    for p in son.rglob('*'):
        if p.is_file():allhash[str(p)]=sha(p)
    dump(ROOT/'data/readonly_baseline.json',allhash)
    inv=dict(records=records,counts=counts,git_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
             unique_files=sum(r['duplicate_of'] is None for r in records),
             dataset_sha={k:hashlib.sha256('\n'.join(f"{r['id']}:{r['sha256']}:{r['metadata_sha256']}" for r in records if r['dataset']==k).encode()).hexdigest() for k in sets},
             exclusions='Historical guard variants, timestep/cache/worker benchmarks and point tracers are not independent formal finite-MB trajectories; no outputs from RBC synthetic tube included.')
    dump(ROOT/'data/trajectory_inventory.json',inv)
    cfg=dict(source=str(src),old_fem=str(old),new_flow=str(new/'frozen_flow/steady_flow_mean_2p0_mmps.vtu'),sonovue=str(son),
             inventory=str(ROOT/'data/trajectory_inventory.json'),output=str(ROOT/'data/local_smoke'),workers=1,
             max_states=300,flow_ids=['FLOW_0P353_MMPS','FLOW_2P0_MMPS'])
    dump(ROOT/'data/local_config.json',cfg)
    remote='/workspace/shear_lift_audit_20260923'
    rc=dict(cfg,source=remote+'/input/src',old_fem=remote+'/input/old_fem',new_flow=remote+'/input/new_flow.vtu',
            sonovue=remote+'/input/sonovue',inventory=remote+'/input/trajectory_inventory.json',output=remote+'/results',workers=7,max_states=None,input_base=remote+'/input')
    dump(ROOT/'data/remote_config.json',rc)
    with tarfile.open(ROOT/'data/remote_input.tar','w') as tar:
        for p in src.rglob('*.py'):tar.add(p,'input/src/'+str(p.relative_to(src)),recursive=False)
        tar.add(old/'frozen_reference','input/old_fem/frozen_reference')
        for p,relative in dependencies:tar.add(p,'input/old_fem/'+relative,recursive=False)
        tar.add(new/'frozen_flow/steady_flow_mean_2p0_mmps.vtu','input/new_flow.vtu')
        tar.add(new/'run/solver.xml','input/new_solver.xml')
        tar.add(new/'reports/physics_validation.json','input/new_physics_validation.json')
        tar.add(son,'input/sonovue')
        tar.add(ROOT/'data/trajectory_inventory.json','input/trajectory_inventory.json')
        for r in records:
            if r['duplicate_of'] is None:
                p=Path(r['local_path']);tar.add(p,'input/'+r['remote_path']);tar.add(p.with_suffix('.json'),'input/'+str(Path(r['remote_path']).with_suffix('.json')))
    print(json.dumps(dict(counts=counts,unique=inv['unique_files'],input_bytes=(ROOT/'data/remote_input.tar').stat().st_size)))

if __name__=='__main__':main()
