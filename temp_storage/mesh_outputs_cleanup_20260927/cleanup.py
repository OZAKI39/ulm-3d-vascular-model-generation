"""Delete only the reviewed obsolete output paths; preserve active mesh data."""
from pathlib import Path
import hashlib
import json
import os
import shutil

C=Path('/home/lzy/projects/formal_3D_flow_solver/FEM_SimVascular')
O=C/'outputs/mesh_and_flow'
F=Path('/home/lzy/projects/ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular/flow_cases/mean-2p0-mmps-A-H0-pressure-v1')
R=Path(__file__).resolve().parent

def sha(p):
    with p.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()

def inventory(root):
    rows=[]
    for p in sorted(root.rglob('*')):
        assert not p.is_symlink(), 'Unexpected symlink: '+str(p)
        if p.is_file():
            s=p.stat()
            rows.append(dict(path=str(p.relative_to(root)),sha256=sha(p),bytes=s.st_size,
                             device=s.st_dev,inode=s.st_ino,hardlinks=s.st_nlink,allocated_bytes=s.st_blocks*512))
    return rows

assert O.is_dir() and not O.is_symlink()
assert not (R/'before.json').exists(), 'Cleanup already started'
remove={
    'vascular_flow':'Unaccepted old transient field: FLOW_SOLVE_FAIL, MASS_BALANCE_FAIL, steady NOT_REACHED',
    'official_fluid_smoke':'Historical official example smoke-test outputs; not the current vascular mesh or field',
    'plot_cache':'Regenerable historical matplotlib cache',
    'api_documentation.json':'Regenerable API documentation snapshot; no active reader',
}
keep={'model','mesh_generation','solver_mesh'}
assert {p.name for p in O.iterdir()}==set(remove)|keep
status=json.loads((C/'reports/mesh_and_flow/final_status.json').read_text())
assert status['status']=='FAIL' and status['steady']=='NOT_REACHED'
assert 'MASS_BALANCE_FAIL' in status['reasons']
assert json.loads((O/'mesh_generation/primary/generation.json').read_text())['status']=='PASS'
mesh_matches={}
for p in sorted((O/'solver_mesh').rglob('*')):
    if p.is_file():
        rel=p.relative_to(O/'solver_mesh');h=sha(p)
        assert h==sha(F/'SV_MESH'/rel),str(rel)
        mesh_matches[str(rel)]=h
field=F/'frozen_flow/steady_flow_mean_2p0_mmps_A_H0.vtu'
field_sha=sha(field)
assert field_sha=='064cbd28f3efa72f426fc946b2f29da21f056c596609095e7283d39070aa55f4'
before=inventory(O)
deleted=[d for d in before if Path(d['path']).parts[0] in remove]
retained=[d for d in before if Path(d['path']).parts[0] in keep]
delete_inodes={}
for d in deleted:delete_inodes.setdefault((d['device'],d['inode']),[]).append(d)
# Every hard link to the reviewed obsolete data is inside the deletion set.
assert all(len(items)==items[0]['hardlinks'] for items in delete_inodes.values())
record=dict(scope=str(O),before=before,removed=deleted,retained=retained,reasons=remove,
            current_H0_mesh_matches=mesh_matches,current_field_sha256=field_sha)
(R/'before.json').write_text(json.dumps(record,indent=2)+'\n')
for name in remove:
    p=O/name
    assert p.parent==O and not p.is_symlink()
    if p.is_dir():shutil.rmtree(p)
    else:p.unlink()
after=inventory(O)
assert {d['path']:d['sha256'] for d in after}=={d['path']:d['sha256'] for d in retained}
assert sha(field)==field_sha
result=dict(all_pass=True,scope=str(O),deleted_entries=list(remove),
            deleted_file_paths=len(deleted),deleted_unique_files=len(delete_inodes),
            reclaimed_allocated_file_bytes=sum(items[0]['allocated_bytes'] for items in delete_inodes.values()),
            remaining_data_bytes=sum(d['bytes'] for d in retained),remaining_data_files=len(retained),
            retained_directories=sorted(keep),retained_hashes_unchanged=True,
            identical_to_current_H0_mesh_files=len(mesh_matches),current_field_sha256=field_sha,
            deleted_data_backed_up=False,server_changed=False)
(R/'result.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
