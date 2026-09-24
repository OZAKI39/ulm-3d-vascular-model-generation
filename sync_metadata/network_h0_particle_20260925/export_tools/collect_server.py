"""Read-only export of retained server source and bounded execution evidence."""
from pathlib import Path
import collections,hashlib,json,os,subprocess,tarfile

DEST=Path('/workspace/archives/github_sync_network_h0_20260924T221559Z')
DEST.mkdir(parents=True,exist_ok=True)
SKIP={'.git','.venv','venv','audit_venv','.env','.cache','external','third_party','__pycache__','.pytest_cache','petsc-3.19.6','petsc_source','site-packages','tmp','frames','build','final_review'}
CODE={'.py','.sh','.cpp','.h','.c','.f90','.f','.toml','.cmake'}
TEXT={'.json','.jsonl','.yaml','.yml','.xml','.csv','.txt','.log','.out','.err','.md','.cfg','.ini','.lua','.patch'}
ARRAY={'.vtu','.vtp','.npz','.bin','.swc'}
MAX=25*1024**2
roots=[p for p in Path('/workspace').iterdir() if p.is_dir() and any(x in p.name for x in ('particle','flow_mean','formal_3D'))]
rows=[];files={};skipped=[]
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
for root in sorted(roots):
    latest=root.name in {'flow_mean_2p0_mmps_A_H0_20260924T181140Z','flow_mean_2p0_mmps_A_H0_TaylorHood_20260924T194503Z','particle_network_flow_mb_validation_v1_20260924T210047Z'}
    for directory,dirs,names in os.walk(root):
        for d in list(dirs):
            if d in SKIP or d.endswith('.egg-info'):
                skipped.append(str(Path(directory)/d));dirs.remove(d)
        for name in sorted(names):
            p=Path(directory)/name
            if not p.is_file() or p.is_symlink():continue
            rel=p.relative_to(root);size=p.stat().st_size;reason=None
            if name in {'.env','credentials','hosts.yml','id_rsa','id_ed25519'} or name.startswith('.env.') or p.suffix in {'.key','.pem','.p12'}:reason='CREDENTIAL_FILE'
            elif p.suffix in CODE or p.suffix in TEXT or name.endswith(('.json.gz','.log.gz','.jsonl.gz')):pass
            elif latest and p.suffix in ARRAY:
                if name.startswith(('result_','stFile_')) and not name.startswith(('result_070','result_071','stFile_071')):reason='INTERMEDIATE_SOLVER_STATE'
            else:reason='DUPLICATE_BINARY_RENDER_OR_DEPLOYMENT_ARCHIVE'
            if size>MAX:reason='NONESSENTIAL_FILE_OVER_25_MIB'
            target=Path('server_evidence')/root.name/rel
            row=dict(source=str(p),target=str(target),size=size,status='EXCLUDED' if reason else 'INCLUDED',reason=reason or 'RETAINED_SERVER_SOURCE_OR_EXECUTION_EVIDENCE')
            if not reason:row['sha256']=sha(p);files[str(target)]=p
            rows.append(row)

# Preserve the actual pinned solver changes and license, while installed third-
# party trees, PETSc/MPI/CUDA shared libraries and build products stay on-server.
source=Path('/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/svMultiPhysics-reuse')
head=subprocess.check_output(['git','-C',str(source),'rev-parse','HEAD'],text=True).strip()
patch=subprocess.check_output(['git','-C',str(source),'diff','--binary','HEAD'])
(DEST/'pinned_solver_changes.patch').write_bytes(patch)
status=subprocess.check_output(['git','-C',str(source),'status','--porcelain=v1','-z']).split(b'\0')
changed=[x[3:].decode() for x in status if x]
for name in changed+['LICENSE','LICENSE.txt']:
    p=source/name
    if not p.is_file() or p.is_symlink():continue
    target='server_evidence/pinned_solver_source/'+name
    rows.append(dict(source=str(p),target=target,size=p.stat().st_size,status='INCLUDED',reason='PINNED_SOLVER_PATCH_SOURCE',sha256=sha(p)));files[target]=p
for name in ['pinned_solver_changes.patch']:
    p=DEST/name;files['server_evidence/pinned_solver_source/'+name]=p
    rows.append(dict(source=str(p),target='server_evidence/pinned_solver_source/'+name,size=p.stat().st_size,status='INCLUDED',reason='PINNED_SOLVER_GIT_DIFF',sha256=sha(p)))
metadata=dict(source=str(source),git_head=head,changed_files=changed,license_source='vendor/svMultiPhysics_stage_q in the repository',executable='/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/sv1_3q/external/build/svmp_gpu_reuse/svMultiPhysics-build/bin/svmultiphysics')
(DEST/'pinned_solver_provenance.json').write_text(json.dumps(metadata,indent=2)+'\n')
p=DEST/'pinned_solver_provenance.json';target='server_evidence/pinned_solver_source/provenance.json';files[target]=p
rows.append(dict(source=str(p),target=target,size=p.stat().st_size,status='INCLUDED',reason='PINNED_SOLVER_PROVENANCE',sha256=sha(p)))
archive=DEST/'server_export.tar.gz'
assert not archive.exists()
with tarfile.open(archive,'w:gz',compresslevel=3) as tar:
    for row in rows:
        if row['status']!='INCLUDED':continue
        p=Path(row['source']);assert sha(p)==row['sha256'],p
        tar.add(p,arcname=row['target'],recursive=False)
manifest=dict(files=rows,skipped_directories=skipped,archive_sha256=sha(archive))
(DEST/'server_inventory.json').write_text(json.dumps(manifest,indent=2)+'\n')
print(json.dumps(dict(selected=sum(r['status']=='INCLUDED' for r in rows),selected_MiB=round(sum(r['size'] for r in rows if r['status']=='INCLUDED')/1024**2,2),archive_MiB=round(archive.stat().st_size/1024**2,2),archive_sha256=manifest['archive_sha256'],excluded_reasons=collections.Counter(r['reason'] for r in rows if r['status']=='EXCLUDED')),indent=2))
