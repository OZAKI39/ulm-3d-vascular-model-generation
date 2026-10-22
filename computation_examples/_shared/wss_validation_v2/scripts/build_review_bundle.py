"""Final evidence inventory and compact reviewer ZIP; refuses an unfinished report."""
from pathlib import Path
import hashlib,json,csv,re,shutil,zipfile,time
V=Path(__file__).resolve().parents[1];C=V.parent
report=(V/'WSS_VALIDATION_REPORT_V2.md').read_text()
assert '执行中版本' not in report and '待实际求解完成后填入' not in report
for stage in [2,3,4]:assert (V/f'stage{stage}/completion.json').exists()
for n in range(1,6):assert list((V/'figures').glob(f'Figure_{n:02d}_*.png'))

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):h.update(chunk)
    return h.hexdigest()

protected=json.loads((V/'evidence/protected_inputs_before.json').read_text())
changes=[p for p,h in protected.items() if not Path(p).exists() or sha(p)!=h]
(V/'evidence/protected_inputs_after.json').write_text(json.dumps(dict(checked_unix=time.time(),checked_files=len(protected),changes=changes,all_unchanged=not changes),indent=2)+'\n')
assert not changes,'Protected original inputs changed; resolve before packaging'
src=V/'evidence/production_source';src.mkdir(exist_ok=True)
for path in [C/'solver_support/src/flow_solver_support/wss.py',C/'solver_support/src/flow_solver_support/wss_case.py',C/'solver_support/src/flow_solver_support/flow_parser.py',C/'rotate_visualization/prepare_surface_data.py']:
    shutil.copy2(path,src/path.name)

artifact_rows=[];case_rows=[]
excerpt_root=V/'evidence/solver_log_excerpts';excerpt_root.mkdir(exist_ok=True)
cases=sorted({p.parent for p in V.glob('stage[234]/**/policy.json') if (p.parent/'SV_MESH').is_dir()})
for case in cases:
    if not (case/'policy.json').exists():continue
    execution=case/'reports/execution.json';reuse=case/'reports/baseline_reuse.json';cancel=case/'reports/intentional_cancellation.json';backend_fail=case/'reports/backend_failure.json'
    d=json.loads(execution.read_text()) if execution.exists() else {}
    status=d.get('status','NOT_EXECUTED')
    if reuse.exists():status='REUSED_QUALIFIED_ORIGINAL_CFD'
    if cancel.exists():status='INTENTIONALLY_STOPPED_NOT_VALIDATION_RESULT'
    if backend_fail.exists():status='FAILED_GPU8_BACKEND_NO_VALID_SOLUTION'
    if (case/'reports/user_cancellation.json').exists():status='SKIPPED_BY_USER'
    case_rows.append(dict(case=case.name,case_relative_path=str(case.relative_to(V)),stage=case.relative_to(V).parts[0],status=status,final_step=d.get('final_step',''),elapsed_s=d.get('elapsed_s',''),source=str(case),reports=str(case/'reports')))
    if (case/'run/solver.log').exists():
        p=case/'run/solver.log';lines=p.read_text(errors='replace').splitlines()
        selected=set(range(min(80,len(lines))))|set(range(max(0,len(lines)-60),len(lines)))
        for i,line in enumerate(lines):
            if re.search(r'^\s*NS\s+\d+-|SV13Q_(BEGIN|END|RECOVERY)|Linear solve (converged|did not converge)|DIVERGED_|WARNING|ERROR|^-ksp_|^-pc_|^-sub_',line):selected.add(i)
        header=f'Case: {case.name}\nRaw log: {p}\nSHA256: {sha(p)}\nSource lines: {len(lines)}; selected lines below retain original line numbers.\n\n'
        (excerpt_root/(str(case.relative_to(V)).replace('/','__')+'.txt')).write_text(header+'\n'.join(f'{i+1}: {lines[i]}' for i in sorted(selected))+'\n')
    for rel in ['run/solver.xml','run/solver.log','run/PETSC_OPTIONS.txt','policy.json','input_hashes.json','initial_state/flow_guess.vtu','SV_MESH/mesh_arrays.npz','SV_MESH/mesh-complete.mesh.vtu','frozen_flow/flow_arrays_si.npz','frozen_flow/manifest.json','wss/data/wall_wss_si.vtp','reports/launch.json','reports/execution.json','reports/flow_quality.json']:
        path=case/rel
        if not path.exists():continue
        server_raw=rel.startswith(('run/','SV_MESH/','frozen_flow/','initial_state/')) or rel in ['policy.json','input_hashes.json','reports/launch.json','reports/execution.json']
        remote='vast4090:/workspace/wss_validation_v2_20260927T1230Z/'+str(path.relative_to(V)) if server_raw and (case/'reports/launch.json').exists() else ''
        artifact_rows.append(dict(case=case.name,path_relative_to_validation=str(path.relative_to(V)),local_path=str(path),server_copy_path=remote,bytes=path.stat().st_size,sha256=sha(path),case_status=status))
for filename,rows in [('ARTIFACT_INDEX.csv',artifact_rows),('CASE_STATUS.csv',case_rows)]:
    with (V/filename).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

files=[]
for name in ['WSS_VALIDATION_REPORT_V2.md','COMMANDS.md','PROGRESS.md','ARTIFACT_INDEX.csv','CASE_STATUS.csv','REVIEW_BUNDLE_README.md']:
    files.append(V/name)
for folder in ['data','figures','scripts','evidence']:
    for p in (V/folder).rglob('*'):
        if p.is_file() and '__pycache__' not in p.parts and p.suffix not in ['.pyc']:
            # NPZ is allowed only for the deliberately small actual-field audit subsets.
            if p.suffix=='.npz' and 'local_wall_snapshots' not in p.parts:continue
            if p.stat().st_size<=5*1024*1024:files.append(p)
for stage in [1,2,3,4]:
    root=V/f'stage{stage}'
    for p in root.glob('*'):
        if p.is_file() and p.suffix in ['.md','.json']:files.append(p)
    for sub in ['data','registration_versions']:
        if (root/sub).exists():files.extend(p for p in (root/sub).rglob('*') if p.is_file())
    for case in [c for c in cases if c.is_relative_to(root)]:
        for p in [case/'policy.json',case/'input_hashes.json',case/'mesh_request.json',case/'mesh_launch.json',case/'run/solver.xml',case/'run/PETSC_OPTIONS.txt',case/'frozen_flow/manifest.json',case/'wss/COMPUTE_VALIDATION.json']:
            if p.exists():files.append(p)
        if (case/'reports').exists():
            files.extend(p for p in (case/'reports').rglob('*') if p.is_file() and p.suffix in ['.json','.csv','.txt','.xml','.py'] and p.stat().st_size<2*1024*1024)
files=sorted(set(files));assert all(p.exists() for p in files)
manifest={str(p.relative_to(V)):dict(bytes=p.stat().st_size,sha256=sha(p)) for p in files}
mf=V/'REVIEW_BUNDLE_MANIFEST.json';mf.write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n');files.append(mf)
archive=V/'WSS_VALIDATION_V2_REVIEW_BUNDLE.zip'
with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=9) as z:
    for p in files:z.write(p,'wss_validation_v2/'+str(p.relative_to(V)))
(V/'WSS_VALIDATION_V2_REVIEW_BUNDLE.sha256').write_text(sha(archive)+'  '+archive.name+'\n')
print(json.dumps(dict(archive=str(archive),bytes=archive.stat().st_size,files=len(files),sha256=sha(archive),protected_files_unchanged=len(protected)),indent=2))
