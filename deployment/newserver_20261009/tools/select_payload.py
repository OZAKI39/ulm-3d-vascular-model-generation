"""Select the current CFD working set, keeping bulky retired evidence in local backup."""
import hashlib, json, os, re, shutil, stat, subprocess
from collections import Counter
from pathlib import Path
from datetime import datetime, timezone

ROOT=Path(__file__).resolve().parents[1]
BACKUP=Path('/home/lzy/projects/server_migration/brava_vast4090_20261009T094459Z')
LOCAL=Path('/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin')
B='/workspace/brava_flow_roi_18mlmin_20260928/'
S='/workspace/formal_3D_flow_solver_FEM_SimVascular_lzy/'
REMOTE='/workspace/brava_cfd'
STAGE=ROOT/'stage'
AUDITS={'speed_event_cause_review','field_reliability_event_review','stabilization_local_conservation_review_20261006','native_p2p1_review_20261006','common_curved_geometry_review_20261006','target_bc_pipe_review_20261006','p2_linear_algebra_review_20261006','p2_linear_algebra_rerun_20261006'}
SKIPPARTS={'.git','__pycache__','.pytest_cache','.mypy_cache','.ruff_cache','.rsync-partial','CMakeFiles'}

def keep_project(rel):
    parts=Path(rel).parts
    if any(p in SKIPPARTS for p in parts) or rel.endswith(('.pyc','.o','.obj','.mp4','.gif','.avi','.pid','.sock')):return False
    if len(parts)==1:return rel.endswith(('.md','.json')) and not rel.startswith(('GITHUB_','SYNC_','EXCLUDED_LARGE','SECOND_'))
    if parts[0] in {'inputs','mesh','code','scripts','gpu_solver_fix','vendor','reports'}:return True
    if parts[0]=='cases':return parts[1]=='balanced_pressure_final'
    if parts[0]=='production_workflow':
        if parts[1]=='runtime':return False
        if parts[1]=='manufactured_pipe':return len(parts)>2 and parts[2]=='R1_GMSH_V2'
        if parts[1]=='gmsh_reference' and len(parts)>2 and parts[2]=='candidate_0':return False
        return True
    if parts[:2]==('population_v1','analysis'):return len(parts)>2 and parts[2] in AUDITS
    if parts[:2]==('population_v1','source'):return True
    if parts[:2]==('population_v1','cfd_wall_audit'):return rel.endswith(('.md','.json','.png','.pdf','.py'))
    if parts[0]=='mesh_particle_validation_v2':
        if len(parts)>2 and parts[1] in {'02_mesh_candidates','03_cfd_runs','04_field_audit'}:
            return parts[2] in {'M1_fixedplc_medium','M1_fixedplc_fine'}
        return len(parts)==2 and rel.endswith(('.md','.json'))
    return False

STACK_MAP={
    'sv1_3l/external/gpu_mpi_fortran/':'deps/openmpi-4.1.6/',
    'sv1_3n/external/solver-inputs/vtk/':'deps/vtk-9.3/',
    'sv1_3n/external/solver-inputs/tetgen/':'deps/tetgen/',
    'sv1_3n/external/petsc325/source_gpu13_configfix/':'sources/petsc-3.25.5/',
    'sv1_3q/external/svMultiPhysics-reuse/':'sources/svMultiPhysics-reuse/',
}
source_rows={}
selected={}
for line in (BACKUP/'inventory/SOURCE_MANIFEST.jsonl').open():
    row=json.loads(line)
    if row['kind'] not in {'file','symlink'}:continue
    path=row['path'];dest=None
    if path.startswith(B):
        rel=path[len(B):]
        if keep_project(rel):dest='project/'+rel
    elif path.startswith(S):
        rel=path[len(S):]
        for prefix,d in STACK_MAP.items():
            if rel.startswith(prefix):
                sub=rel[len(prefix):]
                if any(p in SKIPPARTS or p.startswith('arch-sv13') for p in Path(sub).parts):break
                if sub in {'CTAGS','configure.log','configure.log.bkp','make.log','RDict.db','error.log'}:break
                dest=d+sub;break
        if rel in {'sv1_3q/reports/svmp_reuse_build.json','sv1_3q/configs/policy.json','sv1_3n/use_gpu13_env.sh','sv1_3l/scripts/run_gpu_mpi_fortran.sh'}:
            dest='provenance/old_stack/'+rel
    if dest:
        source_rows[dest]=row
        selected[dest]=dict(source=str(BACKUP/'payload'/path.lstrip('/')),origin='verified_server_backup',expected_sha256=row.get('sha256'),source_path=path)

# The user's current local tree is authoritative for tracked project code/reports.
tracked=subprocess.check_output(['git','-C',str(LOCAL),'ls-files','-z']).decode().split('\0')
for rel in tracked:
    if rel and keep_project(rel) and ((LOCAL/rel).is_file() or (LOCAL/rel).is_symlink()):
        selected['project/'+rel]=dict(source=str(LOCAL/rel),origin='current_local_git_worktree',source_path=str(LOCAL/rel))

# Retain the final saved native state, plus step 10 needed by existing restart regressions.
max_steps={}
for dest in selected:
    m=re.search(r'/(result_|stFile_)(\d+)\.(vtu|bin)$',dest)
    if m:max_steps[(str(Path(dest).parent),m.group(1))]=max(max_steps.get((str(Path(dest).parent),m.group(1)),0),int(m.group(2)))
for dest in list(selected):
    m=re.search(r'/(result_|stFile_)(\d+)\.(vtu|bin)$',dest)
    if m and int(m.group(2)) not in {10,max_steps[(str(Path(dest).parent),m.group(1))]}:del selected[dest]

STAGE.mkdir(exist_ok=True)
manifest=[];conflicts=[];counter=Counter();total=0
for dest,item in sorted(selected.items()):
    src=Path(item['source']);out=STAGE/dest;out.parent.mkdir(parents=True,exist_ok=True)
    rec=dict(destination=dest,**item)
    if src.is_symlink():
        target=os.readlink(src)
        if out.is_symlink():out.unlink()
        out.symlink_to(target)
        rec.update(kind='symlink',target=target)
    else:
        with src.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
        if item.get('expected_sha256'):assert digest==item['expected_sha256'],src
        old=source_rows.get(dest)
        if old and old.get('sha256') and digest!=old['sha256']:
            conflicts.append(dict(destination=dest,selection='current_local_git_worktree',local_sha256=digest,backup_sha256=old['sha256'],backup_path=old['path']))
        shutil.copy2(src,out)
        size=src.stat().st_size;total+=size;counter[dest.split('/')[0]]+=size
        rec.update(kind='file',size_bytes=size,sha256=digest,mode=stat.S_IMODE(src.stat().st_mode))
    manifest.append(rec)
for dest,row in source_rows.items():
    if dest in selected and row['kind']=='symlink':
        pass  # Link targets are audited separately before transfer.
(ROOT/'inventory/SELECTED_MANIFEST.jsonl').write_text(''.join(json.dumps(r)+'\n' for r in manifest))
(ROOT/'inventory/LOCAL_SERVER_DIFFERENCES.json').write_text(json.dumps(conflicts,indent=2)+'\n')
identity=dict(local_head=subprocess.check_output(['git','-C',str(LOCAL),'rev-parse','HEAD'],text=True).strip(),local_status=subprocess.check_output(['git','-C',str(LOCAL),'status','--short'],text=True),backup_root=str(BACKUP),source_manifest_sha256=hashlib.sha256((BACKUP/'inventory/SOURCE_MANIFEST.jsonl').read_bytes()).hexdigest())
(ROOT/'inventory/SOURCE_IDENTITY.json').write_text(json.dumps(identity,indent=2)+'\n')
summary=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),remote_root=REMOTE,files=sum(r['kind']=='file' for r in manifest),symlinks=sum(r['kind']=='symlink' for r in manifest),logical_bytes=total,bytes_by_category=dict(counter),local_server_differences=len(conflicts),selection='Current BraVa CFD code, formal field and BCs, fixed-PLC comparison fields, frozen Gmsh V2 reference, eight current CFD audit/implementation reviews, exact production solver/PETSc source and required MPI/VTK/TetGen dependencies',excluded=['Other projects and acoustic data','Old population state batches and animations','Abandoned mesh/refinement campaigns except current fixed-PLC evidence','Old duplicate CUDA/PETSc installations and generic Python environments','Intermediate repeated time snapshots except final and step10','Build objects and caches','Bulk Git history (local source commit and full selected-file hashes retained)'])
(ROOT/'inventory/SELECTION_SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2),flush=True)
