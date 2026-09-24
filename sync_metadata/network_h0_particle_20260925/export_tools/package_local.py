"""Publish current working files into an isolated checkout; sources stay read-only."""
from pathlib import Path
import collections,csv,hashlib,json,os,shutil,subprocess

A=Path(__file__).resolve().parents[1]
C=json.loads((A/'context.json').read_text());D=Path(C['destination']);W=Path('/home/lzy/projects')
P=W/'ulm_particle_3d_particle0';F=W/'ulm_flow_mean_2p0_mmps/formal_3D_flow_solver/FEM_SimVascular'
O=W/'formal_3D_flow_solver/FEM_SimVascular';V=W/'ulm_3D_vascular';SV='formal_3D_flow_solver/FEM_SimVascular'
M=D/'sync_metadata'/C['snapshot'];M.mkdir(parents=True,exist_ok=True)
SKIP={'.git','.venv','venv','__pycache__','.pytest_cache','.ruff_cache','.mypy_cache','node_modules','jit_cache','.codex_tmp','build','dist','.cache','.env','external','third_party','Ultraliser','external_reference','references','tmp','rotate_visualization_history'}
TEXT={'.json','.jsonl','.yaml','.yml','.xml','.csv','.txt','.log','.out','.err','.md','.html','.sh','.py','.toml','.lua','.f90','.cpp','.h','.c','.cfg','.ini','.mod'}
MAX=25*1024**2
rows=[];selected={};skipped=[]
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def policy(group,rel,p):
    parts=set(Path(rel).parts);n=p.name;size=p.stat().st_size
    if n=='.git':return 'GIT_INTERNAL'
    if n in {'.env','hosts.yml','credentials','credentials.json','id_rsa','id_ed25519'} or n.startswith('.env.') or p.suffix in {'.pem','.key','.p12'}:return 'CREDENTIAL_FILE'
    if group in {'sonovue','required_input'}:return None
    if n.endswith(('.tar','.tar.gz','.tgz','.zip','.7z')):return 'REDUNDANT_ARCHIVE'
    if n.endswith(('.pyc','.pyo','.so','.o','.a','.exe')):return 'REBUILDABLE_BINARY_OR_CACHE'
    if 'frames' in parts:return 'REGENERABLE_RENDER_FRAMES'
    if 'field_diagnostics_revisions' in parts:return 'SUPERSEDED_VISUALIZATION'
    if group=='particle':
        if rel.startswith('outputs/particle8_2a/admission/'):return 'REGENERABLE_ADMISSION_CACHE'
        if '/development_preview/' in '/'+rel:return 'SUPERSEDED_PREVIEW'
        if n in {'trajectory_samples.csv.gz','11_long_lifecycle.csv.gz','C57BL6_RBC_GEOMETRY_VALIDATION_100000.csv'}:return 'DUPLICATE_AGGREGATED_TRACE'
        if 'shear_lift_audit' in parts:
            if '/data/remote_results/states/' in '/'+rel or n=='all_scalar_states.npz' or n.endswith('_eulerian.npz'):return 'REGENERABLE_PER_STATE_AUDIT_CACHE'
            if '/literature/' in '/'+rel and p.suffix in {'.pdf','.html','.txt'}:return 'THIRD_PARTY_LITERATURE'
            if '/data/local_smoke/' in '/'+rel or '/data/local_nearwall_smoke/' in '/'+rel:return 'DUPLICATE_SMOKE'
        if n.endswith('_frames.json') and size>2*1024**2:return 'REPLAY_FRAME_CACHE'
    if group=='flow':
        if '/run/1-procs/' in '/'+rel and (n.startswith('result_') or n.startswith('stFile_')):
            if not n.startswith(('result_070','result_071','stFile_071')):return 'INTERMEDIATE_CFD_STATE'
        if '/streamlines/presentation' in '/'+rel or '/field_diagnostics/animations/' in '/'+rel:return 'OLDER_VISUALIZATION_CURRENT_ROTATE_RETAINED'
        if n=='candidates.npz':return 'STREAMLINE_CANDIDATE_CACHE'
    if group=='old_fem':
        if 'remote_return' in parts or 'previous_local_versions' in parts:return 'DUPLICATE_RETURN'
        if rel.startswith('outputs/') and p.suffix not in TEXT and not rel.startswith('outputs/stage01_7/selected/mesh/'):return 'HISTORICAL_INTERMEDIATE_BINARY'
    if group=='sv_history' and rel.startswith('outputs/') and (p.suffix not in TEXT or size>2*1024**2):return 'HISTORICAL_INTERMEDIATE_BINARY'
    if group=='vascular':
        if rel.startswith('outputs/topbrain_brava_transfer/'):return 'UNRELATED_CONCURRENT_TOPBRAIN_WORK'
        if rel.startswith('outputs/') and p.suffix not in TEXT and not n.endswith('.json.gz'):return 'BULK_UPSTREAM_DATA_REQUIRED_CLOSURE_SELECTED_SEPARATELY'
    if size>MAX:return 'NONESSENTIAL_FILE_OVER_25_MIB'
    return None

def add(group,p,target,rule_rel=None,only_missing=False):
    target=str(target);rel=rule_rel or target
    if p.is_symlink():rows.append(dict(group=group,source=str(p),target=target,status='EXCLUDED',reason='SYMLINK',size=0));return
    reason=policy(group,rel,p)
    row=dict(group=group,source=str(p),target=target,size=p.stat().st_size,status='EXCLUDED' if reason else 'INCLUDED',reason=reason or 'CURRENT_CODE_OR_EVIDENCE')
    if only_missing and ((D/target).exists() or target in selected):row.update(status='ALREADY_PACKAGED',reason='CANONICAL_WORKTREE_TAKES_PRECEDENCE')
    elif not reason:
        if target in selected:
            if sha(p)!=sha(Path(selected[target]['source'])):raise RuntimeError('Conflicting current sources: '+target)
            row.update(status='ALREADY_PACKAGED',reason='IDENTICAL_DEPENDENCY')
        else:selected[target]=row
    rows.append(row)

def walk(group,src,prefix,only_missing=False,extra_skip=()):
    if not src.exists():return
    for directory,dirs,names in os.walk(src,followlinks=False):
        for name in list(dirs):
            if name in SKIP or name in extra_skip or name.endswith('.egg-info'):
                skipped.append(dict(source=str(Path(directory)/name),reason='ENVIRONMENT_VENDOR_CACHE_OR_RAW_DATA'));dirs.remove(name)
        for name in sorted(names):
            p=Path(directory)/name
            if p.is_file() or p.is_symlink():
                rel=p.relative_to(src).as_posix()
                rule_rel=(src.name+'/'+rel) if group=='sv_history' else rel
                add(group,p,(Path(prefix)/rel).as_posix(),rule_rel,only_missing)

# Retire inherited scripts from the published parent using the verified cleanup
# ledger, never by deleting anything in the original workspaces.
cleanup=json.loads(Path('/home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z/audit/local_manifest.json').read_text())
deleted=[]
for row in cleanup['deletions']:
    p=Path(row['path']);target=None
    if p.is_relative_to(P):target=D/p.relative_to(P)
    elif p.is_relative_to(F):target=D/SV/p.relative_to(F)
    elif p.is_relative_to(O):target=D/SV/p.relative_to(O)
    if target and target.exists() and target.is_file():target.unlink();deleted.append(str(target.relative_to(D)))

# The Particle input binding changed since the parent; replace this directory
# exactly so an obsolete old-field VTU cannot masquerade as the active input.
frozen=D/SV/'frozen_reference';assert frozen.is_relative_to(D)
shutil.rmtree(frozen)
walk('particle',P/'particle_3d','particle_3d')
walk('flow',F,SV,extra_skip={'frozen_reference'})
walk('required_input',F/'frozen_reference',SV+'/upstream_stage_q_reference')
walk('required_input',P/SV/'frozen_reference',SV+'/frozen_reference')
walk('required_input',P/SV/'legacy_inputs',SV+'/legacy_inputs')
walk('rotate',O/'rotate_visualization',SV+'/rotate_visualization')
walk('old_fem',W/'formal_3D_flow_solver/FEM','formal_3D_flow_solver/FEM')
walk('sonovue',W/'sonovue_size_distribution_v0','sonovue_size_distribution_v0')
walk('vascular',V,'vascular_network',extra_skip={'vessel_model','archive'})
for folder in ['benchmarks','logs','reports','outputs']:
    walk('sv_history',O/folder,SV+'/'+folder,only_missing=True)

# Exact H0 upstream hash closure, including its single source SWC/ROI, without
# exporting the complete brain dataset. Preserve original byte contracts.
protected=json.loads((V/'reports/a_network_1d0d_boundary_v1/data/protected_input_hashes.json').read_text())
lineage=json.loads((V/'reports/a_network_1d0d_boundary_v1/data/a_network_provenance.json').read_text())['lineage']
required=set(protected)|{lineage['analysis_swc']}
path_map={}
for name in sorted(required):
    p=Path(name)
    if p.is_relative_to(V):target=Path('vascular_network')/p.relative_to(V)
    elif p.is_relative_to(F/'frozen_reference'):target=Path(SV)/'upstream_stage_q_reference'/p.relative_to(F/'frozen_reference')
    elif p.is_relative_to(F):target=Path(SV)/p.relative_to(F)
    elif p.is_relative_to(W/'formal_3D_flow_solver/FEM'):target=Path('formal_3D_flow_solver/FEM')/p.relative_to(W/'formal_3D_flow_solver/FEM')
    else:raise RuntimeError('Unmapped required dependency: '+name)
    if name in protected:assert sha(p)==protected[name],name
    add('required_input',p,target);path_map[name]=str(target)

for row in rows:
    if row['status']!='INCLUDED':continue
    p=Path(row['source']);q=D/row['target'];before=p.stat();h=sha(p)
    q.parent.mkdir(parents=True,exist_ok=True)
    if not q.exists() or sha(q)!=h:shutil.copy2(p,q)
    assert sha(q)==h and p.stat().st_size==before.st_size and p.stat().st_mtime_ns==before.st_mtime_ns,'Concurrent copy change: '+str(p)
    row['sha256']=h

# Preserve the original working-copy navigation and cleanup evidence as records;
# the new branch gets its own relative-path entry page later.
for p,target in [(W/'ACTIVE_VASCULAR_WORKFLOW.md','source_workspace_index.md'),(Path('/home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z/README.md'),'cleanup/README.md'),(Path('/home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z/SUMMARY.json'),'cleanup/SUMMARY.json'),(A/'audit/source_git_before.json','source_git_provenance.json')]:
    q=M/target;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
for name in ['local_plan.json','local_manifest.deleted.jsonl','server_server_manifest.deleted.jsonl','fem_regression_comparison.json']:
    p=Path('/home/lzy/archives/vascular_workflow_cleanup_20260924T214356Z/audit')/name
    if p.exists():shutil.copy2(p,M/'cleanup'/name)
summary=dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=D,text=True).strip(),branch=C['branch'],selected=len(selected),selected_bytes=sum(r['size'] for r in selected.values()),excluded=collections.Counter(r['reason'] for r in rows if r['status']=='EXCLUDED'),skipped_directories=skipped,inherited_retired_paths_removed=sorted(set(deleted)),maximum_new_file_bytes=MAX)
(A/'audit/local_selection.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n')
(M/'local_selection_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
(M/'original_path_map.json').write_text(json.dumps(path_map,ensure_ascii=False,indent=2)+'\n')
with (M/'local_file_inventory.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['status','group','source','target','size','sha256','reason'],extrasaction='ignore');w.writeheader();w.writerows(rows)
print(json.dumps({k:v for k,v in summary.items() if k not in {'skipped_directories','inherited_retired_paths_removed'}},indent=2))
print('INHERITED_RETIRED_REMOVED',len(set(deleted)))
