"""Curated, read-only-source export of the existing FEM / MB / RBC workflow."""
from pathlib import Path
import argparse,os,json,hashlib,csv,shutil,subprocess,collections

OMIT_DIRS={'.git','.venv','venv','__pycache__','.pytest_cache','.ruff_cache','.mypy_cache',
           'node_modules','jit_cache','.codex_tmp','build','dist'}
MAX_NEW_FILE=25*1024**2
TEXT={'.json','.jsonl','.yaml','.yml','.xml','.csv','.txt','.log','.out','.err','.md','.html','.sh','.py','.toml'}

def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
    return h.hexdigest()

def policy(group,rel,p):
    parts=set(Path(rel).parts);n=p.name;size=p.stat().st_size
    if n=='.git':return 'LOCAL_GIT_WORKTREE_METADATA'
    if n in {'.env','hosts.yml','credentials','credentials.json','id_rsa','id_ed25519'} or n.startswith('.env.') or p.suffix in {'.pem','.key','.p12'}:
        return 'CREDENTIAL_OR_LOCAL_SECRET_FILE'
    if group=='sonovue':return None # Its full frozen checksum contract is mandatory.
    if n.endswith(('.tar','.tar.gz','.tgz','.zip','.7z')):return 'REDUNDANT_DEPLOYMENT_OR_ASSET_ARCHIVE'
    if n.endswith(('.pyc','.pyo','.so','.o','.a','.exe')):return 'REBUILDABLE_COMPILED_ARTIFACT'
    if group=='rotate' and (n.endswith('_annotations.json') or 'inspection' in parts):
        return None # Required by the current independent visualization validators.
    if 'frames' in parts or 'inspection' in parts:return 'REGENERABLE_RENDER_FRAMES'
    if 'field_diagnostics_revisions' in parts or 'rotate_visualization_history' in parts:return 'SUPERSEDED_VISUALIZATION_REVISION'
    if group=='particle':
        if rel.startswith('particle_3d/outputs/particle8_2a/admission/'):return 'REGENERABLE_LARGE_ADMISSION_CACHE'
        if '/development_preview/' in '/'+rel:return 'SUPERSEDED_DEVELOPMENT_PREVIEW'
        if rel.endswith('trajectory_samples.csv.gz'):return 'DUPLICATE_TABULAR_EXPORT_CANONICAL_TRAJECTORY_NPZ_RETAINED'
        if n=='11_long_lifecycle.csv.gz':return 'REGENERABLE_LONG_VALIDATION_TRACE_SUMMARY_RETAINED'
        if n=='C57BL6_RBC_GEOMETRY_VALIDATION_100000.csv':return 'REGENERABLE_BULK_GEOMETRY_SAMPLE_SUMMARY_RETAINED'
        if 'shear_lift_audit' in parts:
            if '/data/remote_results/states/' in '/'+rel:return 'REGENERABLE_PER_STATE_LIFT_AUDIT_FORMAL_TRAJECTORIES_RETAINED'
            if n=='all_scalar_states.npz' or n.endswith('_eulerian.npz'):return 'REGENERABLE_FULL_AUDIT_ARRAY_STATISTICS_AND_REPRESENTATIVES_RETAINED'
            if '/literature/' in '/'+rel and p.suffix in {'.pdf','.html','.txt'}:return 'THIRD_PARTY_LITERATURE_DOWNLOAD_CITATION_RECORD_RETAINED'
            if '/data/local_smoke/' in '/'+rel or '/data/local_nearwall_smoke/' in '/'+rel:return 'LOCAL_SMOKE_DUPLICATE_FULL_AUDIT_SUMMARY_RETAINED'
        if n.endswith('_frames.json') and size>2*1024**2:return 'REGENERABLE_REPLAY_FRAME_CACHE'
    if group=='new_flow':
        if '/run/1-procs/' in '/'+rel and (n.startswith('result_') or n.startswith('stFile_')):
            if not n.startswith(('result_070','result_071','stFile_071')):return 'INTERMEDIATE_SOLVER_STATE_STEADY_GATE_AND_FINAL_NATIVE_OUTPUT_RETAINED'
        if '/streamlines/presentation' in '/'+rel:return 'SUPERSEDED_PRESENTATION_LATEST_ROTATE_VISUALIZATION_RETAINED'
        if '/field_diagnostics/animations/' in '/'+rel:return 'SUPERSEDED_ROTATION_LATEST_ROTATE_VISUALIZATION_RETAINED'
        if n=='candidates.npz':return 'REGENERABLE_STREAMLINE_CANDIDATE_CACHE_SELECTED_PATHS_RETAINED'
    if group=='historical_fem':
        if 'remote_return' in parts or 'previous_local_versions' in parts:return 'DUPLICATE_REMOTE_RETURN_OR_HISTORICAL_COPY'
        if rel.startswith('outputs/') and p.suffix not in TEXT and not rel.startswith('outputs/stage01_7/selected/mesh/'):
            return 'HISTORICAL_INTERMEDIATE_BINARY_CURRENT_FROZEN_FLOW_RETAINED'
    if group=='original_sv_history' and rel.startswith('outputs/'):
        if p.suffix not in TEXT or size>2*1024**2:return 'HISTORICAL_SOLVER_INTERMEDIATE_CURATED_FROZEN_PACKAGE_RETAINED'
    if n.endswith('_annotations.json') and size>2*1024**2:return 'REGENERABLE_PER_FRAME_LABEL_LAYOUT'
    if n in {'pre_install_environment.json','final_environment.json'} and size>2*1024**2:
        return 'VERBOSE_SERVER_ENVIRONMENT_DUMP_SUMMARY_AND_LOCKFILES_RETAINED'
    if size>MAX_NEW_FILE:return 'LARGE_REGENERABLE_FILE_OVER_25_MIB'
    return None

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workspace',type=Path,default=Path('/home/lzy/projects'))
    ap.add_argument('--destination',type=Path,required=True);ap.add_argument('--audit',type=Path,required=True)
    ap.add_argument('--copy',action='store_true');args=ap.parse_args();w=args.workspace.resolve();dst=args.destination.resolve()
    args.audit.mkdir(parents=True,exist_ok=True)
    base=set(subprocess.check_output(['git','ls-tree','-r','--name-only','HEAD'],cwd=dst,text=True).splitlines())
    sv='formal_3D_flow_solver/FEM_SimVascular'
    mappings=[('particle',w/'ulm_particle_3d_particle0','',False),
      ('new_flow',w/'ulm_flow_mean_2p0_mmps'/sv/'flow_cases',sv+'/flow_cases',False),
      ('new_flow_code',w/'ulm_flow_mean_2p0_mmps'/sv/'scripts/flow_2mmps',sv+'/scripts/flow_2mmps',False),
      ('new_flow_tests',w/'ulm_flow_mean_2p0_mmps'/sv/'tests/flow_2mmps',sv+'/tests/flow_2mmps',False),
      ('rotate',w/sv/'rotate_visualization',sv+'/rotate_visualization',False),
      ('historical_fem',w/'formal_3D_flow_solver/FEM','formal_3D_flow_solver/FEM',False),
      ('sonovue',w/'sonovue_size_distribution_v0','sonovue_size_distribution_v0',False),
      ('source_reference',w/'ulm_microbubble_traj_gen_2D','source_dependencies/ulm_microbubble_traj_gen_2D',False),
      ('source_reference',w/'vascular_migration_audit','source_dependencies/vascular_migration_audit',False)]
    for folder in ['benchmarks','logs','reports','outputs']:
        mappings.append(('original_sv_history',w/sv/folder,sv+'/'+folder,True))
    rows=[];selected={};skipped_dirs=[]
    for group,src,prefix,only_missing in mappings:
        for directory,dirs,files in os.walk(src,followlinks=False):
            for name in list(dirs):
                if name in OMIT_DIRS or name.endswith('.egg-info'):
                    skipped_dirs.append(dict(source=str(Path(directory)/name),reason='ENVIRONMENT_CACHE_BUILD_OR_GIT_STORAGE'))
                    dirs.remove(name)
            for name in sorted(files):
                p=Path(directory)/name;rel=p.relative_to(src).as_posix();target=(Path(prefix)/rel).as_posix()
                if p.is_symlink():
                    # Scientific source maps contain regular files; external symbolic references stay explicit.
                    rows.append(dict(source=str(p),target=target,group=group,bytes=0,status='EXCLUDED',reason='SYMLINK_REQUIRES_EXPLICIT_REVIEW',link=os.readlink(p)));continue
                if not p.is_file():continue
                if group=='original_sv_history':rule_rel=p.relative_to(w/sv).as_posix()
                else:rule_rel=rel
                reason=policy(group,rule_rel,p)
                if target in base:reason=None # Never break the already published frozen dependency closure.
                row=dict(source=str(p),target=target,group=group,bytes=p.stat().st_size,status='EXCLUDED' if reason else 'INCLUDED',reason=reason or 'CURRENT_WORKFLOW_SOURCE_OR_EVIDENCE')
                if only_missing and ((dst/target).exists() or target in selected):
                    row.update(status='ALREADY_PACKAGED',reason='PUBLISHED_CURATED_COPY_TAKES_PRECEDENCE')
                elif not reason:
                    if target in selected:raise RuntimeError('Unexpected overlapping export target: '+target)
                    selected[target]=row
                rows.append(row)
    if args.copy:
        for row in rows:
            p=Path(row['source'])
            if row['reason'] in {'CREDENTIAL_OR_LOCAL_SECRET_FILE','SYMLINK_REQUIRES_EXPLICIT_REVIEW','LOCAL_GIT_WORKTREE_METADATA'}:continue
            row['source_sha256']=sha(p)
            if row['status']=='INCLUDED':
                q=dst/row['target'];q.parent.mkdir(parents=True,exist_ok=True)
                if not q.exists() or sha(q)!=row['source_sha256']:shutil.copy2(p,q)
                assert sha(q)==row['source_sha256'],row['target']
    counts=collections.defaultdict(lambda:dict(files=0,bytes=0))
    for row in rows:
        key=row['status']+' / '+row['group'];counts[key]['files']+=1;counts[key]['bytes']+=row['bytes']
    summary=dict(copied=args.copy,source_workspace=str(w),destination=str(dst),base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=dst,text=True).strip(),
      maximum_new_file_bytes=MAX_NEW_FILE,counts=dict(counts),skipped_directories=skipped_dirs,
      geometry_scope='Existing v8 source-geometry code and its reviewed ROI fixture retained from published base; full unrelated raw vascular datasets not recopied.',
      frozen_policy='All published base files retained, including frozen SHA256 dependency closure.')
    (args.audit/'snapshot_plan.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2)+'\n')
    (args.audit/'snapshot_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
    fields=['status','group','source','target','bytes','source_sha256','reason']
    with (args.audit/'file_inventory.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
    print(json.dumps(dict(counts),indent=2));print('Skipped directories:',len(skipped_dirs))
    print('Excluded reasons:',dict(collections.Counter(row['reason'] for row in rows if row['status']=='EXCLUDED')))

if __name__=='__main__':main()
