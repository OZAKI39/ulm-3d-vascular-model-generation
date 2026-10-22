#!/usr/bin/env python3
"""Create reviewable project-specific Git snapshots; never edit original inputs."""
import argparse, collections, hashlib, json, os, re, shutil, subprocess
from pathlib import Path

BASE = Path('/home/lzy/projects')
AUDIT = BASE / 'github_sync_audits/flow_projects_20261022'
CE = BASE / 'computation_examples'
BRAVA = CE / 'brava_flow_roi_18mlmin'
BACKUP = BASE / 'server_migration/brava_vast4090_20261009T094459Z/payload/workspace/brava_flow_roi_18mlmin_20260928'
REMOTE = json.loads((AUDIT / 'REMOTE_BRAVA_INVENTORY.json').read_text())
CODE = {'.py','.sh','.cpp','.c','.h','.hpp','.f','.f90','.cmake','.patch','.diff','.toml','.xml','.yaml','.yml'}
SMALL_EVIDENCE = {'.json','.jsonl','.csv','.tsv','.log','.txt','.md','.html','.png','.pdf','.svg','.yaml','.yml','.xml'}
CACHE = {'.git','.venv','__pycache__','.pytest_cache','.ruff_cache','node_modules','CMakeFiles','uv-cache','cupy-cache','.rsync-partial'}

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()

def reason(rel, size, source):
    p=Path(rel); parts=p.parts; s='/'+rel+'/'
    if set(parts)&CACHE or p.name=='.git': return 'cache_or_git_metadata'
    if p.suffix in {'.pyc','.pyo','.o','.a','.so','.deb','.whl'} or '.so.' in p.name: return 'rebuildable_binary_or_dependency'
    if p.name in {'svmultiphysics','svMultiPhysics'}: return 'rebuildable_solver_binary_identity_kept'
    if '/external/SimVascularDistribution/' in s: return 'installed_third_party_distribution_source_recipe_kept'
    if p.suffix in {'.dll','.pyd','.exe'}: return 'installed_third_party_binary'
    if '/vessel_model/' in s and not p.name.startswith('fMOST_0_5_6_0_0_6_0001_02_01.'): return 'unrelated_geometry_library_only_actual_mouse_source_selected'
    if p.name.endswith(':Zone.Identifier'):return 'download_platform_metadata'
    if p.name=='STABILIZATION_TERMS.csv' and (source.parent/'STABILIZATION_TERMS.npz').exists(): return 'duplicate_uncompressed_table_npz_kept'
    if p.name.startswith('NotoSans') and p.suffix in {'.otf','.ttf'}: return 'installable_font'
    if any(x in parts for x in ['frames','render_frames','revisions','before_current_index_update']): return 'superseded_or_intermediate_visualization'
    if p.name in {'candidates.npz','gpu_mesh_input.npz','gpu_derived_field.npz'}: return 'derived_sampling_or_gpu_cache'
    if p.name.endswith('_REVIEW_BUNDLE.zip'): return 'duplicate_review_archive_unpacked_evidence_kept'
    if 'mesh_particle_validation_v2' in parts:
        active = any(x in s for x in ['/M0/','/M1_fixedplc_medium/','/M1_fixedplc_fine/'])
        if not active and p.suffix not in CODE and not (p.suffix in SMALL_EVIDENCE and size<5_000_000):
            return 'superseded_mesh_campaign_bulk_final_reports_kept'
        if p.name in {'FIELD_ARRAYS.npz','CELL_METRICS.csv','WALL_METRICS.csv'}: return 'recomputable_full_domain_diagnostic_array'
    if '/population_v1/runs/' in s and '/realization_0000_contactfix/' not in s and '/states/' in s:
        return 'superseded_population_states_final_accepted_states_kept'
    if '/population_v1/analysis/realization_0000_contactfix/tracks/' in s and p.suffix in {'.npz','.gz'}:
        return 'derived_per_track_copy_canonical_accepted_states_kept'
    if '/stage3/vessel_fine/' in s and p.suffix not in CODE and not (p.suffix in SMALL_EVIDENCE and size<2_000_000):
        return 'unrun_fine_mesh_bulk_skipped_by_user_status_kept'
    if '/failed_attempts/' in s and p.suffix in {'.vtu','.vtp','.npz','.bin'}:
        return 'failed_attempt_bulk_logs_and_identity_kept'
    if re.match(r'(result_\d+\.vtu|stFile_\d+\.bin)$',p.name):
        siblings=list(source.parent.glob('result_*.vtu' if p.suffix=='.vtu' else 'stFile_*.bin')) if source.exists() else []
        if siblings and p.name!=max(x.name for x in siblings) and p.name not in {'result_010.vtu','stFile_010.bin'}:
            return 'intermediate_solver_snapshot_final_and_restart_kept'
    if size>=100*1024*1024: return 'large_recomputable_audit_artifact_local_backup_sha_retained'
    return None

def main():
    ap=argparse.ArgumentParser();ap.add_argument('project',choices=['brava','mouse']);ap.add_argument('--copy-to',type=Path);a=ap.parse_args()
    candidates={}; excluded=[]; links=[]; mappings=[]
    def add_tree(src,dst=None):
        src=Path(src);dst=dst or str(src.relative_to(BASE));mappings.append((src,Path(dst)))
        if src.is_symlink(): links.append((src,dst)); return
        if src.is_file(): candidates[dst]=src; return
        for d,ds,fs in os.walk(src,followlinks=False):
            for n in list(ds):
                p=Path(d)/n
                if p.is_symlink():links.append((p,str(Path(dst)/p.relative_to(src))));ds.remove(n)
                elif n in CACHE: excluded.append({'source':str(p),'reason':'cache_or_git_metadata','directory':True});ds.remove(n)
            for n in fs:
                p=Path(d)/n;dest=str(Path(dst)/p.relative_to(src))
                if p.is_symlink():links.append((p,dest))
                else:candidates.setdefault(dest,p)
    if a.project=='brava':
        add_tree(BRAVA);add_tree(CE/'brava_calibration_equal_split')
        for n in ['src','tests','contracts']:
            add_tree(CE/'_shared/particle_workflow/particle_3d'/n)
        for n in ['pyproject.toml','requirements-particle8.txt','requirements-validation-lock.txt']:
            add_tree(CE/'_shared/particle_workflow/particle_3d'/n)
        for n in ['BRAVA_CONTEXT_HANDOFF.md','brain_report.md']:add_tree(BASE/n,'context/'+n)
        migration=BASE/'server_migration/brava_new_server_20261009'
        for n in ['MIGRATION_AND_VALIDATION_ZH.md','MIGRATION_STATUS.json','DELIVERY_CHECK.json','deploy','tools','server_evidence']:
            add_tree(migration/n,'deployment/newserver_20261009/'+n)
        remote_diffs=[];need=[]
        for row in REMOTE['files']:
            if not row['path'].startswith('project/') or 'size' not in row:continue
            rel=row['path'][8:];dst='computation_examples/brava_flow_roi_18mlmin/'+rel;local=BRAVA/rel
            if local.exists():
                if local.is_file() and sha(local)!=row['sha256']:remote_diffs.append({'path':rel,'local_sha256':sha(local),'server_sha256':row['sha256'],'choice':'current_local_authority'})
                continue
            # Reproducibility evidence on the current server, not retired mesher trees.
            if not (rel.startswith(('gpu_solver_fix/source/','population_v1/analysis/','production_workflow/manufactured_pipe/'))):continue
            p=BACKUP/rel
            why=reason(dst,row['size'],p)
            if why:excluded.append({'source':'new_server:'+row['path'],'destination':dst,'size':row['size'],'sha256':row['sha256'],'reason':why});continue
            if p.is_file() and sha(p)==row['sha256']:candidates[dst]=p
            else:
                p=AUDIT/'remote_selected'/row['path']
                if p.is_file() and sha(p)==row['sha256']:candidates[dst]=p
                else:need.append(row['path'])
        (AUDIT/'REMOTE_LOCAL_COMPARISON.json').write_text(json.dumps(remote_diffs,indent=2))
        if need:
            (AUDIT/'remote_needed.txt').write_text('\n'.join(need)+'\n');raise SystemExit(f'Need remote files: {len(need)}')
    else:
        for p in sorted(CE.glob('mean-*')):add_tree(p)
        for n in ['validation_cases','fem_reference_stage_q']:add_tree(CE/n)
        for n in ['particle_workflow','formal_3D_flow_solver','geometry_network','legacy_flow_workflow','wss_audit','wss_validation_v2','fem_reference_contracts','templates']:
            add_tree(CE/'_shared'/n)
        add_tree(CE/'_shared/active_network_development/ulm_3D_vascular')
        # Shared upstream source required by mouse ROI reconstruction; no human/MeVO outputs.
        up=BASE/'vascular_printing'
        for n in ['vascular_processing','utils','main.py','s1-1_swc_roi_generate_mouse.py','s2_swc_stl_model_generate.py','requirements.txt','requirements-fixture.txt','requirements-vascularmd.txt','cfd_lumen_config.yaml','outputs/sampling/20260825_133201_radius_plus_structure_k5','outputs/model_generate/ultraliser_anchor003274_20260825_133350','outputs/rodent_vasculature/all_run_20260825_133152']:
            if (up/n).exists():add_tree(up/n,'geometry_inputs/'+n)
        raw=up/'vessel_model/T - A high-resolution dataset of mouse brain vasculature/raw_data/analysis_data/analysis_data'
        for sub,ext in [('swc','swc'),('images','tif'),('mask','tif')]:
            p=raw/sub/('fMOST_0_5_6_0_0_6_0001_02_01.'+ext)
            add_tree(p,'geometry_inputs/'+str(p.relative_to(up)))
        for p in (up/'configs').glob('*'):
            if p.is_file() and not any(x in p.name.lower() for x in ['human','brava','mevo','topbrain']):add_tree(p,'geometry_inputs/configs/'+p.name)
    for n in ['distribution_assets','sonovue_size_distribution_v0']:add_tree(CE/'_shared'/n)
    included=[]
    for dst,p in sorted(candidates.items()):
        if not p.is_file():continue
        size=p.stat().st_size;why=reason(dst,size,p)
        if size<1024 and p.read_bytes().startswith(b'version https://git-lfs.github.com/spec/v1\n'):
            why='upstream_test_resource_lfs_pointer_no_local_payload'
        record={'path':dst,'source':str(p),'size':size}
        if why:
            record['reason']=why
            if size>=1_000_000 and why not in {'rebuildable_binary_or_dependency','installed_third_party_distribution_source_recipe_kept'}:record['sha256']=sha(p)
            excluded.append(record);continue
        record['sha256']=sha(p);record['mode']='100755' if p.stat().st_mode&0o111 else '100644';included.append(record)
        if a.copy_to:
            target=a.copy_to/dst;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
    selected={x['path'] for x in included};kept_links=[]
    def map_target(p):
        for src,dst in sorted(mappings,key=lambda x:len(str(x[0])),reverse=True):
            try:return str(dst/p.relative_to(src))
            except ValueError:pass
        try:return str(p.relative_to(BASE))
        except ValueError:return None
    for p,dst in links:
        why=reason(dst,0,p)
        if why:excluded.append({'source':str(p),'path':dst,'reason':why});continue
        target=map_target(p.resolve())
        if target and (target in selected or any(x.startswith(target+'/') for x in selected)):
            link=os.path.relpath(target,Path(dst).parent);rec={'path':dst,'source':str(p),'old_target':os.readlink(p),'target':link}
            kept_links.append(rec)
            if a.copy_to:
                out=a.copy_to/dst;out.parent.mkdir(parents=True,exist_ok=True)
                if not out.exists() and not out.is_symlink():out.symlink_to(link)
        else:excluded.append({'path':dst,'source':str(p),'target':os.readlink(p),'reason':'external_or_excluded_compatibility_target_documented'})
    summary={'project':a.project,'files':len(included),'bytes':sum(x['size'] for x in included),'symlinks':len(kept_links),'excluded_files_or_directories':len(excluded),'excluded_bytes':sum(x.get('size',0) for x in excluded),'policy':'Canonical accepted states retained; no numerical recomputation; ordinary Git files <100 MiB; no LFS.'}
    for name,value in [('INCLUDED',included),('EXCLUDED',excluded),('LINKS',kept_links),('SUMMARY',summary)]:
        (AUDIT/f'{a.project.upper()}_{name}.json').write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
