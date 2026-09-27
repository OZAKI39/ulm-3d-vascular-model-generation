"""Curated, byte-preserving current-workspace snapshot; never modifies sources."""
from pathlib import Path
import os,json,csv,hashlib,shutil,collections,subprocess
DEST=Path(__file__).resolve().parents[3]
B=Path('/home/lzy/projects');META=DEST/'sync_metadata/current_20260927'
SKIP={'.git','.venv','__pycache__','.pytest_cache','.mypy_cache','.ruff_cache','.idea','.vscode','node_modules','.cache','build-wsl'}
EXCLUDED=[];FILES=[];LINKS=[];LFS=[];seen=set()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def omit(p,reason):
    EXCLUDED.append(dict(source=str(p),reason=reason,bytes=p.stat().st_size if p.is_file() else None))
def copy(p,rel=None):
    p=Path(p);rel=Path(rel) if rel else p.relative_to(B)
    if str(rel) in seen:return
    seen.add(str(rel))
    if p.name in SKIP or p.name.endswith(':Zone.Identifier') or p.suffix in ['.pyc','.pyo']:
        omit(p,'runtime/cache/Git metadata or Windows alternate stream');return
    if p.name in ['.gitattributes','.gitignore','.ignore']:
        # Store old attributes as evidence; inherited LFS/EOL filters must not transform scientific blobs.
        if p.is_file():
            d=META/'original_git_rules'/rel;d.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,d)
        omit(p,'repository-local rules archived under sync_metadata/original_git_rules');return
    if 'SimVascularDistribution' in rel.parts:
        omit(p,'official installed mesher/runtime; package URL and SHA retained in mesh reports');return
    if str(rel)=='vascular_printing/third_party/vascularmd/Output':
        omit(p,'upstream demo generated meshes and OpenFOAM output unrelated to current A-H0 geometry');return
    if p.is_symlink():
        target=p.resolve();raw=os.readlink(p)
        # Do not pull unrelated printing-result branches through compatibility links.
        if str(rel).startswith('ulm_3D_vascular/outputs/') and rel.parts[2] not in ['rodent_vasculature','sampling','model_generate']:
            omit(p,'separate printing/old geometry result not on current A-H0 path');return
        if str(rel)=='ulm_3D_vascular/data':
            omit(p,'external human imaging datasets not required for mouse A-H0');return
        if not target.exists():omit(p,'source symlink already broken');return
        try:target_rel=target.relative_to(B)
        except ValueError:omit(p,'external installation outside workspace; retain dependency documentation');return
        dst=DEST/rel;dst.parent.mkdir(parents=True,exist_ok=True)
        link=os.path.relpath(DEST/target_rel,dst.parent)
        if not dst.exists() and not dst.is_symlink():dst.symlink_to(link,target_is_directory=target.is_dir())
        LINKS.append(dict(source=str(p),destination=str(rel),original_target=raw,snapshot_target=link))
        return
    if p.is_dir():
        for f in sorted(p.iterdir()):copy(f,rel/f.name)
        return
    if not p.is_file():return
    if p.stat().st_size < 256 and p.read_bytes().startswith(b'version https://git-lfs.github.com/spec/v1\n'):
        payload=p.read_bytes()
        LFS.append(dict(source=str(p),destination=str(rel),bytes=len(payload),sha256=sha(p),content=payload.decode()))
        omit(p,'upstream LFS placeholder only; exact pointer text preserved in upstream_lfs_pointers.json');return
    # Current primary trajectories and their identity manifests are all retained regardless of suffix.
    if p.stat().st_size>=50*1024**2:
        omit(p,'large reconstructible/source asset >=50 MiB; no current flow/trajectory files meet this limit');return
    if p.name=='candidates.npz' and 'streamlines' in rel.parts:
        omit(p,'41.6 MiB intermediate streamline seed candidates; accepted streamlines and renderer retained');return
    dst=DEST/rel;dst.parent.mkdir(parents=True,exist_ok=True)
    shutil.copy2(p,dst)
    if p.stat().st_mode&0o111:dst.chmod(0o755)
    else:dst.chmod(0o644)
    digest=sha(p);assert sha(dst)==digest
    FILES.append(dict(source=str(p),destination=str(rel),bytes=p.stat().st_size,sha256=digest))

def main():
    for name in ['formal_3D_flow_solver','ulm_flow_mean_2p0_mmps','ulm_particle_formal_p9a5','ulm_3D_vascular','sonovue_size_distribution_v0']:
        copy(B/name)
    for name in ['ACTIVE_VASCULAR_WORKFLOW.md','CURRENT_SERVER_PATHS.md','A_Global_Transient_Evolution.xlsx','A_Global_Transient_Evolution_A_H0.xlsx','FROZEN_SONOVUE_HISTOGRAM.csv']:
        copy(B/name)
    copy(B/'temp_storage/formal_3D_flow_solver_FEM/reports/stage00/source_contract.json')
    copy(B/'temp_storage/README.md')
    for name in ['cleanup_20260927','descriptive_names_20260927','mesh_code_migration_20260927','mesh_outputs_cleanup_20260927','mesh_reports_cleanup_20260927','migration_20260927']:
        copy(B/'temp_storage'/name)
    # Geometry source code shared by N through relative links. No unrelated imaging bundles.
    vp=B/'vascular_printing'
    for p in vp.iterdir():
        if p.is_file() and (p.suffix in ['.py','.yaml','.yml','.toml','.txt','.md','.json'] or p.name.startswith('LICENSE')):copy(p)
    for name in ['utils','config','configs','vascular_processing','tools','tests','test_data','docs','Ultraliser','third_party']:
        copy(vp/name)
    for path in ['outputs/rodent_vasculature/all_run_20260825_133152','outputs/sampling/20260825_133201_radius_plus_structure_k5','outputs/model_generate/ultraliser_anchor003274_20260825_133350']:
        copy(vp/path)
    mouse=vp/'vessel_model/T - A high-resolution dataset of mouse brain vasculature/raw_data/analysis_data/analysis_data'
    for kind,suffix in [('swc','.swc'),('images','.tif'),('mask','.tif')]:
        p=mouse/kind/('fMOST_0_5_6_0_0_6_0001_02_01'+suffix)
        if p.exists():copy(p)
    # Make all retained compatibility directory links resolvable; missing content is explicitly documented.
    for link in LINKS:
        dest=DEST/link['destination'];target=dest.parent/link['snapshot_target']
        if not target.exists():
            src=Path(link['source']).resolve()
            if src.is_file():copy(src)
            elif src.is_dir():
                target.mkdir(parents=True,exist_ok=True)
                (target/'SNAPSHOT_SCOPE.md').write_text('This optional shared directory is not part of the current mouse A/ROI/FEM/particle snapshot. See sync_metadata/current_20260927/EXCLUSIONS.md.\n')
    with (META/'local_file_inventory.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['source','destination','bytes','sha256']);w.writeheader();w.writerows(FILES)
    (META/'excluded_local.json').write_text(json.dumps(EXCLUDED,indent=2)+'\n')
    (META/'symlink_map.json').write_text(json.dumps(LINKS,indent=2)+'\n')
    (META/'upstream_lfs_pointers.json').write_text(json.dumps(LFS,indent=2)+'\n')
    prov={}
    for name in ['formal_3D_flow_solver/FEM_SimVascular','ulm_flow_mean_2p0_mmps','ulm_particle_formal_p9a5','ulm_3D_vascular','formal_3D_flow_solver/FEM_SimVascular/external/flow_solver_source']:
        p=B/name
        prov[name]={k:subprocess.check_output(['git','-C',str(p)]+args,text=True).strip() for k,args in [('commit',['rev-parse','HEAD']),('branch',['rev-parse','--abbrev-ref','HEAD'])]}
        status=subprocess.check_output(['git','-C',str(p),'status','--porcelain=v1','-uno'],text=True)
        prov[name]['tracked_change_count']=len(status.splitlines())
    (META/'source_git_provenance.json').write_text(json.dumps(prov,indent=2)+'\n')
    summary=dict(files=len(FILES),bytes=sum(r['bytes'] for r in FILES),links=len(LINKS),excluded_entries=len(EXCLUDED),upstream_lfs_pointer_records=len(LFS),base_commit='37c40d7c08ca6527c8379c0ffa81436ae708f59f',branch='sync/current-h0-dt1ms-wss-audit-20260927')
    (META/'local_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
