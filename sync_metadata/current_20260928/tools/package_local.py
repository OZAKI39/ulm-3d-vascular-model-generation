"""Read-only collection into a new isolated Git worktree; preserve scientific bytes."""
from pathlib import Path
import collections,csv,hashlib,json,os,shutil,subprocess
D=Path(__file__).resolve().parents[3];M=D/'sync_metadata/current_20260928'
B=Path('/home/lzy/projects');OLD=B/'temp_storage/github_sync_20260927'
SKIP={'.git','.venv','__pycache__','.pytest_cache','.mypy_cache','.ruff_cache','.cache','.vscode','node_modules','SimVascularDistribution','build-wsl'}
FILES=[];LINKS=[];OMIT=[];POINTERS=[];RULES=[]
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(name,value): (M/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def omit(p,rel,reason):
    row=dict(source=str(p),destination=str(rel),reason=reason)
    if p.is_file() and not p.is_symlink():row.update(bytes=p.stat().st_size,sha256=sha(p))
    OMIT.append(row)
def copy(p,rel):
    if p.name in SKIP or p.name.endswith(':Zone.Identifier') or p.suffix in {'.pyc','.pyo','.sock'}:
        omit(p,rel,'Runtime, cache, installed distribution, Git metadata or alternate stream');return
    if p.name in {'.env','.env.local','id_rsa','id_ed25519'}:
        omit(p,rel,'Local credentials/environment outside scientific snapshot');return
    if p.name in {'.gitignore','.gitattributes','.ignore'}:
        if p.is_file():
            dst=M/'original_git_rules'/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
            RULES.append(dict(source=str(p),destination=str(dst.relative_to(D)),bytes=p.stat().st_size,sha256=sha(p)))
        return
    if p.is_symlink():
        raw=os.readlink(p);target=p.resolve()
        if not target.exists():omit(p,rel,'Already broken source link');return
        if target.is_relative_to(OLD):tr=target.relative_to(OLD)
        elif target.is_relative_to(B):tr=target.relative_to(B)
        else:omit(p,rel,'External installed dependency; documented, not vendored');return
        dst=D/rel;dst.parent.mkdir(parents=True,exist_ok=True)
        link=os.path.relpath(D/tr,dst.parent)
        if dst.is_symlink():dst.unlink()
        dst.symlink_to(link,target_is_directory=target.is_dir())
        LINKS.append(dict(source=str(p),destination=str(rel),original_target=raw,snapshot_target=link));return
    if p.is_dir():
        for child in sorted(p.iterdir()):copy(child,rel/child.name)
        return
    if not p.is_file():return
    size=p.stat().st_size
    if size<256 and p.read_bytes().startswith(b'version https://git-lfs.github.com/spec/v1\n'):
        POINTERS.append(dict(source=str(p),destination=str(rel),text=p.read_text(),sha256=sha(p)))
        omit(p,rel,'Upstream LFS placeholder; pointer metadata saved, payload not present');return
    important_gpu=p.name=='gpu_mesh_input.npz' and 'microbubble_' in str(rel)
    if size>=50*1024**2 and not (important_gpu and size<95*1024**2):
        omit(p,rel,'Large non-production asset >=50 MiB; cancelled fine mesh or reconstructible intermediate');return
    if p.name=='candidates.npz' and 'streamlines' in rel.parts:
        omit(p,rel,'Redundant streamline seed candidates; accepted streamlines and renderer retained');return
    dst=D/rel;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dst)
    dst.chmod(0o755 if p.stat().st_mode&0o111 else 0o644)
    h=sha(p);assert sha(dst)==h
    FILES.append(dict(source=str(p),destination=str(rel),bytes=size,sha256=h))
def main():
    # This prior sync worktree is now the actual live network/ROI CFD development root.
    for name in ['ulm_3D_vascular','ulm_flow_mean_2p0_mmps','server_evidence','sync_metadata','temp_storage']:
        copy(OLD/name,Path(name))
    # Shared geometry scope was curated in the previous snapshot; refresh those same paths.
    for root,dirs,files in os.walk(OLD/'vascular_printing',followlinks=False):
        for name in dirs[:]:
            path=Path(root)/name
            if path.is_symlink():copy(path,path.relative_to(OLD));dirs.remove(name)
        for name in files:
            old=Path(root)/name;rel=old.relative_to(OLD);live=B/rel
            copy(live if live.exists() else old,rel)
    for name in ['formal_3D_flow_solver','ulm_particle_formal_p9a5','sonovue_size_distribution_v0']:
        copy(B/name,Path(name))
    for name in ['ACTIVE_VASCULAR_WORKFLOW.md','CURRENT_SERVER_PATHS.md','A_Global_Transient_Evolution.xlsx','A_Global_Transient_Evolution_A_H0.xlsx','FROZEN_SONOVUE_HISTOGRAM.csv']:
        copy(B/name,Path(name))
    # Latest network and flow worktree is authoritative; explicitly record all root conflicts.
    conflicts=[]
    for name in ['ulm_3D_vascular','ulm_flow_mean_2p0_mmps']:
        for r in FILES:
            if not r['destination'].startswith(name+'/'):continue
            p=B/r['destination']
            if p.is_file() and not p.is_symlink() and sha(p)!=r['sha256']:
                conflicts.append(dict(destination=r['destination'],selected=r['source'],other_source=str(p),
                    selected_sha256=r['sha256'],other_sha256=sha(p),reason='Current live network/ROI development is in previous sync worktree'))
    with (M/'local_file_inventory.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['source','destination','bytes','sha256']);w.writeheader();w.writerows(FILES+RULES)
    dump('excluded_local.json',OMIT);dump('symlink_map.json',LINKS);dump('upstream_lfs_pointers.json',POINTERS)
    dump('source_precedence.json',conflicts)
    prov={}
    for root in [OLD,B/'ulm_particle_formal_p9a5',B/'ulm_3D_vascular',B/'formal_3D_flow_solver/FEM_SimVascular']:
        status=subprocess.check_output(['git','-C',str(root),'status','--porcelain=v1','-uno'],text=True)
        prov[str(root)]={'HEAD':subprocess.check_output(['git','-C',str(root),'rev-parse','HEAD'],text=True).strip(),
            'tracked_change_count':len(status.splitlines())}
    dump('source_git_provenance.json',prov)
    (D/'.gitattributes').write_text('# Preserve scientific source/data bytes. No LFS or EOL transformation.\n* -text -filter\n')
    (D/'.gitignore').write_text('**/__pycache__/\n**/.pytest_cache/\n**/.venv/\n**/.cache/\n**/node_modules/\n*.pyc\n*.pyo\n*:Zone.Identifier\n')
    summary=dict(files=len(FILES),bytes=sum(r['bytes'] for r in FILES),symlinks=len(LINKS),
        excluded_entries=len(OMIT),pointer_records=len(POINTERS),source_conflicts=len(conflicts),archived_rules=len(RULES))
    dump('local_summary.json',summary);print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
