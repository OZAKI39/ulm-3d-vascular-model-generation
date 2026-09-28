"""Collect current sources without executing any scientific workflow."""
from pathlib import Path
import collections,csv,hashlib,json,os,shutil,stat,subprocess
D=Path(__file__).resolve().parents[3]
M=D/'sync_metadata/current_20260929'
B=Path('/home/lzy/projects');OLD=B/'temp_storage/github_sync_20260927'
SKIP={'.git','.venv','__pycache__','.pytest_cache','.mypy_cache','.ruff_cache','.cache','mplconfig','.vscode','node_modules','SimVascularDistribution','build-wsl'}
ROWS={};OMIT=[];LINKS=[]
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def dump(name,v): (M/name).write_text(json.dumps(v,ensure_ascii=False,indent=2)+'\n')
def omit(p,rel,reason,hash_file=False):
    row=dict(source=str(p),destination=str(rel),reason=reason)
    if p.is_file() and not p.is_symlink():
        row['bytes']=p.stat().st_size
        if hash_file:row['sha256']=sha(p)
    OMIT.append(row)
def copy(p,rel):
    if 'vessel_fine' in rel.parts and p.name in {'mesh-complete.mesh.vtu','mesh_arrays.npz','volume.vtu'}:
        omit(p,rel,'Previously excluded cancelled fine-CFD mesh; retained locally',True);return
    if str(rel).startswith('vascular_printing/third_party/vascularmd/Output'):
        omit(p,rel,'Third-party example/generated output outside the active flow ROI',p.is_file());return
    if p.name in SKIP or p.name.endswith(':Zone.Identifier') or p.suffix in {'.pyc','.pyo','.sock','.pid'}:
        omit(p,rel,'Runtime, installed environment, Git metadata or cache');return
    if p.name in {'.env','.env.local','id_rsa','id_ed25519','credentials.json'}:
        omit(p,rel,'Local authentication/environment file');return
    if p.name in {'.gitignore','.gitattributes','.ignore'}:
        if p.is_file():
            rel=Path('sync_metadata/current_20260929/original_git_rules')/rel
        else:return
    if p.is_symlink():
        raw=os.readlink(p);target=p.resolve()
        if not target.exists():omit(p,rel,'Source link already unresolved');return
        if target.is_relative_to(OLD):tr=target.relative_to(OLD)
        elif target.is_relative_to(B):tr=target.relative_to(B)
        else:omit(p,rel,'External runtime link');return
        dst=D/rel;dst.parent.mkdir(parents=True,exist_ok=True)
        link=os.path.relpath(D/tr,dst.parent)
        if dst.is_symlink():
            if os.readlink(dst)!=link:dst.unlink();dst.symlink_to(link,target_is_directory=target.is_dir())
        elif not dst.exists():dst.symlink_to(link,target_is_directory=target.is_dir())
        else:
            omit(p,rel,'Keep existing materialized snapshot instead of replacing it with a symlink');return
        LINKS.append(dict(source=str(p),destination=str(rel),original_target=raw,snapshot_target=link));return
    if p.is_dir():
        for child in sorted(p.iterdir()):copy(child,rel/child.name)
        return
    if not p.is_file():return
    size=p.stat().st_size
    if p.suffix.lower() in {'.so','.o','.a','.dll','.exe','.pyd'} or p.name=='svmultiphysics':
        omit(p,rel,'Rebuildable binary; source/build identity retained',True);return
    if size>=95*1024**2:
        omit(p,rel,'Large reconstructible cache/mesh or external dataset; retained at source',True);return
    if size<256 and p.read_bytes().startswith(b'version https://git-lfs.github.com/spec/v1\n'):
        omit(p,rel,'Upstream placeholder, not actual dataset',True);return
    if p.name=='candidates.npz' and 'streamlines' in rel.parts:
        omit(p,rel,'Redundant streamline candidates; accepted paths retained',True);return
    dst=D/rel;dst.parent.mkdir(parents=True,exist_ok=True)
    h=sha(p)
    if not dst.is_file() or dst.is_symlink() or dst.stat().st_size!=size or sha(dst)!=h:
        if dst.is_symlink():dst.unlink()
        shutil.copy2(p,dst)
    dst.chmod(0o755 if p.stat().st_mode&0o111 else 0o644)
    assert sha(dst)==h and sha(p)==h,str(p)
    ROWS[str(rel)]=dict(source=str(p),destination=str(rel),bytes=size,sha256=h)

def main():
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=D,text=True).strip()=='ca0ae424d01717f6a231bbea91405cbea8f095a8'
    # Preserve the proven source precedence of the prior snapshot.
    for name in ['ulm_3D_vascular','ulm_flow_mean_2p0_mmps']:copy(OLD/name,Path(name))
    for name in ['formal_3D_flow_solver','ulm_particle_formal_p9a5','sonovue_size_distribution_v0','brava_flow_roi_18mlmin']:
        copy(B/name,Path(name));print('Collected',name,flush=True)
    # Refresh previously selected shared geometry files, then add current BraVa
    # source/configuration and the actual accepted ROI/core provenance.
    for row in csv.DictReader((D/'sync_metadata/current_20260928/local_file_inventory.csv').open()):
        rel=Path(row['destination'])
        if rel.parts[0]=='vascular_printing' and (B/rel).exists():copy(B/rel,rel)
    V=B/'vascular_printing'
    for p in sorted(V.iterdir()):
        if p.is_file() and not p.is_symlink():copy(p,Path('vascular_printing')/p.name)
    for name in ['config','configs','tools','utils','vascular_processing','tests','docs','third_party/vascularmd']:
        copy(V/name,Path('vascular_printing')/name)
    for name in ['BG001.CNG.swc','BG001_ColorCoded.CNG.swc']:
        copy(V/'vessel_model/T - Brava/swc_files'/name,Path('vascular_printing/vessel_model/T - Brava/swc_files')/name)
    roi=V/'outputs/topbrain_brava_transfer/nn_production/BG001/RMCA'
    for p in sorted(roi.iterdir()):
        if p.is_file():copy(p,p.relative_to(B))
    copy(roi/'refined_roi',(roi/'refined_roi').relative_to(B))
    compact=roi/'compact_manufacturing_roi'
    for p in sorted(compact.iterdir()):
        if p.is_file():copy(p,p.relative_to(B))
    for name in ['candidates/BALANCED','print_fixture_design_all_ports_aligned']:
        copy(compact/name,(compact/name).relative_to(B))
    mold=compact/'final_abs_casting_mold'
    # Geometry/pose provenance is retained. Slicer packages and redundant full
    # casting boxes are unrelated to the solved lumen and are not reuploaded.
    for p in sorted(mold.rglob('*')):
        if p.is_file() and p.suffix.lower() in {'.json','.csv','.md','.txt','.log','.yaml','.yml'}:copy(p,p.relative_to(B))
    for p in sorted(V.glob('outputs/ui*/*')):
        if p.is_file() and p.suffix.lower() in {'.json','.csv','.md','.txt','.log'}:copy(p,p.relative_to(B))
    for name in ['ACTIVE_VASCULAR_WORKFLOW.md','CURRENT_SERVER_PATHS.md','A_Global_Transient_Evolution.xlsx','A_Global_Transient_Evolution_A_H0.xlsx','FROZEN_SONOVUE_HISTOGRAM.csv']:
        copy(B/name,Path(name))
    with (M/'source_files.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['source','destination','bytes','sha256']);w.writeheader();w.writerows(ROWS.values())
    dump('excluded_files.json',OMIT);dump('symlinks.json',LINKS)
    summary=dict(source_files=len(ROWS),logical_bytes=sum(r['bytes'] for r in ROWS.values()),symlinks=len(LINKS),excluded_entries=len(OMIT),scope='Incremental refresh; prior snapshot history retained; no scientific jobs run')
    dump('collection_summary.json',summary);print(json.dumps(summary,indent=2),flush=True)
if __name__=='__main__':main()
