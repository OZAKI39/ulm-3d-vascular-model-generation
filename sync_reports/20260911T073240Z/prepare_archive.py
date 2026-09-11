"""One-off CPU archival inventory/copy; never imports or launches a solver."""
from pathlib import Path
from datetime import datetime, timezone
import collections, hashlib, json, os, shutil, subprocess

ROOT=Path(__file__).resolve().parents[2]
REPORT=Path(__file__).resolve().parent
PROJECTS=Path('/home/lzy/projects')
GIT_DIR=Path(subprocess.check_output(['git','rev-parse','--absolute-git-dir'],cwd=ROOT,text=True).strip())
AUDIT=GIT_DIR/'sync_work'
CONTEXT=json.loads((AUDIT/'context.json').read_text())
REPAIR='rbc_repair_20260910T131105Z'
NATIVE=PROJECTS/f'mirheo_starter/data/single_rbc_repair/{REPAIR}/native/source'

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1048576),b''):h.update(chunk)
    return h.hexdigest()

def write(path,value):
    Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')

def git_state(path):
    r=subprocess.run(['git','status','--short'],cwd=path,capture_output=True,text=True)
    return dict(path=str(path),exit_code=r.returncode,stdout=r.stdout,stderr=r.stderr)

def main():
    if (AUDIT/'copied_inventory.json').exists():raise RuntimeError('ARCHIVE_COPY_ALREADY_DONE')
    sources={};excluded=[];excluded_roots=[]
    def add(path,category,target=None):
        path=Path(path);target=target or path.relative_to(PROJECTS).as_posix()
        if path.is_symlink():
            excluded.append(dict(source_path=str(path),reason='symlink not copied',link=os.readlink(path)));return
        if not path.is_file():raise FileNotFoundError(path)
        previous=sources.get(target)
        if previous and previous['source_path']!=str(path):raise RuntimeError('SOURCE_MAPPING_COLLISION')
        sources[target]=dict(source_path=str(path),target_path=target,category=category)
    def walk(root,category,skip=()):
        root=Path(root)
        if not root.exists():return
        for parent,dirs,names in os.walk(root,followlinks=False):
            for name in list(dirs):
                p=Path(parent)/name
                if name in {'.git','.venv','__pycache__','.pytest_cache','node_modules'} or p in skip or p.is_symlink():
                    dirs.remove(name);excluded_roots.append(dict(source_path=str(p),reason='environment/cache/third-party or separately archived subtree'))
            for name in names:
                p=Path(parent)/name
                if p.suffix in {'.pyc','.so','.a','.o','.exe','.lock','.env'} or name=='.git' or name.startswith('.env'):
                    excluded.append(dict(source_path=str(p),size_bytes=p.stat().st_size,sha256=sha(p),reason='compiled/runtime/environment file'));continue
                add(p,category)

    write(REPORT/'SOURCE_GIT_STATE.json',dict(before=[git_state(PROJECTS/p) for p in ['mirheo_starter','hemocell_starter','cloud_compute','cloud_results']]+[git_state(NATIVE)]))
    walk(PROJECTS/'cloud_compute','cloud orchestration / frozen workers / test evidence')
    walk(PROJECTS/'cloud_results','cloud build and smoke/full-run raw archive and review')
    walk(PROJECTS/'cloud_upload','upload provenance and transfer tools',skip=[PROJECTS/'cloud_upload/20260910T212439Z/source'])
    for relative in ['py_scripts/single_rbc_repair','py_scripts/single_rbc_benchmark','py_scripts/solver_benchmark',
                     f'data/single_rbc_repair/{REPAIR}',f'runs/single_rbc_repair/{REPAIR}',f'test_code/outputs/single_rbc_repair/{REPAIR}']:
        walk(PROJECTS/'mirheo_starter'/relative,'local repaired RBC code / actual A0-A6 evidence',skip=[NATIVE,NATIVE.parent/'build'])
    for pattern in ['*single_rbc*.py','*single_rbc*.cjs','README_single_rbc*.md']:
        for p in (PROJECTS/'mirheo_starter/test_code').glob(pattern):add(p,'RBC checks and documentation')
    for relative in ['README.md','cases','scripts','metadata','logs']:
        p=PROJECTS/'hemocell_starter'/relative
        if p.is_dir():walk(p,'HemoCell source/environment evidence (mostly inherited)')
        else:add(p,'HemoCell documentation')
    # Include all actually modified/untracked native files, plus the previously
    # archived audit subset. Never mistake git diff alone for untracked content.
    changed=subprocess.check_output(['git','diff','--name-only','HEAD','-z'],cwd=NATIVE).split(b'\0')
    changed+=subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'],cwd=NATIVE).split(b'\0')
    for value in changed:
        if value:add(NATIVE/value.decode(),'native local modifications / untracked header')
    prefix=f'mirheo_starter/data/single_rbc_repair/{REPAIR}/native/source/'
    inherited=subprocess.check_output(['git','ls-files','-z',prefix],cwd=ROOT).split(b'\0')
    for value in inherited:
        if value:
            relative=value.decode();add(PROJECTS/relative,'native previously archived audit subset')
    for name in ['LICENSE','CMakeLists.txt','.gitmodules']:
        if (NATIVE/name).exists():add(NATIVE/name,'native license / upstream build metadata')
    build=NATIVE.parent/'build'
    for name in ['CMakeCache.txt','compile_commands.json']:
        if (build/name).exists():add(build/name,'existing native build record',f'{REPORT.relative_to(ROOT).as_posix()}/native_build_evidence/{name}')
    native_state=dict(upstream_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=NATIVE,text=True).strip(),
        upstream_url='https://github.com/cselab/Mirheo.git',submodules=subprocess.check_output(['git','submodule','status'],cwd=NATIVE,text=True),
        modified_and_untracked=[x.decode() for x in changed if x],full_upstream_tree_included=False)
    write(REPORT/'NATIVE_SOURCE_IDENTITY.json',native_state)
    source_patch=subprocess.check_output(['git','diff','--binary','HEAD'],cwd=NATIVE)
    (REPORT/'native_current_tracked_changes.patch').write_bytes(source_patch)
    # Enumerate omitted native source/build files without reading Git metadata.
    for omitted in [NATIVE,build]:
        for parent,dirs,names in os.walk(omitted,followlinks=False):
            dirs[:]=[x for x in dirs if x not in ['.git','__pycache__']]
            for name in names:
                p=Path(parent)/name
                if name=='.git' or p.is_symlink() or any(e['source_path']==str(p) for e in sources.values()):continue
                excluded.append(dict(source_path=str(p),size_bytes=p.stat().st_size,sha256=sha(p),execute_bits=p.stat().st_mode&0o111,
                    reason='upstream tree represented by commit, licenses, original-source audit subset and full local modifications' if omitted==NATIVE else 'compiled build/cache; selected CMake evidence archived separately'))
    copied=[];inherited_count=0
    for target,entry in sorted(sources.items()):
        source=Path(entry['source_path']);before=source.stat();h=sha(source);destination=ROOT/target
        same=destination.is_file() and sha(destination)==h and bool(destination.stat().st_mode&0o111)==bool(before.st_mode&0o111)
        entry.update(size_bytes=before.st_size,sha256=h,execute_bits=before.st_mode&0o111,executable=bool(before.st_mode&0o111),
                     source_mtime_ns=before.st_mtime_ns,copy_type='byte-identical',change='inherited_identical' if same else 'updated' if destination.exists() else 'added')
        if same:inherited_count+=1
        else:
            if destination.is_symlink():raise RuntimeError('TARGET_SYMLINK')
            destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(source,destination)
            assert sha(destination)==h and sha(source)==h
            assert source.stat().st_mtime_ns==before.st_mtime_ns
        copied.append(entry)
    write(AUDIT/'copied_inventory.json',dict(context=CONTEXT,files=copied))
    write(REPORT/'EXCLUDED_FILES.json',dict(files=excluded,excluded_roots=excluded_roots,policy=[
        'No environments, credentials, nested .git, binaries, installation caches or lock files.',
        'Full uploaded source duplicate is represented by transfer manifest and canonical project sources.',
        'All A0-A6 raw scientific outputs and all cloud result archives are included without downsampling.',
        'Full upstream native source/build environments omitted; original source audit subset, modified/untracked sources, patches, license and build identity are present.'],
        impact='No selected scientific observation omitted. Full rebuild environment is not packaged; original experiment missing measurements remain missing.'))
    changed=[e for e in copied if e['change']!='inherited_identical']
    summary=dict(selected_files=len(copied),inherited_identical=inherited_count,new_or_updated_source_files=len(changed),
                 new_or_updated_source_bytes=sum(e['size_bytes'] for e in changed),largest=sorted(changed,key=lambda x:-x['size_bytes'])[:12],
                 files_over_50_MiB=[e for e in changed if e['size_bytes']>50*1024*1024],excluded_files=len(excluded))
    write(REPORT/'LARGE_FILE_REVIEW.json',summary)
    print(json.dumps({k:v for k,v in summary.items() if k!='largest'},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
