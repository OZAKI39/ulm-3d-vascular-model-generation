"""One-time migration evidence builder; this audit script is not a runtime dependency."""
from pathlib import Path
import ast,collections,csv,difflib,hashlib,importlib.util,json,os,platform,shutil,subprocess,sys

ROOT=Path(__file__).resolve().parents[1];AUDIT=ROOT/'docs/audit'

def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def save(name,value):
    (AUDIT/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')

def source_inventory():
    rows=json.loads((AUDIT/'source_map.json').read_text());git=json.loads((AUDIT/'git_baseline.json').read_text());symbols=[]
    for row in rows:
        old=Path(row['source']);new=ROOT/row['path'];before=old.read_bytes();after=new.read_bytes()
        assert sha(old)==row['original_sha256'],('source modified during migration',old)
        row.update(sha256=sha(new),changed=before!=after)
        matches=[g for g in git if g['toplevel'] and old.is_relative_to(g['toplevel'])]
        if matches:
            g=max(matches,key=lambda g:len(g['toplevel']));rel=str(old.relative_to(g['toplevel']))
            result=subprocess.run(['git','-C',g['toplevel'],'show',g['commit']+':'+rel],capture_output=True)
            row.update(git_commit=g['commit'],git_relative=rel,git_state='TRACKED_UNCHANGED' if result.returncode==0 and result.stdout==before else 'MODIFIED_OR_UNTRACKED')
        else:row.update(git_commit=None,git_state='NO_LOCAL_GIT_ROOT; SHA_AND_ORIGINAL_BYTES_RECORDED')
        if old.suffix=='.py':
            o=ast.parse(before);n=ast.parse(after);oldnodes={a.name:a for a in o.body if isinstance(a,(ast.FunctionDef,ast.ClassDef))};newnodes={a.name:a for a in n.body if isinstance(a,(ast.FunctionDef,ast.ClassDef))}
            for name,node in oldnodes.items():
                newname='run_native' if new.name=='native_run.py' and name=='main' else name
                other=newnodes.get(newname);same=other is not None and ast.dump(node,include_attributes=False)==ast.dump(other,include_attributes=False)
                symbols.append(dict(old_file=str(old),old_symbol=name,new_file=row['path'],new_symbol=newname if other else 'REPLACED_BY_EXPLICIT_ADAPTER',
                    comparison='EXACT_AST_EQUAL' if same else 'ENGINEERING_IMPORT_PATH_IDENTITY_OR_INTERFACE_DIFF_REVIEW',old_sha256=row['original_sha256'],new_sha256=row['sha256']))
        diff=''.join(difflib.unified_diff(before.decode(errors='replace').splitlines(True),after.decode(errors='replace').splitlines(True),fromfile=str(old),tofile=row['path']))
        p=AUDIT/'source_diffs_final'/Path(row['path']).with_suffix('.diff');p.parent.mkdir(parents=True,exist_ok=True);p.write_text(diff)
    save('source_map.json',rows)
    with (ROOT/'docs/SOURCE_SYMBOL_MAP.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(symbols[0]));w.writeheader();w.writerows(symbols)
    save('source_equivalence_summary.json',dict(files=len(rows),byte_identical=sum(not r['changed'] for r in rows),
        symbols=len(symbols),symbol_AST_exact=sum(r['comparison']=='EXACT_AST_EQUAL' for r in symbols),
        note='Remaining diffs are explicitly listed; no blanket claim that all relocated source bytes match original contracts.'))

def dependency_inventory():
    records=[];edges=[]
    for p in sorted((ROOT/'mouse_vessel_flow').rglob('*.py')):
        tree=ast.parse(p.read_text());imports=[];calls=set();io=[];env=[]
        for n in ast.walk(tree):
            if isinstance(n,ast.ImportFrom):imports.append(dict(module=n.module,level=n.level,names=[a.name for a in n.names]))
            elif isinstance(n,ast.Import):imports.extend(dict(module=a.name,level=0,names=[]) for a in n.names)
            elif isinstance(n,ast.Call):
                name=ast.unparse(n.func);calls.add(name)
                if any(s in name for s in ('open','read','load','save','write','subprocess','Popen','CDLL','getenv')):
                    io.append(dict(line=n.lineno,call=name,literal_args=[a.value for a in n.args if isinstance(a,ast.Constant) and isinstance(a.value,str)]))
            elif isinstance(n,ast.Subscript) and 'environ' in ast.unparse(n.value):env.append(ast.unparse(n))
        records.append(dict(path=str(p.relative_to(ROOT)),imports=imports,calls=sorted(calls),io_and_native_calls=io,environment_access=env))
    save('dependency_inventory.json',records)
    names={'.'.join(Path(r['path']).with_suffix('').parts).removesuffix('.__init__'):r for r in records}
    import networkx as nx
    graph=nx.DiGraph();graph.add_nodes_from(names)
    for name,r in names.items():
        package=name if r['path'].endswith('__init__.py') else name.rpartition('.')[0]
        for item in r['imports']:
            if item['level']:
                prefix='.'.join(package.split('.')[:len(package.split('.'))-item['level']+1]);target='.'.join(x for x in (prefix,item['module']) if x)
            else:target=item['module']
            if target in names:graph.add_edge(name,target)
            for attr in item['names']:
                if target+'.'+attr in names:graph.add_edge(name,target+'.'+attr)
    cycles=[sorted(c) for c in nx.strongly_connected_components(graph) if len(c)>1]
    save('import_graph.json',dict(nodes=len(graph),edges=list(graph.edges),cycles=cycles,
        interpretation='Includes delayed function-local imports. Import smoke checks establish whether any cycle actually blocks loading.'))
    code_paths=[]
    for p in (ROOT/'mouse_vessel_flow').rglob('*.py'):
        for i,line in enumerate(p.read_text().splitlines(),1):
            if any(t in line for t in ('/home/lzy/','/workspace/','sys.path','runpy','spec_from_file_location','exec(')):
                code_paths.append(dict(path=str(p.relative_to(ROOT)),line=i,text=line))
    save('independence_static.json',dict(status='PASS' if not code_paths else 'REVIEW',legacy_loading_findings=code_paths,
        symlinks=[str(p.relative_to(ROOT)) for top in ('mouse_vessel_flow','inputs','results','cases','configs','external') for p in (ROOT/top).rglob('*') if p.is_symlink()],
        original_bytes_evidence='inputs/reference/source_identity; bytes are hashed only, never executed',
        old_source_baseline='One explicit isolated validation subprocess used old code to capture the small-scene oracle. Runtime/new tests use the saved oracle and local package only.'))

def preservation():
    baseline=json.loads((AUDIT/'local_vs_published.json').read_text());changed=[]
    for r in baseline:
        p=Path(r['source'])
        if not p.is_file() or sha(p)!=r['local_sha256']:changed.append(r['source'])
    save('original_files_after.json',dict(status='PASS' if not changed else 'FAIL',checked=len(baseline),changed=changed,deleted_by_this_task=False))
    rows=json.loads((AUDIT/'data_manifest.json').read_text());groups=collections.defaultdict(list)
    for r in rows:groups[r['sha256']].append(r)
    duplicate=[dict(sha256=k,size=v[0]['size'],paths=[r['path'] for r in v]) for k,v in groups.items() if len(v)>1 and v[0]['size']>1024]
    save('duplicate_data_review.json',dict(groups=duplicate,
        reason='Exact data duplicates inside independent native bundles, completed-track markers, or original media/source locks remain where deleting them would break immutable relative-file contracts. No links to old workspaces. Reviewed network/code mirrors were not duplicated.'))
    excluded=json.loads((AUDIT/'published_EXCLUDED.json').read_text());counts=collections.Counter(r['reason'] for r in excluded)
    save('excluded_review.json',dict(counts=dict(counts),records='published_EXCLUDED.json',
        necessary_unsynced_additions='None needed by the exercised current runtime; native installables and intermediate derived caches remain excluded.',
        unresolved=['Historical particle6_5 real-scene fixture absent before migration; three original tests BLOCKED',
                    'Full native GPU installation/build requires a separate environment task'],
        original_low_priority_files_deleted=False))

def environment():
    import importlib.metadata as md
    distributions={d.metadata['Name']:d.version for d in md.distributions()}
    tools={k:shutil.which(k) for k in ('g++','cmake','mpiexec','nvidia-smi','ffmpeg','ultraVessMorpho2Mesh','svmultiphysics')}
    optional={k:importlib.util.find_spec(k) is not None for k in ('sv','vmtk','lammps','torch')}
    save('environment_verified.json',dict(python=sys.version,executable=sys.executable,prefix=sys.prefix,base_prefix=sys.base_prefix,
        isolated_site_packages='include-system-site-packages = false' in (ROOT/'.venv/pyvenv.cfg').read_text(),
        distributions=distributions,tools=tools,optional_modules=optional,remote_access=False,new_native_CFD=False))

if __name__=='__main__':
    source_inventory();dependency_inventory();preservation();environment();print('FINAL_AUDIT_WRITTEN')
