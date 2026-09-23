"""Bounded read-only discovery. Never imports or executes source-project modules."""
from pathlib import Path
import os, json, hashlib, subprocess, collections, re, stat, ast
import yaml
import xml.etree.ElementTree as ET

S=Path('/home/lzy/projects/ulm_3D_vascular')
A=Path(__file__).resolve().parent
def dump(name,value): (A/name).write_text(json.dumps(value,ensure_ascii=False,indent=2)+'\n')
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def cmd(args,input=None):
    return subprocess.run(args,input=input,capture_output=True,check=True,env={**os.environ,'GIT_OPTIONAL_LOCKS':'0'})

# Correct handling of spaces/non-ASCII paths: NUL-delimited porcelain, no filters.
g=['git','--no-optional-locks','-c','filter.lfs.process=','-c','filter.lfs.clean=','-c','filter.lfs.smudge=',
   '-c','filter.lfs.required=false','-C',str(S)]
raw=cmd(g+['status','--porcelain=v1','-z','--untracked-files=all','--ignore-submodules=all']).stdout
(A/'source_git_status_filters_disabled.z').write_bytes(raw)
rows=[]
parts=iter(raw.split(b'\0'))
for b in parts:
    if not b:continue
    row={'status':b[:2].decode(),'path':b[3:].decode()}
    if 'R' in row['status'] or 'C' in row['status']:row['old_path']=next(parts).decode()
    rows.append(row)
index={}
for b in cmd(g+['ls-files','-s','-z']).stdout.split(b'\0'):
    if b:
        h,p=b.split(b'\t',1);mode,oid,stage=h.decode().split();index[p.decode()]=(mode,oid,stage)
ids=list(dict.fromkeys(index[x['path']][1] for x in rows if x['path'] in index and 'M' in x['status']))
checks=cmd(g+['cat-file','--batch-check=%(objectname) %(objecttype) %(objectsize)'],('\n'.join(ids)+'\n').encode()).stdout.decode().splitlines()
small=[x.split()[0] for x in checks if len(x.split())==3 and int(x.split()[2])<1024]
proc=subprocess.Popen(g+['cat-file','--batch'],stdin=subprocess.PIPE,stdout=subprocess.PIPE)
pointers={}
for oid in small:
    proc.stdin.write((oid+'\n').encode());proc.stdin.flush();hdr=proc.stdout.readline().split();data=proc.stdout.read(int(hdr[2]));proc.stdout.read(1)
    if data.startswith(b'version https://git-lfs.github.com/spec/v1\n'):
        match=re.search(rb'oid sha256:([0-9a-f]{64})\nsize ([0-9]+)',data)
        if match:pointers[oid]=(match[1].decode(),int(match[2]))
proc.stdin.close();proc.wait()
same=[];remaining=[]
for row in rows:
    rel=row['path'];p=S/rel
    if row['status']==' M' and rel in index and index[rel][1] in pointers and p.is_file():
        sha,size=pointers[index[rel][1]]
        if p.stat().st_size==size and digest(p)==sha:
            same.append({**row,'sha256':sha,'size_bytes':size});continue
    remaining.append(row)
dump('source_git_worktree_assessment.json',{'native_git_status':'FAILED_GIT_LFS_UNAVAILABLE','source_is_dirty':bool(remaining),
    'raw_filters_disabled_entries':len(rows),'verified_materialized_lfs_false_positives':len(same),
    'remaining_status_entries':len(remaining),'remaining_status_counts':dict(collections.Counter(x['status'] for x in remaining)),
    'method':'NUL-delimited porcelain with per-command LFS filters disabled; hydrated file size and SHA-256 compared to index LFS pointers. No installation, index refresh, config change, clean/smudge execution or LFS object writes.',
    'actual_or_unresolved_changes':remaining,'lfs_materialized_files_identical_to_index':same,
    'classification':'CONFIRMED_FROM_DATA','parent_status_ignores_submodule_worktrees':True})
(A/'source_git_status_readonly_assessed.txt').write_text('\n'.join(x['status']+' '+x['path'] for x in remaining)+'\n')
sub=[]
for rel in ['Ultraliser','external_reference/LBPM']:
    p=S/rel
    for what,args in [('HEAD',['rev-parse','HEAD']),('remote',['remote','-v']),('status',['status','--short','--untracked-files=no'])]:
        r=subprocess.run(['git','--no-optional-locks','-C',str(p),*args],capture_output=True,text=True,env={**os.environ,'GIT_OPTIONAL_LOCKS':'0'})
        sub.append({'path':rel,'query':what,'exit_code':r.returncode,'stdout':r.stdout,'stderr':r.stderr})
dump('nested_repository_identity.json',sub)

# All authored configuration is small; parse safely and record evidence, never apply it.
configs=[]
for p in sorted((S/'configs').glob('*.yaml'))+[S/'cfd_lumen_config.yaml']:
    data=yaml.safe_load(p.read_text());flat=[]
    def flatten(x,prefix=''):
        if isinstance(x,dict):
            for k,v in x.items():flatten(v,prefix+'.'+str(k) if prefix else str(k))
        else:flat.append({'key':prefix,'value':x})
    flatten(data)
    for r in flat:
        if isinstance(r['value'],str) and (r['key'].startswith('paths.') or r['key'].endswith(('executable','environment_python','runtime_prefix','official_repository','mpi_launcher'))):
            v=r['value']
            r['path_kind']='WINDOWS_PATH' if re.match(r'^[A-Za-z]:',v) else 'POSIX_ABSOLUTE' if v.startswith('/') else 'PROJECT_RELATIVE'
            r['exists_as_current_wsl_path']=False if r['path_kind']=='WINDOWS_PATH' else (Path(v) if v.startswith('/') else S/v).exists()
    configs.append({'path':str(p.relative_to(S)),'sha256':digest(p),'classification':'CONFIRMED_FROM_CONFIG','flattened':flat})
dump('configuration_inventory.json',configs)

# Preserve a bounded map of each important family; no recursive raw filename dump to the UI.
maps={}
for rel in ['utils','configs','docs','tools','patches','test_data','vessel_model','outputs']:
    root=S/rel
    maps[rel]=[str(p.relative_to(S)) for p in sorted(root.glob('*'))]
dump('important_directory_map.json',maps)

formats={'.stl','.obj','.vtk','.vtp','.vtu','.vti','.ply','.msh','.nii','.nii.gz','.dcm','.tif','.tiff','.png','.jpg','.jpeg','.csv','.json','.yaml','.yml','.npy','.npz','.h5','.hdf5','.mat','.swc','.pkl','.graphml','.lua','.lsb','.res'}
data_rows=[];groups=collections.defaultdict(lambda:collections.Counter());examples=collections.defaultdict(list)
for base,ds,fs in os.walk(S):
    ds[:]=[n for n in ds if n not in ('.git','__pycache__','.pytest_cache','.ruff_cache')]
    for n in fs:
        p=Path(base)/n;rel=str(p.relative_to(S));ext='.nii.gz' if n.endswith('.nii.gz') else p.suffix.lower()
        if ext not in formats:continue
        if 'BVLab-Annotation/_internal/' in rel:family='bundled_annotation_runtime';role='Bundled third-party runtime/example data';migration='DO_NOT_MIGRATE'
        elif rel.startswith('vessel_model/'):
            family='/'.join(Path(rel).parts[:2]);role='source dataset candidate; exact purpose needs catalog/header evidence';migration='DATA_ONLY'
        elif rel.startswith('outputs/'):
            family='/'.join(Path(rel).parts[:2]);role='derived scientific output/evidence; per-run classification required';migration='REFERENCE_ONLY' if '/cfd_flow/' in rel else 'DATA_ONLY'
        elif rel.startswith('test_data/'):
            family='test_data';role='minimal regression/geometry input fixture';migration='DATA_ONLY'
        elif rel.startswith(('.codex_tmp/','tmp/')):
            family='scratch';role='temporary visual/debug artifact';migration='DO_NOT_MIGRATE'
        elif rel.startswith(('Ultraliser/','external_reference/')):
            family='/'.join(Path(rel).parts[:2]);role='external geometry tool or old solver reference asset';migration='REFERENCE_ONLY'
        else:family=Path(rel).parts[0];role='configuration, documentation or presentation asset';migration='REFERENCE_ONLY'
        size=p.stat().st_size
        row={'path':str(p),'relative_path':rel,'format':ext,'size_bytes':size,'likely_role':role,'role_confidence':'INFERRED',
            'data_existence_confidence':'CONFIRMED_FROM_DATA','family':family,'referenced_by_code':'See data_family_profiles.json; family-level evidence, not a claim of individual use',
            'referenced_by_config':[],'unit_evidence':'UNIT_UNVERIFIED unless covered by representative headers/config and generating code',
            'coordinate_evidence':'COORDINATE_SYSTEM_UNVERIFIED unless covered by a specific transform record',
            'migration_value':migration}
        for c in configs:
            for x in c['flattened']:
                if isinstance(x['value'],str) and len(x['value'])>10 and (rel==x['value'] or rel.startswith(x['value'].rstrip('/')+'/')):
                    row['referenced_by_config'].append(c['path']+':'+x['key'])
        data_rows.append(row);groups[family][ext]+=1;groups[family]['bytes']+=size
        if len(examples[(family,ext)])<2:examples[(family,ext)].append(rel)
dump('data_file_inventory.json',data_rows)
dump('data_family_summary.json',[{'family':f,'format_counts':dict(v),'examples':{e:xs for (g,e),xs in examples.items() if g==f}} for f,v in sorted(groups.items())])

# Dataset-side physical metadata, read without loading image stacks.
xmlrecords=[]
for p in sorted((S/'vessel_model/T - NNE2/hana_stk').rglob('*.xml')):
    values={}
    for element in ET.parse(p).getroot().iter('Key'):
        k=element.attrib.get('key')
        if k in ('micronsPerPixel_XAxis','micronsPerPixel_YAxis','positionCurrent_ZAxis') and k not in values:values[k]=element.attrib.get('value')
    xmlrecords.append({'path':str(p.relative_to(S)),'sha256':digest(p),'metadata':values,'classification':'CONFIRMED_FROM_DATA',
        'scope':'NNE2 acquisition metadata only; does not establish active fMOST sample calibration'})
dump('nne2_acquisition_metadata.json',xmlrecords)

# Explicit authored-text searches; regex hits are leads, not automatically factual claims.
textfiles=list(S.glob('*.py'))+list((S/'utils').rglob('*.py'))+list((S/'configs').glob('*.yaml'))+list((S/'docs').glob('*.md'))+list((S/'references').glob('*.md'))
categories={'units':r'\b(unit|units|spacing|voxel|resolution|scale|pixel|radius|diameter|velocity|flow_rate|pressure|density|viscosity)\b|_um|_m2_s|_pa',
    'coordinates':r'rotation|translation|transform|affine|origin|lps|ras|spacing_xyz|shape_zyx|flip|swapaxis|transpose',
    'boundaries':r'inlet|outlet|wall|boundary_id|face_id|CellEntityIds|CUT_PORT|TRUE_TERMINAL|port_id',
    'experimental':r'PDMS|microfluid|实验|experimental|measured_flow_direction',
    'quality':r'watertight|non.?manifold|duplicate|degenerate|self.intersect|normal|bounds|remesh|repair|voxel'}
for cat,regex in categories.items():
    rx=re.compile(regex,re.I);hits=[]
    for p in textfiles:
        for i,line in enumerate(p.read_text().splitlines(),1):
            if rx.search(line):hits.append({'file':str(p.relative_to(S)),'line':i,'text':line[:1200]})
    dump('search_'+cat+'.json',hits)
dump('authored_text_source_hashes.json',[{'path':str(p.relative_to(S)),'size_bytes':p.stat().st_size,'sha256':digest(p)} for p in sorted(set(textfiles+[S/'README.md',S/'requirements.txt',S/'.gitignore',S/'.gitattributes']))])
print(json.dumps({'git_raw':len(rows),'lfs_identical':len(same),'remaining_git_changes':len(remaining),'git_change_counts':dict(collections.Counter(x['status'] for x in remaining)),
    'data_files':len(data_rows),'data_families':len(groups),'nne2_xml_files':len(xmlrecords),'configurations':len(configs)},indent=2))
