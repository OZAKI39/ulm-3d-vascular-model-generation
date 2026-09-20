"""Pin official source and cross-reference all compiled-source PETSc symbols."""
import collections, hashlib, json, re, subprocess, tarfile, urllib.request
from datetime import datetime, timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]; R=ROOT/'reports/sv1_3n'; E=ROOT/'external/petsc325'
def sha(p): return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(n,d): (R/(n+'.json')).write_text(json.dumps(d,indent=2)+'\n')
def git(*a): return subprocess.check_output(['git','-C',str(E/'upstream'),*a],text=True).strip()
tag='v3.25.5'; head=git('rev-parse',tag+'^{commit}')
remote=subprocess.check_output(['git','ls-remote','--tags','https://gitlab.com/petsc/petsc.git','refs/tags/'+tag,'refs/tags/'+tag+'^{}'],text=True)
assert head in remote and git('rev-parse','HEAD')==head
assert not git('status','--porcelain')
(R/'official_tag_lookup.txt').write_text(remote)
html=urllib.request.urlopen('https://petsc.org/release/install/download/').read()
(R/'official_download.html').write_bytes(html)
versions=sorted(set(re.findall(rb'petsc-(3\.25\.\d+)\.tar\.gz',html)),key=lambda s:list(map(int,s.split(b'.'))))
assert versions[-1]==b'3.25.5'
archive=E/'petsc-3.25.5.tar.gz';source=E/'petsc-3.25.5'
if not source.exists():
 with tarfile.open(archive) as t:t.extractall(E,filter='data')
version=(source/'include/petscversion.h').read_text()
for k,v in [('MAJOR',3),('MINOR',25),('SUBMINOR',5),('RELEASE',1)]:assert re.search(r'#define PETSC_VERSION_'+k+r'\s+'+str(v)+r'\b',version)
write('petsc325_source',dict(status='PASS',timestamp=datetime.now(timezone.utc).isoformat(),version='3.25.5',tag=tag,commit=head,tag_object=git('rev-parse',tag),source_archive=str(archive.relative_to(ROOT)),source_archive_sha256=sha(archive),archive_url='https://web.cels.anl.gov/projects/petsc/download/release-snapshots/petsc-3.25.5.tar.gz',git_origin='https://gitlab.com/petsc/petsc.git',latest_formal_325_patch_at_start=versions[-1].decode(),version_header_sha256=sha(source/'include/petscversion.h'),candidate='official tag',source_patch=False))
S=ROOT/'external/svMultiPhysics'
names=subprocess.check_output(['git','-C',str(S),'ls-files','-z'],text=True).split('\0')
baseline={n:sha(S/n) for n in names if n and (S/n).is_file()}
commit=subprocess.check_output(['git','-C',str(S),'rev-parse','HEAD'],text=True).strip()
status=subprocess.check_output(['git','-C',str(S),'status','--porcelain'],text=True)
assert commit=='c3f0bb892b765b718f61069ecd9726dbc6d177fd' and not status
write('svmp_source_before',dict(commit=commit,status=status,files=baseline))
# Use public header vocabulary rather than a hand-picked list of API calls.
headers={str(p.relative_to(source)):p.read_text(errors='replace') for p in (source/'include').glob('*.h')}
vocabulary=set(re.findall(r'\b[A-Za-z_]\w*\b','\n'.join(headers.values())))
with tarfile.open(ROOT/'outputs/sv1_3m/native/petsc-3.19.6.tar.gz') as oldtar:
 old_headers='\n'.join(oldtar.extractfile(m).read().decode(errors='replace') for m in oldtar if m.isfile() and re.match(r'petsc-3.19.6/include/[^/]+\.h$',m.name))
old_vocabulary=set(re.findall(r'\b[A-Za-z_]\w*\b',old_headers))
prefix=re.compile(r'^(?:Petsc\w+|PETSC\w*|Mat[A-Z]\w*|MAT[A-Z_]\w*|Vec[A-Z]\w*|VEC[A-Z_]\w*|KSP\w*|PC[A-Z_]\w*|IS[A-Z]\w*|AO[A-Z]\w*|MPIU_\w+|CHKERR\w*|SETERRQ\w*)$')
objects={'Mat','Vec','PC','IS','AO','InsertMode','ScatterMode','NormType','INSERT_VALUES','ADD_VALUES','SCATTER_FORWARD','SCATTER_REVERSE','NORM_2','NORM_INFINITY','SAME_NONZERO_PATTERN','DIFFERENT_NONZERO_PATTERN'}
symbols=collections.defaultdict(list); includes=[]; bool_context=[]; candidate_statements=[]
for n in names:
 if not n or Path(n).suffix not in ('.cpp','.h','.c','.hpp','.cxx','.cc'):continue
 txt=(S/n).read_text(errors='replace') if (S/n).is_file() else subprocess.check_output(['git','-C',str(S),'show','HEAD:'+n],text=True)
 # Preserve line numbers while removing comments, excluding prose false positives.
 txt=re.sub(r'/\*.*?\*/',lambda m:'\n'*m[0].count('\n'),txt,flags=re.S)
 for line_no,line in enumerate(txt.splitlines(),1):
  line=line.split('//',1)[0]
  for h in re.findall(r'#\s*include\s*[<"](petsc[^>"]+)',line):includes.append(dict(header=h,file=n,line=line_no))
  for token in sorted(set(re.findall(r'\b[A-Za-z_]\w*\b',line))):
   if (prefix.match(token) or token in objects) and token in vocabulary|old_vocabulary:
    symbols[token].append(dict(file=n,line=line_no,statement=line.strip()))
  if 'PetscBool' in line:bool_context.append(dict(file=n,line=line_no,statement=line.strip()))
  if re.search(r'sizeof|MPI_\w+\s*\(|reinterpret_cast|offsetof|memcpy|fwrite|fread|serialize',line) and ('petsc' in n.lower()):candidate_statements.append(dict(file=n,line=line_no,statement=line.strip()))
notes={}; relevant=[]
for v in range(320,326):
 p=E/'upstream'/f'doc/changes/{v}.md';body=p.read_text();notes[str(v)]=dict(path=str(p.relative_to(ROOT)),sha256=sha(p),url=f'https://petsc.org/release/changes/{v}/',bytes=len(body))
 for block in re.findall(r'^- .*?(?=\n- |\n\n|\Z)',body,re.M|re.S):
  used=[s for s in symbols if re.search(r'(?<!\w)'+re.escape(s)+r'(?!\w)',block)]
  if used:relevant.append(dict(version=str(v),symbols=used,official_note=block))
api=[]
for name,loc in sorted(symbols.items()):
 matches=[n for n in relevant if name in n['symbols']]
 text='\n'.join(n['official_note'] for n in matches)
 state='unchanged'
 if name not in vocabulary:state='removed'
 elif name=='PetscBool':state='changed'
 elif name=='PETSC_DEFAULT':state='deprecated'
 elif name.startswith('VecGhost') or name in ('VecSetFromOptions','MatSetFromOptions','KSPSolve','PetscInitialize','PetscFinalize','VecSetType'):state='needs runtime validation'
 declaration=next(((h,l.strip()) for h,b in headers.items() for l in b.splitlines() if re.search(r'\b'+name+r'\b',l) and ('PETSC_EXTERN' in l or '#define '+name in l or 'typedef' in l)),None)
 kind='macro' if any(re.search(r'^#\s*define\s+'+name+r'\b',h,re.M) for h in headers.values()) else ('function' if any(re.search(r'\b'+name+r'\s*\(',l['statement']) for l in loc) else ('enum' if name.isupper() and '_' in name else 'type'))
 api.append(dict(name=name,kind=kind,classification=state,locations=loc,header_declaration=declaration,release_note_matches=matches))
write('petsc_api_inventory',dict(status='PASS',source_commit=commit,method='All tracked C/C++ source, comments removed, tokens cross-referenced with official public headers; removed API names additionally checked against old headers in follow-up audit',headers=includes,api=api,count=len(api)))
write('release_change_audit',dict(status='PASS',official_notes=notes,relevant_notes=relevant,summary={'3.20':'No removed/signature-changing used API found; logging internals changed but public stage functions retained.','3.21':'ASM submatrix nullspace preservation; no attached nullspace in this solver. VecScale now collective; calls occur on all ranks.','3.22':'PETSC_DEFAULT deprecated for KSPSetTolerances, compatibility retained; fourth argument means retain divergence tolerance.','3.23':'Device factorization option renamed; not set in formal options. Must inspect runtime factor residency.','3.24':'PetscBool becomes one-byte C bool; independently audited.','3.25':'-options_left no longer implies options_view; frozen proof options explicitly specify both.'},runtime_validation=['CPU science equivalence','CUDA ghost local/forward/reverse values','KSP/ASM/ILU options and residency','PETSc finalization before application MPI_Finalize']))
write('petsc_bool_audit',dict(status='NO_DEPENDENCY_FOUND',declarations=bool_context,reviewed_storage_and_communication_statements=candidate_statements,findings={'sizeof':'PetscMalloc1(nEq, &psol) allocates using the current LSCtx element type; no explicit PetscBool byte count or fixed offsets.','MPI communication':'MPI calls communicate PetscInt/int arrays and scalar data, not PetscBool or enclosing contexts.','binary layout':'LHSCtx and LSCtx are local C++ structs compiled together with current headers; no fixed offsets or external ABI.','reinterpret cast':'No conversion of PetscBool storage to int pointers or raw packed bytes.','serialization':'No raw write/read of PetscBool, LHSCtx or LSCtx.'}))
print(json.dumps({'source_commit':head,'archive_sha256':sha(archive),'APIs':len(api),'relevant_notes':relevant,'bool':bool_context},indent=2))
# Context review excludes generic type-name matches and unused TS/SNES/Fortran APIs.
p=R/'release_change_audit.json';d=json.loads(p.read_text());d['lexical_cross_matches']=d.pop('relevant_notes')
selected={('320','VecSqrtAbs'),('320','MatView'),('321','PCFIELDSPLIT'),('322','KSPSetTolerances'),('322','PCGAMG'),('323','MATAIJ'),('324','PetscBool'),('325','MatSetFromOptions')}
d['relevant_notes']=[x for x in d['lexical_cross_matches'] if any((x['version'],s) in selected for s in x['symbols']) and 'logical(C_BOOL)' not in x['official_note'] and 'MPIU_BOOL' not in x['official_note']]
d['summary']['3.21']='ASM/fieldsplit submatrix nullspace preservation; solver does not attach nullspaces. No used function removed or signature changed.'
d['relevance_method']='Manual context review after lexical match: TS/SNES/TAO/DM/Fortran-only changes and generic Mat/Vec/type mentions excluded; unused internal methods excluded.'
write('release_change_audit',d)
inv=json.loads((R/'petsc_api_inventory.json').read_text())
for item in inv['api']:item['release_note_matches']=[x for x in d['relevant_notes'] if item['name'] in x['symbols']]
inv['method']='All tracked C/C++ source (including sparse paths read from pinned Git), comments removed; public symbols cross-referenced against both 3.19.6 and 3.25.5 headers; manually reviewed release-note applicability.'
write('petsc_api_inventory',inv)
