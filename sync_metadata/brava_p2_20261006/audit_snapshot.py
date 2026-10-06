from pathlib import Path
import os,json,hashlib,re,subprocess,ast
ROOT=Path(__file__).resolve().parents[2];META=Path(__file__).resolve().parent
os.chdir(ROOT)
def sh(*args):return subprocess.check_output(['git',*args],text=True)
def digest(p):
 if p.is_symlink():return hashlib.sha256(os.fsencode(os.readlink(p))).hexdigest()
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
manifest=json.loads((META/'INCLUDED_FILES.json').read_text());checks=[]
def check(n,v):checks.append(dict(check=n,PASS=bool(v)))
for x in manifest['files']:check('file_sha:'+x['path'],digest(ROOT/x['path'])==x['sha256'])
changed=[p for p in sh('diff','--cached','--name-only',manifest['base_GitHub_commit']).splitlines() if p]
check('related_paths_only',all(p=='README.md' or p.startswith(('brava_flow_roi_18mlmin/','ulm_particle_formal_p9a5/particle_3d/src/','ulm_particle_formal_p9a5/particle_3d/tests/','sync_metadata/brava_p2_20261006/')) for p in changed))
check('no_deletions',not sh('diff','--cached','--diff-filter=D','--name-only').strip())
check('root_attributes_unchanged',(ROOT/'.gitattributes').read_bytes()==subprocess.check_output(['git','show','HEAD:.gitattributes']))
check('no_redundant_CSV_in_index',not sh('ls-files','--','brava_flow_roi_18mlmin/population_v1/analysis/stabilization_local_conservation_review_20261006/STABILIZATION_TERMS.csv').strip())
patterns={
 'private_key':rb'-----BEGIN (?:OPENSSH |RSA |EC |DSA |ENCRYPTED )?PRIVATE KEY-----',
 'github_token':rb'(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{50,})',
 'openai_key':rb'\bsk-(?:proj-)?[A-Za-z0-9_-]{32,}',
 'AWS_access_id':rb'\b(?:AKIA|ASIA)[A-Z0-9]{16}\b',
 'credential_in_URL':rb'https?://[^\s/:]+:[^\s/@]{6,}@',
 'bearer_token':rb'(?i)authorization\s*[:=]\s*["\']?Bearer\s+[A-Za-z0-9._-]{24,}',
}
findings=[];text_count=0;syntax=0;jsons=0;outgoing=[]
for name in changed:
 p=ROOT/name
 if not p.exists() and not p.is_symlink():continue
 size=len(os.fsencode(os.readlink(p))) if p.is_symlink() else p.stat().st_size
 outgoing.append(dict(path=name,bytes=size));check('new_file_under_25MB:'+name,size<25_000_000)
 data=os.fsencode(os.readlink(p)) if p.is_symlink() else p.read_bytes()
 if b'\0' in data[:8192]:continue
 try:txt=data.decode('utf-8')
 except UnicodeError:continue
 text_count+=1
 for label,pat in patterns.items():
  if re.search(pat,data):findings.append(dict(path=name,rule=label))
 if p.suffix=='.py':ast.parse(txt,filename=name);syntax+=1
 if p.suffix=='.json':json.loads(txt);jsons+=1
check('credential_pattern_scan',not findings)
# Validate frozen per-round output identity without executing any scientific code.
round_manifests=[]
for d in (ROOT/'brava_flow_roi_18mlmin/population_v1/analysis').glob('*'):
 f=d/'OUTPUT_SHA256.json'
 if not f.exists():continue
 content=json.loads(f.read_text());content=content.get('files',content);verified=0;omitted=[]
 for name,expected in content.items():
  p=d/name
  if not p.exists():
   omitted.append(name);check('declared_omission:'+str(p.relative_to(ROOT)),d.name=='stabilization_local_conservation_review_20261006' and name=='STABILIZATION_TERMS.csv');continue
  # Output manifests hash file contents, following symlinks.
  with p.open('rb') as stream:h=hashlib.file_digest(stream,'sha256').hexdigest()
  check('output_identity:'+str(p.relative_to(ROOT)),h==expected);verified+=1
 round_manifests.append(dict(round=d.name,verified=verified,declared_omissions=omitted))
latest=ROOT/'brava_flow_roi_18mlmin/population_v1/analysis/p2_linear_algebra_rerun_20261006'
check('latest_scientific_gate_preserved',json.loads((latest/'P2_LINEAR_ALGEBRA_DECISION_GATE.json').read_text())['NEXT_RECOMMENDATION']=='STOP_NATIVE_P2P1_ROUTE_RESOURCE_LIMIT')
check('source_project_unchanged',not subprocess.check_output(['git','status','--porcelain'],cwd='/home/lzy/projects/computation_examples/brava_flow_roi_18mlmin',text=True).strip())
result=dict(status='PASS' if all(x['PASS'] for x in checks) else 'FAIL',files_byte_verified=len(manifest['files']),checks=len(checks),passed=sum(x['PASS'] for x in checks),failures=[x for x in checks if not x['PASS']],changed_files=len(outgoing),changed_file_bytes=sum(x['bytes'] for x in outgoing),largest_changed_files=sorted(outgoing,key=lambda x:-x['bytes'])[:15],credential_scan=dict(text_files_scanned=text_count,findings=findings,scope='New/changed UTF-8 files; rule-based scan, not a universal guarantee'),syntax_python_files=syntax,parsed_JSON_files=jsons,round_output_manifests=round_manifests,new_scientific_runs=0,main_or_other_branch_modified=False)
(META/'SYNC_VALIDATION.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['largest_changed_files','round_output_manifests']},indent=2));assert result['status']=='PASS'
