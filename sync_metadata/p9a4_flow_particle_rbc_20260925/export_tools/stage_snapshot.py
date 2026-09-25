from pathlib import Path
import os,json,hashlib,subprocess,shutil,contextlib,io,importlib.util
A=Path(__file__).resolve().parents[1];C=json.loads((A/'context.json').read_text());D=Path(C['destination']);M=D/'sync_metadata'/C['snapshot']
shutil.copy2(Path(__file__),M/'export_tools'/Path(__file__).name)
# Run frozen scientific-byte protection once and include its independent receipt.
spec=importlib.util.spec_from_file_location('inherited',D/'sync_metadata/network_h0_particle_20260925/verify_snapshot.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
s=io.StringIO()
with contextlib.redirect_stdout(s):m.scientific_inputs()
(M/'validation/scientific_input_verification.log').write_text(s.getvalue())
manifest=M/'SNAPSHOT_SHA256.txt'
files=sorted(p for p in D.rglob('*') if p.is_file() and p.name!='.git' and p!=manifest)
for p in files:
 rel=p.relative_to(D)
 assert not p.is_symlink(),str(rel)
 assert not(set(rel.parts)&{'.git','.venv','__pycache__','.pytest_cache'}),str(rel)
 assert '\n' not in str(rel) and '\r' not in str(rel)
 assert p.stat().st_size<50*1024**2,str(rel)
with manifest.open('w') as f:
 for p in files:
  with p.open('rb') as src:h=hashlib.file_digest(src,'sha256').hexdigest()
  f.write(h+'  '+str(p.relative_to(D))+'\n')
files.append(manifest)
listpath=A/'audit/publish_paths.nul';listpath.write_bytes(b'\0'.join(str(p.relative_to(D)).encode() for p in files)+b'\0')
subprocess.run(['git','add','-u'],cwd=D,check=True)
subprocess.run(['git','add','-f','--pathspec-from-file='+str(listpath),'--pathspec-file-nul'],cwd=D,check=True)
indexed=[p for p in subprocess.check_output(['git','ls-files','-z'],cwd=D).decode().split('\0') if p]
assert set(indexed)=={str(p.relative_to(D)) for p in files}
changes=subprocess.check_output(['git','diff','--cached','--name-status','-z'],cwd=D).decode().split('\0')
summary=dict(published_files=len(files),manifest_entries=len(files)-1,published_bytes=sum(p.stat().st_size for p in files),changed_names_and_statuses=changes)
(A/'audit/staged_snapshot.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k!='changed_names_and_statuses'},indent=2))
