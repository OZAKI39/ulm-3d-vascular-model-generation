from pathlib import Path
import sys,subprocess,hashlib,json,shlex
C=Path(__file__).resolve().parent;sys.path.insert(0,str(C));from remote_exec import run,SSH
S=json.loads((C/'TASK_PATHS.json').read_text());R=Path(S['local_result']);label=sys.argv[1] if len(sys.argv)>1 else 'INITIAL'
script='''from pathlib import Path
import hashlib,subprocess,json,shutil
R=Path(%r)
for p in R.rglob('__pycache__'):shutil.rmtree(p)
files=[p for p in sorted(R.rglob('*')) if p.is_file() and p.name!='SHA256SUMS']
# Exclude only the root manifest itself, not historical manifests with the same basename.
files=[p for p in sorted(R.rglob('*')) if p.is_file() and p!=R/'SHA256SUMS']
text=''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+str(p.relative_to(R))+'\\n' for p in files)
(R/'SHA256SUMS').write_text(text)
p=subprocess.run(['sha256sum','--check','--quiet','SHA256SUMS'],cwd=R,capture_output=True,text=True)
assert p.returncode==0,p.stdout+p.stderr
print(json.dumps({'status':'PASS','files':len(files),'bytes':sum(p.stat().st_size for p in files),'manifest_sha256':hashlib.sha256((R/'SHA256SUMS').read_bytes()).hexdigest()}))
'''%S['remote_result']
p=run(script,timeout=150);assert p.returncode==0,p.stdout+p.stderr;remote=json.loads(p.stdout);(C/(label+'_REMOTE_SHA_RECEIPT.json')).write_text(json.dumps(remote,indent=2)+'\n')
p=subprocess.run(['rsync','-az','-e',shlex.join(SSH[:-1]),SSH[-1]+':'+S['remote_result']+'/SHA256SUMS',str(R/'SHA256SUMS')],capture_output=True,text=True);p.check_returncode()
mismatches=[];count=0
for line in (R/'SHA256SUMS').read_text().splitlines():
 h,n=line.split('  ',1);f=R/n
 if not f.is_file() or hashlib.sha256(f.read_bytes()).hexdigest()!=h:mismatches.append(n)
 count+=1
listed={line.split('  ',1)[1] for line in (R/'SHA256SUMS').read_text().splitlines()};extra=[str(p.relative_to(R)) for p in R.rglob('*') if p.is_file() and p!=R/'SHA256SUMS' and str(p.relative_to(R)) not in listed]
local={'status':'PASS' if not mismatches and not extra else 'FAIL','files':count,'mismatches':mismatches,'extra_local_files':extra,'manifest_sha256':hashlib.sha256((R/'SHA256SUMS').read_bytes()).hexdigest()};assert remote['manifest_sha256']==local['manifest_sha256'];(C/(label+'_LOCAL_SHA_RECEIPT.json')).write_text(json.dumps(local,indent=2)+'\n');print(json.dumps({'remote':remote,'local':local},indent=2));assert local['status']=='PASS'
