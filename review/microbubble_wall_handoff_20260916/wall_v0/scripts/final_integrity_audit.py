from pathlib import Path
import hashlib,json,subprocess,sys,shutil
R=Path(__file__).resolve().parents[1];P=json.loads((R/'provenance/TASK_PATHS.json').read_text());W=Path(P['local_work'])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
before=json.loads((R/'provenance/LOCAL_BASELINE_BEFORE.json').read_text());after=[]
for b in before:
 root=Path(b['root']);p=subprocess.run(['sha256sum','--check','--quiet','SHA256SUMS'],cwd=root,capture_output=True,text=True);item=dict(root=str(root),returncode=p.returncode,stdout=p.stdout,stderr=p.stderr,manifest_sha256=sha(root/'SHA256SUMS'));assert p.returncode==0 and item['manifest_sha256']==b['manifest_sha256'],item;after.append(item)
(R/'provenance/LOCAL_BASELINE_AFTER.json').write_text(json.dumps(after,indent=2)+'\n')
ref=Path(P['reference_work'])/'upstream/RigidMultiblobsWall';g=subprocess.run(['git','rev-parse','HEAD'],cwd=ref,capture_output=True,text=True);assert g.stdout.strip()=='8d41e464d7de9b6514a85a741dd227f1219e3a49';status=subprocess.run(['git','diff','--stat'],cwd=ref,capture_output=True,text=True);assert status.returncode==0 and not status.stdout
# Previously accepted tracked modifications, if any, must be considered explicitly.
tracked=subprocess.run(['git','status','--porcelain','--untracked-files=no'],cwd=ref,capture_output=True,text=True);assert not tracked.stdout
(R/'provenance/RMBW_SOURCE_AFTER.json').write_text(json.dumps({'commit':g.stdout.strip(),'tracked_changes':tracked.stdout,'status':'PASS','production_RMBW_source_copy':False},indent=2)+'\n')
manifest=R/'REMOTE_SHA256SUMS';p=subprocess.run(['sha256sum','--check','--quiet',manifest.name],cwd=R,capture_output=True,text=True);assert p.returncode==0,(p.stdout,p.stderr)
(R/'validation/REMOTE_TO_WSL_INTEGRITY.json').write_text(json.dumps({'status':'PASS','remote_manifest_sha256':sha(manifest),'files_verified':len(manifest.read_text().splitlines()),'method':'actual sha256sum --check on downloaded files in original WSL','command_returncode':p.returncode},indent=2)+'\n')
(R/'provenance/control_scripts').mkdir(exist_ok=True)
for path in W.glob('*.py'):shutil.copy2(path,R/'provenance/control_scripts'/path.name)
versions=subprocess.run([sys.executable,'-m','pip','freeze'],capture_output=True,text=True,check=True);(R/'provenance/LOCAL_VALIDATION_REQUIREMENTS.txt').write_text(versions.stdout)
print('INTEGRITY_PASS',len(after),'baselines;',len(manifest.read_text().splitlines()),'remote files')
